
import os

# ===
# LOAD CONFIGURATIONS
# ===

pairsFile = config["sampleInfo"]
outDir = config["outputDir"]
genomes = list(config["genomes"].keys())

# ===
# READ IN SAMPLE INFORMATION (SINGLE-END)
# ===

sampleToR1 = {}
samples = []

with open(pairsFile) as f:
    for lineNum, line in enumerate(f, start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")
    
        sample, r1 = parts
        samples.append(sample)
        sampleToR1[sample] = r1

def getRead1(wildcards):
    return sampleToR1[wildcards.sample]

def getGenomeFasta(wildcards):
    return config["genomes"][wildcards.genome]["fasta"]


# ===
# RULE ALL
# ===

rule all:
    input:
        expand(
            f"{outDir}/rmsScores/{{genome}}/{{sample}}.rms.txt",
            genome=genomes,
            sample=samples
        )


# ===
# PIPELINE
# ===
rule index_genome:
    """
    Generate bowtie2 index files required for mapping
    Input: genome fasta file
    Output: bowtie2 index files
    """
    input:
        fasta = getGenomeFasta
    output:
        index = f"{outDir}/bowtieIndex/{{genome}}.1.bt2"
    params:
        indexPrefix = f"{outDir}/bowtieIndex/{{genome}}"
    conda: "envs/RMS.yaml"
    log: f"{outDir}/logs/index_genome_{{genome}}.log"
    shell:
        """
        bowtie2-build \
        --threads 8 \
        {input.fasta} \
        {params.indexPrefix} \
        2>> {log} 
        """

rule trim_adapters:
    """
    Uses trim galore to detect & trim adapter sequences (READ 1 ONLY, single-end)
    Input: R1 .fastq(.gz)
    Output: trimmed R1 .fq.gz
    """
    input:
        fq1 = getRead1
    output:
        trimmed1 = f"{outDir}/trimmedFastq/{{sample}}_trimmed.fq.gz"
    params:
        outputDir = f"{outDir}/trimmedFastq"
    conda: "envs/RMS.yaml"
    log:
        f"{outDir}/logs/trim_adapters_{{sample}}.log"
    threads: 8
    shell:
        r"""
        trim_galore \
          --cores {threads} \
          --gzip \
          --output_dir {params.outputDir} \
          --basename {wildcards.sample} \
          {input.fq1} 2>> {log}
        """

rule map:
    """
    1. Maps trimmed single-end fastq files to designated genomes.
    2. Converts bowtie output to bed file
    
    Input: trimmed fastq file and bowtie2 index files
    Output: sam file of mapped reads for each genome
    """
    input:
        trimmed = f"{outDir}/trimmedFastq/{{sample}}_trimmed.fq.gz",
        index = f"{outDir}/bowtieIndex/{{genome}}.1.bt2"
    output:
        bed = f"{outDir}/mapped/{{genome}}/{{sample}}.bed"
    params:
        indexPrefix = f"{outDir}/bowtieIndex/{{genome}}"
    conda: "envs/RMS.yaml"
    log:
        f"{outDir}/logs/map_{{sample}}_{{genome}}.log"
    threads: 8
    shell:
        """
        bowtie2 \
        -p {threads} \
        --end-to-end \
        -k 1 \
        -x {params.indexPrefix} \
        -U {input.trimmed} 2>> {log} | \
        samtools view -@ {threads} -b 2>> {log} | \
        samtools sort -@ {threads} 2>> {log} | \
        bedtools bamtobed -i stdin > {output.bed} 2>> {log}
        """

rule count_ends:
    """
    Count the number of 3' read ends at each genomic position
    Handles both strands:
      + strand: 3' end = $3 (end position)
      - strand: 3' end = $2 (start position)
    """
    input:
        bed = f"{outDir}/mapped/{{genome}}/{{sample}}.bed"
    output:
        count3 = f"{outDir}/endCounts/{{genome}}/{{sample}}.3count.txt",
    shell:
        """
        awk '$6 == "+" {{print $1, $3}} $6 == "-" {{print $1, $2}}' {input.bed} | sort | uniq -c | \
        awk '{{print $2, $3, $1}}' | \
        sort -k1,1 -k2,2n > {output.count3}
        """

rule genome_sizes:
    input:
        fasta = getGenomeFasta
    output:
        sizes = f"{outDir}/genomeSizes/{{genome}}.sizes"
    conda: "envs/RMS.yaml"
    log: f"{outDir}/logs/genome_sizes_{{genome}}.log"
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        cut -f1,2 {input.fasta}.fai > {output.sizes}
        """

rule genome_coverage:
    input:
        bed   = f"{outDir}/mapped/{{genome}}/{{sample}}.bed",
        sizes = f"{outDir}/genomeSizes/{{genome}}.sizes"
    output:
        coverage = f"{outDir}/coverage/{{genome}}/{{sample}}.coverage.txt"
    conda: "envs/RMS.yaml"
    log: f"{outDir}/logs/coverage_{{sample}}_{{genome}}.log"
    shell:
        """
        bedtools genomecov \
          -i {input.bed} \
          -g {input.sizes} \
          -d > {output.coverage} 2> {log}
        """


rule compute_rms:
    input:
        count3   = f"{outDir}/endCounts/{{genome}}/{{sample}}.3count.txt",
        coverage = f"{outDir}/coverage/{{genome}}/{{sample}}.coverage.txt"
    output:
        results = f"{outDir}/rmsScores/{{genome}}/{{sample}}.rms.txt"
    conda: "envs/R_base.yaml"
    log: f"{outDir}/logs/compute_rms_{{sample}}_{{genome}}.log"
    script: "scripts/compute_rms.R"