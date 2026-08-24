
sampleSheet = config["sampleSheet"]
outputDir = config["outputDir"]

genome=config["genomeOneCopy"]
genomeIndex=config["genomeIndexOneCopy"]

genomeUnique=config["genomeUnique"]
genomeIndexUnique=config["genomeIndexUnique"]

gtf=config["gtf"]
readLength=config["readLength"]

#====
# Read in sample sheet
#====
with open(sampleSheet) as f:
    samplePaths = [line.strip() for line in f if line.strip()]

# Strip sample names from fastq file (everything before fq.gz)
samples = [
    re.sub(r'\.(fastq|fq)\.gz$', '', os.path.basename(p))
    for p in samplePaths
]

# Dictionary to map sample ID to absolute fastq path
sampleToPath = dict(zip(samples, samplePaths))

# Function to grab sample name from sample path
def getSamplePath(wildcards):
    return sampleToPath[wildcards.sample]

#===
# Rule all
#===



#===
# Index genomes
#===
rule indexOne:
    # Build STAR index files from the genome fasta
    input: 
        fasta = genome,
        gtf = gtf
    output:
        outFile = f"{genomeIndex}/Genome"
    params:
        genomeDir = genomeIndex,
        overhang = readLength - 1
    threads: 8
    conda: "eclip_envs/star.yml"
    shell:
        """
        mkdir -p {params.genomeDir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --genomeDir {params.genomeDir} \
            --genomeFastaFiles {input.fasta}
        """


rule indexUnique:
    # Build STAR index files from the genome fasta
    input: 
        fasta = genomeUnique,
        gtf = gtf
    output:
        outFile = f"{genomeIndexUnique}/Genome"
    params:
        genomeDir = genomeIndexUnique,
        overhang = readLength - 1
    threads: 8
    conda: "eclip_envs/star.yml"
    shell:
        """
        mkdir -p {params.genomeDir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --genomeDir {params.genomeDir} \
            --genomeFastaFiles {input.fasta}
        """

#===
# Mapping
#===
rule mapOne:
    # Map to genome w/ one copy of each snoRNA
    input:
        fastq=getSamplePath,
        index=f"{genomeIndex}/Genome"
    output:
        bam=f"{outputDir}/map_oneCopy/{{sample}}.out.bam"
    params:
        outFilePrefix = f"{outputDir}/map_oneCopy/{{sample}}.out"
        genomeDir = genomeIndex
    log: f"{outputDir}/logs/{{sample}}.log"
    threads: 8
    conda: "eclip_envs/star.yml"
    shell:
        """
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.genomeDir} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outFilePrefix} \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --outSAMattrRGline ID:foo \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Sorted \
            --outStd Log \
            --readFilesIn {input.fastq} \
            --runMode alignReads \
            --runThreadN {threads} > {log} 2>&1
        """

rule mapUnique:
    # map to genome with only unique sequence copies of each snoRNA
    input:
        fastq=getSamplePath,
        index=f"{genomeIndexUnique}/Genome"
    output:
        bam=f"{outputDir}/map_unique/{{sample}}.out.bam"
    params:
        outFilePrefix = f"{outputDir}/map_unique/{{sample}}.out"
        genomeDir = genomeIndexUnique
    log: f"{outputDir}/logs/{{sample}}.log"
    threads: 8
    conda: "eclip_envs/star.yml"
    shell:
        """
        STAR \
            --alignEndsType EndToEnd \
            --genomeDir {params.genomeDir} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression 10 \
            --outFileNamePrefix {params.outFilePrefix} \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterType BySJout \
            --outSAMattrRGline ID:foo \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype BAM Sorted \
            --outStd Log \
            --readFilesIn {input.fastq} \
            --runMode alignReads \
            --runThreadN {threads} > {log} 2>&1
        """
