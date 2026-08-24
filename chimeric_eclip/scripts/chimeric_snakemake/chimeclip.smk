# CourtneyWhilden Dec 2025

import os
import re

###=====
# LOAD CONFIGURATIONS
###=====
sampleSheet = config["IP_samples"]
inputSamples = config["IN_samples"]
outputDir = config["outputDir"]

snoRef = config["snoRNAFasta"]
adapterRef = config["adapterFasta"]
rRNARef = config["rRNAFasta"]
scaRNARef = config["scaRNAFasta"]
tRNARef = config["tRNAFasta"]
snRNARef = config["snRNAFasta"]

#repBase = config["repBaseFasta"]
genome = config["genomeFasta"]
genomeIndex = config["genomeIndexDir"]
gtf = config["genomeGTF"]

scriptDir = config["scriptDir"]

readLength = int(config["readLength"])

###===
# DEFINE SAMPLES AND TARGET DEFINITIONS 
###===

# 1) READ IN IP SAMPLES

# Read in fq.gz files
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

# 2) READ IN INPUT SAMPLES
with open(inputSamples) as f:
    inputSamplePaths = [line.strip() for line in f if line.strip()]

# Strip sample names from fastq file (everything before fq.gz)
inputSamples = [
    re.sub(r'\.(fastq|fq)\.gz$', '', os.path.basename(p))
    for p in inputSamplePaths
]

# Dictionary to map sample ID to absolute fastq path
inputSampleToPath = dict(zip(inputSamples, inputSamplePaths))

# Function to grab sample name from sample path
def getInputSamplePath(wildcards):
    return inputSampleToPath[wildcards.sample]

# Merge IP and input path dictionaries so the trimming rule can serve both.
# If a sample name appears in BOTH lists this will silently overwrite the IP
# path with the input path — make sure sample names are unique across sets.
allSampleToPath = {**sampleToPath, **inputSampleToPath}

def getAllSamplePath(wildcards):
    """Return the raw fastq path for any IP or input sample."""
    return allSampleToPath[wildcards.sample]

# DEFINE TARGETS FOR MAPPING

# Define all target IDs
allTargets = ["rRNA", "tRNA", "scaRNA", "snRNA", "genome"]

# Define repeat target IDs
repeatTargets = ["rRNA", "tRNA", "scaRNA", "snRNA"]

###===
# RULE ALL
###===

rule all:
    input:
        # IP samples → chimera pipeline
        expand(
            f"{outputDir}/chimeraBamLists/{{sample}}.{{target}}.bamList.txt",
            sample=samples,          # IP only
            target=allTargets,
        ),
        # Input samples → CLIP mapping pipeline
        expand(
            f"{outputDir}/mapCLIP/{{sample}}.clip.genomeAligned.out.sam",
            sample=inputSamples,     # input only
        )


rule umi_extract_adapter_trim:
    """
    Extract unique molecular identifiers and trim reads with pre-defined adapter sequences for IP and IN samples.
    """
    input:
        fastq = getAllSamplePath,
        adapters = adapterRef
    output: 
        trimmed = f"{outputDir}/umi_extract_adapter_trim/{{sample}}.trim.fastq.gz"
    threads: 8
    log: f"{outputDir}/logs/{{sample}}.umi_extract.log"
    conda: "envs/umi.yml"
    shell:
        '''
        umi_tools extract \
            --random-seed 1 \
            --bc-pattern NNNNNNNNNN \
            --stdin {input.fastq} \
            --log /dev/null \
            | cutadapt \
            -j {threads} \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 34 \
            -a file:{input.adapters} \
            -O 1 \
            - \
            --report minimal \
            2> {log} \
            | cutadapt \
            -j {threads} \
            --match-read-wildcards \
            --times 1 \
            -e 0.1 \
            --quality-cutoff 6 \
            -m 34 \
            -a file:{input.adapters} \
            -O 5 \
            - \
            --report minimal \
            2>> {log} \
            | cutadapt \
            -j {threads} \
            -u -10 \
            - \
            -o {output.trimmed} \
            --report minimal \
            >> {log} \
            2>> {log}

        ''' 

###===
# CHIMERIC PIPELINE
###===
# # # ALIGN TO SNORNAS
rule index_snoRNA:
    """
    Build bowtie2 index files for snoRNA reference sequences
    """
    input:
        snoFasta = snoRef
    output:
        snoIndex = f"{outputDir}/bowtieIndex/sno.bowtie2.index.1.bt2"
    params:
        prefix = f"{outputDir}/bowtieIndex/sno.bowtie2.index"
    log: f"{outputDir}/logs/build_sno.log"
    conda: "envs/bowtie2.yml"
    shell:
        '''
        bowtie2-build \
        --quiet --offrate 2 \
        {input.snoFasta} {params.prefix} 2> {log}
        '''

rule map_snoRNA:
    """
    Map trimmed fastq files to snoRNA sequences using bowtie2
    """
    input:
        fastqTr = f"{outputDir}/umi_extract_adapter_trim/{{sample}}.trim.fastq.gz",
        snoIndex = f"{outputDir}/bowtieIndex/sno.bowtie2.index.1.bt2"
    output:
        sam = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.sam"
    params:
        indexPrefix = f"{outputDir}/bowtieIndex/sno.bowtie2.index"
    threads: 8
    conda: "envs/bowtie2.yml"
    log: f"{outputDir}/logs/{{sample}}.map_sno.txt"
    shell:
        '''
        bowtie2 \
            -D 20 \
            -R 3 \
            -N 0 \
            -L 16 \
            --local \
            --norc \
            -i S,1,0.50 \
            --score-min L,16,0 \
            --ma 1 \
            --np 0 \
            --mp 2,2 \
            --rdg 5,1 \
            --rfg 5,1 \
            -p {threads} \
            --no-unal \
            -x {params.indexPrefix} \
            -q \
            -U {input.fastqTr} \
            -a \
            -S {output.sam} \
            2> {log}

        '''
    
rule split_snoRNA_sams:
    """
    Split snoRNA mapped sam files to be compatible with
    find putative target script.
    Script from Eric's lab.
    """
    input: 
        sam = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.sam"
    output:
        split = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.AA.bam"
    params:
        script = f"{scriptDir}/split_source_map_sam"
    conda: "envs/chimeras.yml"
    shell:
        '''
        {params.script} \
        {input.sam} \
        {output.split}
        '''

rule find_putative_targets:
    """
    Use snoRNA mapped sam files to identify putative chimeric targets. 
    Script find_putative_target from Eric's lab
    """
    input: 
        req = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.AA.bam"
    output:
        fasta = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.target.fasta"
    conda: "envs/chimeras.yml"
    params: 
        samplePrefix = f"{outputDir}/mapChimeras/{{sample}}",
        script = f"{scriptDir}/find_putative_target"
    shell:
        """
        dnList=(AA AT AG AC AN TA TT TG TC TN GA GT GG GC GN CA CT CG CC CN NA NT NG NC NN)

        for dn in "${{dnList[@]}}"; do
            {params.script} \
                "{params.samplePrefix}.map.to.snoRNA.${{dn}}.bam" \
                "{params.samplePrefix}.map.to.snoRNA.${{dn}}.target.fasta" \
                snoRNA
        done

        cat {params.samplePrefix}.map.to.snoRNA.[ATGCN][ATGCN].target.fasta \
            > {params.samplePrefix}.map.to.snoRNA.target.fasta

        rm {params.samplePrefix}.map.to.snoRNA.[ATGCN][ATGCN].target.fasta || true

        cat {params.samplePrefix}.map.to.snoRNA.[ATGCN][ATGCN].target.csv \
            > {params.samplePrefix}.map.to.snoRNA.target.tmp.csv

        head -1 "{params.samplePrefix}.map.to.snoRNA.target.tmp.csv" > "{params.samplePrefix}.map.to.snoRNA.target.csv"
        tail -n +2 "{params.samplePrefix}.map.to.snoRNA.target.tmp.csv" >> "{params.samplePrefix}.map.to.snoRNA.target.csv"

        rm "{params.samplePrefix}.map.to.snoRNA.target.tmp.csv"
        rm "{params.samplePrefix}.map.to.snoRNA.[ATGCN][ATGCN].target.csv" || true
        rm "{params.samplePrefix}.map.to.snoRNA.[ATGCN][ATGCN].bam" || true
        """

rule map_true_snoRNA_targets:
    input:
        fasta = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.target.fasta",
        snoIndex = f"{outputDir}/bowtieIndex/sno.bowtie2.index.1.bt2"
    output:
        fasta = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.true.target.fasta",
        backMapSam = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.target.back.map.sam" 
    conda: "envs/bowtie2.yml"
    threads: 8
    params:
        samplePrefix = f"{outputDir}/mapChimeras/{{sample}}",
        indexPrefix = f"{outputDir}/bowtieIndex/sno.bowtie2.index"
    shell:
        '''
        bowtie2 \
            -D 20 \
            -R 3 \
            -N 0 \
            -L 16 \
            --local \
            --norc \
            -i S,1,0.50 \
            --score-min L,12,0 \
            --ma 1 \
            --np 0 \
            --mp 2,2 \
            --rdg 5,1 \
            --rfg 5,1 \
            -p {threads} \
            --no-unal \
            -x {params.indexPrefix} \
            --un "{params.samplePrefix}.map.to.snoRNA.true.target.fasta" \
            -f \
            -U "{params.samplePrefix}.map.to.snoRNA.target.fasta" \
            -a \
            -S "{params.samplePrefix}.map.to.snoRNA.target.back.map.sam" \
            2> "{params.samplePrefix}.map.to.snoRNA.target.back.map.log"
        '''

# # # ALIGN TO REPEAT GENOMES (rRNA, snRNA, tRNA, scaRNA)

repeatFastas = {
    "rRNA": rRNARef,
    "tRNA": tRNARef,
    "scaRNA": scaRNARef,
    "snRNA": snRNARef
}

def getRepeatFasta(wildcards):
    # Get fasta path from target name
    return repeatFastas[wildcards.repeatType]

rule index_RNAtarget:
    # Build bowtie index for repeats fasta files
    input:
        fasta = getRepeatFasta
    output:
        index = f"{outputDir}/bowtieIndex/{{repeatType}}.bowtie.index.1.ebwt"
    params: 
        prefix = f"{outputDir}/bowtieIndex/{{repeatType}}.bowtie.index"
    conda: "envs/bowtie.yml"
    log: f"{outputDir}/logs/build_repeats.{{repeatType}}.log"
    shell:
        '''
        bowtie-build \
        --quiet --offrate 2 \
        {input.fasta} {params.prefix} \
        2> {log}
        '''

rule map_RNAtargets:
    # Map true putative targets to repeat genomes
    input: 
        index = f"{outputDir}/bowtieIndex/{{repeatType}}.bowtie.index.1.ebwt",
        fasta = f"{outputDir}/mapChimeras/{{sample}}.map.to.snoRNA.true.target.fasta"
    output:
        sam    = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.repeats.{{repeatType}}.sam",
        unmap  = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{repeatType}}.unmap.fasta",
        logFile = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{repeatType}}.log",
        done   = f"{outputDir}/mapChimeras/{{sample}}.{{repeatType}}.mapRepeats.done"
    params:
        indexPrefix  = f"{outputDir}/bowtieIndex/{{repeatType}}.bowtie.index",            
        samplePrefix = f"{outputDir}/mapChimeras/{{sample}}.{{repeatType}}",
        script = f"{scriptDir}/mapRepeats.sh"
    conda: "envs/bowtie.yml"
    threads: 8
    shell:
        '''
        bowtie \
            -a \
            --best \
            --strata \
            -e 35 \
            -q \
            -l 8 \
            -n 1 \
            -p {threads} \
            -x {params.indexPrefix} \
            --no-unal \
            --norc \
            -f {input.fasta} \
            --sam {output.sam} \
            --un {output.unmap} \
            2> {output.logFile}

            touch {output.done}
        '''

rule compile_unmapped_RNAtarget:
    # Take all unmapped reads and merge them to one file
    # Script compile_unmapped etc. is from ERic's lab.
    input:
        rRNA = f"{outputDir}/mapChimeras/{{sample}}.rRNA.mapRepeats.done",
        tRNA = f"{outputDir}/mapChimeras/{{sample}}.tRNA.mapRepeats.done",
        sca = f"{outputDir}/mapChimeras/{{sample}}.scaRNA.mapRepeats.done",
        snRNA = f"{outputDir}/mapChimeras/{{sample}}.snRNA.mapRepeats.done"
    output:
        fasta = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.RNA.unmap.fasta"
    params:
        prefix = f"{outputDir}/mapChimeras/{{sample}}",
        script = f"{scriptDir}/compile_unmapped_putative_target"
    conda: "envs/chimeras.yml"
    shell:
        '''
        {params.script} \
        {params.prefix} snoRNA rRNA,tRNA,snRNA,scaRNA {output.fasta}

        '''

# # # ALIGN TO GENOME

rule index_genome:
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
    conda: "envs/star.yml"
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

rule map_genome:
    input: 
        genome = f"{genomeIndex}/Genome",
        fasta  = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.RNA.unmap.fasta"
    output:
        sam  = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.genomeAligned.out.sam",     
        done = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.genome.done"                
    params:
        outFilePrefix = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.genome",
        genomeDir = genomeIndex
    log: f"{outputDir}/logs/{{sample}}.map_genome.log"                                         
    threads: 8                                                                       
    conda: "envs/star.yml"
    shell:
        """
        STAR \
            --readFilesType Fastx \
            --alignEndsType EndToEnd \
            --genomeDir {params.genomeDir} \
            --genomeLoad NoSharedMemory \
            --outFileNamePrefix {params.outFilePrefix} \
            --outFilterMatchNminOverLread 0.66 \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:foo \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype SAM \
            --outSAMunmapped None \
            --outStd Log \
            --readFilesIn {input.fasta} \
            --runMode alignReads \
            --runThreadN {threads} > {log} 2>&1

            touch {params.outFilePrefix}.done
        """

# # # IDENTIFY CHIMERAS

def getSamForTarget(wildcards):
    # Define sam file path depending on the target type
    if wildcards.target == "genome":
        # STAR genome-mapping SAM
        return f"{outputDir}/mapChimeras/{wildcards.sample}.snoRNA.genomeAligned.out.sam"
    else:
        # bowtie repeat-mapping SAMs
        return f"{outputDir}/mapChimeras/{wildcards.sample}.snoRNA.repeats.{wildcards.target}.sam"

rule chimeric_sam_to_bam:
    # Convert all mapped sam files to bam
    # sam to bam is from eric's lab
    input:
        sam = getSamForTarget
    output: 
        bam = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.bam"
    params:
        script = f"{scriptDir}/sam_to_bam"
    conda: "envs/chimeras.yml"
    shell:
        """
        {params.script} \
        {input.sam} {output.bam} True
        """

rule identify_chimeric_reads:
    # Identify chimeric reads from genome/repeat mapped bam files
    input:
        bam = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.bam"
    output:
        bam = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.chimeras.bam"
    params:
        script = f"{scriptDir}/identify_chimeric_read_cw",
        countFile = f"{outputDir}/mapChimeras/{{sample}}.{{target}}.chimeras.count.tsv"
    conda: "envs/chimeras.yml"
    shell:
        """
        {params.script} \
        {input.bam} {params.countFile} {output.bam}
        """

rule parse_individual_chimeras:
    # parse individual chimeras is from eric's lab
    input:
        bam = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.chimeras.bam"
    output:
        counts = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.individual.chimeras.count.tsv"
    params:
        script = f"{scriptDir}/parse_individual_chimeras",
        samplePrefix = f"{{sample}}.snoRNA.{{target}}"
    conda: "envs/bam2bw.yml"
    shell:
        """
        #Input samples have so few chimeras that some bams will be empty.
        readCount=$(samtools view -c {input.bam} 2>/dev/null || echo 0)

        if [ "${{readCount}}" -gt 0 ]; then

            {params.script} \
            {input.bam} \
            --uid {params.samplePrefix} \
            --tag {wildcards.target} \
            --minimum_chimeras 1 \
            --count_output {output.counts} \
            --no_bw

        else
            echo "{input.bam} contains no reads"
            touch {output.counts}
        fi
        """

rule get_target_bam_lists:
    """
    Output a list of chimeric bam files for each sample/target combination.
    """
    input: 
        counts = f"{outputDir}/mapChimeras/{{sample}}.snoRNA.{{target}}.individual.chimeras.count.tsv"
    output:
        bamList = f"{outputDir}/chimeraBamLists/{{sample}}.{{target}}.bamList.txt"
    params:
        bamDir = f"{outputDir}/mapChimeras"
    shell:
        """
        find {params.bamDir} \
          -type f \
          -name "{wildcards.sample}.snoRNA.{wildcards.target}.*.chimeras.bam" \
          > {output.bamList}
        """

# CLIP MAPPING OF INPUT SAMPLES

rule clip_map_RNAtargets:
    """
    Sequentially map trimmed reads to RNA targets and output unmapped reads for genome mapping.
    Order: rRNA > snRNA > tRNA > scaRNA
    Reads that fail to map at each step are passed to the next.
    """
    input:
        trimFastq   = f"{outputDir}/umi_extract_adapter_trim/{{sample}}.trim.fastq.gz",
        rRNAIndex   = f"{outputDir}/bowtieIndex/rRNA.bowtie.index.1.ebwt",
        snRNAIndex  = f"{outputDir}/bowtieIndex/snRNA.bowtie.index.1.ebwt",
        tRNAIndex   = f"{outputDir}/bowtieIndex/tRNA.bowtie.index.1.ebwt",
        scaRNAIndex = f"{outputDir}/bowtieIndex/scaRNA.bowtie.index.1.ebwt"
    output:
        rRNAsam   = f"{outputDir}/mapCLIP/{{sample}}.rRNA.sam",
        snRNAsam  = f"{outputDir}/mapCLIP/{{sample}}.snRNA.sam",
        tRNAsam   = f"{outputDir}/mapCLIP/{{sample}}.tRNA.sam",
        scaRNAsam = f"{outputDir}/mapCLIP/{{sample}}.scaRNA.sam",
        unmapped  = f"{outputDir}/mapCLIP/{{sample}}.RNAtargets.unmapped.fasta"
    params:
        rRNAPrefix   = f"{outputDir}/bowtieIndex/rRNA.bowtie.index",
        snRNAPrefix  = f"{outputDir}/bowtieIndex/snRNA.bowtie.index",
        tRNAPrefix   = f"{outputDir}/bowtieIndex/tRNA.bowtie.index",
        scaRNAPrefix = f"{outputDir}/bowtieIndex/scaRNA.bowtie.index",
        # Prefix for intermediate unmapped files — cleaned up at the end
        tmpPrefix    = f"{outputDir}/mapCLIP/{{sample}}.tmp"
    log: f"{outputDir}/logs/{{sample}}.clip_repeat_mapping.log"
    threads: 8
    conda:   "envs/bowtie.yml"
    shell:
        """
        # Step 1: rRNA
        bowtie \
            -a --best --strata -e 35 \
            -q {input.trimFastq} \
            -l 8 -n 1 -p {threads} \
            -x {params.rRNAPrefix} \
            --no-unal --norc \
            --sam {output.rRNAsam} \
            --un  {params.tmpPrefix}.rRNA.unmapped.fastq \
            2>> {log}

        # Step 2: snRNA
        bowtie \
            -a --best --strata -e 35 \
            -q {params.tmpPrefix}.rRNA.unmapped.fastq \
            -l 8 -n 1 -p {threads} \
            -x {params.snRNAPrefix} \
            --no-unal --norc \
            --sam {output.snRNAsam} \
            --un  {params.tmpPrefix}.snRNA.unmapped.fastq \
            2>> {log}

        # Step 3: tRNA
        bowtie \
            -a --best --strata -e 35 \
            -q {params.tmpPrefix}.snRNA.unmapped.fastq \
            -l 8 -n 1 -p {threads} \
            -x {params.tRNAPrefix} \
            --no-unal --norc \
            --sam {output.tRNAsam} \
            --un  {params.tmpPrefix}.tRNA.unmapped.fastq \
            2>> {log}

        # Step 4: filter scaRNA & output unmapped reads
        bowtie \
            -a --best --strata -e 35 \
            -q {params.tmpPrefix}.tRNA.unmapped.fastq \
            -l 8 -n 1 -p {threads} \
            -x {params.scaRNAPrefix} \
            --no-unal --norc \
            --sam {output.scaRNAsam} \
            --un  {output.unmapped} \
            2>> {log}

        # Clean up intermediate unmapped files
        rm -f {params.tmpPrefix}.rRNA.unmapped.fastq \
              {params.tmpPrefix}.snRNA.unmapped.fastq \
              {params.tmpPrefix}.tRNA.unmapped.fastq
        """


rule clip_map_snoRNA:
    input:
        fasta = f"{outputDir}/mapCLIP/{{sample}}.RNAtargets.unmapped.fasta",
        index = f"{outputDir}/bowtieIndex/sno.bowtie2.index.1.bt2"
    output:
        sam   = f"{outputDir}/mapCLIP/{{sample}}.clip.snoRNA.sam",
        unmap = f"{outputDir}/mapCLIP/{{sample}}.clip.snoRNA.unmap.fasta"
    params:
        indexPrefix = f"{outputDir}/bowtieIndex/sno.bowtie2.index"
    threads: 8
    log: f"{outputDir}/logs/{{sample}}.clip_map_snoRNA.log"
    conda: "envs/bowtie2.yml"
    shell:
        '''
        bowtie2 \
            -D 20 \
            -R 3 \
            -N 0 \
            -L 16 \
            --local \
            --norc \
            -i S,1,0.50 \
            --score-min L,12,0 \
            --ma 1 \
            --np 0 \
            --mp 2,2 \
            --rdg 5,1 \
            --rfg 5,1 \
            -p {threads} \
            --no-unal \
            -x {params.indexPrefix} \
            -U {input.fasta} \
            -a \
            -S {output.sam} \
            --un {output.unmap} \
            2> {log}
        '''


rule clip_map_genome:
    input: 
        genome = f"{genomeIndex}/Genome",
        fasta  =f"{outputDir}/mapCLIP/{{sample}}.clip.snoRNA.unmap.fasta"
    output:
        sam  = f"{outputDir}/mapCLIP/{{sample}}.clip.genomeAligned.out.sam",     
        done = f"{outputDir}/mapCLIP/{{sample}}.clip.genome.done"                
    params:
        outFilePrefix = f"{outputDir}/mapCLIP/{{sample}}.clip.genome",
        genomeDir = genomeIndex
    log: f"{outputDir}/logs/{{sample}}.clip_map_genome.log"                                         
    threads: 8                                                                       
    conda: "envs/star.yml"
    shell:
        """
        STAR \
            --readFilesType Fastx \
            --alignEndsType EndToEnd \
            --genomeDir {params.genomeDir} \
            --genomeLoad NoSharedMemory \
            --outFileNamePrefix {params.outFilePrefix} \
            --outFilterMatchNminOverLread 0.66 \
            --outFilterMultimapNmax 1 \
            --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterType BySJout \
            --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:foo \
            --outSAMattributes All \
            --outSAMmode Full \
            --outSAMtype SAM \
            --outSAMunmapped None \
            --outStd Log \
            --readFilesIn {input.fasta} \
            --runMode alignReads \
            --runThreadN {threads} > {log} 2>&1

            touch {params.outFilePrefix}.done
        """