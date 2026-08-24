#!/bin/bash
#SBATCH --array=1-25 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=8 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o logs/%x_%j.o #output file name
#SBATCH -e logs/%x_%j.e #error file name
#SBATCH -p whipple #to run job on whipple lab partition
#SBATCH --mem=150G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

####################################

#Load mamba enviroment
module load Mambaforge
mamba activate RPF_Env_2024

# Pull sample names from Samples.txt
#SAMPLE=EL1
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" config/Samples.txt)

echo 'Begin mapping.sh' ${SAMPLE}

# Paths to reference genomes & annotations
RRNA=/n/lab_storage/whipple_lab/share_root/Lab/Computing_Resources/Genome_STAR/chimeCLIP_rep 
GENOME=/n/lab_storage/whipple_lab/share_root/Lab/Computing_Resources/Genome_STAR/GRCm39_M30
GTF=/n/lab_storage/whipple_lab/share_root/Lab/Computing_Resources/Genome_mmu

fastqDir=/n/lab_storage/whipple_lab/share_root/Lab/data_raw_seq/CMW_EXPT096_riboseq_snord116
outDir=/n/lab_storage/whipple_lab/share_root/Users/cwhilden/PROJECTS/manuscript_2026/riboseq/results

# Adapter trimming
echo 'Begin adapter trimming' 
mkdir -p $outDir/alignments/trimmedfastq

cutadapt -a AGATCGGAAGAGCACACGTCT \
--minimum-length 10 \
$fastqDir/ribo_raw_fastq/${SAMPLE}_CKDL250009974-1A_22J7HJLT4_L1_1.fq.gz \
> $outDir/alignments/trimmedfastq/${SAMPLE}"_trim.fastq" \
2> $outDir/alignments/trimmedfastq/${SAMPLE}"_trimreport.txt"

echo 'Adapter trimming complete'


# Align trimmed reads to rRNA/tRNA/snRNA and keep unmapped reads
echo 'Begin aligning to rRNA' 

# Path to trimmed fastq files
TRIMFASTQ=$outDir/alignments/trimmedfastq/${SAMPLE}"_trim.fastq"
mkdir -p $outDir/alignments/rRNA
# Align using STAR
STAR \
--genomeDir $RRNA \
--runMode alignReads \
--runThreadN 8 \
--outSAMtype BAM Unsorted \
--readFilesIn $TRIMFASTQ \
--outReadsUnmapped Fastx \
--outFileNamePrefix $outDir/alignments/rRNA/${SAMPLE}

echo 'rRNA alignment complete'

# Align remaining reads to genome
echo 'Begin aligning to genome' 

CLEANREADS=$outDir/alignments/rRNA/${SAMPLE}Unmapped.out.mate1
mkdir -p $outDir/alignments/mm39
STAR \
	--sjdbGTFfile $GTF/gencode.vM30.annotation.gtf \
	--runMode alignReads \
	--runThreadN 8 \
	--genomeDir $GENOME \
	--readFilesIn $CLEANREADS \
	--outFileNamePrefix $outDir/alignments/mm39/${SAMPLE}".genomemapped." \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM \
	--outFilterMultimapNmax 1 \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 16 \
	--outFilterMismatchNmax 1

echo 'Genome alignment complete'

# Rename and index bam
mv $outDir/alignments/mm39/${SAMPLE}.genomemapped.Aligned.sortedByCoord.out.bam \
$outDir/alignments/mm39/${SAMPLE}.RPF.bam

echo 'Indexing BAM files sample:' $SAMPLE
samtools index $outDir/alignments/mm39/${SAMPLE}.RPF.bam

####################################################################
# mRNA mapping

# Paths to paired fastq files
FASTQ1=$fastqDir/rna_raw_fastq/${SAMPLE}_1.fq.gz
FASTQ2=$fastqDir/rna_raw_fastq/${SAMPLE}_2.fq.gz

date
echo 'Begin mRNA genome alignment'

STAR \
--runThreadN 8 \
--runMode alignReads \
--genomeDir $GENOME \
--readFilesIn $FASTQ1 $FASTQ2 \
--readFilesCommand zcat \
--outFileNamePrefix $outDir/alignments/mm39/${SAMPLE}.genomemapped \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 70 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbGTFfile $GTF/gencode.vM30.annotation.gtf \
--sjdbOverhang 149 \
--alignEndsType EndToEnd \
--quantMode TranscriptomeSAM


mv $outDir/alignments/mm39/${SAMPLE}.genomemappedAligned.sortedByCoord.out.bam \
$outDir/alignments/mm39/${SAMPLE}.RNA.bam

date
echo 'Finished aligning sample:' $SAMPLE

# Index resulting bam files
date
echo 'Indexing bam file for sample:' $SAMPLE

samtools index \
$outDir/alignments/mm39/${SAMPLE}.RNA.bam

date
echo 'Finished indexing bam file for sample:' $SAMPLE

mamba deactivate

####################################################################
# Generate BigWig files

echo 'Generating BigWig files'

module load Mambaforge
mamba activate bam2bw

bamCoverage \
  -b $outDir/alignments/mm39/${SAMPLE}.RPF.bam \
  -o $outDir/alignments/mm39/${SAMPLE}.RPF.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --extendReads 200\
  --ignoreDuplicates

  bamCoverage \
  -b $outDir/alignments/mm39/${SAMPLE}.RNA.bam \
  -o $outDir/alignments/mm39/${SAMPLE}.RNA.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --extendReads \
  --ignoreDuplicates

# Sort bam files
mkdir -p $outDir/alignments/sortedRpfTrs

samtools \
sort $outDir/alignments/mm39/${SAMPLE}.genomemapped.Aligned.toTranscriptome.out.bam \
-o $outDir/alignments/sortedRpfTrs/${SAMPLE}.RPF.Transcriptome.Sort.bam

samtools index $outDir/alignments/sortedRpfTrs/${SAMPLE}.RPF.Transcriptome.Sort.bam