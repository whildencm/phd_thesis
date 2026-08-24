#!/bin/bash
#SBATCH --array=0-27 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=8 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o logs/%x_%j.o #output file name
#SBATCH -e logs/%x_%j.e #error file name
#SBATCH -p whipple #to run job on whipple lab partition
#SBATCH --mem=100G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

####################################

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" config/Samples.txt)

GTF=/n/lab_storage/whipple_lab/share_root/Lab/Computing_Resources/Genome_mmu

mkdir -p results/counts

module load Mambaforge
mamba activate RPF_Env_2024

# # Filter RPF bam files and count
samtools view -h results/alignments/mm39/${SAMPLE}.RPF.bam | \
awk 'BEGIN {OFS="\t"} /^@/ {print; next} {l=length($10); if (l==32 || l==33 || l==34 || l==35) print}' | \
samtools view -b -o results/alignments/mm39/${SAMPLE}.RPF.sizefiltered.bam

echo 'Indexing BAM files sample:' $SAMPLE
samtools index results/alignments/mm39/${SAMPLE}.RPF.sizefiltered.bam

featureCounts \
-a $GTF/gencode.vM30.annotation.gtf \
-t CDS \
-o results/counts/${SAMPLE}.RPF.featurecounts \
-T 8 \
results/alignments/mm39/${SAMPLE}.RPF.sizefiltered.bam

# mRNA counts
featureCounts \
-a $GTF/gencode.vM30.annotation.gtf \
-p -t CDS \
-o results/counts/${SAMPLE}.mRNA.featurecounts \
-T 8 \
results/alignments/mm39/${SAMPLE}.RNA.bam


####################################################################
# Generate BigWig files

echo 'Generating BigWig files'

module load Mambaforge
mamba activate bam2bw

bamCoverage \
  -b results/alignments/mm39/${SAMPLE}.RPF.sizefiltered.bam \
  -o results/alignments/mm39/${SAMPLE}.RPF.sizefiltered.bw \
  --normalizeUsing RPKM \
  --binSize 10 \
  --extendReads 200\
  --ignoreDuplicates