#!/bin/bash
#SBATCH --array=1 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=4 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o slurm-logs/eCLIP_peakNormalization/%x_%A_%a.o #output file name
#SBATCH -e slurm-logs/eCLIP_peakNormalization/%x_%A_%a.e #error file name
#SBATCH -p whipple
#SBATCH --mem=50G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

module load Mambaforge
mamba activate IDR_py38
#=============
# Positional Arguments
#=============
bam_IP="$1"
bam_IN="$2"
peaks_bed="$3"
outDir="$4"

# Make directories for output
count_ip=$outDir/count_ip
count_in=$outDir/count_in
norm=$outDir/peaksNorm

mkdir -p $count_ip $count_in $norm


# Generate count files 
base_ip=$(basename $bam_IP .bam)
base_in=$(basename $bam_IN .bam)

samtools view -c "$bam_IP" > "$count_ip/${base_ip}.count"

samtools view -c "$bam_IN" > "$count_in/${base_in}.count"


perl scripts/overlap_peakfi_with_bam.pl \
    $bam_IP \
    $bam_IN \
    $peaks_bed \
    "$count_ip/${base_ip}.count" \
    "$count_in/${base_in}.count" \
    $norm/${base_ip}