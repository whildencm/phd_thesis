#!/bin/bash
#SBATCH --array=1 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=12 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o slurm-logs/eCLIP_peakCalling/%x_%A_%a.o #output file name
#SBATCH -e slurm-logs/eCLIP_peakCalling/%x_%A_%a.e #error file name
#SBATCH -p whipple #to run job on whipple lab partition
#SBATCH --mem=100G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

####################################
# Positional args
bam="$1"
outDir="$2"

mkdir -p $outDir

# Load modules
module load Mambaforge
mamba activate clipper3_mm39

base=$(basename $bam ".trim.Aligned.sortedByCoord.out.rmdup.bam")
echo "Running sample " $base "from bam file " $bam

# Run clipper
clipper --processors=12 -b $bam -o ${outDir}/${base}.bed -s mm39