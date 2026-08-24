#!/bin/bash
#SBATCH --array=1-4 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=12 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o slurm-logs/peakCallRepeats/%x_%A_%a.o #output file name
#SBATCH -e slurm-logs/peakCallRepeats/%x_%A_%a.e #error file name
#SBATCH -p whipple #to run job on whipple lab partition
#SBATCH --mem=184G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

####################################
module load Mambaforge
mamba activate /n/whipple_lab/share_root/Users/cwhilden/ENVS/samtools_bedtools

# =========================
# Positional arguments
# =========================
INPUT="$1"    # chimeclip results directory
OUTPUT="$2"   # output directory
FASTA="$3"    # fasta reference directory, same as used for chimeric pipeline
SAMPLE="$4"   # sample ID, IP
inputSample="$5"

if [ "$#" -ne 5 ]; then
  echo "Usage:"
  echo "  sbatch $0 <INPUT> <OUTPUT> <FASTA> <SAMPLE> <inputSample>"
  exit 1
fi

mkdir -p "$OUTPUT"

# Submit a parallel job for each target class
CLASSES=(rRNA snRNA tRNA scaRNA)
CLASS=${CLASSES[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing class ${CLASS} for IP: {}"
##===
# PAD INPUT CLIP SAM FILES
# Input = .sam files for each target class, from chimeric pipeline
# Output = sorted, indexed, padded (+30nt to each coordinate) .bam files.
#====

# Find CLIP sam file corresponding to the target class

clipSAM=${INPUT}/mapCLIP/${inputSample}.${CLASS}.sam
outClipSAM=${OUTPUT}/${CLASS}_clipPaddedBam
mkdir -p $outClipSAM

# Pad the same file and convert to a sorted bam file.

echo "Padding CLIP sam file for class ${CLASS}"

awk '{
    if ($1 ~ /^@/) {
      print $0
    } else {
      $4 = $4 + 30; print
    }
  }' FS='\t' OFS='\t' $clipSAM > ${outClipSAM}/${inputSample}.${CLASS}.padded.sam

  # Convert to bam, sort, then index.
samtools view -b ${outClipSAM}/${inputSample}.${CLASS}.padded.sam -o ${outClipSAM}/${inputSample}.${CLASS}.padded.bam

samtools sort ${outClipSAM}/${inputSample}.${CLASS}.padded.bam -o ${outClipSAM}/${inputSample}.${CLASS}.padded.sort.bam

samtools index ${outClipSAM}/${inputSample}.${CLASS}.padded.sort.bam

# Remove sam and unsorted bam files
rm ${outClipSAM}/${inputSample}.${CLASS}.padded.sam
rm ${outClipSAM}/${inputSample}.${CLASS}.padded.bam

echo "Finished padding CLIP sam file class ${CLASS}"


#===
# FIND AND PAD CHIMERIC BAM FILES
#===

# Find all bams of the input class
shopt -s nullglob
BAMS=( ${INPUT}/mapChimeras/${SAMPLE}.snoRNA.${CLASS}.*.chimeras.bam )
shopt -u nullglob

outChimeraSAM=${OUTPUT}/${CLASS}_chimeraPaddedBam
mkdir -p $outChimeraSAM

# Convert chimeric bam to sam, pad, then convert back to bam
paddedBams=()

for BAM in "${BAMS[@]}"; do
    base=$(basename "$BAM" .bam)

    # Convert to sam
    samtools view -h "$BAM" -o ${outChimeraSAM}/${base}.sam

    # Pad sam file
    awk '{
    if ($1 ~ /^@/) {
      print $0
    } else {
      $4 = $4 + 30; print
    }
    }' FS='\t' OFS='\t' ${outChimeraSAM}/${base}.sam > ${outChimeraSAM}/${base}.padded.sam

    # Convert padded sam to padded bam. Sort & index
    samtools view -b "${outChimeraSAM}/${base}.padded.sam" -o "${outChimeraSAM}/${base}.padded.bam"
    samtools sort "${outChimeraSAM}/${base}.padded.bam" -o "${outChimeraSAM}/${base}.padded.sort.bam"
    samtools index "${outChimeraSAM}/${base}.padded.sort.bam"
    
    # Collect padded bam names into an array
    paddedSortBam="${outChimeraSAM}/${base}.padded.sort.bam"
    paddedBams+=( "$paddedSortBam" )

    rm "${outChimeraSAM}/${base}.sam" "${outChimeraSAM}/${base}.padded.sam" "${outChimeraSAM}/${base}.padded.bam"

done

echo "Have ${#paddedBams[@]} padded chimera BAM(s) for ${CLASS}"

# Get gene IDs from the header of each fasta file, to use as input for clipper
mapfile -t IDS < <(tr -d '\r' < "$FASTA/${CLASS}.fa" \
                   | grep -o '^>[^[:space:]]*' \
                   | sed 's/^>//')

echo "Found ${#IDS[@]} IDs in ${CLASS}.fa"

mamba deactivate

peakDir=${OUTPUT}/${CLASS}_paddedPeaks
mkdir -p $peakDir

mamba activate /n/whipple_lab/share_root/Users/cwhilden/ENVS/clipper_sno

# Run clipper on all padded bams of one class, for each gene ID
parallelJobs="${SLURM_CPUS_PER_TASK:-1}"

for ID in "${IDS[@]}"; do
  echo "Calling peaks for ID=${ID} across ${#paddedBams[@]} BAM(s)…"
  parallel -j "${parallelJobs}" clipper --processors=1 \
    -b {} \
    -g "${ID}" \
    -o "${peakDir}/{/.}.${ID}.bed" \
    -s multimapper_padded \
    ::: "${paddedBams[@]}" > /dev/null
done

mamba deactivate

# Merge adjacent peaks
mamba activate /n/whipple_lab/share_root/Users/cwhilden/ENVS/samtools_bedtools

echo "Merging peaks for ${CLASS}..."

peakBedFiles=${peakDir}/${SAMPLE}*${CLASS}*.bed
mkdir -p $peakDir/merged

for i in $peakBedFiles; do
  base=$(basename $i .bed)
  bedtools sort -i $i | bedtools merge -d 1 -c 4,5,6,7,8 -o distinct,min,distinct,min,max > ${peakDir}/merged/${base}.merged.bed
done