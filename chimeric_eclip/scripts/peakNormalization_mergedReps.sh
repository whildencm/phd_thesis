#!/bin/bash
#SBATCH --array=1 #Number of jobs, corresponding to number of samples
#SBATCH --cpus-per-task=48 #Cores per task
#SBATCH -t 2-00:00 #days-hours:minutes
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -o slurm-logs/%x_%A_%a.o #output file name
#SBATCH -e slurm-logs/%x_%A_%a.e #error file name
#SBATCH --mem=184G #to specify memory alloted for job (whipple partition has 184G per core)
hostname

####################################
# ====
# Positional arguments
# =====
inputDir="$1" # Base output directory of peakCallingRepeats.sh
peakBedDir="$2" # Directory with peak .bed files for chosen target class
outputDir="$3" # Directory where normalized peak output files should go
class="$4"  # Repeat target class
sample="$5" # IP sample ID
inputSample="$6" # Input sample ID

echo "[INFO] Starting job for IP=$sample input=$inputSample class=$class"
echo "[INFO] inputDir=$inputDir"
echo "[INFO] peakBedDir=$peakBedDir"
echo "[INFO] outputDir=$outputDir"

###
# Load modules
###
module load Mambaforge
mamba activate IDR_py38

###
# Prepare output directories
###
mkdir -p $outputDir
clipCountDir="${outputDir}/clipCountFiles"
chimCountDir="${outputDir}/chimeraCountFiles"
normDir="${outputDir}/peaksNorm"

parallelJobs="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$clipCountDir" "$chimCountDir"  "$normDir"

###
# Generate CLIP count file for target class
###
clipBam="${inputDir}/${class}_clipPaddedBam/${inputSample}.${class}.padded.sort.bam"

echo "Counting reads in" ${clipBam}

mkdir -p $clipCountDir

clipCountFile="${clipCountDir}/${inputSample}.${class}.count.txt"
samtools view -c "$clipBam" > "$clipCountFile"

# Fail fast if samtools produces an empty count file
if [[ ! -s "$clipCountFile" ]]; then
  echo "[ERROR] Clip count file is empty: $clipCountFile" >&2
  exit 1
fi

echo "Finished counting reads in" ${clipBam}

###
# Run IDR normalization
###
chimeraBamDir="${inputDir}/${class}_chimeraPaddedBam"

beds=( "$peakBedDir"/*.merged.bed_Intersect )
(( ${#beds[@]} > 0 )) || { echo "[ERROR] No beds found in $peakBedDir" >&2; exit 1; }

echo "[INFO] Found ${#beds[@]} peak BED files"

processOneBed() {
    # Function to build chimeric .bam path from .bed file, 
    # generate chimeric count file,
    # and run IDR normalization

  bed="$1"
  
  base=$(basename "$bed")

  echo "[INFO][$sample] Processing BED: $base"

  core=$(awk -F. '{
    printf ".";
    for (i=2; i<=NF-3; i++) printf "%s%s", (i>2?".":""), $i;
    print ""
  }' <<< "$base")

  chimeraBam="${chimeraBamDir}/${sample}${core}.bam"
  [[ -s "$chimeraBam" ]] || { echo "[SKIP] chimera bam missing: $chimeraBam" >&2; return 0; }

  chimCountFile="${chimCountDir}/${sample}${core}.count.txt"
  if [[ ! -s "$chimCountFile" ]]; then
    echo "[INFO][$sample$core] Counting chimera reads"
    samtools view -c "$chimeraBam" > "$chimCountFile"
  fi

  outPrefix="${normDir}/${sample}${core}"
  
  echo "[INFO][$sample$core] Running overlap_peakfi_with_bam.pl"
  perl scripts/overlap_peakfi_with_bam.pl \
    "$chimeraBam" \
    "$clipBam" \
    "$bed" \
    "$chimCountFile" \
    "$clipCountFile" \
    "$outPrefix"

  echo "[INFO][$sample$core] Finished normalization"

  # Unpad the peak output files
  echo "[INFO][$sample$core] Un-padding normalized peaks"

  # BED6 output
  if [[ -s "$outPrefix" ]]; then
    awk 'BEGIN{OFS="\t"} {
      $2 = $2 - 30;
      $3 = $3 - 30;
      if ($2 < 0) $2 = 0;
      if ($3 < 0) $3 = 0;
      print
    }' "$outPrefix" > "${outPrefix}.depad"
  fi

  # .full output
  if [[ -s "${outPrefix}.full" ]]; then
    awk 'BEGIN{OFS="\t"} {
      $2 = $2 - 30;
      $3 = $3 - 30;
      if ($2 < 0) $2 = 0;
      if ($3 < 0) $3 = 0;
      print
    }' "${outPrefix}.full" > "${outPrefix}.depad.full"
  fi

  echo "[INFO][$sample$core] Finished un-padding"

}

export -f processOneBed
export sample class inputDir peakBedDir outputDir chimeraBamDir clipBam clipCountFile chimCountDir normDir

# run N jobs in parallel (pick based on cores)
echo "[INFO] Launching normalization with $parallelJobs parallel jobs"

parallel -j "${parallelJobs}" processOneBed ::: "${beds[@]}"