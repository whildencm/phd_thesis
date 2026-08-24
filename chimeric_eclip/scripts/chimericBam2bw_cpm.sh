#!/usr/bin/env bash
# For an input .bam and chimeric count summary, export a chimeras per million normalized BW file
# CW2025
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 <totalChimerasFile> <sample> <bamFile> <outFile> <logFile>

Example:
  $0 /path/ChimerasPerClass_AllSamples.csv SSEXPT8_IP_A_S3_R1_001 /path/sample.bam sample.CPM.bw /path/sample.CPM.log
EOF
  exit 2
}

if [ "$#" -ne 5 ]; then
  usage
fi

totalChimerasFile="$1"
sample="$2"
bamFile="$3"
outFile="$4"
logFile="$5"

# Create log directory if needed and initialise log file
mkdir -p "$(dirname "$logFile")"
: > "$logFile"   # truncate/create

# Redirect all stdout and stderr to log file AND terminal
exec > >(tee -a "$logFile") 2>&1

echo "================================================"
echo "  chimericCPM_bw.sh"
echo "  $(date)"
echo "================================================"

# Function to get the total number of chimeras from the summary file
getTotalChimeras() {
  local csvFile="$1"
  local sampleName="$2"

  awk -v samp="$sampleName" -F',' '
    BEGIN { found = 0; col = 0 }
    NR==1 {
      for(i=1;i<=NF;i++){
        h=$i; gsub(/^"|"$/, "", h)
        if(h == samp){ col=i; break }
      }
      if(col==0){ print "ERROR:sample_not_found"; exit 2 }
      next
    }
    {
      sec = $2; gsub(/^"|"$/, "", sec)
      if(sec == "total"){
        val = $col
        gsub(/^"|"$/, "", val)
        gsub(/^[ \t]+|[ \t]+$/, "", val)
        print val
        found = 1
        exit
      }
    }
    END {
      if(found==0) { print "ERROR:total_row_not_found"; exit 3 }
    }
  ' "$csvFile"
}

# Get total chimeras
totalChimeras=$(getTotalChimeras "$totalChimerasFile" "$sample") || true

# Compute CPM scale factor (floating point)
scaleFactor=$(awk -v t="$totalChimeras" 'BEGIN { printf "%.10f", 1000000 / t }')

echo "Sample:              $sample"
echo "Total chimeras:      $totalChimeras"
echo "Scale factor (CPM):  $scaleFactor"
echo "BAM:                 $bamFile"
echo "Output bigWig:       $outFile"
echo "Log file:            $logFile"
echo "------------------------------------------------"

# Run bamCoverage
bamCoverage \
  -b "$bamFile" \
  -o "$outFile" \
  --scaleFactor "$scaleFactor" \
  --normalizeUsing None \
  --binSize 10 \
  --outFileFormat bigwig

echo "------------------------------------------------"
echo "Done: wrote $outFile"
echo "$(date)"
echo "================================================"