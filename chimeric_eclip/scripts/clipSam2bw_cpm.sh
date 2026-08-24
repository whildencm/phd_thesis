#!/usr/bin/env bash
# For converting mapCLIP sam files to bam then to CPM-normalized BW
set -euo pipefail

inSam="$1"
outDir="$2"

mkdir -p "$outDir"

base=$(basename "$inSam" .sam)
echo "Start converting sample $base"

# Convert SAM to BAM and output sorted BAM
samtools view -S -b "$inSam" | samtools sort -o "$outDir/${base}.sort.bam"

# Index sorted BAM
echo "Finished bam sort. Begin indexing"
samtools index "$outDir/${base}.sort.bam"

# Convert sorted BAM to CPM-normalized bigWig
echo "Finished indexing. Begin converting to CPM normalized BW"
bamCoverage \
  -b "$outDir/${base}.sort.bam" \
  -o "$outDir/${base}.clip.cpm.bw" \
  --normalizeUsing CPM \
  --binSize 10 \
  --outFileFormat bigwig

echo "Done: $outDir/${base}.clip.cpm.bw"