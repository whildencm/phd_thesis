#!/bin/bash

sample="$1"
inputDir="$2"
peakSummaryDir="$3"

mkdir -p $peakSummaryDir

outSummary="${peakSummaryDir}/${sample}.repeats.peaks.txt"
saf="${peakSummaryDir}/${sample}.repeats.peaks.saf"
tmp="${saf}.tmp"

mkdir -p "$(dirname "$outSummary")"
: > "$outSummary"
: > "$tmp"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
"filename" \
"sample" \
"targetClass" \
"snoRNA" \
"chromosome" \
"start" \
"end" \
"name" \
"reads in CLIP" \
"reads in INPUT" \
"p-value" \
"chi_val_or_Fisher" \
"Fisher_orChi_square_test" \
"enriched_or_depleted" \
"negative_log10p" \
"log2_fold_change" \
"entropy" \
>> "$outSummary"

printf 'GeneID\tChr\tStart\tEnd\tStrand\n' > "$saf"

found_any=0
for file in "$inputDir/${sample}"*.depad.full; do
    [ -e "$file" ] || continue
    found_any=1
    base=$(basename "$file")

    awk -v base="$base" 'BEGIN{FS=OFS="\t"}
      {
        # Filename pieces: sample.snoRNA.targetClass.target.chimeras.sort.chromosome.norm.bed.depad.full
        split(base, b, /\./)
        sample      = (length(b) >= 1 ? b[1] : "")
        targetClass = (length(b) >= 3 ? b[3] : "")
        target      = (length(b) >= 4 ? b[4] : "")

        chrom = $1
        start = $2
        end   = $3

        # name field has 3–4 colon parts: chrom : start-end : strand [: entropy]
        n = split($4, np, /:/)
        name_region = (n >= 3 ? np[1] ":" np[2] ":" np[3] : $4)
        entropy     = (n >= 4 ? np[n] : "")

        reads_clip   = $5
        reads_input  = $6
        pval         = $7
        chi_or_F_val = $8
        test_type    = $9
        enrich_state = $10
        neglog10p    = $11
        log2fc       = $12

        # Write the full summary row
        print base, sample, targetClass, target, \
              chrom, start, end, name_region, \
              reads_clip, reads_input, pval, chi_or_F_val, test_type, \
              enrich_state, neglog10p, log2fc, entropy

        # Also append the SAF row (GeneID, Chr, Start, End, Strand="+")
        # GeneID is "name" (name_region), Chr is "chromosome" (chrom)
        # Start is "start", End is "end"
        printf "%s\t%s\t%s\t%s\t+\n", name_region, chrom, start, end >> saf_tmp
      }' saf_tmp="$tmp" "$file" >> "$outSummary"
  done

  # De-duplicate SAF rows and finalize
  if (( found_any )); then
    LC_ALL=C sort -u "$tmp" >> "$saf"
  else
    echo "No depadded files for ${sample}; wrote headers only to $outSummary and $saf"
  fi
  rm -f "$tmp"
