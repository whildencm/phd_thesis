#!/bin/bash
# Merge peak replicates CWhilden 2026


# Positional arguments
inputPath="$1"  # peak calling output directory
repA="$2"   # sample names to merge
repB="$3"   # sample names to merge
groupName="$4"  # group name for merged peaks

# Define target classes
classes=(rRNA tRNA snRNA scaRNA)

# Merge replicates & output a 'true merge' and an 'intersection' file
for class in "${classes[@]}"; do
    peakDir="${inputPath}/${class}_paddedPeaks/merged"
    outDir="${inputPath}/${class}_paddedPeaks/${groupName}_ReplicateCoordinates"
    mkdir -p $outDir
    
    # Create a reference from one replicate to pull out target basenames
    peaksRef="${peakDir}/${repA}*.merged.bed"

    # merge each bed file
    for i in $peaksRef; do
        base=$(basename $i)
        newName=${base#${repA}}

        fileA="$peakDir/${repA}${newName}"
        fileB="$peakDir/${repB}${newName}"

        outMg="${outDir}/${groupName}${newName}_FullMerge"

        # require both files to exist and be not empty
        if [[ ! -s "$fileA" || ! -s "$fileB" ]]; then
            echo "[SKIP] $newName  missing/empty: A=$fileA B=$fileB" >&2
            continue
        fi

        # Output "full merged" file
        concatFile="$outDir/ReplicateCoordinates${newName}.concat"
        cat "$fileA" "$fileB" > "$concatFile"

        
        bedtools sort -i "$concatFile" \
        | bedtools merge -d 1 -c 4,5,6,7,8 -o distinct,min,distinct,min,max \
        > $outMg

        # Output "intersection" file
        bedtools intersect -a "$fileA" -b "$fileB" \
            | bedtools sort -i - \
            | bedtools merge -d 1 -c 4,5,6,7,8 -o distinct,min,distinct,min,max \
            > "${outDir}/${groupName}${newName}_Intersect"


        rm $concatFile
    done

done