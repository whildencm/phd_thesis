#!/bin/bash
#!/bin/bash
#SBATCH -n 8             # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH --mail-user=courtneywhilden@g.harvard.edu
#SBATCH --mail-type=END
#SBATCH -t 2-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=100G          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p whipple
#SBATCH -o slurm-logs/eCLIP/%j.o    		# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e slurm-logs/eCLIP/s%j.e    		# File to which STDERR will be written, %j inserts jobid
hostname

#==== Positional arguments

# eCLIP alignment script C Whilden
# Usage: ./script.sh <fastq> <outDir> <genome>

fastq="$1"
outDir="$2"
genome="$3"

# Validate arguments
if [[ -z "$fastq" || -z "$outDir" || -z "$genome" ]]; then
    echo "Usage: $0 <fastq_file> <output_directory> <genome_index>"
    exit 1
fi

mkdir -p "$outDir"

# Extract sample name from file name
base=$(basename "$fastq" ".trim.fastq.gz")

STAR \
--alignEndsType EndToEnd \
--genomeDir "$genome" \
--genomeLoad NoSharedMemory \
--outBAMcompression 10 \
--outFileNamePrefix "${outDir}/${base}" \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 1 \
--outFilterScoreMin 10 \
--outFilterType BySJout \
--outSAMattrRGline ID:foo \
--outSAMattributes All \
--outSAMmode Full \
--outSAMtype BAM Unsorted \
--outStd Log \
--readFilesIn "$fastq" \
--runMode alignReads \
--runThreadN 8
