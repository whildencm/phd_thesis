# RiboMethSeq Pipeline

A Snakemake pipeline for processing RiboMethSeq (RMS) data to quantify 2'-O-methylation
of RNA using single-end sequencing reads.

---

## Overview

This pipeline takes single-end FASTQ files and produces per-position RMS scores for
one or more reference genomes. The steps are:

1. **Index genome** — Build Bowtie2 index from reference FASTA
2. **Trim adapters** — Remove adapter sequences with Trim Galore
3. **Map reads** — Align trimmed reads to the genome with Bowtie2, convert to BED
4. **Count 3' ends** — Tally 3' read ends at each genomic position (strand-aware)
5. **Genome sizes** — Extract chromosome sizes from reference FASTA
6. **Genome coverage** — Compute per-position read depth with bedtools
7. **Compute RMS scores** — Calculate RMS scores from end counts and coverage (R)
