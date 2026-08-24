# C. Whilden 2025
# Goal:
# Export lists of diferentially translated genes for sample subsets

library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(EnhancedVolcano)
library(biomaRt)
library(reshape2)
library(rtracklayer)

# All output files will save here
base_dir <- "/n/lab_storage/whipple_lab/share_root/Users/cwhilden/PROJECTS/manuscript_2026/riboseq"

# Get all Feature Counts files
featurecount_files <- list.files(path = paste0(base_dir,"results/counts"), 
                                 pattern = "\\.featurecounts$", 
                                 full.names = TRUE)

# Create a named list of data frames
count_data <- set_names(featurecount_files, 
                        basename(featurecount_files)) %>% 
  map(~ read_tsv(.x, skip = 1))

# Read in sample IDs and gtf information
sample_id_df <- read.csv(file.path(paste0(base_dir,"config/"), 
                                   "RPF_SampleID.csv"))
sample_id_TE <- read.csv(file.path(paste0(base_dir,"config/", 
                                          "RPF_SampleID_deltaTE.csv")))
sample_names <- sample_id_df$X

coding_genes <- as.data.frame(read.csv(file.path(paste0(base_dir,"config/"), 
                                                 "Protein_coding_genes_GTF.csv")))

gene_map <- read.csv(file.path(paste0(base_dir,"config/"),
                               "Ensembl-to-symbol.csv"), 
                     stringsAsFactors = FALSE)
gene_label_lookup <- setNames(gene_map$Gene, gene_map$Ensembl)


# Generate count matrices from feature counts results -------------------

# Function to extract count matrix by type ("RPF" or "mRNA")
get_counts <- function(type) {
  
  # Validate gene order across all samples
  gene_ids_list <- map(sample_names, function(name) {
    file_key <- paste0(name, ".", type, ".featurecounts")
    if (!file_key %in% names(count_data)) {
      stop(paste("Missing file for:", file_key))
    }
    count_data[[file_key]]$Geneid
  })
  
  # Check that all gene ID vectors are equal
  reference_gene_ids <- gene_ids_list[[1]]
  for (i in seq_along(gene_ids_list)) {
    if (!all(gene_ids_list[[i]] == reference_gene_ids)) {
      stop(paste("Gene order mismatch in sample:", sample_names[i]))
    }
  }
  
  # Extract count column from Feature Counts (column 7)
  count_matrix <- map(sample_names, function(name) {
    file_key <- paste0(name, ".", type, ".featurecounts")
    count_data[[file_key]][[7]]
  }) %>%
    do.call(cbind, .)
  
  # Set row/column names
  colnames(count_matrix) <- sample_names
  rownames(count_matrix) <- gsub("\\.\\d+$", "", reference_gene_ids)
  
  return(as.matrix(count_matrix))
}

# Generate and validate count matrices
riboCount <- get_counts("RPF")
mRNACount <- get_counts("mRNA")

# Validate sample order matches metadata
if (!all(colnames(riboCount) == sample_names)) {
  stop("Sample order mismatch in riboCount matrix.")
}
if (!all(colnames(mRNACount) == sample_names)) {
  stop("Sample order mismatch in mRNACount matrix.")
}

# Filter RPKM matrices for protein coding genes 
rna_pc <- mRNACount[rownames(mRNACount) %in% coding_genes$gene_id, ] %>%
  as.data.frame()

ribo_pc <- riboCount[rownames(riboCount) %in% coding_genes$gene_id, ] %>%
  as.data.frame()

ribo_names <- sample_id_TE %>% filter(SeqType=="RIBO")
colnames(ribo_pc) <- ribo_names$X


# Run DESeq on sample subsets -----------------------------
# Function to generate sub-sample ID data frames as input to the TE calculation
subset_samples <- function(df, age = NULL, sex = NULL, genotype = NULL, litter = NULL) {
  out <- df
  if (!is.null(age)) out <- out[out$Age == age, ]
  if (!is.null(sex)) out <- out[out$Sex == sex, ]
  if (!is.null(genotype)) out <- out[out$Genotype == genotype, ]
  if (!is.null(litter)) out <- out[out$Litter == litter, ]
  rownames(out) <- out$X
  return(out)
}

P0_All <- subset_samples(sample_id_TE, age = "P0")
#P0_F   <- subset_samples(sample_id_TE, age = "P0", sex = "F")
#P0_M   <- subset_samples(sample_id_TE, age = "P0", sex = "M")

P7_All <- subset_samples(sample_id_TE, age = "P7")
#P7_F   <- subset_samples(sample_id_TE, age = "P7", sex = "F")
#P7_M   <- subset_samples(sample_id_TE, age = "P7", sex = "M")


# Function to run DESeq on given sample set
# Exports results matrix and .rnk file for GSEA
deltaTE <- function(samples,
                    rna_pc,
                    ribo_pc,
                    rna_cutoff,
                    ribo_cutoff,
                    gene_label_lookup,
                    base_dir,
                    file_title,
                    design_model) {
  
  # Select samples from rna and ribo count matrices
  rna_samples <- samples %>% filter(SeqType=="RNA")
  rna_subset <- rna_pc[rna_samples$X] 
  ribo_samples <- samples %>% filter(SeqType=="RIBO")
  ribo_subset <- ribo_pc[ribo_samples$X]
  
  # Filter lowly expressed genes
  rna_filtered <- rna_subset[rowSums(rna_subset > rna_cutoff) > ncol(rna_subset) / 2, ]
  ribo_filtered <- ribo_subset[rowSums(ribo_subset > ribo_cutoff) > ncol(ribo_subset) / 2, ]
  
  # Find genes that meet both rna and ribo cutoff filters, and combine into a final df
  common_genes <- intersect(rownames(rna_filtered), rownames(ribo_filtered))
  rna_final <- rna_filtered[common_genes, ]
  ribo_final <- ribo_filtered[common_genes, ]
  
  # DESeq
  
  TE_Mat = DESeqDataSetFromMatrix(countData=cbind(rna_final,ribo_final),
                                  colData=samples,
                                  design= design_model )
  
  # Set factor levels and run DESeq
  TE_Mat$Genotype <- relevel(TE_Mat$Genotype, ref = "WT")
  TE_Mat$SeqType <- relevel(TE_Mat$SeqType, ref = "RNA")
  TE_Mat <- DESeq(TE_Mat)
  
  TE_genotype = results(TE_Mat, name="GenotypeMUT.SeqTypeRIBO")
  
  TE <- as.data.frame(TE_genotype)
  TE$symbol <- gene_label_lookup[rownames(TE)]
  
  write.csv(TE, file.path(base_dir, paste0(file_title,".csv")))
  
  #Output a .rnk file as input to GSEA
  GSEA <- NULL
  GSEA$Gene <- TE$symbol
  GSEA$Rank <- TE$stat
  
  write.table(
    GSEA,
    file = file.path(base_dir, paste0(file_title,".rnk")),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE, row.names=FALSE
  )
  
}

# Run on each sample group
# List of configurations
vars <- list(
  list(name="P0_All", design=~Sex + Litter + Genotype + SeqType + Sex:SeqType + Litter:SeqType + Genotype:SeqType),
  list(name="P7_All", design=~Sex + Litter + Genotype + SeqType + Sex:SeqType + Litter:SeqType + Genotype:SeqType)#,
  #list(name="P0_F",   design=~Litter + Genotype + SeqType + Litter:SeqType + Genotype:SeqType),
  #list(name="P0_M",   design=~Litter + Genotype + SeqType + Litter:SeqType + Genotype:SeqType),
  #list(name="P7_F",   design=~Litter + Genotype + SeqType + Litter:SeqType + Genotype:SeqType),
  #list(name="P7_M",   design=~Litter + Genotype + SeqType + Litter:SeqType + Genotype:SeqType)
)

# Loop through configs
for (var in vars) {
  deltaTE(
    get(var$name),      
    rna_pc,
    ribo_pc,
    rna_cutoff = 10,
    ribo_cutoff = 10,
    gene_label_lookup,
    paste0(base_dir,"results/translationEfficiency"),
    paste0(var$name,"_dTE_Results"),
    design_model = var$design
  )
}


######## RNA only ############

deltaRNA <- function(samples,
                    rna_pc,
                    ribo_pc,
                    rna_cutoff,
                    ribo_cutoff,
                    gene_label_lookup,
                    base_dir,
                    file_title,
                    design_model) {
  
  # Select samples from rna and ribo count matrices
  rna_samples <- samples %>% filter(SeqType=="RNA")
  rna_subset <- rna_pc[rna_samples$X] 
 
  
  # Filter lowly expressed genes
  rna_filtered <- rna_subset[rowSums(rna_subset > rna_cutoff) > ncol(rna_subset) / 2, ]
  
  # DESeq
  
  TE_Mat = DESeqDataSetFromMatrix(countData=cbind(rna_filtered),
                                  colData=rna_samples,
                                  design= design_model )
  
  # Set factor levels and run DESeq
  TE_Mat$Genotype <- relevel(TE_Mat$Genotype, ref = "WT")
  TE_Mat <- DESeq(TE_Mat)
  
  TE_genotype = results(TE_Mat, name="Genotype_MUT_vs_WT")
  
  TE <- as.data.frame(TE_genotype)
  TE$symbol <- gene_label_lookup[rownames(TE)]
  
  write.csv(TE, file.path(base_dir, paste0(file_title,".csv")))
  
  #Output a .rnk file as input to GSEA
  GSEA <- NULL
  GSEA$Gene <- TE$symbol
  GSEA$Rank <- TE$stat
  
  write.table(
    GSEA,
    file = file.path(base_dir, paste0(file_title,".rnk")),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE, row.names=FALSE
  )
  
}

RNAvars <- list(
  list(name="P0_All", design=~Sex + Litter + Genotype + Litter:Genotype + Sex:Litter + Sex:Genotype),
  list(name="P7_All", design=~Sex + Litter + Genotype + Litter:Genotype + Sex:Litter + Sex:Genotype))

# Loop through configs
for (var in RNAvars) {
  deltaRNA(
    get(var$name),      
    rna_pc,
    ribo_pc,
    rna_cutoff = 10,
    ribo_cutoff = 10,
    gene_label_lookup,
    base_dir= paste0(base_dir,"results/translationEfficiency"),
    paste0(var$name,"_dRNA_Results"),
    design_model = var$design
  )
}
