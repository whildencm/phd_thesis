#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse, quietly = TRUE)
  library(tools, quietly = TRUE)
  library(pROC, quietly = TRUE)
})

# =========================
# Import data (CLI args)
# =========================
args <- commandArgs(trailingOnly = TRUE)

usage <- paste0(
  "Usage:\n",
  "  Rscript peakCalling_iN_plots.R <peakDir> <annotFile> <summaryFile> <outDir> <aliasFile>\n\n",
  "Arguments:\n",
  "  peakDir      Directory with output of chimeclip/pipelines/generatePeakSummary.sh\n",
  "  annotFile    Path to annotation CSV of known snoRNA interactions\n",
  "  summaryFile  Path to chimeric count summary CSV (All_ChimeraCounts.csv)\n",
  "  outDir       Output directory\n",
  "  aliasFile    Path to a two-column CSV with headers 'sample' and 'alias'\n"
)

if (length(args) != 5) {
  cat(usage, file = stderr())
  quit(status = 1)
}

peakDir     <- args[[1]]
annotFile   <- args[[2]]
summaryFile <- args[[3]]
outDir      <- args[[4]]
aliasFile   <- args[[5]]

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)


# =========================
# Read alias file
# =========================
aliasTable <- readr::read_csv(aliasFile, show_col_types = FALSE)

if (!all(c("sample", "alias") %in% colnames(aliasTable))) {
  stop(
    "aliasFile must contain exactly the columns 'sample' and 'alias'.\n",
    "Found: ", paste(colnames(aliasTable), collapse = ", ")
  )
}

#===
# Read files in
#===
annot <- read.csv(annotFile) %>%
  mutate(POSITION = suppressWarnings(as.numeric(POSITION)))

dataFiles <- list.files(path = peakDir,
                        pattern = "\\.txt$",
                        full.names = TRUE)

peakNames <- sub("\\.repeats\\.peaks\\.txt$", "", basename(dataFiles))

peakData <- purrr::set_names(dataFiles, peakNames) %>%
  purrr::map(~ readr::read_tsv(.x, show_col_types = FALSE)) %>%
  purrr::map(~ dplyr::filter(.x, enriched_or_depleted == "enriched"))


# Read in chimera count summary file (keep sample columns!)
totalChimeras <- readr::read_csv(summaryFile, show_col_types = FALSE) %>%
  # strip ".sno" suffix from sample columns if present
  rename_with(~ sub("\\.sno$", "", .x), -c(snoRNA, target)) %>%
  # keep snoRNA/target + every other column (sample columns)
  dplyr::select(snoRNA, target, dplyr::everything())

totalChimeras$snoTargetID <- paste0(totalChimeras$snoRNA, "_", totalChimeras$target)

#===
# Peak annotation
#===
annotatePeaks <- function(data,
                          snoAnnot = annot) {

  data$interactionID <- paste0(data$snoRNA, "_", data$chromosome)
  snoAnnot$interactionID <- paste0(snoAnnot$SNO, "_", snoAnnot$TARGET)

  annotLookup <- snoAnnot %>%
    group_by(interactionID) %>%
    summarise(
      knownSite = paste0(sort(unique(as.character(POSITION))), collapse = ";"),
      .groups = "drop"
    )

  data2 <- data %>%
    left_join(annotLookup, by = "interactionID")

  data3 <- data2 %>%
    rowwise() %>%
    mutate(
      peakStatus = {
        if (is.na(knownSite) || knownSite == "") {
          "unknown"
        } else {
          sites <- as.numeric(strsplit(knownSite, ";")[[1]])
          if (any(sites >= start & sites <= end)) "known" else "unknown"
        }
      }
    ) %>%
    ungroup()

  data3$snoTargetID <- paste0(data3$snoRNA, "_", data3$targetClass)

  return(data3)
}

peakDataAnnot <- imap(peakData, ~ annotatePeaks(.x, snoAnnot = annot))

fractionOfChimeras <- function(peakDataAnnot, sample, totalChimeras, outDir) {

  data <- peakDataAnnot[[sample]]

  chimeraCounts <- totalChimeras %>%
    select(snoTargetID, any_of(sample)) %>%
    rename(totalCount = any_of(sample))

  data2 <- left_join(data, chimeraCounts, by = "snoTargetID")

  data3 <- data2 %>%
    mutate(
      clipNumeric = as.numeric(`reads in CLIP`),
      chi_val_or_Fisher = as.numeric(chi_val_or_Fisher),
      fraction = clipNumeric / totalCount
    )

  write.csv(data3, file.path(outDir, paste0(sample, "_Significant_Peak_Results.csv")))

  data3
}

peakDataNorm <- imap(peakDataAnnot, ~ fractionOfChimeras(
  peakDataAnnot = peakDataAnnot,
  sample = .y,
  totalChimeras = totalChimeras,
  outDir = outDir
))

# =========================
# Build samples table from alias file
# =========================
samples <- data.frame(sample = names(peakDataNorm), stringsAsFactors = FALSE) %>%
  left_join(aliasTable, by = "sample")

# Warn about any samples missing from the alias file; fall back to sample name
missing_alias <- samples$sample[is.na(samples$alias)]
if (length(missing_alias) > 0) {
  warning(
    "The following samples have no entry in aliasFile and will keep their ",
    "original name as the alias:\n  ",
    paste(missing_alias, collapse = "\n  ")
  )
  samples$alias[is.na(samples$alias)] <- samples$sample[is.na(samples$alias)]
}

colsToKeep <- c("sample", "targetClass", "snoRNA", "chromosome","start", "end",
                "negative_log10p", "log2_fold_change", "peakStatus",
                "knownSite", "fraction")

# Columns that must be numeric - coerce explicitly to avoid type-clash across samples
numericCols <- c("start", "end", "negative_log10p", "log2_fold_change",
                 "fraction", "sig_fraction")

allBound <- imap_dfr(peakDataNorm, ~ {
  i <- match(.y, samples$sample)
  if (is.na(i)) stop("Sample not found in samples table: ", .y)
  
  .x %>%
    select(any_of(colsToKeep)) %>%
    # Coerce to numeric first so all chunks share the same types before bind_rows
    mutate(across(any_of(numericCols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(sampleAlias = samples$alias[i])
})

write.csv(allBound, file.path(outDir, "All_Samples_Significant_Peak_Results.csv"))

message("Done. Wrote per-sample CSVs and All_Samples_Significant_Peak_Results.csv to: ", outDir)