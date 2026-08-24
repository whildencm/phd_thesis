#!/usr/bin/env Rscript

# C. Whilden 2025
# Goal:
# Export all chimeric counts into a single wide data frame

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

## ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript exportChimeraCounts.R <inputDir> <outputDir>\n",
    "Example:\n",
    "  Rscript exportChimeraCounts.R mapChimeras chimericCountSummary"
  )
}

inputDir <- args[1]
outputDir <- args[2]

## ---- Find chimera count files ----
chimeraFiles <- list.files(
  path = inputDir,
  pattern = "\\.individual\\.chimeras\\.count\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(chimeraFiles) == 0) {
  stop("No .individual.chimeras.count.tsv files found in inputDir")
}

## ---- Read and combine ----
chimeraData <- set_names(chimeraFiles, basename(chimeraFiles)) %>%
  map(~ read_csv(.x, show_col_types = FALSE))

allChimeraDf <- bind_rows(chimeraData) %>%
  filter(
    uID != "uID",
    snoRNA != "snoRNA",
    target != "target"
  ) %>%
  mutate(
    chimeras_count = as.numeric(chimeras_count),
    uID = sub("(\\.sno).*", "\\1", uID)
  )

## ---- Pivot wider (uID = columns) ----
wideDf <- allChimeraDf %>%
  select(snoRNA, target, uID, chimeras_count) %>%
  pivot_wider(
    names_from = uID,
    values_from = chimeras_count,
    values_fill = 0
  )

## ---- Write output ----
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

outFile <- file.path(outputDir, "All_ChimeraCounts.csv")
write.csv(wideDf, outFile, row.names = FALSE)

message("Done.")
message("Wrote file: ", outFile)
