#!/usr/bin/env Rscript
# For each sample, get chimeras per target class and total chimeras


# SETUP
library(tidyverse, quietly=TRUE)

# Read in chimeric summary file
chimeraSummaryFile <- "/n/whipple_lab/share_root/Users/cwhilden/PROJECTS/THESIS/chapter2/summaryFiles/chimericCounts/All_ChimeraCounts.csv"
outFile <- "/n/whipple_lab/share_root/Users/cwhilden/PROJECTS/THESIS/chapter2/summaryFiles/chimericCounts/ChimerasPerClass_AllSamples.csv"

chimeras <- read_csv(chimeraSummaryFile) %>%
  rename_with(~ sub("\\.sno$", "", .x), -c(snoRNA, target)) %>%
  rename_with(~ sub("\\.fq.gz$", "", .x), -c(snoRNA, target))

# Count chimeras per class
byTarget <- chimeras %>%
  group_by(target) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# Calculate total chimeras
byTarget <- byTarget %>%
  bind_rows(
    byTarget %>%
      summarise(
        across(where(is.numeric), ~ sum(.x, na.rm = TRUE))
      ) %>%
      mutate(target = "total")
  )

write.csv(byTarget,outFile)
