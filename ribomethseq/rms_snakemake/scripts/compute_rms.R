log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

# ===
# INPUTS
# ===

count3_file   <- snakemake@input[["count3"]]
coverage_file <- snakemake@input[["coverage"]]
output_file   <- snakemake@output[["results"]]

# ===
# LOAD DATA
# ===

data3 <- read.table(count3_file,   col.names = c("chrom", "position", "count"))  # added chrom
cov   <- read.table(coverage_file, col.names = c("chrom", "position", "depth"))

# ===
# PROCESS END COUNTS
# ===

# Merge by BOTH chrom and position to avoid cross-chromosome contamination
counts3 <- merge(cov, data3, by = c("chrom", "position"), all.x = TRUE)
counts3 <- counts3[order(counts3$chrom, counts3$position), ]
counts3$count[is.na(counts3$count)] <- 0

# ===
# COMPUTE RMS SCORE
# ===

W       <- seq(from = 1, to = 0.5, by = -0.1)
offsets <- 1:6
Wsum    <- sum(W)

compute_rms <- function(end_counts) {
  ni <- end_counts
  n  <- length(ni)

  nj_pos <- outer(offsets, 1:n, function(d, p) p - d)
  nk_pos <- outer(offsets, 1:n, function(d, p) p + d)

  nj_pos[nj_pos < 1] <- NA
  nk_pos[nk_pos > n] <- NA  # use n, not genome_length — windows must not cross chromosome ends

  nj_counts <- matrix(ni[nj_pos], nrow = length(offsets))
  nk_counts <- matrix(ni[nk_pos], nrow = length(offsets))

  Wnj <- colSums(nj_counts * W, na.rm = FALSE)
  Wnk <- colSums(nk_counts * W, na.rm = FALSE)

  rms             <- 1 - (ni / (0.5 * ((Wnj / Wsum) + (Wnk / Wsum))))
  rms[rms < 0]   <- 0
  rms
}

# Apply RMS per chromosome so windows never bleed across boundaries
results_list <- lapply(split(counts3, counts3$chrom), function(chr_df) {
  chr_df       <- chr_df[order(chr_df$position), ]
  chr_df$RMS   <- compute_rms(chr_df$count)
  chr_df
})

results <- do.call(rbind, results_list)

# ===
# WRITE OUTPUT
# ===

results <- data.frame(
  chrom      = results$chrom,
  position   = results$position,
  coverage   = results$depth,
  end_counts = results$count,
  RMS        = results$RMS
)

write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)