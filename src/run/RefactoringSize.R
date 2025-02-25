library(dplyr)
file_path <- "differential_methylation_results.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)

log2_col <- "logFC"
pval_col <- "adj.P.Val"
data$neg_log10_p <- -log10(data[[pval_col]])

# Split into significant and insignificant
sig_data <- data %>% filter(.data[[pval_col]] < 0.05)
insig_data <- data %>% filter(.data[[pval_col]] >= 0.05)
n_sig <- nrow(sig_data)

# Define bins for insignificant points
log2_bins <- seq(-5, 5, length.out = 21)
pval_bins <- seq(0, 3, length.out = 11)


# Add bin labels to insignificant data
insig_data$log2_bin <- cut(insig_data[[log2_col]], breaks = log2_bins, include.lowest = TRUE)
insig_data$pval_bin <- cut(insig_data$neg_log10_p, breaks = pval_bins, include.lowest = TRUE)

# Target total rows: ~30k-45k; assume ~5k significant, subsample ~25k-40k insignificant
target_total <- 35000
target_insig <- target_total - n_sig
n_bins <- length(log2_bins[-1]) * length(pval_bins[-1])
samp_per_bin <- max(50, ceiling(target_insig / n_bins))

# Sample from each bin
set.seed(42)  # For reproducibility
insig_sampled <- insig_data %>%
  group_by(log2_bin, pval_bin) %>%
  sample_n(size = min(n(), samp_per_bin), replace = FALSE) %>%
  ungroup()

# Combine significant and sampled insignificant data
final_data <- bind_rows(sig_data, insig_sampled) %>%
  select(-log2_bin, -pval_bin, -neg_log10_p)
n_final <- nrow(final_data)

# Overwrite the original file
write.csv(final_data, "differential_methylation_results_reduced.csv"
, row.names = FALSE)