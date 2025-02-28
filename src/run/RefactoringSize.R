library(dplyr)
source("src/run/helpers/proportional_resample.R")
file_path <- "differential_methylation_results.csv"
data <- read.csv(file_path, stringsAsFactors = FALSE)

log2_col <- "logFC"
pval_col <- "adj.P.Val"
set.seed(123)
resampled_data <- resample_given_pvalues(data, data[[pval_col]], 35000, sig_threshold=.05)
select(-log2_bin, -pval_bin, -neg_log10_p)

# Overwrite the original file
write.csv(final_data, "differential_methylation_results_reduced.csv"
, row.names = FALSE)