library(minfi)
library(limma)
library(jsonlite)
source("helper.R")
opts <- fromJSON(args[2])

args = commandArgs(trailingOnly = T)

metadata <- read.table(file = args[1], header = T, sep = "\t", check.names = F)

# Convert 'Group' to factor with automatically generated levels
metadata = get_updated_metadata(metadata)

# get m-values
m_values = get_m_values(metadata, opts)

# Create design matrix
design = get_model_matrix(metadata)

colnames(design) <- levels(metadata$condition)

# Fit linear model for M-values
fit <- lmFit(m_values, design)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit)

# Get the coefficient
coeffOpt <- get_coeff(coefficients)

# Extract base folder from the first row of metadata$path
base_folder <- gsub("/Donor.*/.*", "/", metadata$path[1])


# Get differential methylation results
results <- topTable(fit2, coef = coeffOpt, number = Inf, adjust = "fdr")

# Write the results to the CSV file
results <- cbind(CpgId = rownames(results), results)

#refactored results
reduced_results = get_refactored_result(results)

# Get annotation results
annotated_results <- get_annotation_result(reduced_results)
# Save annotated results
 get_compressed_result(annotated_results, file.path(base_folder, "results.tsv.gz"))