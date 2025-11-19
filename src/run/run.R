library(minfi)
library(limma)
library(jsonlite)
source("helper.R")

args = commandArgs(trailingOnly = TRUE)

inputFilePath <- args[1]
outputFilePath <- args[3]
optionsFilePath <- args[2]

options <- fromJSON(optionsFilePath)

metadata <- read.table(file = inputFilePath, header = T, sep = "\t", check.names = F)

# Convert 'Group' to factor with automatically generated levels
metadata = get_updated_metadata(metadata)

# get m-values
m_values = get_m_values(metadata, options)

# Create design matrix
design = get_model_matrix(metadata)

colnames(design) <- levels(metadata$condition)

# Fit linear model for M-values
fit <- lmFit(m_values, design)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit)

# Get the coefficient
coeffOpt <- get_coeff(coefficients)

# Get differential methylation results
results <- topTable(fit2, coef = coeffOpt, number = Inf, adjust = "fdr")

# Binding CpgId to results
results <- cbind(CpgId = rownames(results), results)

#refactored results
reduced_results = get_refactored_result(results)
rm(results)

# Get annotation results
annotated_results <- get_annotation_result(reduced_results)
rm(reduced_results)

# Save annotated results
write.table(annotated_results, outputFilePath, row.names = FALSE, quote = FALSE, sep = "\t")
rm(annotated_results)
gc()