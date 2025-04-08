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

colnames(design) <- levels(metadata$conditions)

# Fit linear model for M-values
fit <- lmFit(m_values, design)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit)

# Get the coefficient
coeffOpt <- get_coeff(coefficients)

# Extract base folder from the first row of metadata$path
base_folder <- gsub("/Donor.*/.*", "/", metadata$path[1])


# Get differential methylation results
coeff <- topTable(fit2, coeffOpt, number = Inf, adjust = "fdr")
results <- file.path(base_folder, "results.csv")

# Write the results to the CSV file
write.csv(coeff, results, row.names = TRUE)

#refactored results
resultdata <- read.csv(results, fileEncoding = "UTF-8")
reduced_results = get_refactored_result(resultdata)
# Write the refactored results to the CSV file
write.csv(reduced_results, file.path(base_folder, "results_reduced.csv"), row.names = TRUE)