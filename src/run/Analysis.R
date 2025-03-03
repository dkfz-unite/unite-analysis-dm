library(minfi)
library(limma)
library(jsonlite)
source("src/run/helpers/helpers.R")
opts <- fromJSON("opts.json")
# Read Inputs
sample_info <- read.csv("sample_sheet.csv")
donor_metadata <- read.csv("donor_metadata.csv")

# Merge IDAT and donor information
sample_sheet <- merge(sample_info, donor_metadata, by = "Sample_ID", all.x = TRUE)

m_values = get_m_values(sample_sheet, opts)
design = get_model_matrix(sample_sheet, opts)

# Fit linear model for M-values
fit <- lmFit(m_values, design)
# Apply empirical Bayes moderation
fit2 <- eBayes(fit)


# Get differential methylation results
results <- topTable(fit2, coef = "SexGroup2", number = Inf, adjust = "fdr")
write.csv(results, "differential_methylation_refactored.csv", row.names = TRUE)