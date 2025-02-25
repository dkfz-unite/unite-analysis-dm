library(minfi)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  

# Read Inputs
sample_info <- read.csv("sample_sheet.csv")
donor_metadata <- read.csv("donor_metadata.csv")

# Merge IDAT and donor information
sample_sheet <- merge(sample_info, donor_metadata, by = "Sample_ID", all.x = TRUE)

# Convert Sex to factor with meaningful labels
sample_sheet$Sex <- factor(sample_sheet$Sex, levels = c("Female", "Male"), labels = c("Group1", "Group2"))

# Read IDAT files
RGset <- read.metharray(sample_sheet$Basename, extended = TRUE)

# Preprocess the data and extract beta and m-values
Mset <- preprocessIllumina(RGset)
beta_values <- getBeta(Mset)
m_values <- getM(Mset)

# Create design matrix
design <- model.matrix(~ Sex, data = sample_sheet)

# Fit linear model for M-values
fit <- lmFit(m_values, design)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit)

# Get differential methylation results
results <- topTable(fit2, coef = "SexGroup2", number = Inf, adjust = "fdr")
write.csv(results, "differential_methylation_refactored.csv", row.names = TRUE)