library(testthat)

# Source the function to be tested
source("src/run/helpers/helpers.R")

test_that("resample_given_pvalues works correctly", {
  # Create a sample data frame
  data <- data.frame(
    value = rnorm(100),
    pvalue = runif(100)
  )
  
  # Create a vector of p-values
  pvalues <- data$pvalue
  
  # Test when target_number is greater than or equal to the number of rows in data
  expect_equal(resample_given_pvalues(data, pvalues, 100), data)
  expect_equal(resample_given_pvalues(data, pvalues, 101), data)
  
  # Test when target_number is less than the number of rows in data
  target_number <- 50
  resampled_data <- resample_given_pvalues(data, pvalues, target_number)
  expect_equal(nrow(resampled_data), target_number)
  
  # Test when pvalues length does not match the number of rows in data
  expect_error(resample_given_pvalues(data, pvalues[1:50], target_number), 
               "Length of pvalues must match number of rows in data")
  
  # Test when sig_threshold is specified
  sig_threshold <- 0.01
  resampled_data <- resample_given_pvalues(data, pvalues, target_number, sig_threshold)
  expect_equal(nrow(resampled_data), target_number)
  
  # Test that significant data is retained
  sig_data <- data[pvalues < sig_threshold,]
  expect_true(all(sig_data$value %in% resampled_data$value))
})

test_that("process_sample_sheet works correctly", {
  # Create a sample data frame
  sample_sheet <- data.frame(
    Sex = c("Male", "Female", "Male", "Female"),
    Age = c(30, 25, 40, 35),
    BMI = c(22, 25, 27, 24)
  )
  
  # Create options list
  opts <- list(
    grouping_variable = "Sex",
    covariates = c("Age", "BMI"),
    are_covariates_factors = c(FALSE, TRUE)
  )
  
  # Process the sample sheet
  processed_sample_sheet <- process_sample_sheet(sample_sheet, opts)
  
  # Test that the grouping variable is converted to a factor
  expect_true(is.factor(processed_sample_sheet$Sex))
  
  # Test that the covariates are converted to factors where specified
  expect_false(is.factor(processed_sample_sheet$Age))
  expect_true(is.factor(processed_sample_sheet$BMI))
})

test_that("get_model_matrix works correctly", {
  # Create a sample data frame
  sample_sheet <- data.frame(
    Sex = factor(c("Male", "Female", "Male", "Female")),
    Age = c(30, 25, 40, 35),
    BMI = (c(22, 25, 27, 24))
  )
  
  # Create options list
  opts <- list(
    grouping_variable = "Sex",
    covariates = c("Age", "BMI"),
    are_covariates_factors = c(FALSE, FALSE)
  )
  
  # Get the model matrix
  design <- get_model_matrix(sample_sheet, opts)
  
  # Test that the design matrix has the correct number of rows and columns
  expect_equal(nrow(design), nrow(sample_sheet))
  expect_equal(ncol(design), length(opts$covariates) + 2)  # Intercept + grouping variable + covariates
})
