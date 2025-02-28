library(testthat)

# Source the function to be tested
source("../src/run/helpers/proportional_resample.R")

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
