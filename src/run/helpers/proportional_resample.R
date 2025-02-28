#' Resample Data Based on P-values
#'
#' This function resamples a given dataset based on p-values, targeting a specific number of rows.
#' It separates the data into significant and non-significant subsets based on a significance threshold,
#' and then randomly samples from the non-significant subset to achieve the target number of rows.
#'
#' @param data A data frame containing the data to be resampled.
#' @param pvalues A numeric vector of p-values corresponding to the rows in the data frame.
#' @param target_number An integer specifying the target number of rows in the resampled data.
#' @param sig_threshold A numeric value specifying the significance threshold for p-values. Default is 0.05.
#' @return A data frame containing the resampled data.
resample_given_pvalues <- function(data, pvalues, target_number, sig_threshold=.05) {

  if (length(pvalues) != nrow(data)) {
    stop("Length of pvalues must match number of rows in data")
  }
  if (target_number >= nrow(data)) {
    return(data)
  }
  
  # split data into significant and non-significant subsets
  sig_data <- data[pvalues < sig_threshold,]
  non_sig_data <- data[pvalues >= sig_threshold,]
  n_sig <- nrow(sig_data)
  n_non_sig <- nrow(non_sig_data)

  # work out the proportion of non-significant data to sample
  target_non_sig <- target_number - n_sig
  # sample randomly from non-significant data this will naturally retain the proportions of high/low p-values
  mask <- sample(1:n_non_sig, target_non_sig)
  
  non_sig_sampled <- non_sig_data[mask,]
  resampled_data <- rbind(sig_data, non_sig_sampled)
  return(resampled_data)
}
