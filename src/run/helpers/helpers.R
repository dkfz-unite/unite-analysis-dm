library(minfi)
opts <- fromJSON("opts.json")
annotation_package <- opts[["annotation_package"]]
library(annotation_package, character.only = TRUE)
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

#' Set factors
#'
#' This function processes the sample sheet by converting the specified grouping variable and covariates (where specified) to factors.
#' The grouping variable and covariates are specified in the options list.
#'
#' @param sample_sheet A data frame containing the sample sheet data.
#' @param opts A list containing options for processing the sample sheet. The list should include:
#'   \itemize{
#'     \item \code{grouping_variable}: A string specifying the name of the grouping variable.
#'     \item \code{covariates}: A character vector specifying the names of the covariates.
#'     \item \code{are_covariates_factors}: A logical vector indicating whether each covariate should be treated as a factor.
#'   }
#' @return The processed sample sheet data frame with the specified grouping variable and covariates converted to factors.

process_sample_sheet <- function(sample_sheet, opts) {
    # Convert the grouping variable to a factor
    sample_sheet[[opts$grouping_variable]] <- as.factor(sample_sheet[[opts$grouping_variable]])
    
    # Convert the covariates that are factors to factors
    for (i in 1:length(opts$covariates)) {
        if (opts$are_covariates_factors[i]) {
            sample_sheet[[opts$covariates[i]]] <- as.factor(sample_sheet[[opts$covariates[i]]])
        }
    }
    return(sample_sheet)
}



#' Get Model Matrix
#'
#' This function creates a design matrix for the DM analysis based on the specified grouping variable and covariates.
#' 
#'
#' @param sample_sheet A data frame containing the sample sheet data.
#' @param opts A list containing options for creating the model matrix. The list should include:
#'   \itemize{
#'     \item \code{grouping_variable}: A string specifying the name of the grouping variable.
#'     \item \code{covariates}: A character vector specifying the names of the covariates.
#'   }
#' @return The design matrix created from the specified grouping variable and covariates.

get_model_matrix <- function(sample_sheet, opts) {
    # create the model formula
    model_string <- paste("~", opts$grouping_variable, sep = "")
    if (length(opts$covariates) > 0) {
        model_string <- paste(model_string, "+", paste(opts$covariates, collapse = " + "), sep = "")
    }
    # Create design matrix
    design <- model.matrix(as.formula(model_string), data = sample_sheet)
    return(design)
}

#' Get M-values
#'
#' This function preprocesses IDAT files based on the specified preprocessing method and extracts M-values.
#' The preprocessing method is specified in the options list.
#'
#' @param sample_sheet A data frame containing the sample sheet data, including the Basename column with paths to the IDAT files.
#' @param opts A list containing options for preprocessing the IDAT files. The list should include:
#'   \itemize{
#'     \item \code{preprocess_method}: A string specifying the preprocessing method to use. Valid options are:
#'       \itemize{
#'         \item \code{"preprocessIllumina"}
#'         \item \code{"preprocessSWAN"}
#'         \item \code{"preprocessQuantile"}
#'         \item \code{"preprocessNoob"}
#'         \item \code{"preprocessRaw"}
#'       }
#'   }
#' @return A matrix of M-values obtained after preprocessing the IDAT files.

get_m_values <- function(sample_sheet, opts) {
    # Read IDAT files
    RGset <- read.metharray(sample_sheet$Basename, extended = TRUE)
    preprocess_method = opts$preprocess_method
    if (preprocess_method == "preprocessIllumina") {
        Mset <- preprocessIllumina(RGset)
    } else if (preprocess_method == "preprocessSWAN") {
        Mset <- preprocessSWAN(RGset)
    } else if (preprocess_method == "preprocessQuantile") {
        Mset <- preprocessQuantile(RGset)
    } else if (preprocess_method == "preprocessNoob") {
        Mset <- preprocessNoob(RGset)
    } else if (preprocess_method == "preprocessRaw") {
        Mset <- preprocessRaw(RGset)
    } else {
        stop("Invalid preprocess method specified")
    }
    # get m_ values
    m_values <- getM(Mset)
    return(m_values)
}