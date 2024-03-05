#' Scale from zero to one
#'
#' Perform scaling of data from zero to one.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#'
#' @details The sample-wise minimum of the features is subtracted from each 
#' feature in e_data, then divided by the difference between the sample-wise 
#' minimum and maximum of the features to get the normalized data. The location 
#' estimates are not applicable for this data and the function returns a NULL 
#' list element as a placeholder. The scale estimates are the sample-wise 
#' feature ranges. All NA values are replaced with zero.
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#'
#' mymetab <- edata_transform(
#'   omicsData = metab_object,
#'   data_scale = "log2"
#' )
#' mymetab <- group_designation(
#'   omicsData = mymetab,
#'   main_effects = "Phenotype"
#' )
#' norm_data <- normalize_zero_one_scaling(
#'   omicsData = mymetab
#' )
#'
#' @author Rachel Richardson
#'
#'
#' @export
#'
normalize_zero_one_scaling <- function(omicsData) {
  
  ## initial checks ##
  
  # Ensure omicsData is an appropriate class
  if (!inherits(omicsData, c(
    "pepData", "proData", "lipidData", "metabData",
    "nmrData"
  ))) {
    # Throw down an error for blatant abuse of the pmartR standards.
    stop(paste("omicsData must be of class 'pepData', 'proData', 'lipidData',",
               "'metabData' or 'nmrData'",
               sep = " "
    ))
  }
  
  # data should be log transformed #
  if (!(get_data_scale(omicsData) %in% c("log2", "log10", "log"))) {
    # Tell the user to read the documentation before running functions willy
    # nilly.
    stop(paste("omicsData$e_data should be log transformed prior to calling",
               "normalize_data. See documentation for the edata_transform",
               "function for more information.",
               sep = " "
    ))
  }
  
  edata_id <- get_edata_cname(omicsData)
  
  # Subset the data prior to normalization -------------------------------------
  
  # subset data using current subset method #
  peps <- all_subset(omicsData$e_data, edata_id)
  
  # Normalize the data ---------------------------------------------------------
  
  # Pull out the name of the normalization function from the input.
  fn_to_use <- zero_one_scale
  
  # Normalize the data according to the method selected.
  norm_results <- fn_to_use(
    e_data = omicsData$e_data,
    edata_id = edata_id
  )
  
  # Check if the normalization will be applied. If return all of the information
  # needed to normalize (which biomolecules to normalize with, normalization
  # method and parameters, ...).
  # Update e_data with the normalized data.
  omicsData$e_data = norm_results$transf_data
  
  # Update the norm_info attribute.
  attributes(omicsData)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "zero_to_one", # added 12/21/17 by KS
    subset_fn = NULL,
    subset_params = NULL,
    norm_fn = "zero_to_one",
    n_features_calc = length(peps),
    prop_features_calc = NULL,
    params = list(
      norm_scale = norm_results$norm_params$scale,
      norm_location = norm_results$norm_params$location,
      bt_scale = norm_results$backtransform_params$scale,
      bt_location = norm_results$backtransform_params$location
    )
  )
  
  # Return the normalized omicsData object with all of its updated attributes!
  return(omicsData)
}


#' Zero to One scaling
#'
#' Re-scales the data to be between 0 and 1
#'
#' @param e_data e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the
#'   number of peptides, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details The sample-wise minimum of the features is subtracted from each 
#' feature in e_data, then divided by the difference between the sample-wise 
#' minimum and maximum of the features to get the normalized data. The location 
#' estimates are not applicable for this data and the function returns a NULL 
#' list element as a placeholder. The scale estimates are the sample-wise 
#' feature ranges. All NA values are replcaed with zero.
#'
#' @return List containing two elements: \code{norm_params} is list with two
#'   elements:
#' \tabular{ll}{
#' scale \tab Range of each sample used in scaling \cr
#' \tab \cr
#' location \tab NULL
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab NULL \cr
#' \tab \cr
#' location \tab NULL
#' \cr
#' }
#'
#' The transformed data is returned as a third
#' list item.
#'
#' @author Lisa Bramer, Kelly Stratton, Rachel Richardson
#'
zero_one_scale <- function(
    e_data,
    edata_id) {
  
  # Determine which column of e_data contains the biomolecule IDs.
  id_col <- which(colnames(e_data) == edata_id)
  
  # Compute the column-wise ranges
  scale_param <- apply(e_data[, -id_col], 2, range, na.rm = TRUE)

  
  # Check the location_param vector for NAs.
  if (any(is.na(unlist(scale_param)))) {
    stop(paste("One or more of the scale parameters (parameter: range) used for zero-to-one scaling ",
               "normalization is NA. Cannot proceed with the normalization.",
               sep = " "
    ))
  }
  

  # normalize the data!!!!!!!!!!!!!!!!!
  range01 <- function(x, scale){
    (x-min(scale))/(max(scale)-min(scale))
  }
  
  e_data_numeric <- e_data[, -id_col]
  
  suppressMessages(
    zeroone_data <- purrr::map_dfc(
      1:ncol(e_data_numeric), 
      function(col) range01(e_data_numeric[,col], scale_param[,col]))
  )
  
  zeroone_data[is.na(zeroone_data)] <- 0
  
  # Replace the old data with the normalized data.
  e_data[, -id_col] <- zeroone_data
  
  # Plop all the information we spent so much time calculating in a list.
  ret_res = list(
    norm_params = list(
      scale = scale_param,
      location = NULL
    ),
    backtransform_params = list(
      scale = NULL,
      location = NULL
    ),
    transf_data = e_data
  )
  
  return(ret_res)
  }
