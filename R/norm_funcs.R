#' Median Center Transformation
#'
#' Normalizes the data via median centering with median of the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param subset_fn character string indicating the subset function to use for normalization.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#' @param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details The sample-wise median of the feature subset specified for normalization is subtracted from each feature in e_data to get the normalized data. The location estimates are the sample-wise medians of the subset data. There are no scale estimates for median centering, though the function returns a NULL list element as a placeholder for a scale estimate. If backtransform is TRUE, the global median of the subset data (across all samples) is added back to the normalized values. Medians are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab NULL \cr
#' \tab \cr
#' location \tab numeric vector of length \code{n} medians for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab NULL \cr
#' \tab \cr
#' location \tab numeric value giving global median across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#'If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, pmartR_groupDF = attr(lipid_object2, "group_DF"))
#' lipids_median <- median_center(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#'
median_center <- function (e_data,
                           edata_id,
                           subset_fn,
                           feature_subset,
                           backtransform = FALSE,
                           apply_norm = FALSE,
                           check.names = TRUE) {
  
  # Determine which column of e_data contains the biomolecule IDs.
  id_col <- which(colnames(e_data) == edata_id)
  
  # Check if all biomolecules will be used.
  if (subset_fn == "all") {
    
    # Compute the column-wise median.
    location_param <- apply(e_data[, -id_col],
                            2, median, na.rm = TRUE)
    scale_param <- NULL
    
    # Runs if any subset function other than "all" is being used.
  } else {
    
    # Determine which rows will be used to calculate the normalization
    # parameters.
    id_rows <- as.character(e_data[, id_col]) %in% as.character(feature_subset)
    
    # Compute the column-wise median.
    location_param <- apply(e_data[id_rows, -id_col],
                            2, median, na.rm = TRUE)
    scale_param <- NULL
    
  }
  
  # Check the location_param vector for NAs.
  if (any(is.na(location_param))) {
    
    stop (paste("One or more of the location parameters used for median",
                "normalization is NA. Cannot proceed with the normalization.",
                sep = " "))
    
  }
  
  # if backtransforming, add back med
  if (backtransform == TRUE) {
    
    # Check which subset function was used and calculate the global median
    # accordingly.
    if (subset_fn == "all") {
      
      # Compute the global median.
      glob_median <- median(as.matrix(e_data[, -id_col]), na.rm = TRUE)
      
    } else {
      
      # Compute the global median (for the subsetted data). So it isn't quite a
      # global median :)
      glob_median <- median(as.matrix(e_data[id_rows, -id_col]), na.rm = TRUE)
      
    }
    
    # Check the glob_median vector for any NAs.
    if (any(is.na(glob_median))) {
      
      stop (paste("One or more of the location parameters used for",
                  "backtransform of the median normalization is NA. Cannot",
                  "proceed with the normalization. Try using backtransform",
                  "== FALSE.",
                  sep = " "))
      
    }
    
  } else {
    
    # Set the backtransform parameters to NULL because they won't be used.
    glob_median <- NULL
    
  }

  # Check the status of the apply_norm argument.
  if (apply_norm == FALSE) {
    
    # Create a list of the parameters used to normalize the data. In this list
    # the backtransform_params will always be NULL because these parameters can
    # only be calculated if the normalization is actually applied to the data.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = NULL,
                                               location = glob_median))
    
    # Runs when apply_norm is TRUE. 
  } else {
    
    # normalize the data!!!!!!!!!!!!!!!!!
    medi_data <- (e_data[, -id_col] - rep(location_param, each = nrow(e_data)))
    
    # if backtransforming, add back the global median
    if (backtransform == TRUE) {
      
      # Add the global median back to the normalized data. This is done because
      # many biologists suffer from a rare fear called negativabundaphobia. This
      # is the fear of seeing a negative abundance value in the normalized data.
      # By adding the median back to the normalized data all abundance values
      # will be positive. This at once relieves the biologist of their fear and
      # retains all of the beautiful features of the normalization process.
      medi_data <- medi_data + glob_median
      
    }
    
    # Replace the old data with the normalized data.
    e_data[, -id_col] <- medi_data
    
    # Plop all the information we spent so much time calculating in a list.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = NULL,
                                               location = glob_median),
                   transf_data = e_data)
    
  }
  
  return (ret_res)

}

#' Mean Center Transformation
#'
#' Normalizes the data via mean centering with mean of the feature subset specified for normalization
#'
#' @param e_data e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param subset_fn character string indicating the subset function to use for normalization.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#' @param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#' @details The sample-wise mean of the feature subset specified for normalization is subtracted from each feature in e_data to get the normalized data. The location estimates are the sample-wise means of the subset data.
#' There are no scale estimates for mean centering, though the function returns a NULL list element as a placeholdfer for a scale estimate.
#' If backtransform is TRUE, the global median of the subset data (across all samples) is added back to the normalized values. Medians are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab NULL \cr
#' \tab \cr
#' location \tab numeric vector of length \code{n} means for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab NULL \cr
#' \tab \cr
#' location \tab numeric value giving global median across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#' If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#'
#' @author Lisa Bramer, Kelly Stratton
#'
mean_center <- function (e_data,
                         edata_id,
                         subset_fn,
                         feature_subset,
                         backtransform = FALSE,
                         apply_norm = FALSE,
                         check.names = TRUE) {
  
  # Determine which column of e_data contains the biomolecule IDs.
  id_col <- which(colnames(e_data) == edata_id)
  
  # Check if all biomolecules will be used.
  if (subset_fn == "all") {
    
    # Compute the column-wise mean.
    location_param <- colMeans(e_data[, -id_col], na.rm = TRUE)
    scale_param <- NULL
    
    # Runs if any subset function other than "all" is being used.
  } else {
    
    # Determine which rows will be used to calculate the normalization
    # parameters.
    id_rows <- as.character(e_data[, id_col]) %in% as.character(feature_subset)
    
    # Compute the column-wise mean.
    location_param <- colMeans(e_data[id_rows, -id_col], na.rm = TRUE)
    scale_param <- NULL
    
  }
  
  # Check the location_param vector for NAs.
  if (any(is.na(location_param))) {
    
    stop (paste("One or more of the location parameters used for mean",
                "normalization is NA. Cannot proceed with the normalization.",
                sep = " "))
    
  }
  
  # if backtransforming, add back med
  if (backtransform == TRUE) {
    
    # Check which subset function was used and calculate the global median
    # accordingly.
    if (subset_fn == "all") {
      
      # Compute the global median.
      glob_median <- median(as.matrix(e_data[, -id_col]), na.rm = TRUE)
      
    } else {
      
      # Compute the global median (for the subsetted data). So it isn't quite a
      # global median :)
      glob_median <- median(as.matrix(e_data[id_rows, -id_col]), na.rm = TRUE)
      
    }
    
    # Check the glob_median vector for any NAs.
    if (any(is.na(glob_median))) {
      
      stop (paste("One or more of the location parameters used for",
                  "backtransform of the mean normalization is NA. Cannot",
                  "proceed with the normalization. Try using backtransform",
                  "== FALSE.",
                  sep = " "))
      
    }
    
  } else {
    
    # Set the backtransform parameters to NULL because they won't be used.
    glob_median <- NULL
    
  }
  
  # Check the status of the apply_norm argument.
  if (apply_norm == FALSE) {
    
    # Create a list of the parameters used to normalize the data. In this list
    # the backtransform_params will always be NULL because these parameters can
    # only be calculated if the normalization is actually applied to the data.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = NULL,
                                               location = glob_median))
    
    # Runs when apply_norm is TRUE. 
  } else {
    
    # normalize the data!!!!!!!!!!!!!!!!!
    mean_data <- (e_data[, -id_col] - rep(location_param, each = nrow(e_data)))
    
    # if backtransforming, mult by weighted-avg-mad and add back med
    if (backtransform == TRUE) {
      
      # Add the global median back to the normalized data. This is done because
      # many biologists suffer from a rare fear called negativabundaphobia. This
      # is the fear of seeing a negative abundance value in the normalized data.
      # By adding the median back to the normalized data all abundance values
      # will be positive. This at once relieves the biologist of their fear and
      # retains all of the beautiful features of the normalization process.
      mean_data <- mean_data + glob_median
      
    }
    
    # Replace the old data with the normalized data.
    e_data[, -id_col] <- mean_data
    
    # Plop all the information we spent so much time calculating in a list.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = NULL,
                                               location = glob_median),
                   transf_data = e_data)
    
  }
  
  return (ret_res)
  
}

#' Z-Score Transformation
#'
#' Normalizes the data via z-score transformation using the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param subset_fn character string indicating the subset function to use for normalization.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#' @param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details Each feature is scaled by subtracting the mean of the feature subset specified for normalization and then dividing the result by the standard deviation (SD) of the feature subset specified for normalization to get the normalized data.
#' The location estimates are the subset means for each sample. The scale estimates are the subset SDs for each sample. If backtransform is TRUE, the normalized feature values are multiplied by a pooled standard deviation (estimated across all samples) and a global mean of the subset data (across all samples) is added back to the normalized values. Means are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab numeric vector of length \code{n} standard deviations for each sample \cr
#' \tab \cr
#' location \tab numeric vector of length \code{n} means for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab numeric value giving the pooled standard deviation across all samples \cr
#' \tab \cr
#' location \tab numeric value giving global mean across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#' If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_edata)
#' peps_subset <- pep_edata$Mass_Tag_ID
#' peps_zscore <- zscore_transform(e_data = pep_edata, edata_id = "Mass_Tag_ID", feature_subset = peps_subset)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton, Bryan Stanfill
#'
zscore_transform <- function (e_data,
                              edata_id,
                              subset_fn,
                              feature_subset,
                              backtransform = FALSE,
                              apply_norm = FALSE,
                              check.names = TRUE) {
  
  # Determine which column of e_data contains the biomolecule IDs.
  id_col <- which(colnames(e_data) == edata_id)
  
  # Check if all biomolecules will be used.
  if (subset_fn == "all") {
    
    # Compute the column-wise mad and median.
    scale_param <- apply(e_data[, -id_col],
                         2, sd, na.rm = TRUE)
    location_param <- colMeans(e_data[, -id_col], na.rm = TRUE)
    
    # Runs if any subset function other than "all" is being used.
  } else {
    
    # Determine which rows will be used to calculate the normalization
    # parameters.
    id_rows <- as.character(e_data[, id_col]) %in% as.character(feature_subset)
    
    # Compute the column-wise mad and median.
    scale_param <- apply(e_data[id_rows, -id_col],
                         2, sd, na.rm = TRUE)
    location_param <- colMeans(e_data[id_rows, -id_col], na.rm = TRUE)
    
  }
  
  # Check the scale_param and location_param vectors for NAs.
  if (any(is.na(location_param)) || any(is.na(scale_param))) {
    
    stop (paste("One or more of the location or scale parameters used for",
                "z-score normalization is NA. Cannot proceed with the",
                "normalization.",
                sep = " "))
    
  }
  
  # if backtransforming, mult by weighted-avg-mad and add back mean
  if (backtransform == TRUE) {
    
    # Check which subset function was used and calculate sample_nonmiss,
    # pooled_sd, and glob_mean accordingly.
    if (subset_fn == "all") {
      
      # Calculate the variation by sample.
      sample_var <- apply(e_data[, -id_col], 2, var, na.rm = TRUE)
      
      # Determine the number of missing values by row.
      sample_nonmiss <- colSums(!is.na(e_data[, -id_col]))
      
      # Compute the global mean.
      glob_mean <- mean(as.matrix(e_data[, -id_col]), na.rm = TRUE)
      
    } else {
      
      # Calculate the variation by sample.
      sample_var <- apply(e_data[id_rows, -id_col], 2, var, na.rm = TRUE)
      
      # Enumerate the number of missing values by row.
      sample_nonmiss <- colSums(!is.na(e_data[id_rows, -id_col]))
      
      # Compute the global mean (for the subsetted data). So it isn't quite a
      # global mean :)
      glob_mean <- mean(as.matrix(e_data[id_rows, -id_col]), na.rm = TRUE)
      
    }
    
    # Calculate the pooled variance.
    pooled_var <- (sum((sample_nonmiss-1) * sample_var)) / 
      (sum(sample_nonmiss) - ncol(e_data[, -id_col]))
    
    # Calculate the pooled standard deviation.
    pooled_sd <- sqrt(pooled_var)
    
    # Check the glob_mean and pooled_sd vectors for any NAs.
    if (any(is.na(pooled_sd)) || any(is.na(glob_mean))) {
      
      stop (paste("One or more of the location or scale parameters used for",
                  "backtransform of the z-score normalization is NA. Cannot",
                  "proceed with the normalization. Try using backtransform",
                  "== FALSE.",
                  sep = " "))
      
    }
    
  } else {
    
    # Set the backtransform parameters to NULL because they won't be used.
    pooled_sd <- NULL
    glob_mean <- NULL
    
  }
  
  # Check the status of the apply_norm argument.
  if (apply_norm == FALSE) {
    
    # Create a list of the parameters used to normalize the data. In this list
    # the backtransform_params will always be NULL because these parameters can
    # only be calculated if the normalization is actually applied to the data.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = pooled_sd,
                                               location = glob_mean))
    
    # Runs when apply_norm is TRUE. 
  } else {
    
    # normalize the data!!!!!!!!!!!!!!!!!
    z_data <- (e_data[, -id_col] - rep(location_param, each = nrow(e_data))) /
      rep(scale_param, each = nrow(e_data))
    
    # if backtransforming, mult by weighted-avg-sd and add back mean
    if (backtransform == TRUE) {
      
      # Add the global mean and multiply the pooled standard deviation to the
      # normalized data. This is done because many biologists suffer from a rare
      # fear called negativabundaphobia. This is the fear of seeing a negative
      # abundance value in the normalized data. By adding the mean back to the
      # normalized data all abundance values will be positive. This at once
      # relieves the biologist of their fear and retains all of the beautiful
      # features of the normalization process.
      z_data <- z_data * pooled_sd + glob_mean
      
    }
    
    # Replace the old data with the normalized data.
    e_data[, -id_col] <- z_data
    
    # Plop all the information we spent so much time calculating in a list.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = pooled_sd,
                                               location = glob_mean),
                   transf_data = e_data)
    
  }
  
  return (ret_res)
  
}

#' Median Absolute Deviation Transformation
#'
#' Calculate normalization parameters for the data via median absolute deviation (MAD) transformation based on the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param subset_fn character string indicating the subset function to use for normalization.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the parameters for backtransforming the data after normalization will be calculated so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#' @param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details Each feature is scaled by subtracting the median of the feature subset specified for normalization and then dividing the result by the median absolute deviation (MAD) of the feature subset specified for normalization to get the normalized data. The location estimates are the subset medians for each sample. The scale estimates are the subset MADs for each sample. Medians are taken ignoring any NA values. If backtransform is TRUE, the normalized feature values are multiplied by a pooled MAD (estimated from all samples) and a global median of the subset data (across all samples) is added back to the normalized values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab numeric vector of length \code{n} median absolute deviations (MAD) for each sample \cr
#' \tab \cr
#' location \tab numeric vector of length \code{n} medians for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab numeric value giving pooled MAD \cr
#' \tab \cr
#' location \tab numeric value giving global median across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#' If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, pmartR_groupDF = attr(lipid_object2, "group_DF"))
#' lipids_mad <- mad_transform(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
#'}
#' @author Lisa Bramer, Kelly Stratton
#'
mad_transform <- function (e_data,
                           edata_id,
                           subset_fn,
                           feature_subset,
                           backtransform = FALSE,
                           apply_norm = FALSE,
                           check.names = TRUE) {
  
  # Determine which column of e_data contains the biomolecule IDs.
  id_col <- which(colnames(e_data) == edata_id)
  
  # Check if all biomolecules will be used.
  if (subset_fn == "all") {
    
    # Compute the column-wise mad and median.
    scale_param <- apply(e_data[, -id_col],
                         2, mad, na.rm = TRUE)
    location_param <- apply(e_data[, -id_col],
                            2, median, na.rm = TRUE)
    
    # Runs if any subset function other than "all" is being used.
  } else {
    
    # Determine which rows will be used to calculate the normalization
    # parameters.
    id_rows <- as.character(e_data[, id_col]) %in% as.character(feature_subset)
    
    # Compute the column-wise mad and median.
    scale_param <- apply(e_data[id_rows, -id_col],
                         2, mad, na.rm = TRUE)
    location_param <- apply(e_data[id_rows, -id_col],
                            2, median, na.rm = TRUE)
    
  }
  
  # Check the scale_param and location_param vectors for NAs.
  if (any(is.na(location_param)) || any(is.na(scale_param))) {
    
    stop (paste("One or more of the location or scale parameters used for MAD",
                "normalization is NA. Cannot proceed with the normalization.",
                sep = " "))
    
  }
  
  # if backtransforming, mult by weighted-avg-mad and add back med
  if (backtransform == TRUE) {
    
    # Check which subset function was used and calculate sample_nonmiss and
    # glob_median accordingly.
    if (subset_fn == "all") {
      
      # Determine the number of missing values by row.
      sample_nonmiss <- colSums(!is.na(e_data[, -id_col]))
      
      # Compute the global median.
      glob_median <- median(as.matrix(e_data[, -id_col]), na.rm = TRUE)
      
    } else {
      
      # Enumerate the number of missing values by row.
      sample_nonmiss <- colSums(!is.na(e_data[id_rows, -id_col]))
      
      # Compute the global median (for the subsetted data). So it isn't quite a
      # global median :)
      glob_median <- median(as.matrix(e_data[id_rows, -id_col]), na.rm = TRUE)
      
    }
    
    # Calculate the pooled mad.
    pooled_mad <- (sum(sample_nonmiss * scale_param)) / sum(sample_nonmiss)
    
    # Check the glob_median and pooled_mad vectors for any NAs.
    if (any(is.na(pooled_mad)) || any(is.na(glob_median))) {
      
      stop (paste("One or more of the location or scale parameters used for",
                  "backtransform of the MAD normalization is NA. Cannot",
                  "proceed with the normalization. Try using backtransform",
                  "== FALSE.",
                  sep = " "))
      
    }
    
  } else {
    
    # Set the backtransform parameters to NULL because they won't be used.
    pooled_mad <- NULL
    glob_median <- NULL
    
  }
  
  # Check the status of the apply_norm argument.
  if (apply_norm == FALSE) {
    
    # Create a list of the parameters used to normalize the data. In this list
    # the backtransform_params will always be NULL because these parameters can
    # only be calculated if the normalization is actually applied to the data.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = pooled_mad,
                                               location = glob_median))
    
    # Runs when apply_norm is TRUE. 
  } else {
    
    # normalize the data!!!!!!!!!!!!!!!!!
    mad_data <- (e_data[, -id_col] - rep(location_param, each = nrow(e_data))) /
      rep(scale_param, each = nrow(e_data))
    
    # if backtransforming, mult by weighted-avg-mad and add back med
    if (backtransform == TRUE) {
      
      # Add the global median and multiply the median absolute deviation to the
      # normalized data. This is done because many biologists suffer from a rare
      # fear called negativabundaphobia. This is the fear of seeing a negative
      # abundance value in the normalized data. By adding the median back to the
      # normalized data all abundance values will be positive. This at once
      # relieves the biologist of their fear and retains all of the beautiful
      # features of the normalization process.
      mad_data <- mad_data * pooled_mad + glob_median
      
    }
    
    # Replace the old data with the normalized data.
    e_data[, -id_col] <- mad_data
    
    # Plop all the information we spent so much time calculating in a list.
    ret_res = list(norm_params = list(scale = scale_param,
                                      location = location_param),
                   backtransform_params = list(scale = pooled_mad,
                                               location = glob_median),
                   transf_data = e_data)
    
  }
  
  return (ret_res)
  
}
