#' Median Absolute Deviation Transformation
#'
#' Calculate normalization parameters for the data via median absolute deviation (MAD) transformation based on the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#'@param backtransform logical argument. If TRUE, the parameters for backtransforming the data after normalization will be calculated so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#'@param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details Each feature is scaled by subtracting the median of the feature subset specified for normalization and then dividing the result by the median absolute deviation (MAD) of the feature subset specified for normalization to get the normalized data. The location estimates are the subset medians for each sample. The scale estimates are the subset MADs for each sample. Medians are taken ignoring any NA values. If backtransform is TRUE, the normalized feature values are multiplied by a pooled MAD (estimated from all samples) and a global median of the subset data (across all samples) is added back to the normalized values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab numeric vector of length \code{n} median absolute deviations (MAD) for each sample
#' \tab \cr
#' location \tab numeric vector of length \code{n} medians for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab numeric value giving pooled MAD
#' \tab \cr
#' location \tab numeric value giving global median across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#' If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, mintR_groupDF = attr(lipid_object2, "group_DF"))
#' lipids_mad <- mad_transform(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
#'}
#' @author Lisa Bramer, Kelly Stratton
#'

mad_transform <- function(e_data, edata_id, feature_subset, backtransform=FALSE, apply_norm = FALSE, check.names=TRUE){

  subset.data <- e_data[as.character(e_data[,edata_id]) %in% as.character(feature_subset),]
  subset.data <- subset.data[, -which(colnames(subset.data)==edata_id)]
  subset.data <- as.matrix(subset.data)

  # get global mad and global median
  location_param <- apply(subset.data, 2, median, na.rm=TRUE)
  scale_param <- apply(subset.data, 2, mad, na.rm=TRUE)


  # if backtransforming, mult by weighted-avg-mad and add back med
  if(backtransform == TRUE){
    # calc a weighted MAD
    sample.nonmiss <- apply(subset.data, 2, function(x) sum(!is.na(x)))
    pooled.mad <- (sum(sample.nonmiss * scale_param))/sum(sample.nonmiss)
    glob.median <- median(subset.data, na.rm=TRUE)
  }else{
    pooled.mad <- NULL
    glob.median <- NULL
  }

  if(apply_norm == FALSE){
  ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = pooled.mad, location = glob.median))
  }else{
    # normalize
    mad.data <- (e_data[,-which(colnames(e_data)==edata_id)] - matrix(rep(location_param, each=nrow(e_data)), nrow=nrow(e_data))) / matrix(rep(scale_param, each=nrow(e_data)), nrow=nrow(e_data))

    # if backtransforming, mult by weighted-avg-mad and add back med
    if(backtransform == TRUE){

      mad.data <- mad.data * pooled.mad + array(median(subset.data, na.rm=TRUE), dim=dim(e_data[,-1]))
    }

    norm_data <- data.frame(e_data[,edata_id], mad.data, check.names=check.names)
    names(norm_data)[1] <- edata_id

    ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = pooled.mad, location = glob.median), transf_data = norm_data)

    if(any(is.na(location_param)) || any(is.na(scale_param))) stop("One or more of the location or scale parameters to use for MAD normalization is NA. Cannot proceed.")

    if(backtransform==TRUE){
      if(any(is.na(pooled.mad)) || any(is.na(glob.median))) stop("One or more of the location or scale parameters to use for backtransform of the MAD normalization is NA. Cannot proceed.")
    }

  }



  return(ret_res)
}

#' Median Center Transformation
#'
#' Normalizes the data via median centering with median of the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#'@param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details The sample-wise median of the feature subset specified for normalization is subtracted from each feature in e_data to get the normalized data. The location estimates are the sample-wise medians of the subset data. There are no scale estimates for median centering, though the function returns a NULL list element as a placeholder for a scale estimate. If backtransform is TRUE, the global median of the subset data (across all samples) is added back to the normalized values. Medians are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab NULL
#' \tab \cr
#' location \tab numeric vector of length \code{n} medians for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab NULL
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
#' lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, mintR_groupDF = attr(lipid_object2, "group_DF"))
#' lipids_median <- median_center(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#'

median_center <- function(e_data, edata_id, feature_subset, backtransform=FALSE, apply_norm = FALSE, check.names=TRUE){

  subset.data <- e_data[as.character(e_data[,edata_id]) %in% as.character(feature_subset),]
  subset.data <- subset.data[, -which(colnames(subset.data)==edata_id)]
  subset.data <- as.matrix(subset.data)

  # get medians per sample to use for normalization
  location_param <- apply(subset.data, 2, median, na.rm=TRUE)
  scale_param <- NULL


  # if backtransforming, subtract subset median from each sample (updated  by KS 9/29/15)
  if(backtransform==TRUE){
    glob.med <- median(subset.data, na.rm=TRUE)
  }else{
    glob.med <- NULL
  }

  if(apply_norm == FALSE){
  ret_res = list(norm_params=list(scale = scale_param, location = location_param), backtransform_params = list(scale = NULL, location = glob.med))
  }else{
    # normalize
    med.diff <- e_data[,-which(colnames(e_data)==edata_id)] - matrix(rep(location_param, each=nrow(e_data)), nrow=nrow(e_data))

    # if backtransforming, subtract subset median from each sample
    if(backtransform==TRUE){
      med.diff <- med.diff + array(median(subset.data, na.rm=TRUE), dim=dim(e_data[,-1]))
    }

    norm_data <- data.frame(e_data[,edata_id], med.diff, check.names=check.names)
    names(norm_data)[1] <- edata_id
    ret_res = list(norm_params=list(scale = scale_param, location = location_param), backtransform_params = list(scale = NULL, location = glob.med), transf_data = norm_data)

    if(any(is.na(location_param))) stop("One or more of the location parameters to use for median normalization is NA. Cannot proceed.")

    if(backtransform==TRUE){
      if(any(is.na(glob.med))) stop("One or more of the location parameters to use for backtransform of the median normalization is NA. Cannot proceed.")
    }

  }
  return(ret_res)
}

#' Mean Center Transformation
#'
#' Normalizes the data via mean centering with mean of the feature subset specified for normalization
#'
#' @param e_data e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#'@param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#' @details The sample-wise mean of the feature subset specified for normalization is subtracted from each feature in e_data to get the normalized data. The location estimates are the sample-wise means of the subset data.
#' There are no scale estimates for mean centering, though the function returns a NULL list element as a placeholdfer for a scale estimate.
#' If backtransform is TRUE, the global median of the subset data (across all samples) is added back to the normalized values. Medians are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab NULL
#' \tab \cr
#' location \tab numeric vector of length \code{n} means for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab NULL
#' \tab \cr
#' location \tab numeric value giving global median across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#'If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#'
#' @author Lisa Bramer, Kelly Stratton
#'

mean_center <- function(e_data, edata_id, feature_subset, backtransform=FALSE, apply_norm = FALSE, check.names=TRUE){

  subset_data <- e_data[as.character(e_data[,edata_id]) %in% as.character(feature_subset),]
  subset_data <- subset_data[, -which(colnames(subset_data)==edata_id)]
  subset_data <- as.matrix(subset_data)

  # get medians per sample to use for normalization
  location_param <- apply(subset_data, 2, mean, na.rm=TRUE)
  scale_param <- NULL



  # if backtransforming, subtract subset median from each sample (updated  by KS 9/29/15)
  if(backtransform==TRUE){
    glob.median <- median(subset_data, na.rm=TRUE)
  }else{
    glob.median <- NULL
  }

if(apply_norm == FALSE){
  ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = NULL, location = glob.median))
}else{
  # normalize
  mean_diff <- e_data[,-which(colnames(e_data)==edata_id)] - matrix(rep(location_param, each=nrow(e_data)), nrow=nrow(e_data))

  # if backtransforming, subtract subset median from each sample #
  if(backtransform==TRUE){
    mean_diff <- mean_diff + array(median(subset_data, na.rm=TRUE), dim=dim(e_data[,-1]))
  }

  norm_data <- data.frame(e_data[,edata_id], mean_diff, check.names=check.names)
  names(norm_data)[1] <- edata_id

  ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = NULL, location = glob.median), transf_data = norm_data)

  if(any(is.na(location_param))) stop("One or more of the location parameters to use for mean normalization is NA. Cannot proceed.")

  if(backtransform==TRUE){
    if(any(is.na(glob.median))) stop("One or more of the location parameters to use for backtransform of the mean normalization is NA. Cannot proceed.")
  }

}

return(ret_res)
}

#' Z-Score Transformation
#'
#' Normalizes the data via z-score transformation using the feature subset specified for normalization
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param feature_subset character vector containing the feature names in the subset to be used for normalization
#' @param backtransform logical argument. If TRUE, the data will be back transformed after normalization so that the values are on a scale similar to their raw values. See details for more information. Defaults to FALSE.
#'@param apply_norm logical argument. If TRUE, the normalization will be applied to the data. Defaults to FALSE.
#'
#' @details Each feature is scaled by subtracting the mean of the feature subset specified for normalization and then dividing the result by the standard deviation (SD) of the feature subset specified for normalization to get the normalized data.
#' The location estimates are the subset means for each sample. The scale estimates are the subset SDs for each sample. If backtransform is TRUE, the normalized feature values are multiplied by a pooled standard deviation (estimated across all samples) and a global mean of the subset data (across all samples) is added back to the normalized values. Means are taken ignoring any NA values.
#'
#' @return List containing two elements: \code{norm_params} is list with two elements:
#' \tabular{ll}{
#' scale \tab numeric vector of length \code{n} standard deviations for each sample
#' \tab \cr
#' location \tab numeric vector of length \code{n} means for each sample
#' \cr
#' }
#'
#' \code{backtransform_params} is a list with two elements:
#' \tabular{ll}{
#' scale \tab numeric value giving the pooled standard deviation across all samples
#' \tab \cr
#' location \tab numeric value giving global mean across all samples
#' \cr
#' }
#' If \code{backtransform} is set to TRUE then each list item under \code{backtransform_params} will be NULL.
#'
#'If \code{apply_norm} is TRUE, the transformed data is returned as a third list item.
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_edata)
#' peps_subset <- pep_edata$Mass_Tag_ID
#' peps_zscore <- zscore_transform(e_data = pep_edata, edata_id = "Mass_Tag_ID", feature_subset = peps_subset)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton, Bryan Stanfill
#'
#'

zscore_transform <- function(e_data, edata_id, feature_subset, backtransform=FALSE, apply_norm = FALSE, check.names=TRUE){

  subset_data <- e_data[as.character(e_data[,edata_id]) %in% as.character(feature_subset),]
  subset_data <- subset_data[, -which(colnames(subset_data)==edata_id)]
  subset_data <- as.matrix(subset_data)

  # get global mad and global median
  location_param <- colMeans(subset_data, na.rm=TRUE)
  scale_param <- apply(subset_data, 2, sd, na.rm=TRUE)


  if(backtransform==TRUE){
    # calc a pooled variance
    sample_var <- apply(subset_data, 2, var, na.rm=TRUE)
    sample_nonmiss <- apply(subset_data, 2, function(x) sum(!is.na(x)))
    pooled_var <- (sum((sample_nonmiss-1) * sample_var))/(sum(sample_nonmiss) - ncol(subset_data))
    pooled_sd <- sqrt(pooled_var)
    glob_mean <- mean(subset_data, na.rm=TRUE)
  }else{
    pooled_sd <- NULL
    glob_mean <- NULL
  }

  if(apply_norm == FALSE){
  ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = pooled_sd, location = glob_mean))
    }else{
      # normalize
      zscore_transf <- (e_data[,-which(colnames(e_data)==edata_id)] - matrix(rep(location_param, each=nrow(e_data)), nrow=nrow(e_data))) / matrix(rep(scale_param, each=nrow(e_data)), nrow=nrow(e_data))

      if(backtransform==TRUE){
        # calc a pooled variance
      zscore_transf <- zscore_transf*pooled_sd + array(mean(subset_data, na.rm=TRUE), dim=dim(e_data[,-1]))
      }
      norm_data <- data.frame(e_data[,edata_id], zscore_transf, check.names=check.names)
      names(norm_data)[1] <- edata_id
      ret_res = list(norm_params = list(scale = scale_param, location = location_param), backtransform_params = list(scale = pooled_sd, location = glob_mean), transf_data = norm_data)

      if(any(is.na(location_param)) || any(is.na(scale_param))) stop("One or more of the location or scale parameters to use for z-score normalization is NA. Cannot proceed.")

      if(backtransform==TRUE){
        if(any(is.na(pooled_sd)) || any(is.na(glob_mean))) stop("One or more of the location or scale parameters to use for backtransform of the z-score normalization is NA. Cannot proceed.")
      }

    }
  return(ret_res)
}
