#'RMD Runs
#'
#'The method computes a robust Mahalanobis distance that can be mapped to a
#'p-value and used to identify outlying samples
#'
#'@param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'  'lipidData', or 'nmrData' usually created by \code{\link{as.pepData}},
#'  \code{\link{as.proData}}, \code{\link{as.metabData}},
#'  \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#'@param ignore_singleton_groups logical indicator of whether to remove
#'  singleton groups or not; defaults to TRUE. If TRUE, rmd_filter results are
#'  returned only for samples in groups of size greater than 1. This is used
#'  when calculating the correlation.
#'  
#'@param metrics A character vector indicating which metrics should be used when
#'       calculating the robust Mahalanobis distance. This vector must contain
#'       between two and five of the following options: "MAD" (Median Absolute
#'       Deviation), "Kurtosis", "Skewness", "Correlation", and
#'       "Proportion_Missing". The default is NULL. When NULL a combination of
#'       metrics will be chosen depending on the class of omicsData.
#'
#'@return a data.frame containing columns for the Sample ID, log2 robust
#'  Mahalanobis distance, p-values, and robust Mahalanobis distance
#'
#'@details The metrics on which the log2 robust Mahalanobis distance is based
#'  can be specified using the \code{metrics} argument. \tabular{ll}{ pepData
#'  \tab For pepData objects, all five of the metrics "MAD", "Kurtosis",
#'  "Skewness", "Correlation", "Proportion_Missing" may be used (this is the
#'  default). \cr proData \tab For proData objects, all five of the metrics
#'  "MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing" may be
#'  used (this is the default). \cr metabData \tab For metabData objects, the
#'  use of "Proportion_Missing" is discouraged due to the general lack of
#'  missing data in metabolomics datasets (the default behavior omits
#'  "Proportion_Missing" from the metrics). \cr lipidData \tab For lipidData
#'  objects, , the use of "Proportion_Missing" is discouraged due to the general
#'  lack of missing data in metabolomics datasets (the default behavior omits
#'  "Proportion_Missing" from the metrics). \cr }
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' metab_object2 <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' metab_object3 <- group_designation(omicsData = metab_object2, main_effects = "Condition")
#' rmd_results <- rmd_filter(omicsData = metab_object3, metrics=c("MAD", "Skewness", "Correlation"))
#' rmd_results <- rmd_filter(omicsData = metab_object2)
#'
#' data(pep_pepData)
#' pep_pepData2 <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' pep_pepData3 <- group_designation(omicsData = pep_pepData2, main_effects = "Condition")
#' rmd_results <- rmd_filter(omicsData = pep_pepData3)
#'}
#'
#'@references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R.,
#'  Pounds, J., and Webb-Robertson, B.J. (2011), \emph{Improved quality control
#'  processing of peptide-centric LC-MS proteomics data}. Bioinformatics.
#'  27(20): 2866-2872.
#'
#'@author Lisa Bramer, Kelly Stratton
#'
#'@export
#'
#'@rdname rmd_filter
#'@name rmd_filter
#'  
rmd_filter <- function(omicsData,
                       ignore_singleton_groups = TRUE,
                       metrics = NULL) {
  
  # Create functions to be used later ------------------------------------------
  
  # This function will be used in two different places later depending on the
  # metrics used in the input. It is here so changes only need to be made in one
  # location of the code.
  mal.fun <- function(x) (t(as.vector(x) - med.mat) %*%
                            solve(as.matrix(cov.mat)) %*%
                            (as.vector(x) - med.mat))
  
  # Run some preliminary checks ------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    # Throw down an error like Zeus with his lightning bolts.
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
    
  }
  
  # check that ignore_singleton_groups is logical #
  if (!is.logical(ignore_singleton_groups)) {
    
    # Stop the illogical user with an error.
    stop ("ignore_singleton_groups must be either TRUE or FALSE")
    
  }
  
  # group_DF attribute is required #
  if (is.null(attr(omicsData, "group_DF"))) {
    
    # Kindly tell the user they are wrong and point them in the right direction
    # to get some much needed help.
    stop (paste("omicsData must contain attribute information for 'group_DF'.",
                "See documentation for group_designation function for more",
                "information.",
                sep = " "))
    
  }
  
  # data should be log transformed #
  if (!attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")) {
    
    # Yet another error message for the user. Poor thing.
    stop (paste("omicsData$e_data should be log transformed prior to calling",
                "rmd_filter. See documentation for edata_transform function",
                "for more information.",
                sep = " "))
    
  }
  
  # Check the number of observations in e_data.
  if (attributes(omicsData)$data_info$num_edata < 50) {
    
    # Warn the user that skimping on the sample size is ALWAYS a bad idea!
    warning (paste("Use the results of the RMD filter with caution due to a",
                   "small number of biomolecules (<50).",
                   sep = " "))
    
  }
  
  # Check if metrics is NULL.
  if (is.null(metrics)) {
    
    # Take a looksie at the omicsData class.
    if (inherits(omicsData, c("pepData", "proData"))) {
      
      # Set default metrics values for peptide and protein data.
      metrics <- c("MAD", "Kurtosis", "Skewness", "Correlation",
                   "Proportion_Missing")
      
    } else if (inherits(omicsData, c("lipidData", "metabData", "nmrData"))) {
      
      # Set default metrics values for lipid and metabolite data.
      metrics <- c("MAD", "Kurtosis", "Skewness", "Correlation")
      
    }
    
    # Check they didn't mess up the metrics vector (if it isn't NULL).
  } else {
    
    # Make sure the entries in metrics are an acceptable type.
    if (!all(metrics %in% c("MAD", "mad", "m", "median absolute deviation",
                            "median_absolute_deviation",
                            "Kurtosis", "kurtosis", "kurt", "k",
                            "Skewness", "skew", "s", "skewness",
                            "Correlation", "corr", "c", "cor", "correlation",
                            "Proportion_Missing", "proportion_missing", "p",
                            "prop_missing", "prop_miss", "proportion_miss",
                            "prop", "proportion", "proportion missing",
                            "prop miss", "proportion miss"))) {
      
      # Stop the user for using unholy elements in the metrics vector.
      stop (paste("One or more of the elements in metrics does not match",
                  "the acceptable values. See details for metrics in the",
                  "rmd_filter documenation.",
                  sep = " "))
      
    }
    
    # Make sure there are not any repeated metrics.
    if (length(unique(metrics)) != length(metrics)) {
      
      # Stop the user for trying to use a metric more than once.
      stop ("One or more of the elements in metrics occurs more than once.")
      
    }
    
    # Make sure the metrics vector as an appropriate number of elements
    if (length(metrics) < 2) {
      
      # Just being Zeus again with the error throwing.
      stop ("The metrics vector must contain at least two elements.")
      
    } else if (length(metrics) > 5) {
      
      # Zeus is nothing compared to pmartR!! BWAHAHAHA!!!
      stop ("The metrics vector cannot contain more than five elements.")
      
    }
    
  }
  
  # Carry out preliminary group calculations -----------------------------------
  
  # Create an object for the group_DF attribute with a name that is almost as
  # long as using the call to the attr function.
  mintR_groupDF = attr(omicsData, "group_DF")
  
  ### Aug 2020: deal with any groups that have only a single sample in them ###
  if (any(table(mintR_groupDF$Group) == 1) && ignore_singleton_groups) {
    
    # which groups have a single sample in them #
    singleton_groups <- names(which(table(mintR_groupDF$Group) == 1))
    singleton_rows <- which(mintR_groupDF$Group %in% singleton_groups)
    
    # filter out samples corresponding to groups of size 1 #
    to_remove <- as.character(mintR_groupDF[singleton_rows,
                                            get_fdata_cname(omicsData)])
    myfilter <- custom_filter(omicsData, f_data_remove = to_remove)
    omicsData <- applyFilt(myfilter, omicsData)
    mintR_groupDF <- attr(omicsData, "group_DF")
    
  }
  
  # Compare the number of samples to the number of metrics used.
  if (dim(omicsData$f_data)[1] < 2*length(metrics)) {
    
    # Warn the user that little sample sizes are the bane of their existence.
    warning (paste("Use the results of the RMD filter with caution due to a",
                   "small number of samples relative to the number of metrics",
                   "being used. Consider reducing the number of metrics being",
                   "used.",
                   sep = " "))
  }
  
  # Compute the metrics used to calculate the RMD ------------------------------
  
  # Extricate the column number from e_data that contains the IDs.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Fish out the column number from f_data that contains sample IDs.
  id_col_f <- which(names(omicsData$f_data) == get_fdata_cname(omicsData))
  
  # Convert metrics to lower case for matching purposes.
  metrics <- tolower(metrics)
  
  # Initialize a vector for counting the number of metrics used.
  metrics_final <- rep(NA, 5)
  
  # Initialize a data frame with the sample ID column from f_data.
  rmd.vals <- data.frame(Sample.ID = omicsData$f_data[, id_col_f])
  
  # Compute the median absolute deviation across the samples (columns).
  if (any(metrics %in% c("mad", "m", "median absolute deviation",
                         "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation",
                               "median_absolute_deviation"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'MAD'.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(omicsData$e_data[, -id_col])[,2]
  }
  
  # Compute the kurtosis across the samples (columns).
  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Kurtosis'.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(omicsData$e_data[, -id_col])[,2]
  }
  
  # Compute the skewness across the samples (columns).
  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Skewness'.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(omicsData$e_data[, -id_col])[,2]
  }
  
  # Compute the correlation across the samples (columns).
  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Correlation'.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = suppressWarnings(run_group_meancor(omicsData,
                                                       mintR_groupDF,
                                                       ignore_singleton_groups)[,2])
  }
  
  # Compute the proportion of missing values across the samples (columns).
  if (any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss",
                         "proportion_miss", "prop", "proportion",
                         "proportion missing", "prop miss", "proportion miss"))) {
    
    # check to see whether there is any missing data
    if (attr(omicsData, 'data_info')$num_miss_obs == 0) {
      
      # if no missing data, we flat out refuse to include prop_missing in the
      # metrics, so give a warning to this effect.
      warning (paste("There are no missing values in e_data, therefore",
                     "Proportion_Missing will not be used as one of the",
                     "metrics for RMD-Runs.",
                     sep = " "))
      
      # If there are missing values compute the proportion of missing values for
      # each sample (column) as usual.
    } else {
      
      ind = which(metrics %in% c("proportion_missing", "p", "prop_missing",
                                 "prop_miss", "proportion_miss", "prop",
                                 "proportion", "proportion missing", "prop miss",
                                 "proportion miss"))
      
      if(length(ind)>1){
        stop("More than one of the entries in metrics matches 'Proportion_Missing'.")}
      
      metrics = metrics[-ind]
      metrics_final[5] = 1
      rmd.vals$Proportion_Missing = run_prop_missing(omicsData$e_data[, -id_col])[,2]
      
      ## proceed, to check the rank of cov.mat ##
      
      ## Conduct Robust PCA ##
      robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1],
                                     k = (ncol(rmd.vals)-1),
                                     mcd = FALSE,
                                     scale = FALSE)
      
      # Check for errors when conducting robust PCA.
      if (inherits(robpca.res, "try-error")) {
        
        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "prop_missing as one of the metrics. Try again, excluding",
                    "Proportion_Missing from the metrics vector.",
                    sep = " "))
      }
      
      ## Calculate Covariance Matrix #
      cov.mat = (robpca.res@loadings %*%
                   diag(robpca.res@eigenvalues) %*%
                   t(robpca.res@loadings))
      
      # Extract the rank of the covariance matrix.
      myrank = qr(cov.mat)$rank
      
      # Check the rank of the covariance matrix compared to the number of
      # metrics used to calculate it.
      if (myrank != sum(metrics_final, na.rm = TRUE)) {
        
        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "prop_missing as one of the metrics. Try again, excluding",
                    "Proportion_Missing from the metrics vector.",
                    sep = " "))
        
      }
      
      ## Calculate Robust Mahalanobis Distance ##
      med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))
      
      rob.dist.vals = try (apply(rmd.vals[,-1], 1, mal.fun))
      
      # Check if there was an error when calculating the the RMD.
      if (inherits(rob.dist.vals, "try-error")) {
        
        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "Proportion_Missing as one of the metrics. Try again,",
                    "excluding Proportion_Missing from the metrics vector.",
                    sep = " "))
        
      }
      
    }
    
  }
  
  # Check the number of metrics used. This value could be less than two if only
  # two metrics were input and one of them was Proportion_Missing. For example,
  # if Proportion_Missing was removed because there were no missing values then
  # there would only be one metric left and a PCA and RMD cannot be calculated
  # for only one metric.
  if (sum(metrics_final, na.rm = TRUE) < 2) {
    
    # Throw an error if only two metrics were given, one of them was
    # Proportion_Missing, and there is no missing data. Give the user a helpful
    # suggestion. We are quite nice :)
    stop (paste("Vector of metrics must contain at least two valid entries.",
                "Try including a metric other than Proportion_Missing.",
                sep = " "))
    
  }
  
  # Extract the names of the metrics that have been used. We cannot simply report
  # the metrics vector used as the input because it is possible that the input
  # metrics vector has had elements removed. For example, if the input data does
  # not have any missing values, and Proportion_Missing was an input metric, then
  # Proportion_Missing will not be calculated. Therefore, this metric will be
  # removed from the metrics_final_txt vector.
  metrics_final_txt <- names(rmd.vals[,-1])
  
  # RMD the heck out of the data -----------------------------------------------
  
  # Check if robust PCA, covariance matrix, and RMD objects have been created 
  if (!exists("robpca.res")) {
    
    ## Conduct Robust PCA ##
    robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1],
                                   k = (ncol(rmd.vals)-1),
                                   mcd = FALSE,
                                   scale = FALSE)
    
    ## Calculate Covariance Matrix #
    cov.mat = (robpca.res@loadings %*%
                 diag(robpca.res@eigenvalues) %*%
                 t(robpca.res@loadings))
    
    ## Calculate Robust Mahalanobis Distance ##
    med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))
    
    rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)
    
  }
  
  log2.dist.vals = log(rob.dist.vals, base = 2)
  
  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))
  
  temp.res = data.frame(Sample.ID = rmd.vals[,1],
                        Log2.md = log2.dist.vals,
                        pvalue = rmd.pvals,
                        rmd.vals[,-1])
  
  names(temp.res)[1] = get_fdata_cname(omicsData)
  
  output = merge(x = mintR_groupDF,
                 y = temp.res,
                 by = get_fdata_cname(omicsData),
                 all = TRUE)
  
  orig_class <- class(output)
  
  class(output) <- c("rmdFilt", orig_class)
  
  # Sum the number of metrics actually calculated.
  attr(output, "df") <- sum(!is.na(metrics_final))
  
  # Report the metrics actually used.
  attr(output, "metrics") <- metrics_final_txt
  
  return(output)
  
}

#' Calculate the Fraction of Missing Data of Sample Runs
#'
#' This function calculates the fraction of missing data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides and \eqn{n} is the number of samples.
#'
#' @return data.frame with two elements: Sample, a character vector giving the sample names; and Prop_missing, a numeric vector giving the fraction of missing values per run
#'
#' @author Lisa Bramer
#'
run_prop_missing <- function(data_only){
  
  # calculate number of missing values per run #
  nummiss <- apply(is.na(data_only), 2, sum)
  
  # calculate the fraction of missing values per run #
  fracmiss <- nummiss/nrow(data_only)
  
  # store data #
  res.final <- data.frame(Sample = names(data_only), Prop_missing = fracmiss, row.names = NULL)
  
  return(res.final)
}

#' Calculate the Median Absolute Deviance (MAD) of Sample Runs
#'
#' This function calculates the median absolute deviance across data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides and \eqn{n} is the number of samples.
#'
#' @details When calculating the MAD within a sample NA values are ignored. If all peptide abundance values are missing within a sample, the MAD is replaced by the overall mean MAD values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the sample names; and MAD, a numeric vector giving the MAD values
#'
#' @author Lisa Bramer
run_mad <- function(data_only){
  
  # calculate MAD #
  mad_val = apply(data_only, 2, function(x) median(abs(x - median(x, na.rm = T)), na.rm = T))
  
  # calculate the number of samples with MAD equal to NA #
  num.miss <- sum(is.na(mad_val))
  
  # if at least one sample has a MAD of NA, replace it with mean MAD value #
  if(num.miss > 0){
    mad_val[is.na(mad_val)] <- mean(mad_val, na.rm = T)
  }
  
  # store data #
  res.final <- data.frame(Sample = names(data_only), MAD = mad_val, row.names = NULL)
  
  return(res.final)
}

#' Calculate the Kurtosis of Sample Runs
#'
#' This function calculates the kurtosis across data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides and \eqn{n} is the number of samples.
#'
#' @details Kurtosis is calculated by method 2 in the \code{e1071} package, which is unbiased under normality. Within a sample NA values are ignorned in the kurtosis calculation. If all peptide abundance values are missing within a sample, the kurtosis is replaced by the overall mean of nonmissing kurtosis values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the sample names; and Kurtosis, a numeric vector giving the kurtosis
#'
#' @author Lisa Bramer
run_kurtosis <- function(data_only){
  
  # calculate kurtosis #
  kurt_res <- apply(data_only, 2, e1071::kurtosis, na.rm = TRUE, type = 2)
  
  # calculate the number of samples with kurtosis equal to NA #
  num.miss <- sum(is.na(kurt_res))
  
  # if at least one sample has a kurtosis of NA, replace it with mean kurtosis #
  if(num.miss > 0){
    kurt_res[is.na(kurt_res)] <- mean(kurt_res, na.rm = T)
  }
  
  # store data #
  res.final <- data.frame(Sample = names(data_only), Kurtosis = kurt_res, row.names = NULL)
  
  return(res.final)
}

#' Calculate the Skewness of Sample Runs
#'
#' This function calculates the skewness across data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides and \eqn{n} is the number of samples.
#'
#' @details Skewness is calculated as a bias-corrected calculation given by method 2 in the \code{e1071} package. Within a sample NA values are ignorned in the skewness calculation. If all peptide abundance values are missing within a sample, the skewness is replaced by the overall mean of nonmissing skewness values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the sample names; and Skewness, a numeric vector giving the skewness values
#'
#' @author Lisa Bramer
run_skewness <- function(data_only){
  
  # calculate skewness #
  skew_res <- apply(data_only, 2, e1071::skewness, na.rm = TRUE, type = 2)
  
  # calculate the number of samples with skewness equal to NA #
  num.miss <- sum(is.na(skew_res))
  
  # if at least one sample has a skewness of NA, replace it with mean skewness #
  if(num.miss > 0){
    skew_res[is.na(skew_res)] <- mean(skew_res, na.rm = T)
  }
  
  # store data #
  res.final <- data.frame(Sample = names(data_only), Skewness = skew_res, row.names = NULL)
  
  return(res.final)
}

#' Calculate the Mean Correlation of a Sample with Respect to Group
#'
#' This function calculates the mean correlation of a sample with all other samples that have the same group membership
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @param mintR_groupDF data.frame created by \code{\link{group_designation}} with columns for sample.id and group.
#' @param use_singletons logical indicator of whether to include singleton groups or not; defaults to FALSE. If FALSE, rmd_filter results are returned only for samples in groups of size greater than 1. This is a pass-through argument from \code{\link{rmd_filter}}.
#'
#' @details Correlation calculations use only complete pairwise observations.
#'
#' @return data.frame with two elements: Sample.ID, a character vector giving the sample names; and Mean_Correlation, a numeric vector giving the mean correlation values
#'
#' @author Lisa Bramer
#'
run_group_meancor <- function(omicsData, mintR_groupDF, ignore_singleton_groups = TRUE){
  
  # if group.data has column for TimeCourse, re-compute group.data including TimeCourse as main effect
  # if mintR_groupDF has TimeCourse variable, re-compute group_data including TimeCourse as main effect
  if(!is.null(attr(mintR_groupDF, "time_course"))){
    if(!is.null(attr(mintR_groupDF, "main_effects"))){
      temp_maineff = c(attr(mintR_groupDF, "main_effects"), attr(mintR_groupDF, "time_course"))}else{
        temp_maineff = attr(mintR_groupDF, "time_course")
      }
    group_data = group_designation(omicsData, temp_maineff)
  }
  
  # create a list of unique groups #
  grps = unique(as.character(mintR_groupDF$Group))
  
  samp_id = attr(omicsData, "cnames")$fdata_cname
  
  # compute all pairwise correlations between samples of the same group #
  nonsingleton_groups <- attributes(mintR_groupDF)$nonsingleton_groups
  singleton_groups <- setdiff(grps, nonsingleton_groups)
  if(length(singleton_groups) > 0 & ignore_singleton_groups == FALSE){
    ## there are singelton groups and we don't want to ignore them in creating the rmd filter ##
    
    # calculate the correlation of the singleton sample to all other samples and 
    # average that value to get the correlation for that particular sample/group
    omicsData_singletons <- omicsData
    omicsData_singletons$f_data$Dummy <- "dummy" # create a dummy grouping variable so that all samples belong to same group
    omicsData_singletons <- group_designation(omicsData_singletons, main_effect = "Dummy")
    
    prwse.grp.cors.all <- cor(omicsData_singletons$e_data[, -which(names(omicsData_singletons$e_data) == get_edata_cname(omicsData_singletons))], use = "pairwise.complete.obs")
    
    prwse.grp.cors.singletons <- list()
    for(i in 1:length(singleton_groups)){
      prwse.grp.cors.singletons[[i]] <- prwse.grp.cors.all
    }
    
    # do the "usual" thing for the nonsingleton groups: 
    # calculate the mean correlation of a sample with all other samples that have the same group membership
    myfilt <- custom_filter(omicsData, f_data_keep = as.character(mintR_groupDF[which(mintR_groupDF$Group %in% nonsingleton_groups), samp_id]))
    omicsData_nonsingletons <- applyFilt(myfilt, omicsData)
    
    # make a list of which columns belong to which groups #
    grp.col.ids = list()
    for(i in 1:length(nonsingleton_groups)){
      
      # pull sample names from group.data in current group #
      nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == nonsingleton_groups[i]), samp_id])
      
      # pull column numbers corresponding to above names #
      grp.col.ids[[i]] = which(names(omicsData_nonsingletons$e_data) %in% nms)
    }
    
    prwse.grp.cors.nonsingletons <- lapply(grp.col.ids, function(x){
      cor(omicsData_nonsingletons$e_data[,x], use = "pairwise.complete.obs")
    })
    
    
    ## combine the singleton and nonsingleton results...make sure the order is correct
    prwse.grp.cors <- c(prwse.grp.cors.nonsingletons, prwse.grp.cors.singletons)
    prws.grp.cors.grpnames <- c(nonsingleton_groups, singleton_groups) 
    
  }else{
    
    # make a list of which columns belong to which groups #
    grp.col.ids = list()
    for(i in 1:length(grps)){
      
      # pull sample names from group.data in current group #
      nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == grps[i]), samp_id])
      
      # pull column numbers corresponding to above names #
      grp.col.ids[[i]] = which(names(omicsData$e_data) %in% nms)
    }
    
    prwse.grp.cors = lapply(grp.col.ids, function(x){
      cor(omicsData$e_data[,x], use = "pairwise.complete.obs")
    })
    # structure of prws.grp.cors: list w/number elements equal to number groups
    # each element is correlation matrix
    # [[1]]
    # Mock1     Mock2     Mock3
    # Mock1 1.0000000 0.9670513 0.9703483
    # Mock2 0.9670513 1.0000000 0.9858240
    # Mock3 0.9703483 0.9858240 1.0000000
    # 
    # [[2]]
    # Infection1 Infection2 Infection3 Infection4 Infection5
    # Infection1  1.0000000  0.9802595  0.9757083  0.9790413  0.9732911
    # Infection2  0.9802595  1.0000000  0.9837869  0.9804636  0.9824165
    # Infection3  0.9757083  0.9837869  1.0000000  0.9786230  0.9801930
    # Infection4  0.9790413  0.9804636  0.9786230  1.0000000  0.9784125
    # Infection5  0.9732911  0.9824165  0.9801930  0.9784125  1.0000000
  }
  
  # turn diagonal into NAs, so we don't include a sample's correlation with itself #
  grp.cors = lapply(prwse.grp.cors, function(x) x*(matrix(1, nrow = nrow(x), ncol = ncol(x)) + diag(NA, nrow(x))))
  
  # compute mean correlation for each sample #
  mean.cor = lapply(grp.cors, function(x) apply(x, 1, mean, na.rm = T))
  
  ## need to adjust the list elements of mean.cor for any singleton groups, 
  ## to just contain the value for the sample in that group
  if(length(singleton_groups) > 0 & ignore_singleton_groups == FALSE){
    mean.cor2 <- mean.cor
    # when I concatenated the pairwise group correlations, I put the singleton groups last #
    for(i in 1:length(singleton_groups)){
      # get the sample name in the current singleton group, and pull that value out of mean.cor #
      cur_singleton <- singleton_groups[i]
      cur_sample <- as.character(omicsData_singletons$f_data[which(mintR_groupDF$Group == cur_singleton), samp_id])
      
      # get the element in mean.cor list that corresponds to this singleton group #
      j <- which(prws.grp.cors.grpnames == cur_singleton)
      mean.cor2[[j]] <- mean.cor[[j]][cur_sample]
      # pull sample names from group.data in current group #
      # nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == nonsingleton_groups[i]), samp_id])
      
    }
  }else{
    mean.cor2 <- mean.cor
  }
  
  # format results #
  # get order to put results in original sample order based on peptide.data #
  temp = match(names(omicsData$e_data)[-1],names(unlist(mean.cor2)))
  res.cor = data.frame(Sample.ID = names(omicsData$e_data)[-1], Mean_Correlation = unlist(mean.cor2)[temp], row.names = NULL)
  
  return(res.cor)
}
