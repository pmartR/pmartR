#' RMD Runs
#'
#' The method computes a robust Mahalanobis distance that can be mapped to a p-value and used to identify outlying samples
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @param ... further arguments
#'
#' @return a data.frame containing columns for the Sample ID, log2 robust Mahalanobis distance, p-values, and robust Mahalanobis distance
#'
#' @details The metrics on which the log2 robust Mahalanobis distance is based can be specified using the \code{metrics} argument.
#' \tabular{ll}{
#' pepData \tab For pepData objects, all five of the metrics "MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing" may be used (this is the default). \cr
#' proData \tab For proData objects, all five of the metrics "MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing" may be used (this is the default). \cr
#' metabData \tab For metabData objects, the use of "Proportion_Missing" is discouraged due to the general lack of missing data in metabolomics datasets (the default behavior omits "Proportion_Missing" from the metrics). \cr
#' lipidData \tab For lipidData objects, , the use of "Proportion_Missing" is discouraged due to the general lack of missing data in metabolomics datasets (the default behavior omits "Proportion_Missing" from the metrics). \cr
#' }
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
#' @references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R., Pounds, J., and Webb-Robertson, B.J. (2011), \emph{Improved quality control processing of peptide-centric LC-MS proteomics data}. Bioinformatics. 27(20): 2866-2872.
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
#'
#'
#'@rdname rmd_filter
#'@name rmd_filter
#'
rmd_filter <- function(omicsData, ...){

  # group_DF attribute is required #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("omicsData must contain attribution information for 'group_DF'. See documentation for group_designation function for more information.")
  }

  # data should be log transformed #
  if(!attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")){
    stop("omicsData$e_data should be log transformed prior to calling rmd_filter. See documentation for edata_transform function for more information.")
  }

  UseMethod("rmd_filter")
}


### Peptide data function ###
#'@export
#'@rdname rmd_filter
#'@name rmd_filter
rmd_filter.pepData <- function(omicsData, metrics=c("MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing")){

  ## some initial checks ##

  # check that omicsData is of class 'pepData', 'proData', 'lipidata', or 'metabData' #
  if(!inherits(omicsData, c("pepData","proData","lipidData","metabData"))) stop("omicsData is not of an appropriate class")

  if(length(metrics) < 2){
    stop("Vector of metrics must contain at least two elements.")
  }else{
    if(length(metrics) > 5){
      stop("Vector of metrics cannot contain more than five elements.")
    }
  }

  mintR_groupDF = attr(omicsData, "group_DF")

  # run_group_meancor already deals with the TimeCourse option #
  edata_id = attr(omicsData, "cnames")$edata_cname
  edata_col_id = which(names(omicsData$e_data) == edata_id)
  dat_only = omicsData$e_data[,-(edata_col_id)]

  samp_id = attr(omicsData, "cnames")$fdata_cname

  ## Calculate metrics for RMD runs ##
  # match the metrics (allow upper/lower case) #
  metrics = tolower(metrics)
  metrics_final = rep(NA, 5)

  rmd.vals = data.frame(Sample.ID = run_prop_missing(dat_only)[,1])

  if(any(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of MAD.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(dat_only)[,2]
  }

  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Kurtosis.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(dat_only)[,2]
  }

  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Skew.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(dat_only)[,2]
  }

  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Correlation.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = run_group_meancor(omicsData, mintR_groupDF)[,2]
  }

  if(any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))){
    ind = which(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Proportion_Missing.")}
    metrics = metrics[-ind]
    metrics_final[5] = 1
    rmd.vals$Proportion_Missing = run_prop_missing(dat_only)[,2]
  }

  if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

  metrics_final_txt <- names(rmd.vals[,-1])

  ## Conduct Robust PCA ##
  robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)

  ## Calculate Covariance Matrix #
  cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)

  ## Calculate Robust Mahalanobis Distance ##
  med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

  mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

  rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)

  log2.dist.vals = log(rob.dist.vals, base = 2)

  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))

  temp.res = data.frame(Sample.ID = rmd.vals[,1], Log2.md = log2.dist.vals, pvalue = rmd.pvals, rmd.vals[,-1])
  names(temp.res)[1] = samp_id

  output = merge(x = mintR_groupDF, y = temp.res, by = samp_id, all = TRUE)


  orig_class <- class(output)

  class(output) <- c("rmdFilt", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data)
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "df") <- sum(!is.na(metrics_final))
  attr(output, "metrics") <- metrics_final_txt

  return(output)
}



### Prodata function ###
#'@export
#'@rdname rmd_filter
#'@name rmd_filter
rmd_filter.proData <- function(omicsData, metrics=c("MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing")){

  ## some initial checks ##

  # check that omicsData is of class 'pepData', 'proData', 'lipidata', or 'metabData' #
  if(!inherits(omicsData, c("pepData","proData","lipidData","metabData"))) stop("omicsData is not of an appropriate class")

  if(length(metrics) < 2){
    stop("Vector of metrics must contain at least two elements.")
  }else{
    if(length(metrics) > 5){
      stop("Vector of metrics cannot contain more than five elements.")
    }
  }

  mintR_groupDF = attr(omicsData, "group_DF")

  # run_group_meancor already deals with the TimeCourse option #
  edata_id = attr(omicsData, "cnames")$edata_cname
  edata_col_id = which(names(omicsData$e_data) == edata_id)
  dat_only = omicsData$e_data[,-(edata_col_id)]

  samp_id = attr(omicsData, "cnames")$fdata_cname

  ## Calculate metrics for RMD runs ##
  # match the metrics (allow upper/lower case) #
  metrics = tolower(metrics)
  metrics_final = rep(NA, 5)

  rmd.vals = data.frame(Sample.ID = run_prop_missing(dat_only)[,1])

  if(any(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of MAD.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(dat_only)[,2]
  }

  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Kurtosis.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(dat_only)[,2]
  }

  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Skew.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(dat_only)[,2]
  }

  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Correlation.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = run_group_meancor(omicsData, mintR_groupDF)[,2]
  }

  if(any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))){
    ind = which(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Proportion_Missing.")}
    metrics = metrics[-ind]
    metrics_final[5] = 1
    rmd.vals$Proportion_Missing = run_prop_missing(dat_only)[,2]
  }

  if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

  metrics_final_txt <- names(rmd.vals[,-1])

  ## Conduct Robust PCA ##
  robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)

  ## Calculate Covariance Matrix #
  cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)

  ## Calculate Robust Mahalanobis Distance ##
  med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

  mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

  rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)

  log2.dist.vals = log(rob.dist.vals, base = 2)

  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))

  temp.res = data.frame(Sample.ID = rmd.vals[,1], Log2.md = log2.dist.vals, pvalue = rmd.pvals, rmd.vals[,-1])
  names(temp.res)[1] = samp_id

  output = merge(x = mintR_groupDF, y = temp.res, by = samp_id, all = TRUE)

  orig_class <- class(output)

  class(output) <- c("rmdFilt", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data)
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "df") <- sum(!is.na(metrics_final))
  attr(output, "metrics") <- metrics_final_txt

  return(output)
}



### Lipiddata function ###
#'@export
#'@rdname rmd_filter
#'@name rmd_filter
rmd_filter.lipidData <- function(omicsData, metrics=c("MAD", "Kurtosis", "Skewness", "Correlation")){

  ## some initial checks ##
  
  # check that omicsData is of class 'pepData', 'proData', 'lipidata', or 'metabData' #
  if(!inherits(omicsData, c("pepData","proData","lipidData","metabData"))) stop("omicsData is not of an appropriate class")

  if(length(metrics) < 2){
    stop("Vector of metrics must contain at least two elements.")
  }else{
    if(length(metrics) > 5){
      stop("Vector of metrics cannot contain more than five elements.")
    }
  }

  mintR_groupDF = attr(omicsData, "group_DF")

  # run_group_meancor already deals with the TimeCourse option #
  edata_id = attr(omicsData, "cnames")$edata_cname
  edata_col_id = which(names(omicsData$e_data) == edata_id)
  dat_only = omicsData$e_data[,-(edata_col_id)]

  samp_id = attr(omicsData, "cnames")$fdata_cname

  ## Calculate metrics for RMD runs ##
  # match the metrics (allow upper/lower case) #
  metrics = tolower(metrics)
  metrics_final = rep(NA, 5)

  rmd.vals = data.frame(Sample.ID = run_prop_missing(dat_only)[,1])

  if(any(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of MAD.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(dat_only)[,2]
  }

  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Kurtosis.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(dat_only)[,2]
  }

  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Skew.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(dat_only)[,2]
  }

  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Correlation.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = run_group_meancor(omicsData, mintR_groupDF)[,2]
  }

  if(any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))){

    # check to see whether there is any missing data
    number_missing = sum(is.na(omicsData$e_data))
    # if no missing data, we flat out refuse to include prop_missing in the metrics, so give a warning to this effect
    if(number_missing == 0){
      warning("There are no missing values in e_data, therefore prop_missing will not be used as one of the metrics for rMd-Runs.")
    }else{

      ind = which(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))
      if(length(ind)>1){stop("More than one of the entries in metric matches to use of Proportion_Missing.")}
      metrics = metrics[-ind]
      metrics_final[5] = 1
      rmd.vals$Proportion_Missing = run_prop_missing(dat_only)[,2]

      ## proceed, to check the rank of cov.mat ##
      # in case the user input 'enough' metrics (2-5) but they were not recognized by the code as valid entries #
      if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

      ## Conduct Robust PCA ##
      robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)
      if(inherits(robpca.res, "try-error")){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }

      ## Calculate Covariance Matrix #
      cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)
      # check the rank #
      myrank = qr(cov.mat)$rank

      if(myrank != length(metrics_final)){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }

      ## Calculate Robust Mahalanobis Distance ##
      med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

      mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

      rob.dist.vals = try(apply(rmd.vals[,-1], 1, mal.fun))
      if(inherits(rob.dist.vals, "try-error")){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }
    }
  }

  # in case the user input 'enough' metrics (2-5) but they were not recognized by the code as valid entries
  if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

  metrics_final_txt <- names(rmd.vals[,-1])

  ## Conduct Robust PCA ##
  robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)

  ## Calculate Covariance Matrix #
  cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)

  ## Calculate Robust Mahalanobis Distance ##
  med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

  mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

  rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)

  log2.dist.vals = log(rob.dist.vals, base = 2)

  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))

  temp.res = data.frame(Sample.ID = rmd.vals[,1], Log2.md = log2.dist.vals, pvalue = rmd.pvals, rmd.vals[,-1])
  names(temp.res)[1] = samp_id

  output = merge(x = mintR_groupDF, y = temp.res, by = samp_id, all = TRUE)

  orig_class <- class(output)

  class(output) <- c("rmdFilt", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data)
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "df") <- sum(!is.na(metrics_final))
  attr(output, "metrics") <- metrics_final_txt

  return(output)
}



### Metabdata function ###
#'@export
#'@rdname rmd_filter
#'@name rmd_filter
rmd_filter.metabData <- function(omicsData, metrics=c("MAD", "Kurtosis", "Skewness", "Correlation")){

  ## some initial checks ##

  # check that omicsData is of class 'pepData', 'proData', 'lipidata', or 'metabData' #
  if(!inherits(omicsData, c("pepData","proData","lipidData","metabData"))) stop("omicsData is not of an appropriate class")

  if(length(metrics) < 2){
    stop("Vector of metrics must contain at least two elements.")
  }else{
    if(length(metrics) > 5){
      stop("Vector of metrics cannot contain more than five elements.")
    }
  }

  mintR_groupDF = attr(omicsData, "group_DF")

  # run_group_meancor already deals with the TimeCourse option #
  edata_id = attr(omicsData, "cnames")$edata_cname
  edata_col_id = which(names(omicsData$e_data) == edata_id)
  dat_only = omicsData$e_data[,-(edata_col_id)]

  samp_id = attr(omicsData, "cnames")$fdata_cname

  ## Calculate metrics for RMD runs ##
  # match the metrics (allow upper/lower case) #
  metrics = tolower(metrics)
  metrics_final = rep(NA, 5)
  metrics_orig = metrics

  rmd.vals = data.frame(Sample.ID = run_prop_missing(dat_only)[,1])

  if(any(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation", "median_absolute_deviation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of MAD.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(dat_only)[,2]
  }

  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Kurtosis.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(dat_only)[,2]
  }

  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Skew.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(dat_only)[,2]
  }

  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){stop("More than one of the entries in metric matches to use of Correlation.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = run_group_meancor(omicsData, mintR_groupDF)[,2]
  }

  if(any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))){

    # check to see whether there is any missing data
    number_missing = sum(is.na(omicsData$e_data))
    # if no missing data, we flat out refuse to include prop_missing in the metrics, so give a warning to this effect
    if(number_missing == 0){
      warning("There are no missing values in e_data, therefore prop_missing will not be used as one of the metrics for rMd-Runs.")
    }else{

      ind = which(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss", "proportion_miss", "prop", "proportion", "proportion missing", "prop miss", "proportion miss"))
      if(length(ind)>1){stop("More than one of the entries in metric matches to use of Proportion_Missing.")}
      metrics = metrics[-ind]
      metrics_final[5] = 1
      rmd.vals$Proportion_Missing = run_prop_missing(dat_only)[,2]

      ## proceed, to check the rank of cov.mat ##
      # in case the user input 'enough' metrics (2-5) but they were not recognized by the code as valid entries #
      if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

      ## Conduct Robust PCA ##
      robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)
      if(inherits(robpca.res, "try-error")){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }

      ## Calculate Covariance Matrix #
      cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)
      # check the rank #
      myrank = qr(cov.mat)$rank

      if(myrank != length(metrics_final)){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }

      ## Calculate Robust Mahalanobis Distance ##
      med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

      mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

      rob.dist.vals = try(apply(rmd.vals[,-1], 1, mal.fun))
      if(inherits(rob.dist.vals, "try-error")){
        stop("There are not enough missing values in e_data to use prop_missing as one of the metrics. Try again, excluding prop_missing.")
      }
    }
  }

  if(sum(metrics_final, na.rm=TRUE)<2){stop("Vector of metrics must contain at least two valid entries.")}

  metrics_final_txt <- names(rmd.vals[,-1])

  ## Conduct Robust PCA ##
  robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1], k = (ncol(rmd.vals)-1), mcd = FALSE, scale = FALSE)

  ## Calculate Covariance Matrix #
  cov.mat = robpca.res@loadings %*% diag(robpca.res@eigenvalues) %*% t(robpca.res@loadings)

  ## Calculate Robust Mahalanobis Distance ##
  med.mat = matrix(apply(rmd.vals[,-1], 2, median), nrow = (ncol(rmd.vals)-1))

  mal.fun <- function(x) t(as.vector(x) - med.mat) %*% solve(as.matrix(cov.mat)) %*% (as.vector(x) - med.mat)

  rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)

  log2.dist.vals = log(rob.dist.vals, base = 2)

  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))

  temp.res = data.frame(Sample.ID = rmd.vals[,1], Log2.md = log2.dist.vals, pvalue = rmd.pvals, rmd.vals[,-1])
  names(temp.res)[1] = samp_id

  output = merge(x = mintR_groupDF, y = temp.res, by = samp_id, all = TRUE)

  orig_class <- class(output)

  class(output) <- c("rmdFilt", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data)
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "df") <- sum(!is.na(metrics_final))
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
#' @details When calculating the MAD within a sample NA values are ignorned. If all peptide abundance values are missing within a sample, the MAD is replaced by the overall mean MAD values for the data.
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
#'
#' @details Correlation calculations use only complete pairwise observations.
#'
#' @return data.frame with two elements: Sample.ID, a character vector giving the sample names; and Mean_Correlation, a numeric vector giving the mean correlation values
#'
#' @author Lisa Bramer
#'


run_group_meancor <- function(omicsData, mintR_groupDF){

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
  # make a list of which columns belong to which groups #
  grp.col.ids = list()
  for(i in 1:length(grps)){

    # pull sample names from group.data in current group #
    nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == grps[i]), samp_id])

    # pull column numbers corresponding to above names #
    grp.col.ids[[i]] = which(names(omicsData$e_data) %in% nms)
  }

  # compute all pairwise correlations between samples of the same group #
  prwse.grp.cors = lapply(grp.col.ids, function(x){
    cor(omicsData$e_data[,x], use = "pairwise.complete.obs")
  })

  # turn diagonal into NAs, so we don't include a sample's correlation with itself #
  grp.cors = lapply(prwse.grp.cors, function(x) x*(matrix(1, nrow = nrow(x), ncol = ncol(x)) + diag(NA, nrow(x))))

  # compute mean correlation for each sample #
  mean.cor = lapply(grp.cors, function(x) apply(x, 1, mean, na.rm = T))

  # format results #
  # get order to put results in original sample order based on peptide.data #
  temp = match(names(omicsData$e_data)[-1],names(unlist(mean.cor)))
  res.cor = data.frame(Sample.ID = names(omicsData$e_data)[-1], Mean_Correlation = unlist(mean.cor)[temp], row.names = NULL)

  return(res.cor)
}
