#' Quantile Normalization
#'
#' Perform quantile normalization 
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#'
#' @details Quantile normalization is an algorithm for normalizing a set of data vectors by giving them the same distribution. It is applied to data on the abundance scale (e.g. not a log scale). It is often used for microarry data. 
#' @return The normalized data is returned in an object of the appropriate S3 class (e.g. pepData), on the same scale as omicsData (e.g. if omicsData contains log2 transformed data, the normalization will be performed on the non-log2 scale and then re-scaled after normalization to be returned on the log2 scale).
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' norm_data <- normalize_quantile(omicsData = lipid_object)
#'}
#'
#' @author Kelly Stratton
#' @references
#'
#' @export
normalize_quantile <- function(omicsData){
  ## initial checks ##
  
  omicsData_orig <- omicsData
  
  # Store data class as it will be referred to often
  dat_class <- class(omicsData)
  
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # data should be on raw scale #
  if(attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")){
    omicsData <- edata_transform(omicsData, "abundance")
  }
  
  # give message if proportion of missing data is > 0
  if(attributes(omicsData)$data_info$prop_missing > 0){
    message(paste("The proportion of missing data is ", round(attributes(omicsData)$data_info$prop_missing, 2), ". Quantile normalization only works with complete data. Consider using SPANS to choose an appropriate normalization method for a dataset that includes missing values.", sep=""))
  }
  
  check_names <- getchecknames(omicsData)
  edata_id <- attr(omicsData, "cnames")$edata_cname
  samp_id <- attr(omicsData, "cnames")$fdata_cname
  
  ## end of initial checks ##
  
  edata <- omicsData$e_data
  edata_cnames <- names(edata)
  edata_rnames <- row.names(edata)
  edata_idcol <- as.character(edata[, edata_id])
  
  # 1. transpose data so it is samples by molecules #
  x <- t(edata[, -1])
  
  # 2. sort each column of x to get x_sort # 
  x_sort <- apply(x, 2, function(c) sort(c, na.last = FALSE))
  x_ord <- apply(x, 2, function(c) order(c, na.last = FALSE))
  
  for(i in 1:ncol(x)){
    x_ord[, i] <- match(x_sort[, i], x[, i])
  }
  
  # 3. take means across rows of x_sort and assign the mean to each element in the row to get x_sort_prime #
  row_means <- rowMeans(x_sort, na.rm = TRUE)
  x_sort_prime <- array(row_means, dim = dim(x_sort))
  
  # 4. get x_norm by rearranging each column of x_sort_prime to have same ordering as the original x #
  x_norm <- x
  temp_x_sort_prime <- rep(NA, nrow(x_norm))
  for(i in 1:ncol(x_norm)){
    temp_x_sort_prime[x_ord[, i]] <- x_sort_prime[, i]
    x_norm[, i] <- temp_x_sort_prime
    x_norm[, i] <- as.numeric(x_norm[, i])
  }
  
  # re-assemble omicsData #
  edata_norm <- data.frame(cbind(edata_idcol, t(x_norm)))
  names(edata_norm) <- edata_cnames
  
  edata_norm[, edata_id] <- as.character(edata_norm[, edata_id])
  edata_norm[, -which(names(edata_norm) == edata_id)] <- apply(edata_norm[, -which(names(edata_norm) == edata_id)], 1, as.numeric)
  
  
  # assign attributes # 
  omicsData$e_data <- edata_norm
  attributes(omicsData)$data_info$norm_info$is_normalized <- TRUE
  attributes(omicsData)$data_info$norm_info$norm_type <- "quantile" # new attribute as of 12/21/17
  # none of these other attributes are applicable for quantile normalization #
  # attributes(omicsData)$data_info$norm_info$subset_fn = subset_fn
  # attributes(omicsData)$data_info$norm_info$subset_params = params
  # attributes(omicsData)$data_info$norm_info$norm_fn = norm_fn
  # attributes(omicsData)$data_info$norm_info$n_features_calc = length(peps)
  # attributes(omicsData)$data_info$norm_info$prop_features_calc = length(peps)/nrow(omicsData$e_data)
  # attributes(omicsData)$data_info$norm_info$params$norm_scale = norm_results$norm_params$scale
  # attributes(omicsData)$data_info$norm_info$params$norm_location = norm_results$norm_params$location
  # attributes(omicsData)$data_info$norm_info$params$bt_scale = norm_results$backtransform_params$scale
  # attributes(omicsData)$data_info$norm_info$params$bt_location = norm_results$backtransform_params$location
  
  if(attributes(omicsData_orig)$data_info$data_scale != "abundance"){
    omicsData <- edata_transform(omicsData, attributes(omicsData_orig)$data_info$data_scale)
  }
  
  res <- omicsData
  
  return(res)
  
}