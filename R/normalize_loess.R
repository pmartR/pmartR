#' Loess Normalization
#'
#' Perform Loess normalization 
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#' @param min_prop numeric threshold between 0 and 1 giving the minimum value for the proportion of features subset (rows of \code{e_data}); default value is 0.4.
#'
#' @details Intensity-dependent normalization. We use the robust scatter plot smoother ‘lowess’, implemented in the statistical software package R (17), to perform a local A-dependent normalization log2R/G → log2R/G – c(A) = log2R/[k(A)G] where c(A) is the lowess fit to the MA-plot. The lowess scatter plot smoother performs robust locally linear fits. In particular, it will not be affected by a small percentage of differentially expressed genes, which will appear as outliers in the MA-plot. The user-defined parameter f is the fraction of the data used for smoothing at each point; the larger the f value, the smoother the fit. We typically use f  = 40%. It is applied to data on the abundance scale (e.g. not a log scale). It is often used for microarry data. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC100354/)
#'  
#' @references Strand, M. (2003). \emph{Contrast normalization of oligonucleotide arrays}. Journal of Computational Biology. 10(1): 95-102. Laurent, G., Leslie C., Benjamin, M. B., Rafael, A. I. (2004). \emph{affy—analysis of Affymetrix GeneChip data at the probe level}. Bioinformatics. 20(3): 307–315.  
#' 
#' @return The normalized data is returned in an object of the appropriate S3 class (e.g. pepData), on the same scale as omicsData (e.g. if omicsData contains log2 transformed data, the normalization will be performed on the non-log2 scale and then re-scaled after normalization to be returned on the log2 scale).
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' norm_data <- normalize_lowess(omicsData = lipid_object, subset_fn = "all", norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
#'}
#'
#' @author Kelly Stratton
#' @references
#'
#' @export
normalize_loess <- function(omicsData, min_prop = 0.4){
  ## initial checks ##
  
  omicsData_orig <- omicsData # keep the "orig" version to grab attributes...will need this later, especially to make sure data goes back onto the scale it came into the function on
  
  # Store data class as it will be referred to often
  dat_class <- class(omicsData)
  
  # check that omicsData is of the appropriate class
  if(!(inherits(dat_class, c("proData","pepData","lipidData", "metabData")))) stop("omicsData is not an object of appropriate class")
  
  # data should be on raw scale #
  if(attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")){
    omicsData <- edata_transform(omicsData, "abundance")
  }
  
  # check that min_prop is greater than 0 but less than or equal to 1
  if(!is.null(min_prop)){
    if(min_prop <= 0 | min_prop > 1) stop("min_prop must be greater than zero but less than or equal to 1")
  }
  
  # give message if proportion of missing data is > 0.1
  if(attributes(omicsData)$data_info$prop_missing > 0.05){
    message(paste("The proportion of missing data is ", round(attributes(omicsData)$data_info$prop_missing, 2), ". We do not recommend using quantile normalization when the proportion of missing data exceeds 0.05 or 0.1, as it may result in heavily skewed normalized data that does not accurately respresent the data observed.", sep=""))
  }
  
  check_names <- getchecknames(omicsData)
  edata_id <- attr(omicsData, "cnames")$edata_cname
  samp_id <- attr(omicsData, "cnames")$fdata_cname
  
  ## end of initial checks ##
  
  edata <- omicsData$e_data
  edata_cnames <- names(edata)
  edata_rnames <- row.names(edata)
  edata_idcol <- as.character(edata[, edata_id])
  
  # apply normalization #
  
  # re-assemble omicsData #
  edata_norm <- data.frame(cbind(edata_idcol, t(x_norm)))
  names(edata_norm) <- edata_cnames
  
  edata_norm[, edata_id] <- as.character(edata_norm[, edata_id])
  edata_norm[, -which(names(edata_norm) == edata_id)] <- apply(edata_norm[, -which(names(edata_norm) == edata_id)], 1, as.numeric)
  
  
  # assign attributes # 
  omicsData$e_data <- edata_norm
  attributes(omicsData)$data_info$data_norm <- TRUE
  attributes(omicsData)$norm_info$norm_type <- "quantile" # new attribute as of 12/21/17
  # none of these other attributes are applicable for quantile normalization #
  # attributes(omicsData)$norm_info$subset_fn = subset_fn
  # attributes(omicsData)$norm_info$subset_params = params
  # attributes(omicsData)$norm_info$norm_fn = norm_fn
  # attributes(omicsData)$norm_info$n_features_calc = length(peps)
  # attributes(omicsData)$norm_info$prop_features_calc = length(peps)/nrow(omicsData$e_data)
  # attributes(omicsData)$norm_info$params$norm_scale = norm_results$norm_params$scale
  # attributes(omicsData)$norm_info$params$norm_location = norm_results$norm_params$location
  # attributes(omicsData)$norm_info$params$bt_scale = norm_results$backtransform_params$scale
  # attributes(omicsData)$norm_info$params$bt_location = norm_results$backtransform_params$location
  
  if(attributes(omicsData_orig)$data_info$data_scale != "abundance"){
    omicsData <- edata_transform(omicsData, attributes(omicsData_orig)$data_info$data_scale)
  }
  
  res <- omicsData
  
  return(res)
  
}