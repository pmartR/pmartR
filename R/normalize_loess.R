#' Loess Normalization
#'
#' Perform Loess normalization 
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#' @param method character string specifying which variant of the cyclic loess method to use. Options are "fast" (default), "affy", or "pairs"
#' @param span span of loess smoothing window, between 0 and 1.
#' 
#' @details A wrapper for the normalizeCyclicLoess function from the limma package. 
#' 
#' @return The normalized data is returned in an object of the appropriate S3 class (e.g. pepData), on the same scale as omicsData (e.g. if omicsData contains log2 transformed data, the normalization will be performed on the non-log2 scale and then re-scaled after normalization to be returned on the log2 scale).
#' 
#' @examples
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' pep_object = edata_transform(pep_object, "log2")
#' result = normalize_loess(pep_object)
#'
#'}
#'
#'@seealso \code{\link{normalizeCyclicLoess}} in the \code{limma} package
#'
#' @references Bolstad, B. M., Irizarry R. A., Astrand, M., and Speed, T. P. (2003). \emph{A comparison of normalization methods for high density oligonucleotide array data based on bias and variance.} Bioinformatics 19, 185-193. Ballman, KV Grill, DE, Oberg, AL and Therneau, TM (2004). \emph{Faster cyclic loess: normalizing RNA arrays via linear models.} Bioinformatics 20, 2778-2786. 
#'
#' @export
normalize_loess<- function(omicsData, method = "fast", span = .4){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData", "nmrData"))) stop("omicsData is not an object of appropriate class")
  
  # data should be on log scale #
  if(!(attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log"))) stop("omicsData must be on log scale")

  #check that norm attribute is False
  if(get_data_norm(omicsData)) stop("data has already been normalized")
  # if(attr(omicsData, "data_info")$norm_info$is_normalized == TRUE) stop("data has already been normalized")
  
  #extract cname
  edata_cname = attr(omicsData, "cnames")$edata_cname
  
  #remove edata_cname col from e_data
  e_data = omicsData$e_data
  e_data = e_data[, -which(names(e_data) == edata_cname)]
  
  #apply normalizeCyclicLoess function to e_data
  result = limma::normalizeCyclicLoess(e_data, span = span, method = method)
  
  #update e_data with normalized result
  omicsData$e_data[,-which(names(omicsData$e_data) == edata_cname)] = result
  
  #update norm attribute in omicsData
  attr(omicsData, "data_info")$norm_info$is_normalized = TRUE
  
  return(omicsData)
}