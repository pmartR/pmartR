#' Creates R markdown document report
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' result = analysis_log(metab_object)
#'}
#'
#' @export
#' 

analysis_log<- function(omicsData){
  library(rmarkdown)
  
  # check that omicsData is of correct class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'.")
  
  data <- omicsData
  classes <- class(data)
  
  params <- list(data=data, classes=classes)
  
  render("vignettes/analysis_log.Rmd", params = params)
  
}  