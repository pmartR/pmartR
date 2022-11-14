#'Creates R markdown document report
#'
#'@param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'  'lipidData', 'nmrData', or 'seqData', created by
#'  \code{\link{as.pepData}}, \code{\link{as.proData}},
#'  \code{\link{as.metabData}}, \code{\link{as.lipidData}},
#'  \code{\link{as.nmrData}}, or \code{\link{as.seqData}}, respectively.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(omicsData = metab_object)
#' result <- analysis_log(omicsData = metab_object)
#'}
#'
#' @export
#'
analysis_log <- function(omicsData) {
  library(rmarkdown)
  
  # check that omicsData is of correct class #
  if (!inherits(omicsData,
                c(
                  "pepData",
                  "proData",
                  "metabData",
                  "lipidData",
                  "nmrData",
                  "seqData"
                )))
    stop(
      "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', or 'seqData'."
    )
  
  data <- omicsData
  classes <- class(data)
  
  params <- list(data = data, classes = classes)
  
  render("vignettes/analysis_log.Rmd", params = params)
  
}  