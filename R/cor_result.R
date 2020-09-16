#' Correlation matrix of biomolecule data
#'
#' This function returns an object of class corRes (correlation Result)
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData' created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#'
#' @details The Pearson correlation between samples is calculated based on biomolecules that are observed in both samples.  the \eqn{n \times n} correlation matrix of normalized data. See \code{\link{cor}} for further details.
#'
#' @return An \eqn{n \times n} matrix giving the correlation between samples.
#'
#' @examples
#' dontrun{
#' library(pmartR)
#' data(pep_object)
#' my_correlation <- cor_result(omicsData = pep_object)
#'}
#' @author Kelly Stratton, Lisa Bramer
#'
#' @export

cor_result <- function(omicsData){

  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")

  edata_id = attr(omicsData, "cnames")$edata_cname
  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata) == edata_id)]

  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata) == edata_id)]

  CM = cor(edata, use = "pairwise.complete.obs")

  rownames(CM) <- names(edata)
  output <- CM

  orig_class <- class(output)

  class(output) <- c("corRes", orig_class)

  attr(output, "sample_names") <- names(edata)
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
 
   if(inherits(omicsData, "isobaricpepData")){
    attr(output, "isobaric_norm") <- ifelse(attr(omicsData,"isobaric_info")$norm_info$is_normalized == TRUE, TRUE, FALSE)
    attr(output, "is_normalized") <- ifelse(attr(omicsData,"data_info")$norm_info$is_normalized == TRUE, TRUE, FALSE)
  }
  else if(inherits(omicsData, "pepData")){
    attr(output, "is_normalized") <- ifelse(attr(omicsData,"data_info")$norm_info$is_normalized == TRUE, TRUE, FALSE)
  }


  return(output)
}
