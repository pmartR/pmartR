#' Correlation matrix of biomolecule data
#'
#' This function returns an object of class corRes (correlation Result)
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData' created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#'
#' @details The Pearson correlation between samples is calculated based on
#'   biomolecules that are observed in both samples.  the \eqn{n \times n}
#'   correlation matrix of normalized data. See \code{\link{cor}} for further
#'   details.
#'
#' @return An \eqn{n \times n} matrix giving the correlation between samples.
#'
#' @examples
#' \dontrun{
#' library(pmartR)
#' 
#' data(pep_object)
#' 
#' my_correlation <- cor_result(omicsData = pep_object)
#' }
#' 
#' @author Kelly Stratton, Lisa Bramer
#'
#' @export
#' 
cor_result <- function(omicsData){

  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  # Check the data scale. Data must be on one of the log scales.
  if (get_data_scale(omicsData) == "abundance") {
    
    # Welcome to the pit of despair!
    stop ("Data must be on the log scale.")
    
  }

  edata_cname <- attr(omicsData, "cnames")$edata_cname
  id_col <- which(colnames(omicsData$e_data) == edata_cname)

  output <- cor(omicsData$e_data[, -id_col],
                use = "pairwise.complete.obs")

  rownames(output) <- names(omicsData$e_data[, -id_col])

  orig_class <- class(output)

  class(output) <- c("corRes", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data[, -id_col])
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "is_normalized") <- attr(omicsData, "data_info")$norm_info$is_normalized
  
  if (inherits(omicsData, "isobaricpepData")) {
    
    attr(output, "isobaric_norm") <- attr(omicsData,"isobaric_info")$norm_info$is_normalized
    attr(output, "is_normalized") <- attr(omicsData,"data_info")$norm_info$is_normalized
    
  } else if (inherits(omicsData, "nmrData")) {
    
    attr(output, "nmr_norm") <- attr(omicsData, "nmr_info")$norm_info$is_normalized
    attr(output, "is_normalized") <- attr(omicsData,"data_info")$norm_info$is_normalized
    
  }
  
  return (output)
  
}
