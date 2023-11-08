#' Compute Correlation matrix of biomolecule data
#'
#' This function returns an object of class corRes (correlation Result)
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData',
#'   'proData', 'nmrData', or 'seqData', created by
#'   \code{\link{as.lipidData}}, \code{\link{as.metabData}},
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.nmrData}}, or \code{\link{as.seqData}}, respectively.
#'
#' @details The pairwise correlations between samples are calculated based on
#'   biomolecules that are observed in both samples. For seqData objects,
#'   Spearman correlation is used. For all other data types, Pearson correlation
#'   is used and data must be log transformed. See \code{\link{cor}} for further details.
#'
#' @return An \eqn{n \times n} matrix of class corRes giving the correlation
#'   between samples.
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#'
#' mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' my_correlation <- cor_result(omicsData = mymetab)
#'
#' \dontrun{
#' myseq_correlation <- cor_result(omicsData = rnaseq_object)
#' }
#'
#' @author Kelly Stratton, Lisa Bramer
#'
#' @seealso \code{\link{edata_transform}}
#'
#' @export
#'
cor_result <- function(omicsData) {
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData", "seqData")))
    stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', or 'seqData'")

  # Check the data scale. Data must be on one of the log scales.
  if (get_data_scale(omicsData) == "abundance") {
    # Welcome to the pit of despair!
    stop("Data must be on the log scale.")
  }

  # Check the data scale. Data must be on one of the log scales.
  if (inherits(omicsData, "seqData") && get_data_scale(omicsData) != "counts") {
    # Welcome to the pit of despair!
    stop("seqData must be untransformed prior to running cor_result.")
  }

  edata_cname <- get_edata_cname(omicsData)
  id_col <- which(colnames(omicsData$e_data) == edata_cname)

  output <- cor(omicsData$e_data[, -id_col],
    use = "pairwise.complete.obs",
    method = ifelse(inherits(omicsData, "seqData"),
      "spearman", "pearson"
    )
  )

  rownames(output) <- names(omicsData$e_data[, -id_col])

  orig_class <- class(output)

  class(output) <- c("corRes", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data[, -id_col])
  attr(output, "group_DF") <- get_group_DF(omicsData)
  attr(output, "is_normalized") <- get_data_norm(omicsData)
  attr(output, "cor_method") <- ifelse(inherits(omicsData, "seqData"),
    "spearman", "pearson"
  )

  if (inherits(omicsData, "isobaricpepData")) {
    # attr(output, "isobaric_norm") <- attr(omicsData,"isobaric_info")$norm_info$is_normalized
    # attr(output, "is_normalized") <- attr(omicsData,"data_info")$norm_info$is_normalized
    attr(output, "isobaric_norm") <- get_isobaric_norm(omicsData)
  } else if (inherits(omicsData, "nmrData")) {
    # attr(output, "nmr_norm") <- attr(omicsData, "nmr_info")$norm_info$is_normalized
    # attr(output, "is_normalized") <- attr(omicsData,"data_info")$norm_info$is_normalized
    attr(output, "nmr_norm") <- get_nmr_norm(omicsData)
  }

  return(output)
}
