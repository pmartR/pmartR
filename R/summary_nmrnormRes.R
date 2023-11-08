#' Summary of nmrnormRes Object
#'
#' For creating a summary of an S3 object of type 'nmrnormRes'
#'
#' @param object object of type nmrnormRes, created by
#'   \code{\link{normalize_nmr}}
#' @param ... further arguments passed to or from other methods.
#'
#' @return data frame object
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mynmr <- edata_transform(
#'   omicsData = nmr_identified_object,
#'   data_scale = "log2"
#' )
#' nmr_norm <- normalize_nmr(
#'   omicsData = mynmr, apply_norm = FALSE,
#'   sample_property_cname = "Concentration"
#' )
#' mysummary <- summary(nmr_norm)
#'
#' @export
#' @rdname summary-nmrnormRes
#' @name summary-nmrnormRes
#'
summary.nmrnormRes <- function(object, ...) {
  nmrnormRes_object <- object

  # check for an nmrnormRes object #
  if (!inherits(nmrnormRes_object, "nmrnormRes")) stop("object must be of class 'nmrnormRes'")

  # extract attributes from nmrnomrRes_object
  sample_property_cname <- attr(nmrnormRes_object, "nmr_info")$sample_property_cname
  metabolite_name <- attr(nmrnormRes_object, "nmr_info")$metabolite_name
  if (is.null(sample_property_cname)) {
    normalized_using <- metabolite_name
  } else {
    if (is.null(metabolite_name)) {
      normalized_using <- sample_property_cname
    }
  }
  # value_col_ind = which(names(nmrnormRes_object) == "value")

  res_median = median(nmrnormRes_object$value, na.rm = TRUE)
  res_sd = sd(nmrnormRes_object$value, na.rm = TRUE)

  final_res = as.data.frame(cbind(normalized_using, res_median, res_sd))
  row.names(final_res) = NULL

  if (is.null(sample_property_cname)) {
    # normalized_using <- metabolite_name
    names(final_res) = c("Metabolite", "Median", "SD")
  } else {
    if (is.null(metabolite_name)) {
      # normalized_using <- sample_property_name
      names(final_res) = c("Sample Property", "Median", "SD")
    }
  }


  return(final_res)
}
