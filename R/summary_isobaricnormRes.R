#' Summary for isobaricnormRes Object
#'
#' For creating a summary of an S3 object of type 'isobaricnormRes'
#'
#' @param object object of type isobaricnormRes, created by
#'   \code{\link{normalize_isobaric}}
#' @param ... further arguments passed to or from other methods.
#'
#' @return data frame object
#'
#' @examples
#' library(pmartRdata)
#' myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
#' myiso_norm <- normalize_isobaric(
#'   omicsData = myiso, exp_cname = "Plex",
#'   apply_norm = FALSE,
#'   refpool_cname = "Virus",
#'   refpool_notation = "Pool"
#' )
#' mysummary <- summary(myiso_norm)
#'
#' @export
#' @rdname summary-isobaricnormRes
#' @name summary-pmartR
#'
summary.isobaricnormRes <- function(object, ...) {
  isobaricnormRes_object <- object

  # check for an isobaricnormRes object #
  if (!inherits(isobaricnormRes_object, "isobaricnormRes")) {
    stop("object must be of class 'isobaricnormRes'")
  }

  # extract attributes from isobaricnomrRes_object
  exp_cname <- attr(isobaricnormRes_object, "isobaric_info")$exp_cname
  fdata_cname <- attr(isobaricnormRes_object, "cnames")$fdata_cname

  # Fish out all unique experiments.
  xprmnts <- unique(isobaricnormRes_object$f_data[, exp_cname])

  # Assemble a data frame with the experiment levels/values.
  dfxprmnts <- data.frame(xprmnts)
  names(dfxprmnts) <- exp_cname

  # Create a list that will hold the samples belonging to each experiment.
  smpls <- vector(
    mode = "list",
    length = length(xprmnts)
  )

  # Generate a vector for the median for the samples in each experiment.
  medi <- vector(
    mode = "numeric",
    length = length(xprmnts)
  )

  # Produce a vector for the standard deviation of the samples in each
  # experiment.
  stdev <- vector(
    mode = "numeric",
    length = length(xprmnts)
  )

  # Loop through each level of experiment and extract all sample names
  # corresponding to each experiment.
  for (e in 1:length(xprmnts)) {
    # Grab the row indices of f_data for the eth experiment.
    idx <- which(isobaricnormRes_object$f_data[, exp_cname] == xprmnts[[e]])

    # Seize the sample names corresponding to the eth experiment.
    smpls[[e]] <- isobaricnormRes_object$f_data[idx, fdata_cname]

    # Compute the median for the samples in the eth experiment.
    medi[[e]] <- median(isobaricnormRes_object$e_data[, smpls[[e]]],
      na.rm = TRUE
    )

    # Calculate the standard deviation for the samples in the eth experiment.
    stdev[[e]] <- sd(isobaricnormRes_object$e_data[, smpls[[e]]],
      na.rm = TRUE
    )
  }

  # Unite the median and standard deviation vectors in a data frame.
  return(data.frame(dfxprmnts,
    Median = medi,
    SD = stdev
  ))
}
