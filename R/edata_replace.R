#' Replace Values Equal to x with y
#'
#' This function finds all values of x in the e_data element of omicsData and
#' replaces them with y
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData' created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @param x value to be replaced, usually numeric or NA
#'
#' @param y replacement value, usually numeric or NA
#'
#' @param threshold Positive numeric value.  Observed values below this
#'  threshold will be replaced by `y` (in addition to all `x` values).
#' 
#' @details This function is often used to replace any 0 values in peptide,
#'   protein, metabolite, or lipid data with NA's. For omicsData on the
#'   abundance scale, when the omicsData object is created, any 0's in e_data
#'   are automatically converted to NA's. For omicsData on the count scale (e.g.
#'   seqData objects), when the omicsData object is created, any NA's in e_data
#'   are automatically converted to 0's.
#'
#' @return data object of the same class as omicsData
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mymetab <- edata_replace(omicsData = metab_object, x = 0, y = NA)
#'
#' @author Kelly Stratton
#'
#' @export
#'
edata_replace <- function(omicsData, x, y, threshold = NULL) {
  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c(
    "pepData", "proData", "metabData",
    "lipidData", "nmrData" # , "seqData"
  ))) {
    # Bestow an error on the user because the input does not have the correct
    # class.
    stop(paste("omicsData must be of class 'pepData', 'proData',",
      "'metabData', 'lipidData', or 'nmrData'. ",
      sep = " "
    ))
  }

  # Acquire the index of the edata_cname column.
  id_col <- which(names(omicsData$e_data) == attr(
    omicsData,
    "cnames"
  )$edata_cname)

  # get count of the number of values replaced #
  if (is.na(x)) {
    num_replaced <- sum(is.na(omicsData$e_data[, -id_col]),
      na.rm = TRUE
    )

    replace_inds <- which(is.na(omicsData$e_data[, -id_col]),
      arr.ind = TRUE
    )
  } else {
    num_replaced <- sum(omicsData$e_data[, -id_col] == x,
      na.rm = TRUE
    )

    replace_inds <- which(omicsData$e_data[, -id_col] == x,
      arr.ind = TRUE
    )
  }

  # also count the number of values below the threshold and get indices to 
  # replace
  if (!is.null(threshold)) {
    cond_thresh <- omicsData$e_data[, -id_col] < threshold
    cond_match <- if(is.na(x)) {
      !is.na(omicsData$e_data[, -id_col])
    } else {
      !omicsData$e_data[, -id_col] == x
    }

    cond <- cond_thresh & cond_match

    num_thresholded <- sum(cond,
      na.rm = TRUE
    )

    threshold_inds = which(cond,
      arr.ind = TRUE
    )
  }

  # Update the e_data data frame in the omicsData object.
  omicsData$e_data[, -id_col][replace_inds] <- y

  # replace values below the threshold with y
  if (!is.null(threshold)) {
    omicsData$e_data[, -id_col][threshold_inds] <- y
  }

  # Update the data_info attribute of the omicsData object.
  attr(omicsData, 'data_info') <- set_data_info(
    e_data = omicsData$e_data,
    edata_cname = get_edata_cname(omicsData),
    data_scale_orig = get_data_scale_orig(omicsData),
    data_scale = get_data_scale(omicsData),
    data_types = get_data_info(omicsData)$data_types,
    norm_info = get_data_info(omicsData)$norm_info,
    is_normalized = get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = get_data_info(omicsData)$batch_info,
    is_bc = get_data_info(omicsData)$batch_info$is_bc
  )

  # Report the number of replaced elements in e_data
  message(paste(num_replaced,
    "instances of",
    x,
    "have been replaced with",
    y,
    sep = " "
  ))

  # Report the number of values below the threshold that were replaced
  if (!is.null(threshold)) {
    message(paste(num_thresholded,
      "values below the threshold",
      threshold,
      "have been replaced with",
      y,
      sep = " "
    ))
  }

  # Return the updated omicsData object.
  return(omicsData)
}

#' Replace x with y for a single vector
#'
#' @param one_vector numeric vector
#' @param x value to be replaced
#' @param y replacement value
#'
#' @return numeric vector
#'
#' @author Kelly Stratton
#'
vector_replace <- function(one_vector, x, y) {
  # find indices where the value is x #
  if (is.na(x)) {
    inds <- is.na(one_vector)
  } else {
    inds <- which(one_vector == x)
  }

  # Check if any values in the input vector match the value to be replaced. If
  # zero values in the input vector match x then the which function will return
  # integer(0) or the sum of is.na(x) will be zero.
  if (length(inds) == 0 || sum(inds) == 0) {
    # Return the one_vector unchanged.
    return(one_vector)
  } else {
    # Replace x with y.
    one_vector[inds] <- y

    # Return one_vector with the updated values.
    return(one_vector)
  }
}
