#' Creates an object of class naRes (NA Result)
#'
#' This function takes in an omicsData object, and outputs a list of two data
#' frames, one containing the number of missing values by sample, and the other
#' containing the number of missing values by molecule
#'
#' @param omicsData an object of class "pepData", "proData", "metabData",
#'   "lipidData", "nmrData", or "seqData", created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, or
#'   \code{\link{as.seqData}}, respectively.
#'
#' @return S3 object of class naRes, which is a list of two data frames, one
#'   containing the number of missing values per sample, and the other
#'   containing the number of missing values per molecule. For count data,
#'   zeroes represent missing values; for abundance data, NA's represent missing
#'   values. This object can be used with 'plot' and 'summary' methods to
#'   examine the missing values in the dataset.
#'
#' @examples
#' library(pmartRdata)
#' result1 = missingval_result(omicsData = lipid_neg_object)
#' result2 = missingval_result(omicsData = metab_object)
#'
#' @export
#'
missingval_result <- function(omicsData) {
  # Add a check

  # check for correct class
  if (!inherits(omicsData, c(
    "pepData", "proData", "lipidData",
    "metabData", "nmrData", "seqData"
  ))) {
    stop(paste("omicsData must have class of the following, 'pepData',",
      "'proData', 'lipidData', 'metabData', 'nmrData', 'seqData'",
      sep = " "
    ))
  }


  # pulling cname attr
  # edata_cname<- attr(omicsData, "cnames")$edata_cname
  edata_cname <- get_edata_cname(omicsData)
  edata_cname_id <- which(names(omicsData$e_data) == edata_cname)
  # fdata_cname<- attr(omicsData, "cnames")$fdata_cname
  fdata_cname <- get_fdata_cname(omicsData)

  # Count the number of NA or zeros values per column.
  if(inherits(omicsData, "seqData")){
    res_per_col<- colSums((omicsData$e_data[, -edata_cname_id]) == 0)
    res_per_col_non <- colSums((omicsData$e_data[, -edata_cname_id]) != 0)
    res_by_sample<- data.frame(
      "sample_names" = names(omicsData$e_data[, -edata_cname_id]),
      "num_zeros" = as.numeric(res_per_col),
      "num_nonzeros" = as.numeric(res_per_col_non)
    )
  } else {
    res_per_col<- colSums(is.na(omicsData$e_data[, -edata_cname_id]))
    res_per_col_non <- colSums(!is.na(omicsData$e_data[, -edata_cname_id]))
    res_by_sample<- data.frame(
      "sample_names" = names(omicsData$e_data[, -edata_cname_id]),
      "num_NA" = as.numeric(res_per_col),
      "num_non_NA" = as.numeric(res_per_col_non)
    )
  }

  names(res_by_sample)[1] <- fdata_cname

  # Merge res_by_sample with f_data to get additional columns of f_data. For
  # example, the Group and VizSampNames columns. Group is used to color the plot
  # and VizSampNames is used to display shorter sample names.
  res_by_sample <- merge(res_by_sample, omicsData$f_data, by = fdata_cname)

  # Check if the group designation function has been run. The group_DF info will
  # be used to add the "Group" column to res_by_sample. This column may contain
  # the same data as another column in f_data but it will have a different name
  # from the f_data column.
  if (!is.null(attr(omicsData, "group_DF"))) {
    res_by_sample <- merge(res_by_sample, attr(omicsData, "group_DF"))
  }

  # Count the number of NA values per row.

  if (inherits(omicsData, "seqData")) {
    res_per_row <- rowSums(omicsData$e_data[, -edata_cname_id] == 0)
    res_per_row_non <- rowSums(omicsData$e_data[, -edata_cname_id] != 0)
    
    res_by_molecule <- data.frame("molecule"= omicsData$e_data[, edata_cname_id],
                                 "num_zeros"= as.numeric(res_per_row),
                                 "num_nonzeros" = as.numeric(res_per_row_non))
    names(res_by_molecule)[1] <- edata_cname

    result <- list(
      "zeros.by.sample" = res_by_sample,
      "zeros.by.molecule" = res_by_molecule
    )
    class(result) <- "naRes"

    attr(result, "cnames") <- list(
      "edata_cname" = edata_cname,
      "fdata_cname" = fdata_cname
    )
  } else {
    res_per_row <- rowSums(is.na(omicsData$e_data[, -edata_cname_id]))
    res_per_row_non <- rowSums(!is.na(omicsData$e_data[, -edata_cname_id]))
    res_by_molecule <- data.frame("molecule"= omicsData$e_data[, edata_cname_id],
                                 "num_NA"= as.numeric(res_per_row),
                                 "num_non_NA" = as.numeric(res_per_row_non))
    names(res_by_molecule)[1] <- edata_cname

    result <- list(
      "na.by.sample" = res_by_sample,
      "na.by.molecule" = res_by_molecule
    )
    class(result) <- "naRes"

    attr(result, "cnames") <- list(
      "edata_cname" = edata_cname,
      "fdata_cname" = fdata_cname
    )
  }

  return(result)
}
