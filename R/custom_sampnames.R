#' Creates custom sample names to be used in plots
#'
#' This helper function creates custom sample names for plot data object
#' functions
#'
#' @param omicsData an object of the class 'pepData', 'proData',
#'   'metabData', 'lipidData', 'nmrData', or 'seqData', created by \code{as.pepData},
#'   \code{as.proData}, \code{as.metabData}, \code{as.lipidData}, \code{as.nmrData},
#'   or \code{as.seqData}, respectively.
#' @param firstn an integer specifying the first n characters to keep as the
#'   sample name. This argument is optional.
#' @param from an integer specifying the start of the range of characters to
#'   keep as the sample name. This argument is optional. If this argument is
#'   specified, 'to' must also be specified.
#' @param to an integer specifying the end of the range of characters to keep as
#'   the sample name. This argument is optional. If this argument is specified,
#'   'from' must also be specified.
#' @param delim character delimiter by which to separate sample name components.
#'   This argument is optional. If this argument is specified, 'components' must
#'   also be specified.
#' @param components integer vector specifying which components separated by the
#'   delimiter should be kept as the custom sample name. This argument is
#'   optional. If this argument is specified, 'delim' must also be specified.
#' @param pattern character string specifying the regex pattern to use to
#'   extract substrings from the sample names
#' @param ... extra arguments passed to regexpr if pattern is specified
#'
#' @return Object of same class as omicsData, with added column in f_data named
#'   'VizSampNames'.
#'
#' @details This function can be used to create custom (and shorter) sample
#'   names to be used when plotting so that axis labels are not so long that
#'   they interfere with the graph itself. To use the custom sample names when
#'   plotting, specify the optional argument 'use_VizSampNames = TRUE'.
#'
#' @examples
#' library(pmartRdata)
#'
#' mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' plot(mypep)
#'
#' # specify new names using firstn argument
#' results <- custom_sampnames(omicsData = mypep, firstn = 9)
#' plot(results, use_VizSampNames = TRUE)
#'
#' # specify new names using from and to arguments
#' results <- custom_sampnames(omicsData = mypep, from = 1, to = 9)
#' plot(results, use_VizSampNames = TRUE)
#'
#' # specify new names using delim and components arguments
#' results <- custom_sampnames(omicsData = mypep, delim = "_", components = c(1, 2))
#' plot(results, use_VizSampNames = TRUE)
#'
#' ## specify new names using pattern arguments (regex)
#'
#' # match everything after "Sample_"
#' pattern1 <- "[0-9]+_[0-9A-Za-z]+_[A-Z]"
#'
#' results <- custom_sampnames(omicsData = mypep, pattern = pattern1)
#' plot(results, use_VizSampNames = TRUE)
#'
#' # match "Sample_" and the number after it
#' pattern2 <- "^Sample_[0-9]+"
#'
#' results <- custom_sampnames(omicsData = mypep, pattern = pattern2)
#' plot(results, use_VizSampNames = TRUE)
#'
#' @export
#'
custom_sampnames <- function(omicsData, firstn = NULL, from = NULL,
                             to = NULL, delim = NULL, components = NULL,
                             pattern = NULL, ...) {
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c(
    "pepData", "proData", "metabData",
    "isobaricpepData", "lipidData", "nmrData",
    "seqData"
  ))) {
    # Throw an error that the input for omicsData is not the appropriate class.
    stop(paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
      "'isobaricpepData', 'lipidData', 'nmrData', or 'seqData'",
      sep = ' '
    ))
  }

  if (sum(c(is.null(from), is.null(to))) == 1) {
    stop("Please specify both 'from' and 'to' if either argument is used.")
  }

  if (sum(c(is.null(delim), is.null(components))) == 1) {
    stop("Please specify both 'delim' and 'components' if either argument is used.")
  }

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c(
    "pepData", "proData", "metabData",
    "lipidData", "nmrData", "seqData"
  ))) {
    # Bestow an error on the user because the input does not have the correct
    # class.
    stop(paste("omicsData must be of class 'pepData', 'proData',",
      "'metabData', 'lipidData', 'nmrData', or 'seqData'. ",
      sep = " "
    ))
  }

  # extract sample names from omicsData object
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  names = as.list(as.character(omicsData$f_data[[fdata_cname]]))

  output = subset_names(names,
    firstn = firstn, from = from, to = to,
    delim = delim, components = components,
    pattern = pattern, ...
  )

  omicsData$f_data[["VizSampNames"]] = output

  return(omicsData)
}

subset_names <- function(x, firstn, from, to, delim, components, pattern, ...) {
  if (!is.null(firstn)) {
    if (!is.null(from) | !is.null(to) | !is.null(delim) | !is.null(components) | !is.null(pattern)) stop("only firstn argument is needed")
    if (!is.numeric(firstn)) stop("firstn must be a non-negative numeric value less than the number of characters in the smallest sample name")

    output = lapply(x, function(x) {
      temp = strsplit(x, split = "")[[1]]
      if (firstn > length(temp)) warning(paste("there are less than", firstn, "characters in a sample name, nothing was truncated", sep = " "))
      temp <- temp[1:firstn]
      temp2 = paste(temp[!is.na(temp)], collapse = "")
      return(temp2)
    })
    output = unlist(output)
  } else if ((!is.null(from)) & (!is.null(to))) {
    if (!is.null(firstn) | !is.null(delim) | !is.null(components) | !is.null(pattern)) stop("only from and to arguments are needed")
    if (!is.numeric(from) | !is.numeric(to)) stop("'from' and 'to' must be non-negative numeric values")
    if (to < from) stop("'to' must be less than 'from'")

    output = lapply(x, function(x) {
      temp = strsplit(x, split = "")
      if (!(from %in% 1:length(temp[[1]])) & !(to %in% 1:length(temp[[1]]))) stop(paste(from, to, "are not in the range of the length of at least one sample name", sep = " "))
      temp2 = paste(temp[[1]][from:min(to, length(temp[[1]]))], collapse = "")
      return(temp2)
    })

    output = unlist(output)
  } else if ((!is.null(delim)) & (!is.null(components))) {
    if (!is.null(firstn) | !is.null(from) | !is.null(to) | !is.null(pattern)) stop("only delim and components arguments are needed")

    output = lapply(x, function(x) {
      temp = strsplit(x, split = delim)[[1]]
      if (length(components) > length(temp)) stop(paste("the length of 'components vector must be less than", length(temp), sep = " "))
      if (!any(components %in% 1:length(temp))) stop("none of the indices specified in 'components' match indices of the split sample name")
      temp <- temp[components]
      temp2 = paste(temp[!is.na(temp)], collapse = delim)
      return(temp2)
    })

    output = unlist(output)
  } else if (!is.null(pattern)) {
    matches = regexpr(pattern, x, ...)

    output = character(0)

    for (i in 1:length(matches)) {
      tmp_chr = substr(x[i], matches[i], matches[i] + attributes(matches)$match.length[i] - 1)
      output = c(output, tmp_chr)
    }
  }

  # Check if all short sample names are unique.
  if (length(unique(output)) != length(x)) {
    display_df = data.frame("Names" = x, "Extracted" = output)

    # Holy non-unique sample names, Batman!
    stop(
      "The input used does not produce a unique sample name for each sample.\n",
      paste(capture.output(print(display_df)), collapse = "\n")
    )
  }

  return(output)
}
