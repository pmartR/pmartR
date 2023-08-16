#' Convert Data of class MSnSet to pmartR pepData Class
#'
#' Converts an object of class MSnSet to an object of the class 'pepData'
#'
#' @param msnset_object an object of class MSnSet, which stores quantification
#'   data and meta data. Creating an MSnSet object is described in the MSNbase
#'   package io vignette.
#' @param edata_cname character string specifying the name of the column
#'   containing the peptide identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the protein identifiers (or other mapping variable) in
#'   \code{e_meta} (if applicable).
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param data_scale Scale of the data provided in \code{e_data}. Acceptable
#'   values are 'log2', 'log10', 'log', and 'abundance', which indicate data is
#'   log base 2, base 10, natural log transformed, and raw abundance,
#'   respectively.
#' @param check.names deprecated
#'
#' @return pepData object
#'
#' @details The MSnbase package is available via Bioconductor
#'
#' @references Gatto L, Lilley K (2012). “MSnbase - an R/Bioconductor package
#'   for isobaric tagged mass spectrometry data visualization, processing and
#'   quantitation.” Bioinformatics, 28, 288-289.
#'
#' @references Gatto L, Gibb S, Rainer J (2020). “MSnbase, efficient and elegant
#'   R-based processing and visualisation of raw mass spectrometry data.”
#'   bioRxiv.
#'
#' @examples
#' \dontrun{
#' library(MSnbase)
#' data("msnset")
#' result = MSnSet2pepData(msnset,
#'   data_scale = "log2",
#'   edata_cname = "UniqueID",
#'   fdata_cname = "SampleID",
#'   emeta_cname = "UniqueID"
#' )
#' }
#'
#' @export
#'
MSnSet2pepData <- function(msnset_object, data_scale, edata_cname = "UniqueID",
                           fdata_cname = "SampleID", emeta_cname = "UniqueID",
                           check.names = NULL) {
  if (!missing(check.names))
    warning("check.names parameter is deprecated")

  # check that msnset_object is of correct class
  if (!inherits(msnset_object, "MSnSet"))
    stop("msnset_object must be of class 'MSnSet'")

  # check that data_scale is one of the acceptable options #
  if (!(data_scale %in% c('log2', 'log10', 'log', 'count', 'abundance')))
    stop(paste(data_scale, " is not a valid option for 'data_scale'", sep = ""))

  msnset_edata <- msnset_object@assayData$exprs
  if (any(dim(msnset_edata) == 0))
    stop("msnset_object@assayData must not have empty rows or columns ")
  msnset_edata <- as.data.frame(msnset_edata)
  msnset_edata <- cbind(row.names(msnset_edata), msnset_edata)
  row.names(msnset_edata) <- NULL
  names(msnset_edata)[1] <- "UniqueID"

  msnset_fdata <- msnset_object@phenoData@data
  if (any(dim(msnset_fdata) == 0))
    stop("msnset_object@phenoData must not have empty rows or columns ")
  msnset_fdata <- as.data.frame(msnset_fdata)
  msnset_fdata <- cbind(row.names(msnset_fdata), msnset_fdata)
  row.names(msnset_fdata) <- NULL
  names(msnset_fdata)[1] <- "SampleID"

  msnset_emeta <- msnset_object@featureData@data
  if (any(dim(msnset_emeta) == 0))
    stop("msnset_object@featureData must not have empty rows or columns ")
  msnset_emeta <- as.data.frame(msnset_emeta)
  msnset_emeta <- cbind(row.names(msnset_emeta), msnset_emeta)
  row.names(msnset_emeta) <- NULL
  names(msnset_emeta)[1] <- "UniqueID"

  res <- as.pepData(
    e_data = msnset_edata, f_data = msnset_fdata,
    e_meta = msnset_emeta, edata_cname = edata_cname,
    fdata_cname = fdata_cname, emeta_cname = emeta_cname,
    data_scale = data_scale
  )

  return(res)
}
