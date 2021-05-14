#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of peptide-level data to an
#' object of the class 'pepData'. Objects of the class 'pepData' are lists with
#' two obligatory components \code{e_data} and \code{f_data}. An optional list
#' component \code{e_meta} is used if analysis or visualization at other levels
#' (e.g. protein) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of peptides observed and \eqn{n} is the number of
#'   samples (an additional peptide identifier/name column should also be
#'   present anywhere in the data.frame). Each row corresponds to data for each
#'   peptide. One column specifying a unique identifier for each peptide (row)
#'   must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with at least \eqn{p} rows. Each row
#'   corresponds to a peptide with one column giving peptide names (must be
#'   named the same as the column in \code{e_data}) and other columns giving
#'   meta information (e.g. mappings of peptides to proteins).
#' @param edata_cname character string specifying the name of the column
#'   containing the peptide identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the protein identifiers (or other mapping variable) in
#'   \code{e_meta} (if applicable). Defaults to NULL. If \code{e_meta} is NULL,
#'   then either do not specify \code{emeta_cname} or specify it as NULL. If
#'   \code{e_meta} is NULL, then specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in
#'   \code{f_data} containing the identifiers for the biological samples if the
#'   observations represent technical replicates.  This column is used to
#'   collapse the data when \code{combine_techreps} is called on this object.
#'   Defaults to NULL (no technical replicates).
#'
#' @param data_scale_orig A character string indicating what scale the data are
#'   on. Acceptable values are "abundance", "log", "log2", and "log10".
#'
#' @param ... further arguments
#'
#' @details Objects of class 'pepData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'   'log10', 'log', and 'abundance', which indicate data is log base 2, base
#'   10, natural log transformed, and raw abundance, respectively. Default is
#'   'abundance'. \cr \tab \cr is_normalized \tab A logical argument, specifying
#'   whether the data has been normalized or not. Default value is FALSE. \cr
#'   \tab \cr norm_info \tab Default value is an empty list, which will be
#'   populated with a single named element \code{is_normalized = is_normalized}.
#'   When a normalization is applied to the data, this becomes populated with a
#'   list containing the normalization function, normalization subset and subset
#'   parameters, the location and scale parameters used to normalize the data,
#'   and the location and scale parameters used to backtransform the data (if
#'   applicable). \cr \tab \cr data_types \tab Character string describing the
#'   type of data (e.g.'Positive ion'). Default value is NULL. \cr \tab \cr
#'   check.names \tab Logical defaults to TRUE. Indicates whether 'check.names'
#'   attribute of returned omicsData object is TRUE or FALSE. \cr } Computed
#'   values included in the \code{data_info} attribute are as follows:
#'   \tabular{ll}{ num_edata \tab The number of unique \code{edata_cname}
#'   entries.\cr \tab \cr num_miss_obs \tab The number of missing
#'   observations.\cr \tab \cr num_emeta \tab The number of unique
#'   \code{emeta_cname} entries. \cr \tab \cr prop_missing \tab The proportion
#'   of \code{e_data} values that are NA. \cr \tab \cr num_samps \tab The number
#'   of samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'   \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_edata")
#' data("pep_fdata")
#' data("pep_emeta")
#' mypepData <- as.pepData(e_data = pep_edata,
#'                         e_meta = pep_emeta,
#'                         f_data = pep_fdata,
#'                         edata_cname = "Mass_Tag_ID",
#'                         fdata_cname = "SampleID",
#'                         emeta_cname = "Mass_Tag_ID")
#' }
#' 
#' @author Kelly Stratton, Lisa Bramer
#' @seealso \code{\link{as.proData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.metabData}}
#'
#' @export
#' 
as.pepData <- function (e_data, f_data, e_meta = NULL,
                        edata_cname, fdata_cname, emeta_cname = NULL,
                        techrep_cname = NULL, data_scale_orig = NULL, ...) {
  .as.pepData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
              techrep_cname, data_scale_orig, ...)
}

## peptide data ##
.as.pepData <- function (e_data, f_data, e_meta = NULL, edata_cname, 
                         fdata_cname, emeta_cname = NULL,
                         techrep_cname = NULL,
                         data_scale_orig = NULL,
                         data_scale = data_scale_orig,
                         is_normalized = FALSE, norm_info = list(),
                         data_types = NULL, check.names = TRUE) {
  
  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('peptides', 'as.pepData')
  
  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale_orig = data_scale_orig,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)
  
  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname
  
  # set column name attributes #
  attr(res, "cnames") = list(edata_cname = edata_cname,
                             emeta_cname = emeta_cname,
                             fdata_cname = fdata_cname,
                             techrep_cname = techrep_cname)
  
  # Compute the data_info attributes.
  attr(res, "data_info") <- set_data_info(e_data = res$e_data,
                                          edata_cname = edata_cname,
                                          data_scale_orig = data_scale_orig,
                                          data_scale = data_scale,
                                          data_types = data_types,
                                          norm_info = norm_info,
                                          is_normalized = is_normalized)
  
  # set check.names attribute #
  attr(res, "check.names") = check.names 
  
  # set meta data attributes #
  attr(res, "meta_info") <- set_meta_info(e_meta = res$e_meta,
                                          emeta_cname = emeta_cname)
  
  # set group dataframe attribute to NULL, will be filled in after running
  # group_designation function #
  attr(res, "group_DF") = NULL
  
  # Initialize the filters attribute with a list. This list will be populated
  # with filter class objects as filters are applied. Multiple filters can be
  # implemented on one data set.
  attr(res, "filters") <- list()
  
  # set class of list #
  class(res) = "pepData"
  
  return(res)
  
}
