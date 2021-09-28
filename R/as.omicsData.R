#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of isobaric peptide-level data
#' to an object of the class 'isobaricpepData'. Objects of the class
#' 'isobaricpepData' are lists with two obligatory components \code{e_data} and
#' \code{f_data}. An optional list component \code{e_meta} is used if analysis
#' or visualization at other levels (e.g. protein) is also desired.
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
#' @param ... further arguments
#'
#' @details The class 'isobaricpepData' is meant to deal with peptide data
#'   generated on instruments where a reference pool for normalization is
#'   available (e.g. TMT, iTRAQ).
#'
#'  If your data has already undergone normalization to the reference pool, you
#'  should speficy \code{isobaric_norm = T}.
#'
#'  Objects of class 'isobaricpepData' contain some attributes that are
#'  referenced by downstream functions. These attributes can be changed from
#'  their default value by manual specification. A list of these attributes as
#'  well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'  Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'  'log10', 'log', and 'abundance', which indicate data is log base 2, base 10,
#'  natural log transformed, and raw abundance, respectively. Default is
#'  'abundance'. \cr \tab \cr is_normalized \tab A logical argument, specifying
#'  whether the data has been normalized or not. Default value is FALSE. \cr
#'  \tab \cr isobaric_norm \tab A logical argument, specifying whether the data
#'  has been normalized to the approporiate reference pool sample or not.
#'  Default value is FALSE \cr \tab \cr norm_info \tab Default value is an empty
#'  list, which will be populated with a single named element
#'  \code{is_normalized = is_normalized}. When a normalization is applied to the
#'  data, this becomes populated with a list containing the normalization
#'  function, normalization subset and subset parameters, the location and scale
#'  parameters used to normalize the data, and the location and scale parameters
#'  used to backtransform the data (if applicable). \cr \tab \cr data_types \tab
#'  Character string describing the type of data (e.g.'Positive ion'). Default
#'  value is NULL. \cr \tab \cr check.names \tab Logical defaults to TRUE.
#'  Indicates whether 'check.names' attribute of returned omicsData object is
#'  TRUE or FALSE. \cr } Computed values included in the \code{data_info}
#'  attribute are as follows: \tabular{ll}{ num_edata \tab The number of unique
#'  \code{edata_cname} entries.\cr \tab \cr num_miss_obs \tab The number of
#'  missing observations.\cr \tab \cr num_emeta \tab The number of unique
#'  \code{emeta_cname} entries. \cr \tab \cr prop_missing \tab The proportion of
#'  \code{e_data} values that are NA. \cr \tab \cr num_samps \tab The number of
#'  samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'  \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#'  \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("isobaric_edata")
#' data("isobaric_fdata")
#' data("isobaric_emeta")
#' mypepData <- as.isobaricpepData(e_data = isobaric_edata,
#'                                 e_meta = isobaric_emeta,
#'                                 f_data = isobaric_fdata,
#'                                 edata_cname = "Peptide",
#'                                 fdata_cname = "Sample",
#'                                 emeta_cname = "Protein")
#' }
#'
#' @author Lisa Bramer
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{normalize_isobaric}}
#'
#'@export
#'
as.isobaricpepData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                                fdata_cname, emeta_cname = NULL,
                                techrep_cname = NULL, ...) {
  .as.isobaricpepData(e_data, f_data, e_meta, edata_cname, fdata_cname,
                      emeta_cname, techrep_cname, ...)
}

## peptide data ##
.as.isobaricpepData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                                 fdata_cname, emeta_cname = NULL,
                                 techrep_cname = NULL,
                                 data_scale = "abundance",
                                 is_normalized = FALSE, isobaric_norm = FALSE,
                                 norm_info = list(),  data_types = NULL,
                                 check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('peptides', 'as.isobaricpepData')

  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

  # Set isobaric specific attribute.
  attr(res, "isobaric_info") <- set_isobaric_info(exp_cname = NA,
                                                  channel_cname = NA,
                                                  refpool_channel = NA,
                                                  refpool_cname = NA,
                                                  refpool_notation = NA,
                                                  norm_info = norm_info,
                                                  isobaric_norm = isobaric_norm)

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
  class(res) = c("isobaricpepData", "pepData")

  return(res)

}

#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of lipid-level data to an
#' object of the class 'lipidData'. Objects of the class 'lipidData' are lists
#' with two obligatory components \code{e_data} and \code{f_data}. An optional
#' list component \code{e_meta} is used if analysis or visualization at other
#' levels (e.g. lipid identification) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of lipids observed and \eqn{n} is the number of
#'   samples. Each row corresponds to data for each lipid. One column specifying
#'   a unique identifier for each lipid (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to a lipid with one column giving lipid names (must be named the same as
#'   the column in \code{e_data}) and other columns giving meta information.
#' @param edata_cname character string specifying the name of the column
#'   containing the lipid identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the mapped identifiers in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then either do not specify
#'   \code{emeta_cname} or specify it as NULL. If \code{e_meta} is NULL, then
#'   specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in
#'   \code{f_data} containing the identifiers for the biological samples if the
#'   observations represent technical replicates.  This column is used to
#'   collapse the data when \code{combine_techreps} is called on this object.
#'   Defaults to NULL (no technical replicates).
#'
#' @param ... further arguments
#'
#' @details Objects of class 'lipidData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'   'log10', 'log', and 'abundance', which indicate data is log base 2, base
#'   10, natural log transformed, and raw abundance, respectively. Default
#'   values is 'abundance'. \cr \tab \cr is_normalized \tab A logical argument,
#'   specifying whether the data has been normalized or not. Default value is
#'   FALSE. \cr \tab \cr norm_info \tab Default value is an empty list, which
#'   will be populated with a single named element \code{is_normalized =
#'   is_normalized}. When a normalization is applied to the data, this becomes
#'   populated with a list containing the normalization function, normalization
#'   subset and subset parameters, the location and scale parameters used to
#'   normalize the data, and the location and scale parameters used to
#'   backtransform the data (if applicable). \cr \tab \cr data_types \tab
#'   Character string describing the type of data (e.g. 'Positive ion'). Default
#'   value is NULL. \cr \tab \cr check.names \tab Logical defaults to TRUE.
#'   Indicates whether 'check.names' attribute of returned omicsData object is
#'   TRUE or FALSE. \cr } Computed values included in the \code{data_info}
#'   attribute are as follows: \tabular{ll}{ \code{num_edata} \tab The number of
#'   unique \code{edata_cname} entries.\cr \tab \cr \code{num_miss_obs} \tab The
#'   number of missing observations.\cr \tab \cr \code{num_emeta} \tab The
#'   number of unique \code{emeta_cname} entries. \cr \tab \cr
#'   \code{prop_missing} \tab The proportion of \code{e_data} values that are
#'   NA. \cr \tab \cr \code{num_samps} \tab The number of samples that make up
#'   the columns of \code{e_data}.\cr \tab \cr \code{meta_info} \tab A logical
#'   argument, specifying where the \code{e_meta} is provided.\cr \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("lipid_edata")
#' data("lipid_fdata")
#' mylipidData <- as.lipidData(e_data = lipid_edata,
#'                             f_data = lipid_fdata,
#'                             edata_cname = "LipidCommonName",
#'                             fdata_cname = "Sample_Name")
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#' @seealso \code{\link{as.metabData}}
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.proData}}
#'
#' @export
#'
as.lipidData <- function (e_data, f_data, e_meta = NULL,
                          edata_cname, fdata_cname, emeta_cname = NULL,
                          techrep_cname = NULL, ...) {
  .as.lipidData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
                techrep_cname, ...)
}

## lipid data ##
.as.lipidData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                           fdata_cname, emeta_cname = NULL,
                           techrep_cname = NULL,
                           data_scale = "abundance",
                           is_normalized = FALSE, norm_info = list(),
                           data_types = NULL, check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('lipids', 'as.lipidData')

  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

  #set check.names attribute #
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
  class(res) = "lipidData"

  return(res)

}

#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of metabolomic-level data to an
#' object of the class 'metabData'. Objects of the class 'metabData' are lists
#' with two obligatory components \code{e_data} and \code{f_data}. An optional
#' list component \code{e_meta} is used if analysis or visualization at other
#' levels (e.g. metabolite identification) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of metabolites observed and \eqn{n} is the number of
#'   samples. Each row corresponds to data for each metabolite. One column
#'   specifying a unique identifier for each metabolite (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to a metabolite with one column giving metabolite names (must be named the
#'   same as the column in \code{e_data}) and other columns giving meta
#'   information.
#' @param edata_cname character string specifying the name of the column
#'   containing the metabolite identifiers in \code{e_data} and \code{e_meta}
#'   (if applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the mapped identifiers in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then either do not specify
#'   \code{emeta_cname} or specify it as NULL. If \code{e_meta} is NULL, then
#'   specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in
#'   \code{f_data} containing the identifiers for the biological samples if the
#'   observations represent technical replicates.  This column is used to
#'   collapse the data when \code{combine_techreps} is called on this object.
#'   Defaults to NULL (no technical replicates).
#'
#' @param ... further arguments
#'
#' @details Objects of class 'metabData' contain some attributes that are
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
#'   \tab A logical argument, specifying where the \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("metab_edata")
#' data("metab_fdata")
#' mymetabData <- as.metabData(e_data = metab_edata,
#'                             f_data = metab_fdata,
#'                             edata_cname = "Metabolite",
#'                             fdata_cname = "SampleID")
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.proData}}
#'
#' @export
#'
as.metabData <- function (e_data, f_data, e_meta = NULL,
                          edata_cname, fdata_cname, emeta_cname = NULL,
                          techrep_cname = NULL, ...) {
  .as.metabData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
                techrep_cname, ...)
}

## metabolite data ##
.as.metabData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                           fdata_cname, emeta_cname = NULL,
                           techrep_cname = NULL,
                           data_scale = "abundance",
                           is_normalized = FALSE, norm_info = list(),
                           data_types = NULL, check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('metabolites', 'as.metabData')

  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

  #set check.names attribute #
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
  class(res) = "metabData"

  return(res)

}

#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of NMR-generated
#' metabolomic-level data to an object of the class 'nmrData'. Objects of the
#' class 'nmrData' are lists with two obligatory components \code{e_data} and
#' \code{f_data}. An optional list component \code{e_meta} is used if analysis
#' or visualization at other levels (e.g. metabolite identification) is also
#' desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of metabolites observed and \eqn{n} is the number of
#'   samples. Each row corresponds to data for each metabolite. One column
#'   specifying a unique identifier for each metabolite (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to a metabolite with one column giving metabolite names (must be named the
#'   same as the column in \code{e_data}) and other columns giving meta
#'   information.
#' @param edata_cname character string specifying the name of the column
#'   containing the metabolite identifiers in \code{e_data} and \code{e_meta}
#'   (if applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the mapped identifiers in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then either do not specify
#'   \code{emeta_cname} or specify it as NULL. If \code{e_meta} is NULL, then
#'   specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in
#'   \code{f_data} containing the identifiers for the biological samples if the
#'   observations represent technical replicates.  This column is used to
#'   collapse the data when \code{combine_techreps} is called on this object.
#'   Defaults to NULL (no technical replicates).
#'
#' @param ... further arguments
#'
#' @details Objects of class 'nmrData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'   'log10', 'log', and 'abundance', which indicate data is log base 2, base
#'   10, natural log transformed, and raw abundance, respectively. Default is
#'   'abundance'. \cr \tab \cr is_normalized \tab A logical argument, specifying
#'   whether the data has been normalized or not. Default value is FALSE. \cr
#'   \tab \cr nmr_norm \tab A logical argument, specifying whether the data has
#'   been normalized either to a spiked in metabolite or to a property taking
#'   sample-specific values \cr #' \tab \cr norm_info \tab Default value is an
#'   empty list, which will be populated with a single named element
#'   \code{is_normalized = is_normalized}. When a normalization is applied to
#'   the data, this becomes populated with a list containing the normalization
#'   function, normalization subset and subset parameters, the location and
#'   scale parameters used to normalize the data, and the location and scale
#'   parameters used to backtransform the data (if applicable). \cr \tab \cr
#'   data_types \tab Character string describing the type of data (e.g.'binned'
#'   or 'identified', for NMR data). Default value is NULL. \cr \tab \cr
#'   check.names \tab Logical defaults to TRUE. Indicates whether 'check.names'
#'   attribute of returned omicsData object is TRUE or FALSE. \cr } Computed
#'   values included in the \code{data_info} attribute are as follows:
#'   \tabular{ll}{ num_edata \tab The number of unique \code{edata_cname}
#'   entries.\cr \tab \cr num_miss_obs \tab The number of missing
#'   observations.\cr \tab \cr num_emeta \tab The number of unique
#'   \code{emeta_cname} entries. \cr \tab \cr prop_missing \tab The proportion
#'   of \code{e_data} values that are NA. \cr \tab \cr num_samps \tab The number
#'   of samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'   \tab A logical argument, specifying where the \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("nmr_edata_identified")
#' data("nmr_fdata_identified")
#' mynmrData <- as.nmrData(e_data = nmr_edata_identified,
#'                         f_data = nmr_fdata_identified,
#'                         edata_cname = "Metabolite",
#'                         fdata_cname = "SampleID",
#'                         data_type = "identified")
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#' @seealso \code{\link{as.metabData}}
#' @seealso \code{\link{as.isobaricpepData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.proData}}
#'
#' @export
#'
as.nmrData <- function (e_data, f_data, e_meta = NULL,
                        edata_cname, fdata_cname, emeta_cname = NULL,
                        techrep_cname = NULL, ...) {
  .as.nmrData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
              techrep_cname, ...)
}

## metabolite data ##
.as.nmrData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                         fdata_cname, emeta_cname = NULL,
                         techrep_cname = NULL,
                         data_scale = "abundance",
                         is_normalized = FALSE, nmr_norm = FALSE,
                         norm_info = list(), data_types = NULL,
                         check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('metabolites', 'as.nmrData')

  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

  # Set nmr specific attribute.
  attr(res, "nmr_info") <- set_nmr_info(metabolite_name = NA,
                                        sample_property_cname = NA,
                                        norm_info = norm_info,
                                        nmr_norm = nmr_norm,
                                        backtransform = NA)

  #set check.names attribute #
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
  class(res) = c("nmrData")

  return(res)

}

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
                        techrep_cname = NULL, ...) {
  .as.pepData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
              techrep_cname, ...)
}

## peptide data ##
.as.pepData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                         fdata_cname, emeta_cname = NULL,
                         techrep_cname = NULL,
                         data_scale = "abundance",
                         is_normalized = FALSE, norm_info = list(),
                         data_types = NULL, check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

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
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of protein-level data to an
#' object of the class 'proData'. Objects of the class 'proData' are lists with
#' two obligatory components \code{e_data} and \code{f_data}. An optional list
#' component \code{e_meta} is used if analysis or visualization at other levels
#' (e.g. gene) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of proteins observed and \eqn{n} is the number of
#'   samples. Each row corresponds to data for each protein. One column
#'   specifying a unique identifier for each protein (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to a protein with one column giving protein names (must be named the same
#'   as the column in \code{e_data}) and other columns giving meta information
#'   (e.g. mappings of proteins to genes).
#' @param edata_cname character string specifying the name of the column
#'   containing the protein identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the gene identifiers (or other mapping variable) in
#'   \code{e_meta} (if applicable). Defaults to NULL. If \code{e_meta} is NULL,
#'   then either do not specify \code{emeta_cname} or specify it as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in
#'   \code{f_data} containing the identifiers for the biological samples if the
#'   observations represent technical replicates.  This column is used to
#'   collapse the data when \code{combine_techreps} is called on this object.
#'   Defaults to NULL (no technical replicates).
#'
#' @param ... further arguments
#'
#' @details Objects of class 'proData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'   'log10', 'log', and 'abundance', which indicate data is log base 2, base
#'   10, natural log transformed, and raw abundance, respectively. Default
#'   values is 'abundance'. \cr \tab \cr is_normalized \tab A logical argument,
#'   specifying whether the data has been normalized or not. Default value is
#'   FALSE. \cr \tab \cr norm_info \tab Default value is an empty list, which
#'   will be populated with a single named element \code{is_normalized =
#'   is_normalized}. When a normalization is applied to the data, this becomes
#'   populated with a list containing the normalization function, normalization
#'   subset and subset parameters, the location and scale parameters used to
#'   normalize the data, and the location and scale parameters used to
#'   backtransform the data (if applicable). \cr \tab \cr data_types \tab
#'   Character string describing the type of data (e.g.'Positive ion'). Default
#'   value is NULL. \cr \tab\cr check.names \tab Logical defaults to TRUE.
#'   Indicates whether 'check.names' attribute of returned omicsData object is
#'   TRUE or FALSE. \cr } Computed values included in the \code{data_info}
#'   attribute are as follows: \tabular{ll}{ num_edata \tab The number of unique
#'   \code{edata_cname} entries.\cr \tab \cr num_miss_obs \tab The number of
#'   missing observations.\cr \tab \cr num_emeta \tab The number of unique
#'   \code{emeta_cname} entries. \cr \tab \cr prop_missing \tab The proportion
#'   of \code{e_data} values that are NA. \cr \tab \cr num_samps \tab The number
#'   of samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'   \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pro_edata")
#' data("pro_fdata")
#' myproData <- as.proData(e_data = pro_edata,
#'                         f_data = pro_fdata,
#'                         edata_cname = "Reference",
#'                         fdata_cname = "SampleID",
#'                         is_normalized = TRUE)
#' }
#'
#' @author Kelly Stratton, Lisa Bramer
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.metabData}}
#'
#' @export
#'
as.proData <- function (e_data, f_data, e_meta = NULL,
                        edata_cname, fdata_cname, emeta_cname = NULL,
                        techrep_cname = NULL, ...) {
  .as.proData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
              techrep_cname, ...)
}

## protein data ##
.as.proData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                         fdata_cname, emeta_cname = NULL,
                         techrep_cname = NULL,
                         data_scale = "abundance",
                         is_normalized = FALSE, norm_info = list(),
                         data_types = NULL, check.names = TRUE) {

  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale

  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('proteins', 'as.proData')

  # Perform pre analysis checks. Return updated data frames if they all receive
  # a gold star.
  res <- pre_flight(e_data = e_data,
                    f_data = f_data,
                    e_meta = e_meta,
                    edata_cname = edata_cname,
                    fdata_cname = fdata_cname,
                    emeta_cname = emeta_cname,
                    techrep_cname = techrep_cname,
                    data_scale = data_scale,
                    is_normalized = is_normalized,
                    norm_info = norm_info,
                    data_types = data_types,
                    check.names = check.names,
                    dType = dType)

  # Set the (possibly new) emeta_cname.
  emeta_cname <- res$emeta_cname

  # Remove the emeta_cname element from the res list. That way only the e_data,
  # f_data, and e_meta (when applicable) data frames will be part of the output.
  # emeta_cname will no longer be an element of and omicsData object.
  res <- res %>%
    purrr::list_modify("emeta_cname" = NULL)

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

  #set check.names attribute #
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

  # Set the protein quantitation attribute to NA. This will be updated to one of
  # rollup, qrollup, rrollup, or zrollup if/when a pepData object is rolled up
  # to a proData object.
  attr(res, "pro_quant_info") <- list(method = NA)

  # set class of list #
  class(res) = "proData"

  return(res)

}
