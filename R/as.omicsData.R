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
                                 is_bc = FALSE, batch_info = list(),
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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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
                           is_bc = FALSE, batch_info = list(),
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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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
                           is_bc = FALSE, batch_info = list(),
                           data_types = NULL, check.names = TRUE
                           ) {

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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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
                         is_bc = FALSE, batch_info = list(),
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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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
                         is_bc = FALSE, batch_info = list(),
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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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
                         is_bc = FALSE, batch_info = list(),
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
                    dType = dType,
                    is_bc = is_bc,
                    batch_info = batch_info)

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
                                          is_normalized = is_normalized,
                                          batch_info = batch_info,
                                          is_bc = is_bc)

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

#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of RNA-seq transcript data to
#' an object of the class 'seqData'. Objects of the class 'seqData' are lists
#' with two obligatory components \code{e_data} and \code{f_data}. An optional
#' list component \code{e_meta} is used if analysis or visualization at other
#' levels (e.g. gene, protein, pathway) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of RNA transcripts observed and \eqn{n} is the number
#'   of samples (an additional transcript identifier/name column should also be
#'   present anywhere in the data.frame). Each row corresponds to data for each
#'   transcript One column specifying a unique identifier for each transcript
#'   (row) must be present. All counts are required to be raw for processing.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample. For library size normalization, this can be provided
#'   here or calculated from columns in e_data.
#' @param e_meta an optional data.frame with at least \eqn{p} rows. Each row
#'   corresponds to a transcript with one column giving transcript names (must be
#'   named the same as the column in \code{e_data}) and other columns giving
#'   meta information (e.g. mappings of transcripts to genes or proteins).
#' @param edata_cname character string specifying the name of the column
#'   containing the transcript identifiers in \code{e_data} and \code{e_meta} 
#'   (if applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the gene identifiers (or other mapping variable) in
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
#' @details Objects of class 'seqData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Only 'counts' is valid for 
#'   'seqData'. \cr \tab \cr is_normalized \tab A logical argument, specifying
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
#' data("seq_edata")
#' data("seq_fdata")
#' data("seq_emeta")
#' mypepData <- as.seqData(e_data = seq_edata,
#'                         e_meta = seq_emeta,
#'                         f_data = seq_fdata,
#'                         edata_cname = "Seq_Tag_ID",
#'                         fdata_cname = "SampleID",
#'                         emeta_cname = "Seq_Tag_ID")
#' }
#'
#' @author Kelly Stratton, Lisa Bramer
#' @seealso \code{\link{as.proData}}
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.metabData}}
#' @seealso \code{\link{as.nmrData}}
#'
#' @export
#'
as.seqData <- function (e_data, f_data, e_meta = NULL,
                        edata_cname, fdata_cname, emeta_cname = NULL,
                        techrep_cname = NULL, ...) {
  .as.seqData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname,
              techrep_cname, ...)
}

## RNA-seq data ##
.as.seqData <- function (e_data, f_data, e_meta = NULL, edata_cname,
                         fdata_cname, emeta_cname = NULL,
                         techrep_cname = NULL,
                         data_scale = "counts",
                         is_normalized = FALSE, norm_info = list(),
                         data_types = NULL, check.names = TRUE) {
  
  # Set the original data scale to the input data scale.
  data_scale_orig <- data_scale
  
  # Define the dType variable. This is used for customizing the warnings and
  # errors according to the data type (peptide, protein, lipid, ...).
  dType <- c('RNA transcripts', 'as.seqData')
  
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
  
  ### Should we add these to pre_flight? 
  # Analyses must have raw counts
  
  nums <- e_data[which(colnames(e_data) != edata_cname)]
  notint <- rowSums(nums%%1) !=  0
  if(any(notint)){
    warning("Non-integers detected. Analyses supported by pmartR for RNA-seq data require raw counts.")
  }
  
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
  class(res) = "seqData"
  
  return(res)
  
}


# pre_flight checks all of the data frames, id column names, and other input
# values to ensure they are all the correct class. It also checks the e_data,
# f_data, and e_meta data frames to make sure they are the right dimension.
#
# pre_flight returns e_data, f_data, e_meta, and emeta_cname. Each of these
# objects could have been modified in this function.
#
pre_flight <- function (e_data,
                        f_data,
                        e_meta,
                        edata_cname,
                        fdata_cname,
                        emeta_cname,
                        techrep_cname,
                        data_scale,
                        is_normalized,
                        norm_info,
                        data_types,
                        check.names,
                        dType,
                        is_bc,
                        batch_info) {

  # Verify classes/values of arguments -----------------------------------------
  
  # Check if e_data is a tibble or data.table. If it is convert to a data frame.
  if (inherits(e_data, "tbl_df") ||
      inherits(e_data, "tbl") ||
      inherits(e_data, "data.table")) {

    e_data <- data.frame(e_data, check.names = check.names)

  }
  
  # Make sure e_data is a data.frame.
  if (!inherits(e_data, "data.frame")) {
    
    stop ("e_data must be of class 'data.frame'")
    
  }
  
  # Check if f_data is a tibble or data.table. If it is convert to a data frame.
  if (inherits(f_data, "tbl_df") ||
      inherits(f_data, "tbl") ||
      inherits(f_data, "data.table")) {

    f_data <- data.frame(f_data, check.names = check.names)

  }
  
  # Make sure f_data is a data.frame.
  if (!inherits(f_data, "data.frame")) {
    
    stop ("f_data must be of class 'data.frame'")
    
  }
  
  # Determine if e_meta is present.
  if (!is.null(e_meta)) {
    
    # Check if e_meta is a tibble or data.table. If it is convert to a data
    # frame.
    if (inherits(e_meta, "tbl_df") ||
        inherits(e_meta, "tbl") ||
        inherits(e_meta, "data.table")) {

      e_meta <- data.frame(e_meta, check.names = check.names)

    }
    
    # Test that e_meta is a data frame.
    if (!inherits(e_meta, "data.frame")) {
      
      # Lay down an error because e_meta is not a data frame.
      stop ("e_meta must be of class 'data.frame'")
      
    }
    
  }
  
  # Investigate the columns of edata and ensure they are the correct class and
  # do not contain any erroneous values.
  str_col(edata = e_data,
          edata_cname = edata_cname)
  
  # check remaining input params are of correct class # added by KGS 5/1/2020
  ## cnames are character strings
  if (!is.null(edata_cname)) {
    if (!inherits(edata_cname, "character")) {
      stop ("edata_cname must be of the class 'character'")}
  }
  if (!is.null(fdata_cname)){
    if (!inherits(fdata_cname, "character")) {
      stop ("fdata_cname must be of the class 'character'")}
  }
  if (!is.null(emeta_cname)){
    if (!inherits(emeta_cname, "character")) {
      stop ("emeta_cname must be of the class 'character'")}
  }
  
  # Inspect the is_normalized argument.
  if (!is.null(is_normalized)) {
    
    # Make sure it is logical (not irrational:)).
    if(!inherits(is_normalized, "logical")) {
      
      # BAM!! pmart 1 user 0.
      stop ("is_normalized must be of the class 'logical'")
      
    }
    
  }

  # Inspect the is_bc argument
  if (!is.null(is_bc)) {

    # Make sure it is logical
    if (!inherits(is_bc, "logical")) {

      stop ("is_bc must be of the class 'logical'")

    }

  }

  # Examine the norm_info argument. Ensure it is a list.
  if (!inherits(norm_info, "list")) {
    
    # Throw an error at the user.
    stop ("norm_info must be of the class 'list'")
    
  }

  # Examine the batch_info argument. Ensure it is a list.
  if (!inherits(batch_info, "list")) {

    # Throw an error at the user
    stop ("batch_info must be of the class 'list'")
  }

  # Check the data_types argument.
  if (!is.null(data_types)) {
    
    # Confirm it is a character string.
    if (!inherits(data_types, "character")) {
      
      # Deliver a devastating blow to the users morale.
      stop ("data_types must be of the class 'character'")
      
    }
    
  }
  
  # Investigate the check.names argument.
  if (!is.null(check.names)) {
    
    # Establish its class.
    if (!inherits(check.names, "logical")) {
      
      # Shell out an error that check.names must be logical.
      stop ("check.names must be of the class 'logical'")
      
    }
    
  }
  
  
  # Make sure data_scale is one of the acceptable strings
  if(dType[[2]] == "as.seqData"){
    
    if (data_scale != c("counts")) {
      
      # Throw an error because data_scale is not an acceptable form.
      stop ("data_scale must be 'counts' for as.seqData")
      
    }
    
  } else {
    if (!(data_scale %in% c("abundance", "log", "log2", "log10"))) {
      
      # Throw an error because data_scale is not an acceptable form.
      stop (paste0("data_scale must be one of the following for ",
                  dType[[2]], ":",
                  " 'abundance', 'log', 'log2', or 'log10'."))
      
    }
  }
  
  # Check if techrep_cname is null or not.
  if (!is.null(techrep_cname)) {
    
    # Check that techrep_cname is a character string.
    if (!inherits(techrep_cname, "character") || length(techrep_cname) == 0) {
      
      # Use an error to let the user know there is little hope.
      stop (paste("techrep_cname must be a character string specifying a",
                  "column in f_data",
                  sep = ' '))
      
    }
    
  }
  
  # Check column names in e_data, f_data, and e_meta ---------------------------
  
  # Ensure the ID column exists in e_data.
  if (!(edata_cname %in% names(e_data))) {
    
    # Return an error if the ID column doesn't exist.
    stop (paste("The",
                dType[[1]],
                "ID column",
                edata_cname,
                "is not found in e_data. See details of",
                dType[[2]],
                "for specifying column names.",
                sep = " "))
    
  }
  
  # Check if techrep_cname is null or not.
  if (!is.null(techrep_cname)) {
    
    # Check that techrep_cname is in f_data and is not the same as fdata_cname.
    if (!(techrep_cname %in%
          colnames(f_data[, -which(names(f_data) == fdata_cname)]))) {
      
      stop(paste("Specified technical replicate column was not found in",
                 "f_data or was the same as fdata_cname",
                 sep = ' '))
      
    }
    
    # Check that the tech rep column does not have a unique value in each row.
    if (length(unique(f_data[,techrep_cname])) == nrow(f_data)) {
      
      stop (paste("Specified technical replicate column had a unique value for",
                  "each row.  Values should specify groups of technical",
                  "replicates belonging to a biological sample.",
                  sep = ' '))
      
    }
    
  }
  
  # Check if e_meta is NULL and emeta_cname is non-NULL #
  if (is.null(e_meta) && !is.null(emeta_cname)) {
    
    # Set emeta_cname to null and state that it will not be used.
    emeta_cname <- NULL
    message("emeta_cname set to NULL, no e_meta object was provided.")
    
  }
  
  # Verify e_meta is not null.
  if (!is.null(e_meta)) {
    
    # Confirm the peptide ID column exists in e_meta.
    if (!(edata_cname %in% names(e_meta))) {
      
      # Return an error if the ID column does not exist in e_meta.
      stop (paste("The",
                  dType[[1]],
                  "ID column",
                  edata_cname,
                  "is not found in e_meta. The column name containing the",
                  dType[[1]],
                  "IDs must match for e_data and e_meta. See details of",
                  dType[[2]],
                  "for specifying column names.",
                  sep = " "))
      
    }
    
    # Check if the emeta_cnames argument is null.
    if (is.null(emeta_cname)) {
      
      stop ("Since e_meta is non-NULL, emeta_cname must also be non-NULL.")
      
    } else {
      
      # If emeta_cname is not null ensure this column is present in e_meta.
      if (!(emeta_cname %in% names(e_meta))) {
        
        stop (paste("Mapping variable column",
                    emeta_cname,
                    "not found in e_meta. See details of",
                    dType[[2]],
                    "for specifying column names.",
                    sep = " "))
        
      }
      
    }
    
  }
  
  # Verify that the Sample column name is in the f_data column names #
  if (!(fdata_cname %in% names(f_data))) {
    
    # If this sample column name is not present return an error.
    stop (paste("Sample column",
                fdata_cname,
                "not found in f_data. See details of",
                dType[[2]],
                "for specifying column names.",
                sep = " "))
    
  }
  
  # Make sure the word 'Group' does not appear in f_data column names.
  if ("Group" %in% names(f_data)) {
    
    # Find the column number where the name "Group" occurs. This index will be
    # used to change the name to "group" because the group_designation function
    # creates a column named "Group". Some functions merge data frames with
    # f_data and having two columns named "Group" causes issues.
    group_idx <- which(names(f_data) == "Group")
    
    # Change the name to "group" so the merge function doesn't get confused.
    names(f_data)[group_idx] <- "group"
    
    # Let the user know they have overstepped their bounds and must be put in
    # their appropriate place.
    message(paste("A column in f_data is named 'Group'. This name is reserved",
                  "for use in the group_designation funtion. The column name",
                  "has been changed to 'group'.",
                  sep = " "))
    
  }
  
  # Ensure the data frames agree with each other -------------------------------
  
  # check that all samples in e_data (column names of e_data) are present in
  # f_data (rows of f_data in the fdata_cname column) #
  edat_sampid = which(names(e_data) == edata_cname) # fixed by KS 10/13/2016
  samps.miss = sum(!(names(e_data[,-edat_sampid]) %in% f_data[,fdata_cname]))
  if( samps.miss > 0) stop (paste( samps.miss,
                                   " samples from e_data not found in f_data",
                                   sep = ""))
  
  # check for any extra samples in f_data than in e_data - necessary to remove
  # before group_designation function #
  if (any(!(f_data[, fdata_cname] %in% names(e_data)))){
    
    # Remove rows found in f_data that are not also in e_data.
    f_data <- f_data[-which(!(f_data[, fdata_cname] %in% names(e_data))), ]
    
    # Throw down a warning that the extra rows in f_data were removed.
    warning (paste("Extra samples were found in f_data that were not in",
                   "e_data. These have been removed from f_data.",
                   sep = ' '))
    
  }
  
  # if e_meta is provided, remove any extra features that are not also found in
  # e_data.
  if(!is.null(e_meta)){
    if(any(!(e_meta[,edata_cname] %in% e_data[,edata_cname]))){
      
      # Remove any rows in e_meta corresponding to IDs that are not also found
      # in e_data.
      e_meta <- e_meta[-which(!(e_meta[, edata_cname] %in%
                                  e_data[, edata_cname])),]
      
      # Slam the user with a warning that the e_meta data frame was modified.
      warning (paste("Extra",
                     dType[[1]],
                     "were found in e_meta that were not in e_data.",
                     "These have been removed from e_meta.",
                     sep = " "))
    }
  }
  
  # Execute checks on e_data ---------------------------------------------------
  
  # Ensure the rows in e_data are unique.
  if (nrow(e_data) != length(unique(e_data[, edata_cname]))) {
    
    # Rewrite e_data with only the unique rows of the data frame.
    e_data <- unique(e_data)
    
    # Check if the unique data frame has non unique IDs.
    if (nrow(e_data) != length(unique(e_data[, edata_cname]))) {
      
      # Return an error if some IDs are repeated in e_data.
      stop ("The 'edata_cname' identifier is non-unique.")
      
    }
    
  }
  
  # Verify the data scale and if there are zeros in edata.
  if (data_scale == 'abundance' && any(na.omit(e_data == 0))) {
    
    # Exchange 0 for NA in edata.
    e_data <- replace_zeros(edata = e_data,
                            edata_cname = edata_cname)
    
  }
  
  # Perform checks on f_data ---------------------------------------------------
  
  # check that f_data has at least 2 columns #
  if (ncol(f_data) < 2) stop ("f_data must contain at least 2 columns")
  
  # Check if techrep_cname is null or not.
  if (!is.null(techrep_cname)) {
    
    # Check that the tech rep column does not have a unique value in each row.
    if (length(unique(f_data[,techrep_cname])) == nrow(f_data)) {
      
      stop (paste("Specified technical replicate column had a unique value for",
                  "each row.  Values should specify groups of technical",
                  "replicates belonging to a biological sample.",
                  sep = ' '))
      
    }
    
  }
  
  # Conduct checks on e_meta ---------------------------------------------------
  
  # If e_meta is provided, check that all peptides, proteins, ... in e_data also
  # occur in e_meta.
  if (!is.null(e_meta)){
    
    # If e_data has more rows than e_meta return an error.
    if (sum(!(e_data[,edata_cname] %in% e_meta[,edata_cname])) > 0) {
      
      stop (paste("Not all",
                  dType[[1]],
                  "in e_data are present in e_meta.",
                  sep = " "))
      
    }
    
  }
  
  # if e_meta is provided check that there are no duplicates of edata_cname and 
  # emeta_cname combinations in the e_meta
  if(!is.null(e_meta)){
    # identify which columns are edata and emeta cnames in the emeta dataset
    id_edat_cname = which(colnames(e_meta) == edata_cname)
    id_emet_cname = which(colnames(e_meta) == emeta_cname)
    
    # find the subset of just those 2 columns
    sub_emeta <- unique(data.frame(e_meta[,id_edat_cname],e_meta[,id_emet_cname]))
    
    # are the two number of rows equal to each other?
    if(nrow(sub_emeta) != nrow(e_meta)){
      stop (paste("Not all e_data cname and e_meta cname combinations are unique"))
    }
  }

  # Return the (possibly updated) data frames and cnames.
  return (list(e_data = e_data,
               f_data = f_data,
               e_meta = e_meta,
               emeta_cname = emeta_cname))
  
}

# Check to see if Excel or MATLAB is ruining our lives with their silly number
# conventions.
#
# str_edata checks the structure of each column in e_data. If a number is
# divided by zero in Excel (#DIV/0!) or if a number is +- infinity (#NUM!) the
# column(s) containing these values will be read in as character vectors instead
# of numeric. The same is true with MATLAB (so I have heard) when it exports a
# file containing NaNs. When R reads in this file the column(s) containing NaN
# will be read in as character vectors.
#
# edata - The e_data data frame.
# edata_cname - The name of the identification column in e_data.
#
# str_col only returns an error if one occurs. Otherwise it does not return
# anything.
#
str_col <- function (edata,
                     edata_cname) {
  
  # Determine the index of the id column
  id_col <- which(names(edata) == edata_cname)
  
  # Determine the number of columns in edata without the id column.
  n_col <- length(edata)
  
  # Create a vector that will hold the column indices for non numeric columns.
  non_numeric <- NULL
  
  # Create a vector that will hold the column indices for columns that contain
  # infinite values.
  too_vast <- NULL
  
  # Create a counter that will fill in the vector containing the column number
  # for the non-numeric columns.
  v <- 0
  
  # Create a counter that will fill in the vector containing the column number
  # for the columns containing infinite values.
  a <- 0
  
  # Loop through each column of edata (minus the id column). The call to seq_len
  # will loop through each index of edata except the one belonging to the column
  # that contains the IDs. For example, if the first column contains the IDs
  # the for loop will start at 2.
  for (e in seq_len(n_col)[-id_col]) {
    
    # Check if the current column is not numeric.
    if (!inherits(edata[, e],
                  c('integer', 'double', 'numeric'))) {
      
      # Update the counter used to fill in the non_numeric vector.
      v <- v + 1
      
      # Append the column number to the non_numeric vector.
      non_numeric[[v]] <- e
      
      # If the column is numeric the following code will run.
    } else {
      
      # Check if the column contains infinite values.
      if (any(is.infinite(edata[, e]))) {
        
        # Update the counter used to fill in the too_vast vector.
        a <- a + 1
        
        # Include the column number in the too_vast vector.
        too_vast[[a]] <- e
        
      }
      
    }
    
  }
  
  # Confirm whether or not there are any non-numeric columns
  if (length(non_numeric) > 0) {
    
    # Check the length of non_numeric.
    if (length(non_numeric) == 1) {
      
      # Use proper English grammar when there is one naughty column.
      grammar <- c('Column', 'contains')
      
    } else {
      
      # Use proper English grammar when there are multiple naughty columns.
      grammar <- c('Columns', 'contain')
      
    }
    
    # Forcefully tell the user their data is not acceptable. In other words, we
    # don't want their crap data because our life is crazy enough already.
    stop (paste(grammar[[1]],
                knitr::combine_words(non_numeric),
                'of e_data',
                grammar[[2]],
                'non-numeric values.',
                sep = ' '))
    
  }
  
  # Confirm whether or not there are any columns containing infinity.
  if (length(too_vast) > 0) {
    
    # Check the length of too_vast.
    if (length(too_vast) == 1) {
      
      # Use proper English grammar when there is one boundless column.
      grammar <- c('Column', 'contains')
      
    } else {
      
      # Use proper English grammar when there are multiple boundless columns.
      grammar <- c('Columns', 'contain')
      
    }
    
    # Throw down an error letting the user know it is impossible to have an
    # infinite number of something in their sample. (How would it all fit?)
    stop (paste(grammar[[1]],
                knitr::combine_words(too_vast),
                'of e_data',
                grammar[[2]],
                'infinite values.',
                sep = ' '))
    
  }
  
}

#' Replace 0 with NA
#'
#' This function finds all instances of 0 in e_data and replaces them with NA.
#'
#' @param e_data A \eqn{p \times n + 1} data frame of expression data, where
#'        \eqn{p} is the number of xxx observed and \eqn{n} is the
#'        number of samples.
#'
#' @param edata_cname A character string specifying the name of the ID column in
#'        the e_data data frame.
#'
#' @details This function is used in the as.pepData, as.proData, as.lipidData,
#'          as.metabData, as.isobaricpepData, and as.nmrData functions to
#'          replace any 0 values with NAs.
#'
#' @return An updated e_data data frame where all instances of 0 have been
#'         replaced with NA.
#'
replace_zeros <- function(edata,
                          edata_cname) {
  
  # Acquire the index of the edata_cname column.
  id_col <- which(names(edata) == edata_cname)
  
  # Enumerate the number of zeros to be replaced with NA
  n_zeros <- sum(edata[, -id_col] == 0,
                 na.rm = TRUE)
  
  # Loop through each column in e_data.
  for (e in 1:ncol(edata)) {
    
    # Check if the current column is NOT the ID column.
    if (e != id_col) {
      
      # Find indices where the value is 0.
      inds <- which(edata[, e] == 0)
      
      # Check if any values in the eth column of edata match 0. If there are no
      # matches then which() will return integer(0) -- which has length 0.
      if (length(inds) == 0) {
        
        # Move to the next column in the data frame.
        next
        
        # If the length of inds is greater than 0 enter the else statement.
      } else {
        
        # Replace 0 with NA.
        edata[inds, e] <- NA
        
      }
      
      # If e is equal to the index of the ID column enter the else statement.
    } else {
      
      # Move to the next column in the data frame.
      next
      
    }
    
  }
  
  # Report the number of replaced elements in e_data
  message(paste(n_zeros,
                "instances of",
                0,
                "have been replaced with",
                NA,
                sep = " "))
  
  # Return the updated edata object.
  return (edata)
  
}
