#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of peptide-level data to an object of the class 'pepData'. Objects of the class 'pepData' are lists with two obligatory components \code{e_data} and \code{f_data}. An optional list component \code{e_meta} is used if analysis or visualization at other levels (e.g. protein) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where \eqn{p} is the number of peptides observed and \eqn{n} is the number of samples (an additional peptide identifier/name column should also be present anywhere in the data.frame). Each row corresponds to data for each peptide. One column specifying a unique identifier for each peptide (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a sample with one column giving the unique sample identifiers found in e_data column names and other columns providing qualitative and/or quantitative traits of each sample.
#' @param e_meta an optional data.frame with at least \eqn{p} rows. Each row corresponds to a peptide with one column giving peptide names (must be named the same as the column in \code{e_data}) and other columns giving meta information (e.g. mappings of peptides to proteins).
#' @param edata_cname character string specifying the name of the column containing the peptide identifiers in \code{e_data} and \code{e_meta} (if applicable).
#' @param emeta_cname character string specifying the name of the column containing the protein identifiers (or other mapping variable) in \code{e_meta} (if applicable). Defaults to NULL. If \code{e_meta} is NULL, then either do not specify \code{emeta_cname} or specify it as NULL. If \code{e_meta} is NULL, then specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column containing the sample identifiers in \code{f_data}.
#' @param ... further arguments
#'
#' @details Objects of class 'pepData' contain some attributes that are referenced by downstream functions. These attributes can be changed from their default value by manual specification. A list of these attributes as well as their default values are as follows:
#' \tabular{ll}{
#' data_scale \tab Scale of the data provided in \code{e_data}. Acceptable values are 'log2', 'log10', 'log', and 'abundance', which indicate data is log base 2, base 10, natural log transformed, and raw abundance, respectively. Default is 'abundance'. \cr
#' \tab \cr
#' data_norm \tab A logical argument, specifying whether the data has been normalized or not. Default value is FALSE. \cr
#' \tab \cr
#' norm_info \tab Default value is NULL. When a normalization is applied to the data, this becomes populated with a list containing the normalization function, normalization subset and subset parameters, the location and scale parameters used to normalize the data, and the location and scale parameters used to backtransform the data (if applicable). \cr
#' \tab \cr
#' data_types \tab Character string describing the type of data (e.g.'Positive ion'). Default value is NULL. \cr
#' }
#' Computed values included in the \code{data_info} attribute are as follows:
#' \tabular{ll}{
#' num_edata \tab The number of unique \code{edata_cname} entries.\cr
#' \tab \cr
#' num_miss_obs \tab The number of missing observations.\cr
#' \tab \cr
#' num_emeta \tab The number of unique \code{emeta_cname} entries. \cr
#' \tab \cr
#' prop_missing \tab The proportion of \code{e_data} values that are NA. \cr
#' \tab \cr
#' num_samps \tab The number of samples that make up the columns of \code{e_data}.\cr
#' \tab \cr
#' meta_info \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#' \tab \cr
#' }
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_edata")
#' data("pep_fdata")
#' data("pep_emeta")
#' mypepData <- as.pepData(e_data = pep_edata, e_meta = pep_emeta, f_data = pep_fdata, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID")
#'}
#' @author Kelly Stratton, Lisa Bramer
#' @seealso \code{\link{as.proData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.metabData}}
#'
#' @export
as.pepData <- function(e_data, f_data, e_meta = NULL, edata_cname, fdata_cname, emeta_cname = NULL, ...){
  .as.pepData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname, ...)
}

## peptide data ##
.as.pepData <- function(e_data, f_data, e_meta = NULL, edata_cname, fdata_cname,
                        emeta_cname = NULL, data_scale = "abundance",
                        data_norm = FALSE, norm_info = NULL, data_types=NULL){

  # initial checks #

  # check that e_data and f_data are data.frames #
  if(class(e_data) != "data.frame") stop("e_data must be of the class 'data.frame'")
  if(class(f_data) != "data.frame") stop("f_data must be of the class 'data.frame'")

  # check to see if e_meta is NULL, if not check that it is a data.frame #
  if(!is.null(e_meta)){ if(class(e_meta) != "data.frame") stop("e_meta must be of the class 'data.frame'")}

  # check that the peptide column exists in e_data and e_meta (if applicable) #
  if(!(edata_cname %in% names(e_data))) stop(paste("Peptide column ", edata_cname," not found in e_data. See details of as.pepData for specifying column names.", sep = ""))
  if(!is.null(e_meta)){
    if(!(edata_cname %in% names(e_meta))) stop(paste("Peptide column ", edata_cname," not found in e_meta. Column names for peptide names must match for e_data and e_meta. See details of as.pepData for specifying column names.", sep = ""))
  }

  # if e_meta is NULL set emeta_cname to NULL #
  if(is.null(e_meta)) emeta_cname = NULL

  # if e_meta is not NULL check that emeta_cname column is found #
  if(!is.null(e_meta)){
    if(!is.null(emeta_cname)){
      if(!(emeta_cname %in% names(e_meta))) stop(paste("Protein column ", emeta_cname, " not found in e_meta. See details of as.pepData for specifying column names. Blarg.", sep = "") )
    }
  }

  # check that the Sample column name is in f_data column names #
  if(!(fdata_cname %in% names(f_data))) stop(paste("Sample column ", fdata_cname, " not found in f_data. See details of as.pepData for specifying column names.", sep = ""))

  # check that all samples in e_data are present in f_data #
  edat_sampid = which(names(e_data) == edata_cname) # fixed by KS 10/13/2016
  samps.miss = sum(!(names(e_data[,-edat_sampid]) %in% f_data[,fdata_cname]))
  if( samps.miss > 0) stop(paste( samps.miss, " samples from e_data not found in f_data", sep = ""))

  # check for any extra samples in f_data than in e_data - necessary to remove before group_designation function #
  if(any(!(f_data[,fdata_cname] %in% names(e_data)))){
    f_data <- f_data[-which(!(f_data[,fdata_cname] %in% names(e_data))),]
    warning("Extra samples were found in f_data that were not in e_data. These have been removed from f_data.")
  }

  # check that f_data has at least 2 columns #
  if(ncol(f_data) < 2) stop("f_data must contain at least 2 columns")

  # if e_meta is provided, check that all peptides in e_data occur in e_meta #
  if(!is.null(e_meta)){
    if(sum(!(e_data[,edata_cname] %in% e_meta[,edata_cname])) > 0 ) stop("Not all peptides in e_data are present in e_meta")
  }

  # if e_meta is provided, remove any extra features that were provided #
  if(!is.null(e_meta)){
    if(any(!(e_meta[,edata_cname] %in% e_data[,edata_cname]))){
      e_meta <- e_meta[-which(!(e_meta[,edata_cname] %in% e_data[,edata_cname])),]
      warning("Extra peptides were found in e_meta that were not in e_data. These have been removed from e_meta.")
    }
  }

  # check that data_scale is one of the acceptable options #
  if(!(data_scale %in% c('log2', 'log10', 'log', 'count', 'abundance'))) stop(paste(data_scale, " is not a valid option for 'data_scale'. See details of as.pepData for specifics.", sep=""))

  # if e_meta is NULL, set emeta_cname to NULL #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }

  # check that e_data has unique rows #
  if(nrow(e_data) == length(unique(e_data[, edata_cname]))){
    e_data <- e_data
  }else{
    e_data_unique <- unique(e_data)
    if(nrow(e_data_unique) == length(unique(e_data_unique[, edata_cname]))){
      e_data <- e_data_unique
    }else{
      stop("The 'edata_cname' identifier is non-unique.")
    }
  }

  # store results #
  res = list(e_data = e_data, f_data = f_data, e_meta = e_meta)

  # set column name attributes #
  attr(res, "cnames") = list(edata_cname = edata_cname, emeta_cname = emeta_cname, fdata_cname = fdata_cname)

  # count missing values in e_data #
  num_miss_obs = sum(is.na(e_data[,-which(names(e_data)==edata_cname)]))
  prop_missing = mean(is.na(e_data[,-which(names(e_data)==edata_cname)]))

  # number of unique peptides #
  num_edata = length(unique(e_data[, edata_cname]))

  # number of samples #
  num_samps = ncol(e_data) - 1

  if(!is.null(e_meta)){
    # number of unique proteins that map to a peptide in e_data #
    if(!is.null(emeta_cname)){
      num_emeta = length(unique(e_meta[which(as.character(e_meta[, edata_cname]) %in% as.character(e_data[, edata_cname])), emeta_cname]))
    }else{num_emeta = NULL}
  }else{
    num_emeta = NULL
  }

  # set data information attributes #
  attr(res, "data_info") = list(data_scale = data_scale, data_norm = data_norm, norm_info = norm_info, num_edata = num_edata, num_miss_obs = num_miss_obs, num_emeta = num_emeta, prop_missing = prop_missing, num_samps = num_samps, data_types = data_types)

  # set meta data attributes #
  if(!is.null(e_meta)){
    attr(res, "meta_info") = TRUE
  }else{ attr(res, "meta_info") = FALSE}


  # set group dataframe attribute to NULL, will be filled in after running group_designation function #
  attr(res, "group_DF") = NULL

  # set filters attributes #
  if(!is.null(attr(res, "filters"))){
    attr(res, "filters") = attr(res, "filters")
  }else{ attr(res, "filters") = list()}

  # set class of list #
  class(res) = "pepData"

  return(res)

}
