#' Convert Data to Appropriate pmartR Class
#'
#' Converts a list object or several data.frames of NMR-generated metabolomic-level data to an object of the class 'nmrData'. Objects of the class 'nmrData' are lists with two obligatory components \code{e_data} and \code{f_data}. An optional list component \code{e_meta} is used if analysis or visualization at other levels (e.g. metabolite identification) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where \eqn{p} is the number of metabolites observed and \eqn{n} is the number of samples. Each row corresponds to data for each metabolite. One column specifying a unique identifier for each metabolite (row) must be present.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a sample with one column giving the unique sample identifiers found in e_data column names and other columns providing qualitative and/or quantitative traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds to a metabolite with one column giving metabolite names (must be named the same as the column in \code{e_data}) and other columns giving meta information.
#' @param edata_cname character string specifying the name of the column containing the metabolite identifiers in \code{e_data} and \code{e_meta} (if applicable).
#' @param emeta_cname character string specifying the name of the column containing the mapped identifiers in \code{e_meta} (if applicable). Defaults to NULL. If \code{e_meta} is NULL, then either do not specify \code{emeta_cname} or specify it as NULL. If \code{e_meta} is NULL, then specify \code{emeta_cname} as NULL.
#' @param fdata_cname character string specifying the name of the column containing the sample identifiers in \code{f_data}.
#' @param techrep_cname character string specifying the name of the column in \code{f_data} containing the identifiers for the biological samples if the observations represent technical replicates.  This column is used to collapse the data when \code{combine_techreps} is called on this object.  Defaults to NULL (no technical replicates). 
#' @param ... further arguments
#'
#' @details Objects of class 'nmrData' contain some attributes that are referenced by downstream functions. These attributes can be changed from their default value by manual specification. A list of these attributes as well as their default values are as follows:
#' \tabular{ll}{
#' data_scale \tab Scale of the data provided in \code{e_data}. Acceptable values are 'log2', 'log10', 'log', and 'abundance', which indicate data is log base 2, base 10, natural log transformed, and raw abundance, respectively. Default is 'abundance'. \cr
#' \tab \cr
#' is_normalized \tab A logical argument, specifying whether the data has been normalized or not. Default value is FALSE. \cr
#' \tab \cr
#' norm_info \tab Default value is an empty list, which will be populated with a single named element \code{is_normalized = is_normalized}. When a normalization is applied to the data, this becomes populated with a list containing the normalization function, normalization subset and subset parameters, the location and scale parameters used to normalize the data, and the location and scale parameters used to backtransform the data (if applicable). \cr
#' \tab \cr
#' data_types \tab Character string describing the type of data (e.g.'binned' or 'identified', for NMR data). Default value is NULL. \cr
#' \tab \cr
#' check.names \tab Logical defaults to TRUE. Indicates whether 'check.names' attribute of returned omicsData object is TRUE or FALSE. \cr
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
#' meta_info \tab A logical argument, specifying where the \code{e_meta} is provided.\cr
#' \tab \cr
#' }
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("nmr_edata")
#' data("nmr_fdata")
#' mynmrData <- as.nmrData(e_data = nmr_edata, f_data = nmr_fdata, edata_cname = "Metabolite", fdata_cname = "SampleID")
#'}
#'
#' @author Lisa Bramer, Kelly Stratton
#' @seealso \code{\link{as.metabData}}
#' @seealso \code{\link{as.isobaricpepData}}
#' @seealso \code{\link{as.lipidData}}
#' @seealso \code{\link{as.pepData}}
#' @seealso \code{\link{as.proData}}
#'
#' @export
as.nmrData <- function(e_data, f_data, e_meta = NULL, edata_cname, fdata_cname, emeta_cname = NULL, techrep_cname = NULL, ...){
  .as.nmrData(e_data, f_data, e_meta, edata_cname, fdata_cname, emeta_cname, techrep_cname, ...)
}

## metabolite data ##
.as.nmrData <- function(e_data, f_data, e_meta = NULL, edata_cname, fdata_cname,
                          emeta_cname = NULL, techrep_cname = NULL, data_scale = "abundance",
                          is_normalized = FALSE, norm_info=list(),
                          data_types=NULL, check.names = TRUE 
                          ){

  # initial checks #

  # check that e_data and f_data are data.frames #
  if(!inherits(e_data, "data.frame")) stop("e_data must be of the class 'data.frame'")
  if(!inherits(f_data, "data.frame")) stop("f_data must be of the class 'data.frame'")

  # check to see if e_meta is NULL, if not check that it is a data.frame #
  if(!is.null(e_meta)){ if(!inherits(e_meta, "data.frame")) stop("e_meta must be of the class 'data.frame'")}
  
  
  # check that remaining input params are of correct class # added by KGS 5/1/2020
  ## cnames are character strings
  if(!is.null(edata_cname)){ if(!inherits(edata_cname, "character")) stop("edata_cname must be of the class 'character'")}
  if(!is.null(fdata_cname)){ if(!inherits(fdata_cname, "character")) stop("fdata_cname must be of the class 'character'")}
  if(!is.null(emeta_cname)){ if(!inherits(emeta_cname, "character")) stop("emeta_cname must be of the class 'character'")}
  
  ## is_normalized is logical
  if(!is.null(is_normalized)){ if(!inherits(is_normalized, "logical")) stop("is_normalized must be of the class 'logical'")}
  
  ## norm_info is list
  if(!inherits(norm_info, "list")) stop("norm_info must be of the class 'list'")
  
  ## data_types is character if not missing
  if(!is.null(edata_cname)){ if(!inherits(edata_cname, "character")) stop("edata_cname must be of the class 'character'")}
  
  ## check.names is logical
  if(!is.null(check.names)){ if(!inherits(check.names, "logical")) stop("check.names must be of the class 'logical'")}

  
  # check that the metabolite column exists in e_data and e_meta (if applicable) #
  if(!(edata_cname %in% names(e_data))) stop(paste("Metabolite column ", edata_cname," not found in e_data. See details of as.nmrData for specifying column names.", sep = ""))

  if(!is.null(e_meta)){
    if(!(edata_cname %in% names(e_meta))) stop(paste("Metabolite column ", edata_cname," not found in e_meta. Column names for metabolite names must match for e_data and e_meta. See details of as.nmrData for specifying column names.", sep = ""))
  }

  # if e_meta is NULL and emeta_cname is non-NULL #
  if(is.null(e_meta) & !is.null(emeta_cname)){
    message("emeta_cname will not be used, no e_meta was provided") 
    emeta_cname = NULL
  } 
  
  # if e_meta is non-NULL and emeta_cname is NULL #
  if(!is.null(e_meta) & is.null(emeta_cname)) stop("if e_meta is non-NULL, emeta_cname must also be non-NULL")

  # if e_meta is not NULL check that identification columns are found #
  if(!is.null(e_meta)){
    if(!is.null(emeta_cname)){
      if(!(emeta_cname %in% names(e_meta))) stop(paste("Metabolite identification column ", emeta_cname, " not found in e_meta. See details of as.nmrData for specifying column names.", sep = "") )
    }
  }

  # check that the Sample column name is in f_data column names #
  if(!(fdata_cname %in% names(f_data))) stop(paste("Sample column ", fdata_cname, " not found in f_data. See details of as.nmrData for specifying column names.", sep = ""))

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

  # if e_meta is provided, check that all metabolites in e_data occur in e_meta #
  if(!is.null(e_meta)){
    if(sum(!(e_data[,edata_cname] %in% e_meta[,edata_cname])) > 0 ) stop("Not all metabolites in e_data are present in e_meta")
  }

  # if e_meta is provided, remove any extra features that were provided #
  if(!is.null(e_meta)){
    if(any(!(e_meta[,edata_cname] %in% e_data[,edata_cname]))){
      e_meta <- e_meta[-which(!(e_meta[,edata_cname] %in% e_data[,edata_cname])),]
      warning("Extra metabolites were found in e_meta that were not in e_data. These have been removed from e_meta.")
    }
  }

  # check that data_scale is one of the acceptable options #
  data_scale <- tolower(data_scale) # added by KGS 5/1/2020
  if(!(data_scale %in% c('log2', 'log10', 'log', 'abundance'))) stop(paste(data_scale, " is not a valid option for 'data_scale'. See details of as.nmrData for specifics.", sep=""))
  
  

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
  
  # check that technical replicate identifier column specifies at least one biological sample with 2 or more technical replicates.
  if(!is.null(techrep_cname)){
    if(!inherits(techrep_cname, "character") | length(techrep_cname) == 0) stop("techrep_cname must be a character string specifying a column in f_data")
    if(!(techrep_cname %in% colnames(f_data[,-which(names(f_data) == fdata_cname)]))) stop("Specified technical replicate column was not found in f_data or was the same as fdata_cname")
    if(length(unique(f_data[,techrep_cname])) == nrow(f_data)) stop("Specified technical replicate column had a unique value for each row.  Values should specify groups of technical replicates belonging to a biological sample")
  }

  # store results #
  res = list(e_data = e_data, f_data = f_data, e_meta = e_meta)

  # set column name attributes #
  attr(res, "cnames") = list(edata_cname = edata_cname, emeta_cname = emeta_cname, fdata_cname = fdata_cname, techrep_cname = techrep_cname)
  
  # count missing values in e_data #
  num_miss_obs = sum(is.na(e_data[,-which(names(e_data)==edata_cname)]))
  prop_missing = mean(is.na(e_data[,-which(names(e_data)==edata_cname)]))

  # number of unique metabolites #
  num_edata = length(unique(e_data[, edata_cname]))

  # number of samples #
  num_samps = ncol(e_data) - 1

  if(!is.null(e_meta)){
    # number of metabolite identifiers that map to a metabolite in e_data #
    if(!is.null(emeta_cname)){
      num_emeta = length(unique(e_meta[which(as.character(e_meta[, edata_cname]) %in% as.character(e_data[, edata_cname])), emeta_cname]))
    }else{num_emeta = NULL}
  }else{
    num_emeta = NULL
  }
  
  # set data information attributes #
  norm_info$is_normalized = is_normalized
  attr(res, "data_info") = list(data_scale = data_scale, norm_info = norm_info, num_edata = num_edata, num_miss_obs = num_miss_obs,
                                num_emeta = num_emeta, prop_missing = prop_missing, num_samps = num_samps, data_types = data_types)
  #set check.names attribute #
  attr(res, "check.names") = check.names 
  
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
  class(res) = c("nmrData") 

  return(res)

}


