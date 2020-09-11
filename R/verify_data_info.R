#' Verifies some basic attributes of an omicsData object.  Intended use is in testing scripts and internally within some functions.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData' created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @details This function will invisibly return TRUE or FALSE depending on whether or not the objects contents are consistent with its attributes.  If FALSE, there will also be a print output of which values are inconsistent.
#'
#' @return TRUE if the object is consistent for certain parts of the object. FALSE otherwise. 
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(pep_edata)
#' bool_val <- verify_data_info(pep_edata)
#'}
#'

verify_data_info <- function(omicsData){
  
  # save typing
  edata <- omicsData$e_data
  fdata <- omicsData$f_data
  emeta <- omicsData$e_meta
  edata_ind <- which(colnames(edata) == attributes(omicsData)$cnames$edata_cname)
  fdata_ind <- which(colnames(fdata) == attributes(omicsData)$cnames$fdata_cname)
  group_DF <- attributes(omicsData)$group_DF
  
  if(!is.null(emeta)){
    emeta_ind <- which(colnames(emeta) == attributes(omicsData)$cnames$edata_cname)
  }
  
  n = ncol(edata) - 1
  m = nrow(edata)
  #
  
  if(any(is.null(edata), is.null(fdata))) stop("One or both of edata/fdata are missing")
  if(any(c(length(edata_ind), length(fdata_ind)) == 0)) stop("One or more of your e_data or f_data had no column corresponding to the provided column names")
  if(m < 2) stop("Your e_data object did not have at least 2 columns (ID and data)")
  if(n < 1) stop("Your e_data object did not have at least 1 row/feature/molecule")
    
  # edata and fdata have the same sample ID's
  sample_IDs <- all(colnames(edata)[-edata_ind] %in% as.character(fdata[,fdata_ind]))
  
  # number of columns in edata - 1 is the same as the num samps attribute and the rows of fdata
  n_samps <- all((ncol(edata) - 1) == c(attributes(omicsData)$data_info$num_samps, nrow(fdata))) 
  
  # num missing in e_data is the same as attribute
  n_missing <- sum(is.na(edata[-edata_ind])) == attributes(omicsData)$data_info$num_miss_obs
  
  # prop missing in edata is the same as attribute
  prop_missing <- sum(is.na(edata[,-edata_ind]))/(m*n) == attributes(omicsData)$data_info$prop_missing
  
  # number of unique biomolecules same as attribute
  num_edata <- length(unique(edata[,edata_ind])) == attributes(omicsData)$data_info$num_edata
  
  # if group DF exists, it should have the same column as fdata
  group_DF_names <- if(!is.null(group_DF)) all(group_DF[,which(colnames(group_DF) == attributes(omicsData)$cnames$fdata_cname)] %in% fdata[,fdata_ind]) else NA
  
  # emeta_ID == edata_ID
  edata_emeta_ID_col <- if(!is.null(emeta)) length(setdiff(edata[,edata_ind], emeta[,emeta_ind])) == 0 else NA
  
  # store T/F values in a named list
  TF_list <- list(sample_IDs = sample_IDs, n_samps = n_samps, n_missing = n_missing, prop_missing = prop_missing, num_edata = num_edata,
                  group_DF_names = group_DF_names, edata_emeta_ID_col = edata_emeta_ID_col)
  
  # if there are bad attributes, alert user to which ones
  if(!all(unlist(TF_list), na.rm = TRUE)){
    bad_attrs <- which(TF_list == FALSE)
    warning(paste0("The following items had inconsistencies: ", paste(names(TF_list[bad_attrs]), collapse = " | ")))
  }
  
  # return T if data is okay, F otherwise
  return(invisible(all(unlist(TF_list), na.rm = TRUE)))
}