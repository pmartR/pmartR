#' Molecule filter object
#'
#' This function returns a moleculeFilt object for use with \code{\link{applyFilt}}
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param min_num an integer value specifying the minimum number of times each feature must be observed across all samples. Default value is 2.
#'
#' @details Attribute of molecule_filt object is "total_poss_obs", the number of total possible observations for each feature (same as the number of samples)
#'
#' @return Object of class moleculeFilt (also a data.frame) that contains the molecule identifier and the number of samples for which the molecule was measured (not NA)
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' to_filter <- molecule_filter(omicsData = pep_object)
#' summary(to_filter, min_num = 2)
#'}
#'
#' @author Kelly Stratton
#'
#' @export

molecule_filter <- function(omicsData){
  ## some initial checks ##
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
               "'lipidData', or 'nmrData'",
               sep = ' '))
    
  }
  
  # Extricate the column number of the ID column.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Compute the number of non-missing values for each row.
  num_obs <- rowSums(!is.na(omicsData$e_data[, -id_col]))
  
  # Create a data frame with the ID column and the number of non-missing values.
  output <- data.frame(omicsData$e_data[, id_col], num_obs)
  names(output) <- c(get_edata_cname(omicsData), "Num_Observations")
  
  # Extract the 'data.frame' class from the the output data frame.
  orig_class <- class(output)
  
  # Create the moleculeFilt class and attach the data.frame class to it as well.
  class(output) <- c("moleculeFilt", orig_class)
  
  # Fabricate an attribute that has the total number of samples (columns in 
  # e_data minus the ID column). This will be used to ensure someone doesn't try
  # to filter e_data using a threshold larger than the number of samples.
  attr(output, "num_samps") <- get_data_info(omicsData)$num_samps
  
  # Return the completed object!!!
  return(output)
  
}


# feature_filter <- function(omicsData, min_num=2){
#   ## some initial checks ##
#
#   # check that omicsData is of appropriate class #
#   if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
#
#   edata_id = attr(omicsData, "cnames")$edata_cname
#   emeta_id = attr(omicsData, "cnames")$emeta_cname
#   samp_id = attr(omicsData, "cnames")$fdata_cname
#
#
#   # check that min_num is numeric and >=1 #
#   if(!inherits(min_num, "numeric") | min_num < 1) stop("min_num must be an integer greater than or equal to 1")
#   # check that min_num is less than the number of samples #
#   if(min_num > summary(omicsData)$num_samps) stop("min_num cannont be greater than the number of samples")
#
#
#   nonmiss <- !is.na(omicsData$e_data[, -which(names(omicsData$e_data) == edata_id)])
#
#   # get row sums of nonmiss = number of times each feature is observed #
#   num_obs <- rowSums(nonmiss)
#
#   # get indices for which ones don't meet the min requirement #
#   inds <- which(num_obs < min_num)
#
#   if(length(inds) < 1){
#     filter.edata <- NULL
#   }else{
#     filter.edata <- omicsData$e_data[, which(names(omicsData$e_data) == edata_id)][inds]
#   }
#
#   return(list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL))
# }
