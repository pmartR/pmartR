#' Molecule filter object
#'
#' This function returns a moleculeFilt object for use with \code{\link{applyFilt}}
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
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
  if(!class(omicsData) %in% c("pepData", "proData", "metabData", "lipidData")) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")

  edata_id = attr(omicsData, "cnames")$edata_cname
  emeta_id = attr(omicsData, "cnames")$emeta_cname
  samp_id = attr(omicsData, "cnames")$fdata_cname

#   # check that min_num is numeric and >=1 #
#   if(class(min_num) != "numeric"| min_num < 1) stop("min_num must be an integer greater than or equal to 1")
#   # check that min_num is less than the number of samples #
#   if(min_num > summary(omicsData)$num_samps) stop("min_num cannont be greater than the number of samples")


  nonmiss <- !is.na(omicsData$e_data[, -which(names(omicsData$e_data) == edata_id)])

  # get row sums of nonmiss = number of times each feature is observed #
  num_obs <- rowSums(nonmiss)

  # output #
  output <- data.frame(omicsData$e_data[, edata_id], num_obs)
  names(output) <- c(edata_id, "Num_Observations")

  orig_class <- class(output)
  class(output) <- c("moleculeFilt", orig_class)
  attr(output, "total_poss_obs") <- length(samp_id)
  #attr(output, "min_num") <- min_num

  return(output)

#   # get indices for which ones don't meet the min requirement #
#   inds <- which(num_obs < min_num)
#
#   if(length(inds) < 1){
#     filter.edata <- NULL
#   }else{
#     filter.edata <- omicsData$e_data[, which(names(omicsData$e_data) == edata_id)][inds]
#   }

#  return(list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL))
}


# feature_filter <- function(omicsData, min_num=2){
#   ## some initial checks ##
#
#   # check that omicsData is of appropriate class #
#   if(!class(omicsData) %in% c("pepData", "proData", "metabData", "lipidData")) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
#
#   edata_id = attr(omicsData, "cnames")$edata_cname
#   emeta_id = attr(omicsData, "cnames")$emeta_cname
#   samp_id = attr(omicsData, "cnames")$fdata_cname
#
#
#   # check that min_num is numeric and >=1 #
#   if(class(min_num) != "numeric"| min_num < 1) stop("min_num must be an integer greater than or equal to 1")
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
