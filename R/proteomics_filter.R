#' Proteomics filter object
#'
#' This function counts the number of peptides that map to each protein and/or the number of proteins to which each individual peptide maps.
#'
#' @param omicsData an object of class "pepData", the a result of \code{\link{as.pepData}}. The e_meta component of omicsData must be nonempty.
#'
#' @details min_num_peps is an attribute of the proteomicsFilt object
#'
#' @return a list of two character vectors, the first, \code{peptides_filt}, giving degenerate peptide names. The second, \code{proteins_filt}, gives the names of proteins which no longer have peptides mapping to them in the dataset.
#'
#' @examples
#' dontrun{
#' library(MSomicsQC)
#' data("pep_pepData")
#' my_filter <- proteomics_filter(omicsData = pep_pepData)
#' summary(my_filter, min_num_peps = 3)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export

proteomics_filter <- function(omicsData){

  ## initial checks ##

  # check that omicsData is of class 'pepData' #
  if(class(omicsData) != "pepData") stop("omicsData must be of class 'pepData'")

  # check that e_meta is not NULL #
  if(is.null(omicsData$e_meta)) stop("e_meta must be non-NULL")

  # get peptide and protein column names #
  pep_id = attr(omicsData, "cnames")$edata_cname
  pro_id = attr(omicsData, "cnames")$emeta_cname

  # check that protein column is not NULL #
  #if(is.null(pro_id)) stop("Cannot find degenerate peptides; protein column name is NULL.")

  # check that peptide and protein column names are non-null #
  if(is.null(pep_id)) stop("Peptide column name is NULL")
  if(is.null(pro_id)) stop("Protein column name is NULL")



  ## degenerate peptide portion ##

  # count the number of rows for each peptide #
  bypep = dplyr::group_by_(omicsData$e_meta, pep_id)

  count_bypep = dplyr::summarise_(.data = bypep, n = paste("length(",pro_id,")", sep = ""))
  count_bypep <- data.frame(count_bypep)

#   # pull peptides with more than one row #
#   degen_peps = as.character(data.frame(count_bypep[which(count_bypep$n > 1), ])[, pep_id])
#
#   #####################################
#   ## identify any proteins that now will not have peptides mapping to them ##
#   ## find rows in e_meta that correspond to peptides to be filtered ##
#   pepfilt.ids = which(omicsData$e_meta[,pep_id] %in% degen_peps)
#
#   ## find the proteins that are in the filter list but are not in the unfiltered lists ##
#   add_prots = as.character(setdiff(omicsData$e_meta[pepfilt.ids, pro_id], omicsData$e_meta[-pepfilt.ids, pro_id]))
#
#   if(length(add_prots)==0){add_prots = NULL}
#   if(length(degen_peps)==0){degen_peps = NULL}
#   return(list(peptides_filt = degen_peps, proteins_filt = add_prots))


  ## protein filter portion ##



  # group by proteins #
  bypro = dplyr::group_by_(omicsData$e_meta, pro_id)

  # count the number of rows for each protein #
  count_bypro = dplyr::summarise_(bypro, n = paste("length(",pro_id,")",sep=""))
  count_bypro <- data.frame(count_bypro)

#   # pull proteins with less than min_num_peps #
#   pro_filt = as.character(data.frame(count_bypro[which(count_bypro$n < min_num_peps), ])[, pro_id])
#
#
#   # determine which peptides no longer have a protein to map to  #
#   ## find rows in peptide.info that correspond to proteins to be filtered ##
#   protfilt.ids = which(omicsData$e_meta[,pro_id] %in% pro_filt)
#
#   ## find the peptides that are in the filter list but are not in the unfiltered lists ##
#   pep_filt = as.character(setdiff(omicsData$e_meta[protfilt.ids, pep_id], omicsData$e_meta[-protfilt.ids, pep_id]))
#
#   if(length(pep_filt)==0){pep_filt = NULL}
#   if(length(pro_filt)==0){pro_filt = NULL}
#
#   return(list(proteins_filt = pro_filt, peptides_filt = pep_filt))

  output <- list(counts_by_pep = count_bypep, counts_by_pro = count_bypro)

  orig_class <- class(output)
  class(output) <- c("proteomicsFilt", orig_class)

  attr(output, "e_meta") <- omicsData$e_meta

  # attr(output, "min_num_peps") <- min_num_peps

  return(output)
}
