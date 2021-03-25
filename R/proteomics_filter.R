#' Proteomics filter object
#'
#' This function counts the number of peptides that map to each protein and/or
#' the number of proteins to which each individual peptide maps.
#'
#' @param omicsData an object of class "pepData", the a result of
#'        \code{\link{as.pepData}}. The e_meta component of omicsData must be
#'        nonempty.
#'
#' @return A list with two elements. The first element is a data frame of counts 
#'         for each unique peptide. The second element is also a data frame.
#'         This data frame contains the counts for the number of peptides that
#'         map to each unique protein.
#'
#' @examples
#' dontrun{
#' library(pmartR)
#' data("pep_object")
#' my_filter <- proteomics_filter(omicsData = pep_object)
#' summary(my_filter, min_num_peps = 3)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton
#' 
#' @importFrom dplyr group_by tally
#'
#' @export

proteomics_filter <- function(omicsData){

  # Preliminary checks and setup -----------------------------------------------

  # check that omicsData is of class 'pepData' #
  if (!inherits(omicsData, "pepData")) {
    
    # Let the user know that if the data does not contain peptides or proteins
    # you CANNOT APPLY A FILTER USING PEPTIDES AND PROTEINS. What is going on
    # upstairs!?!
    stop ("omicsData must be of class 'pepData'") 
    
  }

  # check that e_meta is not NULL #
  if(is.null(omicsData$e_meta)) stop("e_meta must be non-NULL")

  # get peptide and protein column names #
  pep_id = attr(omicsData, "cnames")$edata_cname
  pro_id = attr(omicsData, "cnames")$emeta_cname

  # check that peptide and protein column names are non-null #
  if(is.null(pep_id)) stop("Peptide column name is NULL")
  if(is.null(pro_id)) stop("Protein column name is NULL")
  
  # Count peptides and proteins ------------------------------------------------

  # Count the number of rows for each peptide in e_meta.
  pepCount <- omicsData$e_meta %>%
    # Group the e_meta data frame by the peptide identifiers in the ID column.
    group_by(.data[[pep_id]]) %>%
    # Count the number of times each peptide identifier occurs (the number of 
    # rows an identifier appears in).
    tally() %>%
    # Convert pepCount to a data frame from a tibble.
    data.frame()

  # Count the number of peptides associated with each protein.
  proCount <- omicsData$e_meta %>%
    # Group the e_meta data frame by the protein identifiers in the protein ID
    # column.
    group_by(.data[[pro_id]]) %>%
    # Count the number of times each protein identifier occurs (the number of 
    # rows an identifier appears in).
    tally() %>%
    # Convert proCount to a data frame from a tibble.
    data.frame()

  # Generate a list containing the data frames for the two counts.
  output <- list(counts_by_pep = pepCount,
                 counts_by_pro = proCount)

  # Preserve the list class.
  orig_class <- class(output)
  
  # Incorporate the proteomicsFilt class along with the list class.
  class(output) <- c("proteomicsFilt", orig_class)

  # Return the list of data frames containing peptide and protein counts.
  # We can count!!
  return(output)
  
}
