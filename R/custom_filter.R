#' Custom Filter
#'
#' This function creates a customFilt S3 object based on user-specified items to
#' filter out of the dataset
#'
#' @param omicsData an object of class "pepData", "proData", "metabData",
#' "lipidData", or "nmrData, created by \code{\link{as.pepData}},
#' \code{\link{as.proData}}, \code{\link{as.metabData}},
#' \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#' 
#' @param e_data_remove character vector specifying the names of the e_data
#' identifiers to remove from the data. This argument can only be specified with
#' other 'remove' arguments.
#' 
#' @param f_data_remove character vector specifying the names of f_data
#' identifiers to remove from the data. This argument can only be specified with
#' other 'remove' arguments.
#' 
#' @param e_meta_remove character vector specifying the names of the e_meta
#' identifiers to remove from the data. This argument can only be specified with
#' other 'remove' arguments.
#' 
#' @param e_data_keep character vector specifying the names of the e_data
#' identifiers to keep from the data. This argument can only be specified with
#' other 'keep' arguments.
#' 
#' @param f_data_keep character vector specifying the names of f_data
#' identifiers to keep from the data. This argument can only be specified with
#' other 'keep' arguments.
#' 
#' @param e_meta_keep character vector specifying the names of the e_meta
#' identifiers to keep from the data. This argument can only be specified with
#' other 'keep' arguments.
#'
#' @return An S3 object of class 'customFilt', which is a list with 3 elements:
#' e_data_remove, f_data_remove, and e_meta_remove.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("metab_object")
#' to_filter <- custom_filter(metab_object, e_data_remove = "fumaric acid",
#'                            f_data_remove = "Infection1")
#' summary(to_filter)
#' to_filter2 <- custom_filter(metab_object, e_data_remove = "fumaric acid")
#' summary(to_filter2)
#'}
#'
#' @author Kelly Stratton
#'
#' @export
#' 
custom_filter <- function (omicsData,
                           e_data_remove = NULL,
                           f_data_remove = NULL,
                           e_meta_remove = NULL,
                           e_data_keep = NULL,
                           f_data_keep = NULL,
                           e_meta_keep = NULL ) {
  
  # Run some preliminary checks ------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    # Follow the instructions foul creature!!!
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = " "))
    
  }
  
  # check that not all remove and keep arguments are NULL.
  if (is.null(c(e_data_remove, f_data_remove, e_meta_remove,
                e_data_keep, f_data_keep, e_meta_keep))) {
    
    # Stop the user because they are trying to take the long way to the exact
    # same data set.
    stop ("No items have been identified for filtering.")
    
  }
  
  #check that both keep and remove arguments are not non-NULL
  if (!is.null(c(e_data_remove, f_data_remove, e_meta_remove)) && 
      !is.null(c(e_data_keep, f_data_keep, e_meta_keep))) {
    
    stop (paste("Cannot have both remove arguments and keep arguments",
                "be non-NULL. Create separate filter objects for the",
                "remove arguments and keep arguments.",
                sep = " "))
    
  }
  
  # Extricate the names of the different ID columns.
  edata_id = attr(omicsData, "cnames")$edata_cname
  emeta_id = attr(omicsData, "cnames")$emeta_cname
  samp_id = attr(omicsData, "cnames")$fdata_cname
  
  # Create filter_object for remove arguments ----------------------------------
  
  if (!is.null(c(e_data_remove, f_data_remove, e_meta_remove))) {
    
    # checks for e_data_remove #
    if (!is.null(e_data_remove)) {
      
      # check that e_data_remove are all in omicsData #
      if (!(all(e_data_remove %in% omicsData$e_data[, edata_id]))) {
        
        # Throw an error because the user tried to filter nonexistent
        # biomolecules.
        stop ("Not all of the items in e_data_remove are found in e_data.")
        
      }
      
      # check that e_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_data[, edata_id] %in% e_data_remove)) {
        
        # Throw an error because the user tried to filter every biomolecule.
        stop("e_data_remove specifies all the items in the data.")
        
      }
      
    }
      
    # checks for f_data_remove #
    if (!is.null(f_data_remove)) {
      
      # check that f_data_remove are all in omicsData #
      if (!(all(f_data_remove %in% omicsData$f_data[, samp_id]))) {
        
        # Throw an error because the user tried to filter nonexistent samples.
        stop ("Not all of the items in f_data_remove are found in f_data.")
        
      }
      
      # check that f_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$f_data[, samp_id] %in% f_data_remove)) {
        
        # Throw an error because the user tried to filter every sample.
        stop ("f_data_remove specifies all the items in f_data.")
        
      }
      
    }
    
    # checks for e_meta_remove #
    if (!is.null(e_meta_remove)) {
      
      # check that e_meta_remove are all in omicsData #
      if (!(all(e_meta_remove %in% omicsData$e_meta[, emeta_id]))) {
        
        # Throw an error because the user tried to filter nonexistent mapping
        # variables.
        stop ("Not all of the items in e_meta_remove are found in e_meta.")
        
      }
      
      # check that e_meta_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_meta[, emeta_id] %in% e_meta_remove)) {
        
        # Throw an error because the user tried to filter every single mapping
        # variable. They must be tired of their work and they want to be done
        # as soon as possible. Not having any data to analyze is a very
        # effective way of finishing quickly.
        stop ("e_meta_remove specifies all the items in e_meta.")
        
      }
    }
    
    # Fashion a filter object with the remove elements.
    filter_object <- list(e_data_remove = e_data_remove,
                          f_data_remove = f_data_remove,
                          e_meta_remove = e_meta_remove)
    
  }
  
  # Create filter_object for keep arguments ----------------------------------
  
  if (!is.null(c(e_data_keep, f_data_keep, e_meta_keep))) {
    
    # checks for e_data_keep #
    if(!is.null(e_data_keep)){
      
      # check that e_data_keep are all in omicsData #
      if (!(all(e_data_keep %in% omicsData$e_data[, edata_id]))) {
        
        # Stop the greedy user from trying to keep things that aren't theirs!!
        stop ("Not all of the items in e_data_keep are found in e_data.")
        
      }
      
      # check that e_data_keep doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_data[, edata_id] %in% e_data_keep)) {
        
        # Let the greedy user know that they are keeping everything. Why run the
        # filter in the first place!?!
        stop ("e_data_keep specifies all the items in e_data.")
        
      }
      
    }
    
    # checks for f_data_keep #
    if (!is.null(f_data_keep)) {
      
      # check that f_data_keep are all in omicsData #
      if (!(all(f_data_keep %in% omicsData$f_data[, samp_id]))) {
        
        # Stop the greedy user from trying to keep samples that aren't theirs!
        stop ("Not all of the items in f_data_keep are found in f_data.")
        
      }
      
      # check that f_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$f_data[, samp_id] %in% f_data_keep)) {
        
        # Stop the greedy user from keeping all of the samples. Why filter?!
        stop ("f_data_keep specifies all the items in f_data.")
        
      }
      
    }
    
    # checks for e_meta_remove #
    if (!is.null(e_meta_keep)) {
      
      # check that e_meta_keep are all in omicsData #
      if (!(all(e_meta_keep %in% omicsData$e_meta[, emeta_id]))) {
        
        # Stop the greedy user from trying to keep mapping variables that do not
        # belong to them.
        stop ("Not all of the items in e_meta_keep are found in e_meta.")
        
      }
      
      # check that e_meta_keep doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_meta[, emeta_id] %in% e_meta_keep)) {
        
        # Stop the greedy user from keeping ALL of the mapping variables.
        stop ("e_meta_keep specifies all the items in e_meta")
        
      }
      
    }
    
    # Manufacture a filter object with the keep elements.
    filter_object <- list(e_data_keep = e_data_keep,
                          f_data_keep = f_data_keep,
                          e_meta_keep = e_meta_keep)

  } 
  
  
  # Set filter_object class and attributes -------------------------------------
  
  class(filter_object) <- c("customFilt", "list")
  
  # Save the counts of the biomolecules, samples, and mapping variables.
  attr(filter_object,
       "num_samples") <- length(unique(omicsData$f_data[, samp_id]))
  attr(filter_object,
       "num_edata") <-  length(unique(omicsData$e_data[, edata_id]))
  attr(filter_object,
       "num_emeta") <- if (!is.null(emeta_id)) {
   
    length(unique(omicsData$e_meta[, emeta_id])) 
    
  }
  
  # Save the ID column names. These attributes are used in the
  # summary.customFilt method.
  attr(filter_object, "cnames")$edata_cname = edata_id
  attr(filter_object, "cnames")$emeta_cname = emeta_id
  attr(filter_object, "cnames")$fdata_cname = samp_id
  
  # Save the entire omicsData object as an attribute. This is used in the
  # summary.customFilt method.
  attr(filter_object, "omicsData") = omicsData # added 12/5/2017 by KS #
  
  # Return the customated filter object. Good on us!!!
  return (filter_object)
  
}
