# pre_flight checks all of the data frames, id column names, and other input
# values to ensure they are all the correct class. It also checks the e_data,
# f_data, and e_meta data frames to make sure they are the right dimension.
#
# preflight returns e_data, f_data, e_meta, and emeta_cname. Each of these
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
                        dType) {
  
  # Verify classes/values of arguments -----------------------------------------

  # Check that e_data and f_data are data.frames.
  if (!inherits(e_data, "data.frame")) {
    
    stop ("e_data must be of the class 'data.frame'")
    
  }
  
  if (!inherits(f_data, "data.frame")) {
    
    stop ("f_data must be of the class 'data.frame'") 
    
  }

  # Determine if e_meta is present.
  if (!is.null(e_meta)) {

    # Test that e_meta is a data frame.
    if (!inherits(e_meta, "data.frame")) {
      
      # Lay down an error because e_meta is not a data frame.
      stop ("e_meta must be of the class 'data.frame'")
      
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
  
  # Examine the norm_info argument. Ensure it is a list.
  if (!inherits(norm_info, "list")) {
    
    # Throw an error at the user.
    stop ("norm_info must be of the class 'list'")
    
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
  
  # Make sure data_scale is one of the acceptable strings.
  if (!(data_scale %in% c("abundance", "log", "log2", "log10"))) {
    
    # Throw an error because data_scale is not an acceptable form.
    stop (paste("data_scale must be one of the following:",
                "'abundance', 'log', 'log2', or 'log10'.",
                sep = " "))
    
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
