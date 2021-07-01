#' Normalize an object of class isobaricpepData
#'
#' The samples are normalized to their corresponding reference pool sample
#'
#' @param omicsData an object of the class 'isobaricpepData'
#' @param apply_norm logical, indicates whether normalization should be applied
#'   to omicsData$e_data
#' @param exp_cname character string specifying the name of the column
#'   containing the experiment/plate information in \code{f_data}.
#' @param channel_cname optional character string specifying the name of the
#'   column containing the instrument channel a sample was run on in
#'   \code{f_data}. This argument is optional, see Details for how to specify
#'   information regarding reference pool samples. If using this argument, the
#'   'refpool_channel' argument must also be specified; in this case,
#'   'refpool_cname' and 'refpool_notation' should not be specified.
#' @param refpool_channel optional character string specifying which channel
#'   contained the reference pool sample, only used when this remains the same
#'   from experiment to experiment. This argument is optional, see Details for
#'   how to specify information regarding reference pool samples. If using this
#'   argument, the 'channel_cname' argument must also be specified; in this
#'   case, 'refpool_cname' and 'refpool_notation' should not be specified.
#' @param refpool_cname optional character string specifying the name of the
#'   column containing information about which samples are reference samples in
#'   \code{f_data}. This argument is optional, see Details for how to specify
#'   information regarding reference pool samples. If using this argument, the
#'   'refpool_notation' argument must also be specified; in this case,
#'   'channel_cname' and 'refpool_channel' should not be specified.
#' @param refpool_notation optional character string specifying the value in the
#'   refpool_channel column which denotes that a sample is a reference sample.
#'   This argument is optional, see Details for how to specify information
#'   regarding reference pool samples. If using this argument, the
#'   'refpool_cname' argument must also be specified; in this case,
#'   'channel_cname' and 'refpool_channel' should not be specified.
#' 
#' @details There are two ways to specify the information needed for identifying
#'   reference samples which should be used for normalization: \enumerate{ \item
#'   specify \code{channel_cname} and \code{refpool_channel}. This should be
#'   used when the reference sample for each experiment/plate was always located
#'   in the same channel. Here \code{channel_cname} gives the column name for
#'   the column in \code{f_data} which gives information about which channel
#'   each sample was run on, and \code{refpool_channel} is a character string
#'   specifying the value in \code{channel_colname} that corresponds to the
#'   reference sample channel. \item specify \code{refpool_cname} and
#'   \code{refpool_notation}. This should be used when the reference sample is
#'   not in a consistent channel across experiments/plates. Here,
#'   \code{refpool_cname} gives the name of the column in \code{f_data} which
#'   indicates whether a sample is a reference or not, and
#'   \code{refpool_notation} is a character string giving the value used to
#'   denote a reference sample in that column. } In both cases you must specify
#'   \code{exp_cname} which gives the column name for the column in
#'   \code{f_data} containing information about which experiment/plate a sample
#'   was run on.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(isobaric_object)
#'
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set",
#'                                    apply_norm = TRUE,
#'                                    channel_cname = "iTRAQ.Channel",
#'                                    refpool_channel = "116")
#'
#' # alternate specification: #
#' data(isobaric_object)
#'
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set",
#'                                    apply_norm = TRUE,
#'                                    refpool_cname = "Reference",
#'                                    refpool_notation = "Yes")
#'
#' }
#'
#' @export
#' 
normalize_isobaric <- function (omicsData, exp_cname = NULL, apply_norm = FALSE,
                                channel_cname = NULL, refpool_channel = NULL,
                                refpool_cname = NULL, refpool_notation = NULL) {
  
  # initial checks #
  
  #check that omicsData is of correct class
  if (!inherits(omicsData, "isobaricpepData")) {
    
    stop ("omicsData must be of the class 'isobaricpepData'")
    
  }
  
  # check that the data has not already been isobaric normalized #
  if (attr(omicsData, "isobaric_info")$norm_info$is_normalized == TRUE) {
    
    stop (paste("omicsData is already normalized with respect to the isobaric",
                "labels, per attributes assigned upon data object creation",
                sep = " "))
    
  }
  
  #check that omicsData$e_data is log transformed
  if (!(get_data_scale(omicsData) %in% c('log2', 'log10', 'log'))) {
    
    stop ("omicsData$e_data must be log transformed")
    
  }
  
  # check that exp_cname is in f_data #
  if (!(exp_cname %in% names(omicsData$f_data))) {
    
    stop (paste("Experiment column", exp_cname, "is not found in f_data.",
                sep = " "))
    
  }
  
  #check that apply_norm is of class logical
  if (!is.logical(apply_norm)) stop ("apply_norm must be of class 'logical'")
  
  # check that channel_cname is in f_data, if not NULL #
  if (!is.null(channel_cname)) {
    
    if (!(channel_cname %in% names(omicsData$f_data))) {
      
      stop (paste("Channel column", channel_cname,
                  "is not found in f_data. See details of as.isobaricpepData",
                  "for specifying column names.",
                  sep = " "))
      
    }
    
  }
  
  # check that refpool_cname is in f_data, if not NULL #
  if (!is.null(refpool_cname)) {
    
    if (!(refpool_cname %in% names(omicsData$f_data))) {
      
      stop(paste("Reference pool column", refpool_cname,
                 "is not found in f_data. See details of as.isobaricpepData",
                 "for specifying column names.",
                 sep = " "))
      
    }
    
  }
  
  # make sure the reference pool info is specified appropriately #
  # possibility 1: specify refpool_cname and refpool_notation #
  poss1 = !is.null(refpool_cname) & !is.null(refpool_notation)
  # possibility 2: specify refpool_channel and channel_cname#
  poss2 = !is.null(refpool_channel) & !is.null(channel_cname) 
  
  # throw an error if neither or both of these are true #
  if ((poss1 + poss2) != 1) {
    
    stop (paste("Reference samples information was not correctly specified.",
                "See Details and Examples for more information.",
                sep = " "))
    
  }
  
  # Pluck out the name of the column in e_data with the peptide IDs.
  edata_cname <- get_edata_cname(omicsData)
  
  # Fish out the name of the column in f_data containing the sample names.
  fdata_cname <- get_fdata_cname(omicsData)
  
  ### ### ### ### ###
  # Convert columns of f_data to character vectors because all inputs for
  # distinguishing reference samples are character strings.
  ### ### ### ### ###
  
  # Convert the column in f_data containing the sample names into a character
  # vector.
  omicsData$f_data[ , fdata_cname] <- as.character(
    omicsData$f_data[ , fdata_cname]
  )
  
  # Prepare possibility 1 info -------------------------------------------------
  
  # if possibility 1 is used, check that refpool_cname is a value seen in each
  # experiment.
  if (poss1 == TRUE) {
    
    # Make sure refpool_notation is a character string.
    if (!is.character(refpool_notation)) {
      
      # Throw an error for being out of character!
      stop ("refpool_notation must be a character string")
      
    }
    
    # Mutate the column indicating whether a sample is a reference sample or not
    # into a character vector.
    omicsData$f_data[, refpool_cname] <- as.character(
      omicsData$f_data[, refpool_cname]
    )
    
    # Divide the sample names by experiment.
    idx = split(as.character(omicsData$f_data[ , refpool_cname]),
                omicsData$f_data[ , exp_cname])
    
    # Check that there is a reference sample in each experiment.
    temp_check = lapply(idx, function(x) refpool_notation %in% x)
    if (sum(unlist(temp_check)) != length(temp_check)) {
      
      stop(paste("refpool_notation =", refpool_notation,
                 "is not in every experiment. See Details and Examples for",
                 "more information.",
                 sep = " "))
      
    }
    
  }
  
  # Prepare possibility 2 info -------------------------------------------------
  
  # if possibility 2 is used, check that refpool_channel is a value seen in each
  # experiment.
  if (poss2 == TRUE) {
    
    # Ensure refpool_channel is a character string.
    if (!is.character(refpool_channel)) {
      
      # Throw an error for being out of character!
      stop ("refpool_channel must be a character string")
      
    }
    
    # Metamorphose the column in f_data containing the channels into a character
    # vector.
    omicsData$f_data[, channel_cname] <- as.character(
      omicsData$f_data[, channel_cname]
    )
    
    # Divide the sample names by experiment.
    idx = split(omicsData$f_data[ , channel_cname],
                omicsData$f_data[ , exp_cname])
    
    # Ensure there is a reference sample in each experiment.
    temp_check = lapply(idx, function(x) refpool_channel %in% x)
    if (sum(unlist(temp_check)) != length(temp_check)) {
      
      stop (paste("refpool_channel =", refpool_channel,
                  "is not in every experiment. See Details and Examples for",
                  "more information.",
                  sep = " "))
      
    }
    
  }
  
  # Set reference column IDs and sample names ----------------------------------
  
  # Set the reference column and reference sample names according to the user's
  # input.
  if (!is.null(refpool_channel) && !is.null(channel_cname)) {
    
    # Set the name of the reference column.
    ref_col <- channel_cname
    
    # Set the name of the reference samples
    ref_name <- refpool_channel
    
  } else if (!is.null(refpool_cname) && !is.null(refpool_notation)) {
    
    # Set the name of the reference column.
    ref_col <- refpool_cname
    
    # Set the name of the reference samples
    ref_name <- refpool_notation
    
  }
  
  # Carry out normalization/create isobaricnormRes object ----------------------
  
  #case where apply_norm is TRUE
  if (apply_norm == TRUE) {
    
    # Normalize the data to the isobaric reference samples.
    omicsData <- isonorm(omicsData = omicsData,
                         exp_cname = exp_cname,
                         ref_col = ref_col,
                         ref_name = ref_name)
    
    # Include information for the isobaric_info attribute.
    attr(omicsData, "isobaric_info") <- set_isobaric_info(
      exp_cname = exp_cname, 
      channel_cname = channel_cname, 
      refpool_channel = refpool_channel, 
      refpool_cname = refpool_cname, 
      refpool_notation = refpool_notation,
      norm_info = list(),
      isobaric_norm = TRUE
    )
    
    # Update the data_info attribute because the reference samples have been
    # removed from e_data and f_data.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale_orig = get_data_scale_orig(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = get_data_info(omicsData)$data_types,
      norm_info = get_data_info(omicsData)$norm_info,
      is_normalized = get_data_info(omicsData)$norm_info$is_normalized
    )
    
    # Return the normalized omicsData object along with its updated attributes.
    return (omicsData)
    
    # Runs when apply_norm is false.
  } else {
    
    # Obtain the sample names corresponding to the reference samples.
    rfrnc_nms <- omicsData$f_data[
      omicsData$f_data[, ref_col] == ref_name, fdata_cname
    ]
    
    # Subset the columns of e_data pertaining to the reference samples.
    edata <- omicsData$e_data[, which(
      names(omicsData$e_data) %in% c(edata_cname, rfrnc_nms)
    )]
    
    # Subset the rows of f_data corresponding to the reference samples.
    fdata <- omicsData$f_data[which(
      omicsData$f_data[, fdata_cname] %in% rfrnc_nms
    ), ]
    
    # Create a list containing the columns and rows from e_data and f_data
    # that correspond to the reference samples.
    result <- list(e_data = edata,
                   f_data = fdata)
    
    # Set the class of the reference samples.
    class(result) = "isobaricnormRes"
    
    # Add helpful attributes to the isobaricnormRes object.
    attr(result, "cnames") <- list(edata_cname = edata_cname,
                                   fdata_cname = fdata_cname)
    attr(result, "isobaric_info") <- set_isobaric_info(
      exp_cname = exp_cname, 
      channel_cname = channel_cname, 
      refpool_channel = refpool_channel, 
      refpool_cname = refpool_cname, 
      refpool_notation = refpool_notation,
      norm_info = list(),
      isobaric_norm = FALSE
    )
    
    # Return the isobaricnormRes object!!!
    return(result)
    
  }
  
}

# Carries out the normalization for isobaric data. Within each experiment each
# sample is normalized to the corresponding reference sample.
#
# @param omicsData An isobaricpepData object.
#
# @param exp_cname A character string specifying the name of the column
#   containing the experiment/plate information in \code{f_data}.
# @param ref_col A character string indicating which column in f_data contains
#   the information identifying reference samples.
#
# @param ref_name A character string identifying the reference sample.
#
# @return An isobaricpepData object. The data in e_data has been normalized to
#   the reference samples and all columns corresponding to the reference samples
#   have been removed. The rows in f_data corresponding to the reference samples
#   have also been removed.
#   
isonorm <- function (omicsData, exp_cname,
                     ref_col, ref_name) {
  
  # Find the name of the sample column in f_data
  f_name <- get_fdata_cname(omicsData)
  
  # Fish out all experiment levels.
  xprmnt <- unique(omicsData$f_data[, exp_cname])
  
  # Create a list that will hold the sample names for each experiment. It is 
  # named no_ref even though at this point it has the names of the reference
  # samples. Later in the function the reference sample names will be removed
  # from this list.
  no_ref <- vector(mode = "list",
                   length = length(xprmnt))
  
  # Forge a vector that houses the sample names for the reference samples in
  # each experiment.
  ref <- vector(length = length(xprmnt))
  
  # Loop through each experiment and extract all samples belonging to that
  # experiment and normalizing the data in an isobaric fashion.
  for (e in 1:length(xprmnt)) {
    
    # Convert the sample names to a character string. This will remove any
    # factors if the sample column is a factor.
    no_ref[[e]] <- as.character(
      omicsData$f_data[omicsData$f_data[, exp_cname] == xprmnt[[e]], f_name] 
    )
    
    # Convert the sample names for the reference samples to a character string.
    # The reference sample name is extracted by experiment and name of the
    # reference sample (either channel name or reference pool name).
    ref[[e]] <- as.character(
      omicsData$f_data[(omicsData$f_data[, exp_cname] == xprmnt[[e]] &
                          omicsData$f_data[, ref_col] == ref_name), f_name]
    )
    
    # Pluck out the indices for the samples in the eth experiment.
    idx_no_ref <- which(names(omicsData$e_data) %in% no_ref[[e]])
    
    # Go fishin for the index of the reference sample for the eth experiment.
    idx_ref <- which(names(omicsData$e_data) %in% ref[[e]])
    
    # Separate the reference sample from the others for the eth experiment.
    idx_no_ref <- idx_no_ref[idx_no_ref != idx_ref]
    
    # Isobaricate the heck out of the data!
    # Divide the non-reference samples by the reference sample. Since we are on
    # the log scale we will subtract the reference sample from the non-reference
    # samples.
    omicsData$e_data[, idx_no_ref] <- (omicsData$e_data[, idx_no_ref] -
                                         omicsData$e_data[, idx_ref])
    
  }
  
  # Remove all columns and rows from e_data and f_data that correspond to the
  # reference samples.
  omicsData$e_data <- omicsData$e_data[, -which(
    names(omicsData$e_data) %in% ref
  )]
  omicsData$f_data <- omicsData$f_data[-which(
    omicsData$f_data[, f_name] %in% ref
  ), ]
  
  # Return the normalized omicsData object! Woot!!
  return (omicsData)
  
}
