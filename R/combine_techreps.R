#' Combine technical replicates of an omicsData object
#'
#' For each biomolecule, aggregates the technical replicates within biological
#' samples using a specified aggregation method
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData',
#'   'proData', or 'nmrData', usually created by \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, or \code{\link{as.nmrData}}, respectively.  The
#'   parameter techrep_cnames must have been specified when creating this
#'   object.
#' @param combine_fn a character string specifying the function used to
#'   aggregate across technical replicates, currently only supports 'sum' and 
#'   'mean'. Defaults to 'sum' for seqData and mean for all other omicsData.
#' @param bio_sample_names a character string specifying a column in
#'   \code{f_data} which contains names by which to label aggregated samples in
#'   \code{omicsData$e_data} and \code{omicsData$f_data} OR a character vector
#'   with number of elements equal to the number of biological samples.  If a
#'   column name is specified, it should have a one-to-one correspondence with
#'   the technical replicate ID column in \code{f_data}.  Defaults to NULL, in
#'   which case default names are used according to the technical replicate ID
#'   column.
#'
#' @details Loss of information after aggregation \tabular{ll}{ f_data: \tab
#'   If there are columns of f_data that have more than 1 value per biological
#'   sample, then for each biological sample, only the first value in that
#'   column will be retained.  Technical replicate specific information will be
#'   lost.\cr group information:  \tab If a grouping structure has been set
#'   using a main effect from f_data that has more than 1 level within any given
#'   biological sample, that grouping structure will be removed.  Call
#'   \code{group_designation} again on the aggregated data to assign a grouping
#'   structure. \cr sample names:  \tab Identifiers for each biological sample
#'   will replace the identifiers for technical replicates as column names in
#'   e_data as well as the identifier column \code{attr(omicsData,
#'   'fdata_cname')} in f_data. \cr }
#'
#' @return An object with the same class as omicsData that has been aggregated
#'   to the biological sample level
#'
#' @examples
#' library(pmartR)
#' library(pmartRdata)
#'
#' data(techrep_pep_object)
#'
#' pep_object_averaged = combine_techreps(techrep_pep_object)
#'
#' @author Daniel Claborne
#'
#' @export
#' 
combine_techreps <- function (omicsData, combine_fn = NULL,
                              bio_sample_names = NULL) {
  
  # check that omicsData is of pmartR S3 class#
  if(!inherits(omicsData, 
               c("pepData", "proData", "lipidData", "metabData", "nmrData", "seqData"))
     ) stop("omicsData must be of class 'pepData', 'proData', 'lipidData', 'metabData', 'nmrData', or 'seqData'")
  
  if(is.null(combine_fn)){
    combine_fn <- ifelse(inherits(omicsData, "seqData"), "sum", "mean")
  }
  
  f_data = omicsData$f_data
  e_data = omicsData$e_data
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  edata_cname = attr(omicsData, "cnames")$edata_cname
  techrep_cname = attr(omicsData, "cnames")$techrep_cname
  
  # object had a techrep column specified
  if(is.null(attr(omicsData, "cnames")$techrep_cname)) stop("This object did not have technical replicates specified.  Specify the argument techrep_cname when creating your omicsData object or edit the attribute attr(omicsData, 'cnames')$techrep_cname.")
  
  # isobaric object has been reference sample normalized
  if(inherits(omicsData, "isobaricpepData")){
    if(!attr(res, "isobaric_info")$norm_info$is_normalized) stop("isobaricpepData objects must have been normalized to the appropriate reference pool samples before combining technical replicates")
  }
  
  # legit display name column or vector of display names
  if(!is.null(bio_sample_names)){
    if(!inherits(bio_sample_names, "character") | length(techrep_cname) == 0) stop("bio_sample_names must be a character string specifying a column in f_data")
    if(length(bio_sample_names) == 1){
      if(!(bio_sample_names %in% colnames(f_data[,-which(names(f_data) == fdata_cname)]))) stop("Specified display name column was not found in f_data or was the same as fdata_cname")
      one_to_one <- f_data[c(bio_sample_names, techrep_cname)] %>% unique() %>% nrow() 
      unique_bio_sample_names <- length(unique(f_data[,bio_sample_names]))
      if(any(c(one_to_one, unique_bio_sample_names) != length(unique(f_data[,techrep_cname])))) stop("Specified display name column did not have a one-to-one correspondence with the techrep ID column")
    } else if(length(bio_sample_names) > 1){
      if(length(unique(bio_sample_names)) != length(unique(f_data[,techrep_cname]))) stop("character vector of sample names does not have the same number of names as the number of biological samples")
    }
  }

  # create a list with the biological sample identifiers as names and the corresponding technical replicate names (column names of edata) as values
  # used below in constructing new e_data and also stored as an attribute
  bio_sample_list = list()
  
  for(el in as.character(unique(f_data[,which(names(f_data) == techrep_cname)]))){
    bio_sample_list[[el]] = unique(
      f_data %>% 
        dplyr::filter(!!rlang::sym(techrep_cname) == el) %>%
        {.[,fdata_cname]}) %>% as.character()
  }
  #
  
  if(combine_fn == "sum" && 
     length(unique(sapply(bio_sample_list, length))) != 1){
    stop("Differing number of technical replicates per sample; sum is an invalid combine option.")
  }
  
  ### Do certain columns in f_data have multiple values per biological sample?.....
  
  # get number of distinct levels for each column per biological sample.
  distinct <- f_data %>% 
    dplyr::group_by(!!rlang::sym(techrep_cname)) %>% 
    dplyr::summarise_all(dplyr::n_distinct) %>% 
    dplyr::select(-dplyr::one_of(c(techrep_cname, fdata_cname))) %>% 
    as.data.frame() 
  
  # names of columns with multiple values per bio sample, if any.
  names_mult <- lapply(distinct, function(x){any(x > 1)}) %>% 
    {which(.==TRUE)} %>% 
    names()
  
  ### ....if so, throw a warning saying that some information will be lost in the new f_data.
  if(length(names_mult) > 0) warning(paste0("The following columns in f_data have multiple values across technical replicates in a biological sample:  [", paste(names_mult, collapse = ", "), "] The first value within each biological sample will be retained, the rest discarded."))
  
  # create new, collapsed f_data object
  new_fdata <- f_data %>% 
    dplyr::group_by(!!rlang::sym(techrep_cname)) %>%
    dplyr::slice(1) %>% 
    dplyr::select(!!rlang::sym(techrep_cname), dplyr::everything(), -dplyr::one_of(fdata_cname)) %>%
    as.data.frame()
  
  # create new, collapsed e_data object, averaged over technical replicates
  new_edata <- e_data[get_edata_cname(omicsData)]
  for(el in names(bio_sample_list)){
    edata_subsample <- e_data %>% dplyr::select(dplyr::one_of(bio_sample_list[[el]]))
    if(combine_fn == "mean"){
      new_edata[el] <- rowMeans(edata_subsample, na.rm = TRUE)
    }else if(combine_fn == "sum"){
      new_edata[el] <- rowSums(edata_subsample, na.rm = TRUE)
    }
    # other combine methods coming soon! #
  }
  
  # NaN to NA
  new_edata[is.na(new_edata)] <- NA
  if(inherits(omicsData, "seqData")) new_edata[is.na(new_edata)] <- 0
  
  # Assign column names to e_data + 
  # Assign new ID column to f_data
  if(!is.null(bio_sample_names)){
    # if a column was specified, reassign that column as the id column
    if(length(bio_sample_names) == 1){
      new_fdata <- dplyr::select(new_fdata, bio_sample_names, dplyr::everything())
      attr(omicsData, "cnames")$fdata_cname = bio_sample_names
      bio_sample_names <- new_fdata[,bio_sample_names] 
    }
    else if(length(bio_sample_names) > 1){
      bsn <- bio_sample_names
      names(bsn) <- f_data[,techrep_cname]
      new_fdata[techrep_cname] <- bsn[new_fdata[[techrep_cname]]]
      names(bio_sample_list) <- bsn[names(bio_sample_list)]
      bio_sample_names <- bsn[colnames(new_edata)[-which(colnames(new_edata) == edata_cname)]]
      attr(omicsData, "cnames")$fdata_cname <- techrep_cname
    }
    colnames(new_edata)[-which(colnames(new_edata) == edata_cname)] <- bio_sample_names
  }
  else attr(omicsData, "cnames")$fdata_cname <- techrep_cname
  
  # Check for bad grouping structure and make new grouping DF
  if(!is.null(attr(omicsData, "group_DF"))){ ## needed if done before grouping
    
    # gives number of unique main effect levels in group_DF for a given group of technical replictes...
    multiple_groups <- get_group_DF(omicsData) %>%
      dplyr::left_join(f_data[c(fdata_cname, techrep_cname)], by = fdata_cname) %>%
      dplyr::group_by(!!rlang::sym(techrep_cname)) %>%
      dplyr::summarise_all(dplyr::n_distinct) %>%
      dplyr::select("Group")
    
    # ... if any group of technical replicates spans multiple main effects, throw a warning and discard the grouping structure.
    if(any(multiple_groups > 1)){
      warning("A grouping structure was present that assigned multiple groups within a biological sample, this grouping structure will be discarded, run group_designation() again to set a new grouping structure")
      new_group_DF <- NULL
    }
    # otherwise collapse the grouping structure around the newly created f_data
    else{
      new_group_DF <- get_group_DF(omicsData) %>%
        dplyr::left_join(f_data[c(fdata_cname, techrep_cname)], by = fdata_cname) %>% 
        dplyr::group_by(!!rlang::sym(techrep_cname)) %>%
        dplyr::slice(1) %>%
        dplyr::select(!!rlang::sym(techrep_cname), dplyr::everything(), -dplyr::one_of(fdata_cname)) %>%
        as.data.frame()
        
        colnames(new_group_DF)[which(colnames(new_group_DF) == techrep_cname)] <- get_fdata_cname(omicsData) # this attribute will always have been reset at this point
        if(!is.null(bio_sample_names)) new_group_DF[,get_fdata_cname(omicsData)] <- bio_sample_names # display names is always a vector of values at this point
        
    }
    
  }else new_group_DF <- NULL

  # store new data and reset attributes
  omicsData$e_data <- new_edata
  omicsData$f_data <- new_fdata
  id_col <- which(names(omicsData$e_data)==edata_cname)
  attr(omicsData, "cnames")$techrep_cname <- NULL
  attr(omicsData, "group_DF") <- new_group_DF 
  attr(omicsData, "data_info")$num_samps = ncol(omicsData$e_data) - 1
  attr(omicsData, "data_info")$num_edata = length(unique(omicsData$e_data[, edata_cname]))
  if(inherits(omicsData, "seqData")){
    attr(omicsData, "data_info")$num_zero_obs <- sum(omicsData$e_data[,-id_col] == 0)
    attr(omicsData, "data_info")$prop_zeros  <- mean(omicsData$e_data[,-id_col] == 0)
  } else {
    attr(omicsData, "data_info")$num_miss_obs = sum(is.na(omicsData$e_data[,-id_col]))
    attr(omicsData, "data_info")$prop_missing = mean(is.na(omicsData$e_data[,-id_col]))
  }

  # technical replicate specific attributes
  attributes(omicsData)$tech_rep_info$tech_reps_by_sample = bio_sample_list
  attributes(omicsData)$tech_rep_info$combine_method = combine_fn
  
  return(omicsData)
  
}
