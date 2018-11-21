#' Combine technical replicates of an omicsData object
#'
#' For each biomolecule, aggregates the technical replicates within biological samples using a specified aggregation method
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData', or 'proData', usually created by \code{\link{as.lipidData}}, \code{\link{as.metabData}}, \code{\link{as.pepData}}, or \code{\link{as.proData}}, respectively.  The parameter techrep_cnames must have been specified when creating this object.
#' @param combine_fn a character string specifying the function used to aggregate across technical replicates, currently only supports "mean".
#' 
#' @details Loss of information after aggregation
#' \tabular{ll}{
#' f_data: \tab   If there are columns of f_data that have more than 1 value per biological sample, then for each biological sample, only the first value in that column will be retained.  Technical replicate specific information will be lost.\cr
#' group information:  \tab If a grouping structure has been set using a main effect from f_data that has more than 1 level within any given biological sample, that grouping structure will be removed.  Call \code{group_designation} again on the aggregated data to assign a grouping structure.
#' }
#' 
#' @return An object with the same class as omicsData that has been aggregated to the biological sample level
#' 
#' @examples 
#' # UNDER CONSTRUCTION #
#' 
#' @author Daniel Claborne
#' 
#' @export
combine_techreps <- function(omicsData, combine_fn = "mean"){
  if(is.null(attr(omicsData, "cnames")$techrep_cname)) stop("This object did not have technical replicates specified.  Specify the argument techrep_cname when creating your omicsData object.")
  if(inherits(omicsData, "isobaricpepData")){
    if(!attr(res, "data_info")$isobaric_norm) stop("isobaricpepData objects must have been normalized to the appropriate reference pool samples before combining technical replicates")
  }
  
  f_data = omicsData$f_data
  e_data = omicsData$e_data
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  edata_cname = attr(omicsData, "cnames")$edata_cname
  techrep_cname = attr(omicsData, "cnames")$techrep_cname
  
  # create a list with the biological sample identifiers as names and the corresponding technical replicate names (column names of edata) as values
  # used below in constructing new e_data and also stored as an attribute
  bio_sample_list = list()
  
  for(el in as.character(unique(f_data[,which(names(f_data) == techrep_cname)]))){
    bio_sample_list[[el]] = unique(f_data %>% dplyr::filter(!!rlang::sym(techrep_cname) == el) %>% {.[,fdata_cname]}) %>% as.character()
  }
  #
  
  ### Do certain columns in f_data have multiple values per biological sample?.....
  
  # distinct per biological sample
  distinct <- fdata %>% 
    dplyr::group_by(!!rlang::sym(techrep_cname)) %>% 
    dplyr::summarise_all(n_distinct) %>% 
    dplyr::select(-dplyr::one_of(c(techrep_cname, fdata_cname))) %>% 
    as.data.frame() 
  
  # names of columns with multiple values per bio sample.
  names_mult <- lapply(distinct, function(x){any(x > 1)}) %>% 
    {which(.==TRUE)} %>% 
    names()
  
  #### ....if so, throw a warning saying that some information will be lost in the new f_data.
  if(length(multiple) > 0) warning(paste0("The following columns in f_data have multiple values across technical replicates in a biological sample:  [", paste(names_mult, collapse = ", "), "] The first value within each biological sample will be retained, the rest discarded."))
  
  # create new, collapsed f_data object
  new_fdata <- f_data %>% 
    dplyr::group_by(!!rlang::sym(techrep_cname)) %>%
    dplyr::slice(1) %>% 
    dplyr::select(!!rlang::sym(techrep_cname), dplyr::everything(), -dplyr::one_of(fdata_cname)) %>%
    as.data.frame()
  
  # create new, collapsed e_data object, averaged over technical replicates
  new_edata <- e_data[which(names(e_data) == attr(omicsData, "cnames")$edata_cname)]
  for(el in names(bio_sample_list)){
    edata_subsample <- e_data %>% dplyr::select(one_of(bio_sample_list[[el]]))
    if(combine_fn == "mean") new_edata[el] = rowMeans(edata_subsample, na.rm = TRUE)  
  }
  
  # NaN to NA
  new_edata[is.na(new_edata)] = NA
  
  # Check for bad grouping structure and make new grouping DF
  if(!is.null(attr(omicsData, "group_DF"))){
    
    # gives number of unique main effect levels in group_DF for a given group of technical replictes...
    multiple_groups <- attr(omicsData, "group_DF") %>% 
      dplyr::left_join(fdata[c(fdata_cname, techrep_cname)], by = fdata_cname) %>%
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
      new_group_DF <- attr(omicsData, "group_DF") %>%
        dplyr::left_join(fdata[c(fdata_cname, techrep_cname)], by = fdata_cname) %>% 
        dplyr::group_by(!!rlang::sym(techrep_cname)) %>%
        dplyr::slice(1) %>%
        dplyr::select(!!rlang::sym(techrep_cname), dplyr::everything(), -dplyr::one_of(fdata_cname)) %>%
        as.data.frame()
    }
  }else new_group_DF <- NULL
  
  omicsData$e_data <- new_edata
  omicsData$f_data <- new_fdata
  attr(omicsData, "group_DF") <- new_group_DF 
  attr(omicsData, "cnames")$fdata_cname = techrep_cname
  attr(omicsData, "data_info")$num_samps = ncol(omicsData$e_data) - 1
  attr(omicsData, "data_info")$num_edata = length(unique(omicsData$e_data[, edata_cname]))
  attr(omicsData, "data_info")$num_miss_obs = sum(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
  attr(omicsData, "data_info")$prop_missing = mean(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
  
  # technical replicate specific attributes
  attributes(omicsData)$tech_rep_info$tech_reps_by_sample = bio_sample_list
  attributes(omicsData)$tech_rep_info$combine_method = combine_fn
  
  return(omicsData)
  
}

