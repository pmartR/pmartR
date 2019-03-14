#' Produce a basic summary of a pmartR filter S3 Object
#'
#' This function will provide basic summary statistics for filter objects from the pmartR package.
#'
#' @param filter_object S3 object of class 'moleculeFilt', 'proteomicsFilt', 'imdanovaFilt', 'rmdFilt', 'cvFilt', or 'customFilt', created by \code{\link{molecule_filter}}, \code{\link{proteomics_filter}}, \code{\link{imdanova_filter}}, \code{\link{rmd_filter}}, \code{\link{cv_filter}}, or \code{\link{custom_filter}} respectively.
#'
#' @return a summary table for the pmartR filter object. If assigned to a variable, the elements of the summary table are saved in a list format.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' myfilter <- molecule_filter(pep_object)
#' summary(myfilter)
#' summary(myfilter, min_num = 2)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton, Thomas Johansen
#'
#' @export


#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
#'@param min_num an integer value specifying the minimum number of times each feature must be observed across all samples. Default value is 2.
summary.moleculeFilt <- function(filter_object, min_num=NULL){
  
  if(!is.null(min_num)) {
    # check that min_num is not a vector #
    if(length(min_num) > 1) stop("min_num must be of length 1")
    # check that min_num is numeric >= 0 #
    if(!is.numeric(min_num) | min_num < 0) stop("min_num must be an integer >= 0")
    # check that min_num is an integer #
    if(min_num %% 1 != 0) stop("min_num must be an integer >= 0")
    # check that min_num is less than the max number of observations #
    if(min_num > max(filter_object$Num_Observations)) stop("min_num cannot be greater than the number of samples")
  }
  
  # return the numeric version of plot, the threshold used, the number that would be tested and the number that would not be tested
  
  # how many peptides appear in the dataset once, twice, 3 times, etc.
  cut_data <- table(cut(filter_object$Num_Observations, breaks = -1:max(filter_object$Num_Observations)))
  cumcounts <- cumsum(cut_data)
  pep_observation_counts <- data.frame(num_observations=0:(length(cumcounts)-1), frequency_counts=cumcounts)
  
  if(!is.null(min_num)){
    # get number molecules tested
    num_not_tested <- pep_observation_counts$frequency_counts[pep_observation_counts$num_observations==(min_num-1)]
    
    # get number molecules not tested
    num_tested <- pep_observation_counts$frequency_counts[nrow(pep_observation_counts)]-num_not_tested
    
  }else{
    num_tested = "NULL"
    num_not_tested = "NULL"
    min_num = "NULL"
  }
  
  res <- list(pep_observation_counts=pep_observation_counts[-1,], min_num=min_num, num_not_filtered=num_tested, num_filtered=num_not_tested)
  
  
  # create output #
  cat("\nSummary of Molecule Filter \n----------------------------------")
  catmat <- data.frame(c("Minimum Number:"=min_num, "Filtered:"=res$num_filtered, "Not Filtered:"=res$num_not_filtered, " "="","Molecule Observation Counts:"=""))
  colnames(catmat) <- NULL
  cat(capture.output(print(catmat, row.names = TRUE)), sep="\n")
  cat("\n")
  cat(capture.output(print(res$pep_observation_counts, row.names = FALSE)), sep = "\n")
  
  return(invisible(res))
}




#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
#'@param min_num_peps an optional integer value between 1 and the maximum number of peptides that map to a protein in the data. The value specifies the minimum number of peptides that must map to a protein. Any protein with less than \code{min_num_peps} mapping to it will be returned as a protein that should be filtered. Default value is NULL.
#'@param degen_peps logical indicator of whether to filter out degenerate peptides (TRUE) or not (FALSE). Default value is FALSE.
summary.proteomicsFilt <- function(filter_object, min_num_peps=NULL, degen_peps=FALSE){
  
  # error checks for min_num_peps, if not NULL #
  if(!is.null(min_num_peps)) {
    # check that min_num_peps is not a vector #
    if(length(min_num_peps) > 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is numeric and >=1 #
    if(!is.numeric(min_num_peps) | min_num_peps < 1) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if(min_num_peps %% 1 != 0) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if(length(min_num_peps) != 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if(min_num_peps > max(max(filter_object$counts_by_pep$n),max(filter_object$counts_by_pro$n))) stop("min_num_peps cannot be greater than the total number of peptides")
  }
  # check that degen_peps is logical #
  if(!is.logical(degen_peps)) stop("degen_peps must be either TRUE or FALSE")
  
  
  pep_sum <- summary(filter_object$counts_by_pep$n)
  pro_sum <- summary(filter_object$counts_by_pro$n)
  
  if(!is.null(min_num_peps)){
    count_bypro <- filter_object$counts_by_pro
    count_bypep <- filter_object$counts_by_pep
    
    pro_id <- names(count_bypro)[names(count_bypro) != "n"]
    pep_id <- names(count_bypep)[names(count_bypep) != "n"]
    pro_filt <- as.character(data.frame(count_bypro[which(count_bypro$n < min_num_peps), ])[, pro_id])
    
    # determine which peptides no longer have a protein to map to  #
    ## find rows in peptide.info that correspond to proteins to be filtered ##
    protfilt.ids <- which(attr(filter_object, "e_meta")[,pro_id] %in% pro_filt)
    
    
    ## find the peptides that are in the filter list but are not in the unfiltered lists ##
    pep_filt <- as.character(setdiff(attr(filter_object, "e_meta")[protfilt.ids, pep_id], attr(filter_object, "e_meta")[-protfilt.ids, pep_id]))
    
    
    if(length(pep_filt)==0) {pep_filt = NULL}
    if(length(pro_filt)==0) {pro_filt = NULL}
    
    filter_object_new2 <- list(proteins_filt = pro_filt, peptides_filt = pep_filt)
    
  }else{
    filter_object_new2 <- list(peptides_filt = c(), proteins_filt = c())
  }
  
  if(degen_peps){
    count_bypro <- filter_object$counts_by_pro
    count_bypep <- filter_object$counts_by_pep
    
    pro_id <- names(count_bypro)[names(count_bypro) != "n"]
    pep_id <- names(count_bypep)[names(count_bypep) != "n"]
    degen_peptides <- as.character(data.frame(count_bypep[which(count_bypep$n > 1), ])[, pep_id])
    
    ## identify any proteins that now will not have peptides mapping to them ##
    ## find rows in e_meta that correspond to peptides to be filtered ##
    pepfilt.ids <- which(attr(filter_object, "e_meta")[,pep_id] %in% degen_peptides)
    
    ## find the proteins that are in the filter list but are not in the unfiltered lists ##
    add_prots <- as.character(setdiff(attr(filter_object, "e_meta")[pepfilt.ids, pro_id], attr(filter_object, "e_meta")[-pepfilt.ids, pro_id]))
    
    if(length(add_prots)==0) {add_prots = NULL}
    if(length(degen_peptides)==0) {degen_peptides = NULL}
    filter_object_new1 <- list(peptides_filt = degen_peptides, proteins_filt = add_prots)
  } else {
    filter_object_new1 <- list(peptides_filt = c(), proteins_filt = c())
  }
  
  ## consolidate filter_object_new1 and filter_object_new2 ##
  filter_object_new <- list(proteins_filt = unique(c(filter_object_new1$proteins_filt, filter_object_new2$proteins_filt)), peptides_filt = unique(c(filter_object_new1$peptides_filt, filter_object_new2$peptides_filt)))
  
  num_filtered <- lapply(filter_object_new, length)
  # return the 5 number summary for both parts of the filter, give total number of peps and/or prots filtered
  res <- list(num_per_pep=pep_sum, num_per_pro=pro_sum, num_pep_filtered=num_filtered$peptides_filt, num_pro_filtered=num_filtered$proteins_filt)
  
  
  # create output #
  catmat <- data.frame(c(pep_sum, " ", res$num_pep_filtered, nrow(filter_object$counts_by_pep)-res$num_pep_filtered),
                       c(pro_sum, " ", res$num_pro_filtered, nrow(filter_object$counts_by_pro)-res$num_pro_filtered))
  colnames(catmat) <- c("Obs. Per Peptide", "Obs. Per Protein")
  rownames(catmat) <- c(names(pep_sum), " ", "Filtered", "Not Filtered")
  
  cat("\nSummary of Proteomics Filter \n\n")
  cat(capture.output(catmat), sep="\n")
  cat("\n")
  
  return(invisible(res))
}




#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
#'@param min_nonmiss_gtest the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Suggested value is 3.
#'@param min_nonmiss_anova the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Suggested value is 2.
summary.imdanovaFilt <- function(filter_object, min_nonmiss_anova=NULL, min_nonmiss_gtest=NULL){
  
  ## initial checks ##
  
  # it is fine if both min_nonmiss_anova and min_nonmiss_gtest are NULL in the summary function #
  
  # check that if they aren't NULL, min_nonmiss_anova and min_nonmiss_gtest are numeric, >=2 and >=3, respectively, and neither are bigger than the minimum group size (group_sizes in an attribute of the filter_object, see below) #
  if(!is.null(min_nonmiss_anova)) {
    # check that min_nonmiss_anova is not a vector #
    if(length(min_nonmiss_anova) > 1) stop("min_nonmiss_anova must be of length 1")
    # check that min_nonmiss_anova is numeric >= 2 #
    if(!is.numeric(min_nonmiss_anova) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is an integer #
    if(min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is less than the minimum group size #
    if(min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_anova cannot be greater than the minimum group size")
  }
  if(!is.null(min_nonmiss_gtest)) {
    # check that min_nonmiss_gtest is not a vector #
    if(length(min_nonmiss_gtest) > 1) stop("min_nonmiss_gtest must be of length 1")
    # check that min_nonmiss_gtest is numeric >= 3 #
    if(!is.numeric(min_nonmiss_gtest) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is an integer #
    if(min_nonmiss_gtest %% 1 != 0) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is less than the minimum group size #
    if(min_nonmiss_gtest > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_gtest cannot be greater than the minimum group size")
  }
  
  ## end of initial checks ##
  
  
  my_names <- as.character(attr(filter_object, "group_sizes")$Group)
  
  # if min_nonmiss_anova is not provided, the vector of zeros will indicate nothing needs removing
  inds_rm_anova <- rep(0, nrow(filter_object))
  
  # if min_nonmiss_gtest is not provided, the vector of zeros will indicate nothing needs removing
  inds_rm_gtest <- rep(0, nrow(filter_object))
  
  
  # if min_nonmiss_anova is provided
  if(!is.null(min_nonmiss_anova)){
    # need at least n=min_nonmiss_anova per group in at least 2 groups in order to keep the peptide
    temp_anova <- (filter_object[,which(names(filter_object) %in% my_names)] >= min_nonmiss_anova)
    
    # sum the number of groups that meet the nonmissing per group requirement
    # these are the rows to remove since they do not have at least 2 groups not meeting nonmissing requirements
    inds_rm_anova <- rowSums(temp_anova) < 2
  }
  
  
  # if min_nonmiss_gtest is provided
  if(!is.null(min_nonmiss_gtest)){
    # need at least min_nonmiss_gtest peptide IDs in one group #
    temp_gtest <- (filter_object[,which(names(filter_object) %in% my_names)] >= min_nonmiss_gtest)
    
    # sum the number of groups that meet the nonmissing per group requirement
    inds_rm_gtest <- rowSums(temp_gtest) == 0 # remove these rows since no groups meet the requirement
  }
  
  
  #if-statement for the case where both min_nonmiss_anova and min_nonmiss_gtest are non-NULL. Note: num_not_tested will be the count of (inds_rm_anova + inds_rm_gtest) that sum to 2. That is the intersection of both filters. 
  if(!is.null(min_nonmiss_anova) & !is.null(min_nonmiss_gtest)){
    # get number molecules not tested
    num_filtered <- sum((inds_rm_anova + inds_rm_gtest) > 1)
    
    # get number molecules tested
    num_not_filtered <- nrow(filter_object) - num_filtered
    
  }
  
  else if(!is.null(min_nonmiss_anova) | !is.null(min_nonmiss_gtest)){
    # get number molecules not tested
    num_filtered <- sum((inds_rm_anova + inds_rm_gtest) > 0)
    
    # get number molecules tested
    num_not_filtered <- nrow(filter_object) - num_filtered
    
  }
  else{
    num_not_filtered = NULL
    num_filtered = NULL
  }
  
  res <- list(pep_observation_counts = nrow(filter_object), num_filtered = num_filtered, num_not_filtered = num_not_filtered)
  
  # create output #
  catmat <- data.frame(sapply(res, function(x) ifelse(is.null(x),"NULL",x)))
  rownames(catmat) <- c("Total Observations: ", "Filtered: ", "Not Filtered: ")
  colnames(catmat) <- NULL
  
  cat("\nSummary of IMD Filter\n")
  cat(capture.output(catmat), sep="\n")
  
  return(invisible(res))
}




#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
#'@param pvalue_threshold A threshold for the Robust Mahalanobis Distance (RMD) p-value. All samples below the threshold will be filtered out. Default value is NULL.
summary.rmdFilt <- function(filter_object, pvalue_threshold = NULL){
  
  # check that pvalue_threshold is numeric [0,1] #
  if(!is.null(pvalue_threshold)) {
    if(!is.numeric(pvalue_threshold)) stop("pvalue_threshold must be numeric between 0 and 1")
    if(pvalue_threshold < 0 | pvalue_threshold > 1) stop("pvalue_threshold must be numeric between 0 and 1")
  }
  
  # get metrics used #
  samp_id <- names(attr(filter_object, "group_DF"))[1]
  metrics <- attributes(filter_object)$metrics
  
  # get samples filtered out, if pvalue_threshold is specified #
  filt <- NULL
  if(!is.null(pvalue_threshold)) {
    filt <- as.character(filter_object[filter_object$pvalue < pvalue_threshold, samp_id])
    if(length(filt)==0) filt <- NULL
  }
  
  # for display #
  if(is.null(filt)) {filt <- "NULL"}
  
  res <- list(pvalue = summary(filter_object$pvalue), metrics = metrics, filtered_samples = filt)
  
  # create output #
  cat(c("\nSummary of RMD Filter\n----------------------\nP-values:\n", capture.output(res$pvalue)), sep = "\n")
  cat(c("\nMetrics Used:", paste(res$metrics, collapse = ", "), "\n"))
  cat(c("\nFiltered Samples:", paste(res$filtered_samples, collapse = ", "), "\n\n"))
  
  return(invisible(res))
}




#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
#'@param cv_threshold numeric value greater than 1 and less than the value given by filter_object$CV_pooled.  CV values above cv_threshold are filtered out. Default value is NULL.
summary.cvFilt <- function(filter_object, cv_threshold = NULL){
  
  # checks for cv_threshold if not null #
  if(!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if(!is.numeric(cv_threshold)) stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if(length(cv_threshold)>1) stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if(cv_threshold <= 1 | cv_threshold >= max(filter_object$CV_pooled, na.rm = TRUE)) stop("cv_threshold must be greater than 1 and less than the maximum CV_pooled value")
  }
  
  # get rid of NAs #
  new_object <- filter_object[!is.na(filter_object$CV_pooled),]
  
  # get summary of CVs #
  CVs <- summary(new_object$CV_pooled)
  
  # get total NAs, total non-NAs #
  tot_NAs <- attributes(filter_object)$tot_nas
  tot_non_NAs <- nrow(filter_object) - tot_NAs
  
  # get biomolecules to filter if cv_threshold is not NULL #
  filt <- NULL
  if(!is.null(cv_threshold)) {
    filt <- as.character(new_object[new_object$CV_pooled > cv_threshold, 1])
    if(length(filt)==0) filt <- NULL
  }
  
  res <- list(CVs = CVs, tot_NAs = tot_NAs, tot_non_NAs = tot_non_NAs, filtered_biomolecules = length(filt))
  
  # create output #
  cat(c("\nSummary of Coefficient of Variation (CV) Filter\n----------------------\nCVs:\n", capture.output(res$CVs)), sep = "\n")
  cat(c("\nTotal NAs:", res$tot_NAs))
  cat(c("\nTotal Non-NAs:", res$tot_non_NAs, "\n"))
  cat(c("\nNumber Filtered Biomolecules:", paste(res$filtered_biomolecules, collapse = ", "), "\n\n"))
  
  return(invisible(res))
}




#'@export
#'@rdname summary-pmartR-filter
#'@name summary-pmartR-filter
summary.customFilt <- function(filter_object){
  
  # get omicsData object #
  omicsData <- attr(filter_object, "omicsData")
  summary_orig <- summary(omicsData)
  
  # get names #
  edata_id <- attr(filter_object, "cnames")$edata_cname
  emeta_id <- attr(filter_object, "cnames")$emeta_cname
  samp_id <- attr(filter_object, "cnames")$fdata_cname
  
  # apply the filter #
  filtered_data <- applyFilt(filter_object, omicsData)
  summary_filt <- summary(filtered_data)
  
  #if filter_object contains removes
  if(!is.null(filter_object$e_data_remove)|!is.null(filter_object$f_data_remove)|!is.null(filter_object$e_meta_remove))
  {
    # samples #
    num_samples <- attributes(filter_object)$num_samples
    if(!is.null(filter_object$f_data_remove)) {
      samps_filt <- length(filter_object$f_data_remove)
      samps_left <- num_samples - samps_filt
    } else {
      samps_filt <- 0
      samps_left <- num_samples
    }
    
    # e_data #
    num_edata <- attributes(filter_object)$num_edata
    if(!is.null(filter_object$e_data_remove)) {
      edata_filt <- length(filter_object$e_data_remove)
      edata_left <- num_edata - edata_filt
    } else {
      edata_filt <- num_edata - nrow(filtered_data$e_data)
      edata_left <- nrow(filtered_data$e_data)
    }
    
    # e_meta #
    if(!is.null(filtered_data$e_meta)){
      num_emeta <- attributes(filter_object)$num_emeta
      if(!is.null(filter_object$e_meta_remove)) {
        emeta_filt <- length(filter_object$e_meta_remove)
        emeta_left <- num_emeta - emeta_filt
      } else {
        emeta_filt <- num_emeta - length(unique(filtered_data$e_meta[, emeta_id] ))
        emeta_left <- length(unique(filtered_data$e_meta[, emeta_id] ))
      }
    }
    # Display #
    
    samp_id <- paste(samp_id, "s (f_data)", sep="")
    edata_id <- paste(edata_id, "s (e_data)", sep="")
    display_emeta_id <- paste(emeta_id, "s (e_meta)", sep="")
    
    ## construct data frame ##
    if(is.null(emeta_id)) {
      disp <- data.frame(Filtered = c(samps_filt, edata_filt), Remaining = c(samps_left, edata_left), Total = c(num_samples, num_edata))
      rownames(disp) <- c(samp_id, edata_id)
    } else {
      disp <- data.frame(Filtered = c(samps_filt, edata_filt, emeta_filt), Remaining = c(samps_left, edata_left, emeta_left), Total = c(num_samples, num_edata, num_emeta))
      rownames(disp) <- c(samp_id, edata_id, display_emeta_id)
    }
    
  }
  
  #if filter_object contains keeps
  if(!is.null(filter_object$e_data_keep)|!is.null(filter_object$f_data_keep)|!is.null(filter_object$e_meta_keep))
  {
    # samples #
    num_samples <- attributes(filter_object)$num_samples
    if(!is.null(filter_object$f_data_keep)) {
      samps_keep <- length(filter_object$f_data_keep)
      samps_discard <- num_samples - samps_keep
    } else {
      samps_keep <- num_samples
      samps_discard <- 0
    }
    
    # e_data #
    num_edata <- attributes(filter_object)$num_edata
    if(!is.null(filter_object$e_data_keep)) {
      edata_keep <- length(filter_object$e_data_keep)
      edata_discard <- num_edata - edata_keep
    } else {
      edata_keep <- nrow(filtered_data$e_data)
      edata_discard <- num_edata - nrow(filtered_data$e_data)
    }
    
    # e_meta #
    num_emeta <- attributes(filter_object)$num_emeta
    if(!is.null(filter_object$e_meta_keep)) {
      emeta_keep <- length(filter_object$e_meta_keep)
      emeta_discard <- num_emeta - emeta_keep
    } else {
      emeta_keep <- length(unique(filtered_data$e_meta[, emeta_id] ))
      emeta_discard <- num_emeta - length(unique(filtered_data$e_meta[, emeta_id] ))
    }
    
    # Display #
    
    samp_id <- paste(samp_id, "s (f_data)", sep="")
    edata_id <- paste(edata_id, "s (e_data)", sep="")
    display_emeta_id <- paste(emeta_id, "s (e_meta)", sep="")
    
    ## construct data frame ##
    if(is.null(emeta_id)) {
      disp <- data.frame(Kept = c(samps_keep, edata_keep), Discarded = c(samps_discard, edata_discard), Total = c(num_samples, num_edata))
      rownames(disp) <- c(samp_id, edata_id)
    } else {
      disp <- data.frame(Kept = c(samps_keep, edata_keep, emeta_keep), Discarded = c(samps_discard, edata_discard, emeta_discard), Total = c(num_samples, num_edata, num_emeta))
      rownames(disp) <- c(samp_id, edata_id, display_emeta_id)
    }
    
  }
  
  
  ## Display output ##
  cat("\nSummary of Custom Filter\n\n")
  cat(capture.output(disp), sep = "\n")
  cat("\n")
  
  return(invisible(disp))
}




