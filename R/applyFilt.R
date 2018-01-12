#' Apply a S3 filter  object to an MSomics S3 object
#'
#' This function takes a filter object of class 'cvFilt', 'rmdFilt', 'moleculeFilt', 'proteomicsFilt', 'imdanovaFilt', or 'customFilt' and applies the filter to a dataset of class \code{pepData}, \code{proData}, \code{lipidData}, or \code{metabData}.
#'
#' @param filter_object an object of the class 'cvFilt', 'proteomicsFilt', 'rmdFilt', 'moleculeFilt', 'imdanovaFilt', or 'customFilt' created by \code{cv_filter}, \code{proteomics_filter}, \code{rmd_filter}, \code{molecule_filter}, or \code{imdanova_filter}, respectively.
#' @param omicsData an object of the class \code{pepData}, \code{proData}, \code{lipidData}, or \code{metabData} usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param ... further arguments
#'
#' @return An object of the class \code{pepData}, \code{proData}, \code{lipidData}, or \code{metabData},  with specified cname_ids, edata_cnames, and emeta_cnames filtered out of the appropriate datasets.
#'
#' @details Various further arguments can be specified depending on the class of the \code{filter_object} being applied.
#' For a \code{filter_object} of type 'moleculeFilt':
#' \tabular{ll}{
#' \code{min_num} \tab an integer value specifying the minimum number of times each biomolecule must be observed across all samples in order to retain the biomolecule. Default value is 2. \cr
#' }
#' For a \code{filter_object} of type 'cvFilt':
#' \tabular{ll}{
#' \code{cv_threshold} \tab an integer value specifying the maximum coefficient of variation (CV) threshold for the biomolecules. Biomolecules with CV greater than this threshold will be filtered. Default is 150. \cr
#' }
#' For a \code{filter_object} of type 'rmdFilt':
#' \tabular{ll}{
#' \code{pvalue_threshold} \tab numeric value between 0 and 1, specifying the p-value, below which samples will be removed from the dataset. Default is 0.001. \cr
#' }
#' For a \code{filter_object} of type 'proteomicsFilt':
#' \tabular{ll}{
#' \code{min_num_peps} \tab an optional integer value between 1 and the maximum number of peptides that map to a protein in omicsData. The value specifies the minimum number of peptides that must map to a protein. Any protein with less than \code{min_num_peps} mapping to it will be returned as a protein that should be filtered. Default value is NULL. \cr
#' \code{degen_peps} \tab logical indicator of whether to filter out degenerate peptides (TRUE) or not (FALSE). Default value is FALSE.\cr
#' }
#' For a \code{filter_object} of type 'imdanovaFilt':
#' \tabular{ll}{
#' \code{min_nonmiss_anova} \tab integer value specifying the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Default value is 2. \cr
#' \code{min_nonmiss_gtest} \tab integer value specifying the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Default value is 3.\cr
#' }
#' There are no futher arguments for a \code{filter_object} of type 'customFilt'.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' to_filter <- molecule_filter(omicsData = pep_object)
#' pep_object2 <- apply(filter_object = to_filter, omicsData = pep_object, min_num = 2)
#' print(str(attributes(pep_object2)$filters))
#' to_filter2 <- imdanova_filter(omicsData = pep_object2)
#' pep_object3 <- apply(filter_object = to_filter2, omicsData = pep_object2, min_num_peps = 4)
#' print(str(attributes(pep_object3)$filters))
#' }
#'
#' @seealso \code{\link{molecule_filter}}
#' @seealso \code{\link{imdanova_filter}}
#' @seealso \code{\link{rmd_filter}}
#' @seealso \code{\link{cv_filter}}
#' @seealso \code{\link{proteomics_filter}}
#' @seealso \code{\link{custom_filter}}
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
applyFilt <- function(filter_object, omicsData, ...){

  # check that omicsData is of mintR S3 class#
  if(!(class(omicsData) %in% c("pepData", "proData", "lipidData", "metabData"))) stop("omicsData must be of class 'pepData', 'proData', 'lipidData', or 'metabData'")

  # check that filter_object is of an appropriate class #
  if(!(any(class(filter_object) %in% c("cvFilt", "proteomicsFilt", "moleculeFilt", "rmdFilt", "imdanovaFilt", "customFilt")))) stop("filter_object must be of class 'cvFilt', 'proteomicsFilt', 'moleculeFilt', 'rmdFilt', 'imdanovaFilt', or 'customFilt.")

  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname

  # generate warnings if data type is not present but user asks to filter #
  if(!is.null(filter_object$filter_edata) & is.null(edata_cname)) warning("e_data identifier column specified in filter_object is not present in omicsData$e_data. Specified e_data features cannot be filtered from data.")

  if(!is.null(filter_object$filter_emeta) & is.null(emeta_cname)) warning("e_meta identifier column specified in filter_object is not present in omicsData$e_meta. Specified e_meta features cannot be filtered from data.")

  if(!is.null(filter_object$filter_samples) & is.null(samp_cname)) warning("Sample identifier column specified in filter_object is not present in omicsData. Specified samples cannot be filtered from data.")


  UseMethod("applyFilt")
}

# function for moleculeFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.moleculeFilt <- function(filter_object, omicsData, min_num=2){

  # check to see whether a moleculeFilt has already been run on omicsData #
  if("moleculeFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    min_num_prev <- attributes(omicsData)$filters$moleculeFilt$threshold

    stop(paste("A molecule filter has already been run on this dataset, using a 'min_num' of ", min_num_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous moleculeFilt, so go ahead and run it like normal #


    # check that min_num is numeric and >=1 #
    if(!(class(min_num) %in% c("numeric","integer")) | min_num < 1) stop("min_num must be an integer greater than or equal to 1")
    # check that min_num is an integer #
    if(min_num %% 1 != 0) stop("min_num must be an integer greater than or equal to 1")
    # check that min_num is less than the number of samples #
    if(min_num > max(filter_object$Num_Observations)) stop("min_num cannot be greater than the number of samples")
    # check that min_num is of length 1 #
    if(length(min_num) != 1) stop("min_num must be of length 1")

    edata_cname <- attributes(omicsData)$cnames$edata_cname
    emeta_cname <- attributes(omicsData)$cnames$emeta_cname
    
    num_obs <- filter_object$Num_Observations
    #min_num <- attr(filter_object, "min_num")

    # get indices for which ones don't meet the min requirement #
    inds <- which(num_obs < min_num)

    if(length(inds) < 1){
      filter.edata <- NULL
    }
    
    else{
      filter.edata <- omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)][inds]
    }
    
    #checking if filter specifies all of omicsData$e_data
    if(all(omicsData$e_data[[edata_cname]] %in% filter.edata)) {stop("filter_object specifies all samples in omicsData")}
    
    filter_object_new = list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL)

    # call the function that does the filter application
    results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    results <- omicsData
    results$e_data <- results_pieces$temp.pep2
    results$f_data <- results_pieces$temp.samp2
    results$e_meta <- results_pieces$temp.meta1

    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }

    # set attributes for which filters were run
    attr(results, "filters")$moleculeFilt <- list(report_text = "", threshold = c(), filtered = c())
    attr(results, "filters")$moleculeFilt$report_text <- paste("A molecule filter was applied to the data, removing ", edata_cname, "s ", "that were present in fewer than ", min_num, " samples. A total of ", length(filter.edata), " ", edata_cname, "s ", "were filtered out of the dataset by this filter.", sep="")
    attr(results, "filters")$moleculeFilt$threshold <- min_num
    attr(results, "filters")$moleculeFilt$filtered <- filter.edata

  }

  return(results)
}

# function for cvFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.cvFilt <- function(filter_object, omicsData, cv_threshold = 150){

  # check to see whether a cvFilt has already been run on omicsData #
  if("cvFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    cv_threshold_prev <- attributes(omicsData)$filters$cvFilt$threshold

    stop(paste("A CV filter has already been run on this dataset, using a 'cv_threshold' of ", cv_threshold_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous cvFilt, so go ahead and run it like normal #

    # check that cv_threshold is greater than 0 #
    if(cv_threshold <= 0) stop("cv_threshold must be greater than 0.")

    samp_cname <- attributes(omicsData)$cnames$fdata_cname
    edata_cname <- attributes(omicsData)$cnames$edata_cname
    emeta_cname <- attributes(omicsData)$cnames$emeta_cname

    # determine which peptides have a CV greater than the threshold #
    p.ids = which(filter_object$CV_pooled > cv_threshold)
  
    # return peptide names to be filtered #
    if(length(p.ids) > 0){p_filt = as.character(filter_object[p.ids,1])}else{p_filt = NULL}
    
    #checking if filter specifies all peptides
    if(all(omicsData$e_data[[edata_cname]] %in% p_filt)) {stop("filter_object specifies all peptides in omicsData")}
    
    filter_object_new = list(edata_filt = p_filt, emeta_filt = NULL, samples_filt = NULL)

    # call the function that does the filter application
    results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    results <- omicsData
    results$e_data <- results_pieces$temp.pep2
    results$f_data <- results_pieces$temp.samp2
    results$e_meta <- results_pieces$temp.meta1

    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }

    # set attributes for which filters were run

    attr(results, "filters")$cvFilt <- list(report_text = "", threshold = c(), filtered = c())
    attr(results, "filters")$cvFilt$report_text <- paste("A coefficient of variation (CV) filter was applied to the data, removing ", edata_cname, "s ", "with a CV greater than ", cv_threshold, ". A total of ", length(p.ids), " ", edata_cname, "s were filtered out of the dataset by this filter.", sep="")
    attr(results, "filters")$cvFilt$threshold <- cv_threshold
    attr(results, "filters")$cvFilt$filtered <- p.ids

  }

  return(results)
}



# function for rmdFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.rmdFilt <- function(filter_object, omicsData, pvalue_threshold=0.001){


  # check to see whether a rmdFilt has already been run on omicsData #
  if("rmdFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    threshold_prev <- attributes(omicsData)$filters$rmdFilt$threshold

    stop(paste("A rMd filter has already been run on this dataset, using a 'pvalue_threshold' of ", threshold_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous rmdFilt, so go ahead and run it like normal #

    # check that filter_object is of class rmdFilt #
    if(!"rmdFilt" %in% class(filter_object)) stop("filter_object must be of the class 'rmdFilt'. See rmd_filter for details.")

    # check that pvalue_threshold is between 0 and 1 #
    if(pvalue_threshold < 0 | pvalue_threshold > 1) stop("pvalue_threshold must be between 0 and 1.")

    samp_cname <- attributes(omicsData)$cnames$fdata_cname
    edata_cname <- attributes(omicsData)$cnames$edata_cname
    emeta_cname <- attributes(omicsData)$cnames$emeta_cname

    # determine which samples have a pvalue less than the threshold #
    samp.ids = which(filter_object$pvalue < pvalue_threshold)
    
    # return sample names to be filtered #
    if(length(samp.ids) > 0){samp_filt = as.character(filter_object[samp.ids,1])}else{samp_filt = NULL}
    
    #checking that filter_object does not specify all samples
    if(all(as.character(omicsData$f_data[[samp_cname]]) %in% samp_filt)) {stop("samples_filt specifies all samples")}
    
    filter_object_new = list(edata_filt = NULL, emeta_filt = NULL, samples_filt = samp_filt)

    # call the function that does the filter application
    results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    results <- omicsData
    results$e_data <- results_pieces$temp.pep2
    results$f_data <- results_pieces$temp.samp2
    results$e_meta <- results_pieces$temp.meta1

    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }

    # set attributes for which filters were run
    attr(results, "filters")$rmdFilt <- list(report_text = "", threshold = c(), filtered = c())

    attr(results, "filters")$rmdFilt$report_text <- paste("A robust Mahalanobis distance (rMd) filter was applied to the data, removing ", samp_cname, "s ", "with an rMd-associated p-value less than ", pvalue_threshold, ". A total of ", length(samp_filt), " ", samp_cname, "s were filtered out of the dataset by this filter.", sep="")

    attr(results, "filters")$rmdFilt$threshold <- pvalue_threshold
    attr(results, "filters")$rmdFilt$filtered <- samp_filt

  }

  return(results)
}



# function for proteomicsFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
#' @export
applyFilt.proteomicsFilt <- function(filter_object, omicsData, min_num_peps=NULL, degen_peps=FALSE){
  # #' @describeIn MSomics_filter MSomics_filter for proteomicsFilt S3 object


  # check to see whether a "proteomicsFilt" has already been run on omicsData #
  if("proteomicsFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    threshold_prev <- attributes(omicsData)$filters$proteomicsFilt$threshold$min_num_peps
    degen_peps_prev <- attributes(omicsData)$filters$proteomicsFilt$threshold$degen_peps

    stop(paste("A proteomics filter has already been run on this dataset, using a 'min_num_peps' of ", threshold_prev, " and 'degen_peps' of ", degen_peps_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous rmdFilt, so go ahead and run it like normal #

    # error checks for min_num_peps, if not NULL #
    if(!is.null(min_num_peps)) {
      # check that min_num_peps is numeric and >=1 #
      if(class(min_num_peps) != "numeric"| min_num_peps < 1) stop("min_num_peps must be an integer greater than or equal to 1")
      # check that min_num_peps is an integer #
      if(min_num_peps %% 1 != 0) stop("min_num_peps must be an integer greater than or equal to 1")
      # check that min_num_peps is of length 1 #
      if(length(min_num_peps) != 1) stop("min_num_peps must be of length 1")
      # check that min_num_peps is less than the total number of peptides #
      if(min_num_peps > nrow(omicsData$e_data)) stop("min_num_peps cannot be greater than the total number of peptides")
    }
    # check that degen_peps is logical #
    if(class(degen_peps) != "logical") stop("degen_peps must be either TRUE or FALSE")


    pep_id = attr(omicsData, "cnames")$edata_cname
    pro_id = attr(omicsData, "cnames")$emeta_cname
    samp_cname <- attributes(omicsData)$cnames$fdata_cname

    ## degenerate peptides portion ##
    # pull peptides with more than one row #
    if(degen_peps == TRUE){
      count_bypep <- filter_object$counts_by_pep

      degen_peptides = as.character(data.frame(count_bypep[which(count_bypep$n > 1), ])[, pep_id])

      ## identify any proteins that now will not have peptides mapping to them ##
      ## find rows in e_meta that correspond to peptides to be filtered ##
      pepfilt.ids = which(omicsData$e_meta[,pep_id] %in% degen_peptides)

      ## find the proteins that are in the filter list but are not in the unfiltered lists ##
      add_prots = as.character(setdiff(omicsData$e_meta[pepfilt.ids, pro_id], omicsData$e_meta[-pepfilt.ids, pro_id]))

      if(length(add_prots)==0){add_prots = NULL}
      if(length(degen_peptides)==0){degen_peptides = NULL}
      filter_object_new1 <- list(edata_filt = degen_peptides, emeta_filt = add_prots)
    }else{
      filter_object_new1 <- list(edata_filt = c(), emeta_filt = c())
    }


    ## protein filter portion ##
    # pull proteins with less than min_num_peps #
    if(!is.null(min_num_peps)){
      count_bypro <- filter_object$counts_by_pro

      pro_filt = as.character(data.frame(count_bypro[which(count_bypro$n < min_num_peps), ])[, pro_id])

      # determine which peptides no longer have a protein to map to  #
      ## find rows in peptide.info that correspond to proteins to be filtered ##
      protfilt.ids = which(omicsData$e_meta[,pro_id] %in% pro_filt)

      ## find the peptides that are in the filter list but are not in the unfiltered lists ##
      pep_filt = as.character(setdiff(omicsData$e_meta[protfilt.ids, pep_id], omicsData$e_meta[-protfilt.ids, pep_id]))

      if(length(pep_filt)==0){pep_filt = NULL}
      if(length(pro_filt)==0){pro_filt = NULL}

      filter_object_new2 <- list(emeta_filt = pro_filt, edata_filt = pep_filt)
    }else{
      filter_object_new2 <- list(edata_filt = c(), emeta_filt = c())
    }


    ## consolidate filter_object_new1 and filter_object_new2 ##
    filter_object_new <- list(emeta_filt = unique(c(filter_object_new1$emeta_filt, filter_object_new2$emeta_filt)), edata_filt = unique(c(filter_object_new1$edata_filt, filter_object_new2$edata_filt)))
    
    #checking that filter_object_new does not specify all of e_data(peps) or all of e_meta(protiens) in omicsData
    if(all(omicsData$e_meta[[pro_id]] %in% filter_object_new$emeta_filt)) {stop("filter_object specifies all proteins in e_meta")}
    if(all(omicsData$e_data[[pep_id]] %in% filter_object_new$edata_filt)) {stop("filter_object specifies all peps in e_data")}

    # call the function that does the filter application
    results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    results <- omicsData
    results$e_data <- results_pieces$temp.pep2
    results$f_data <- results_pieces$temp.samp2
    results$e_meta <- results_pieces$temp.meta1

    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }

    # set attributes for which filters were run
    attr(results, "filters")$proteomicsFilt <- list(report_text = "", threshold = c(), filtered = list())
    attr(results, "filters")$proteomicsFilt$threshold <- data.frame(min_num_peps = min_num_peps, degen_peps = as.character(degen_peps))
    attr(results, "filters")$proteomicsFilt$filtered <- filter_object_new
    if(degen_peps == TRUE & is.null(min_num_peps)){
      attr(results, "filters")$proteomicsFilt$report_text <- paste("A degenerate peptide filter was applied to the data, identifying ", pep_id, "s that map to more than one ", pro_id, ". The filter identified ", length(filter_object_new1$peptides_filt), " such ", pep_id, "s in the data. Associated with these ", pep_id, "s were ", length(filter_object_new1$proteins_filt), " ", pro_id, "s that no longer had ", pep_id, "s mapping to them.", sep = "")

    }else if(degen_peps == FALSE & !is.null(min_num_peps)){
      attr(results, "filters")$proteomicsFilt$report_text <- paste("A protein filter was applied to the data, identifying ", pro_id, "s that do not have at least ", min_num_peps, " ", pep_id, "s mapping to them. This filter identified ", length(filter_object_new2$proteins_filt), " ", pro_id, "s from the data. Associated with these ", pro_id, "s were ", length(filter_object_new2$peptides_filt), " ", pep_id, "s that no longer mapped to a ", pro_id, ".", sep="")

    }else if(degen_peps == TRUE & !is.null(min_num_peps)) {
      attr(results, "filters")$proteomicsFilt$report_text <- paste("A degenerate peptide filter was applied to the data, identifying ", pep_id, "s that map to more than one ", pro_id, ". The filter identified ", length(filter_object_new1$peptides_filt), " such ", pep_id, "s in the data. Associated with these ", pep_id, "s were ", length(filter_object_new1$proteins_filt), " ", pro_id, "s that no longer had ", pep_id, "s mapping to them. Additionally, a protein filter was applied to the data, identifying ", pro_id, "s that do not have at least ", min_num_peps, " ", pep_id, "s mapping to them. This filter identified ", length(filter_object_new2$proteins_filt), " ", pro_id, "s from the data. Associated with these ", pro_id, "s were ", length(filter_object_new2$peptides_filt), " ", pep_id, "s that no longer mapped to a ", pro_id, ". Together, the two filters removed ", length(filter_object_new$peptides_filt), " ", pep_id, "s and ", length(filter_object_new$proteins_filt), " ", pro_id, "s from the data.", sep="")

    }

  }

  return(results)


}

# function for imdanovaFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.imdanovaFilt <- function(filter_object, omicsData, min_nonmiss_anova=NULL, min_nonmiss_gtest=NULL){
  # #' @details If filter_method="combined" is specified, then both the \code{anova_filter} and \code{gtest_filter} are applied to the data, and the intersection of features from the two filters is the set returned. For ANOVA, features that do not have at least \code{min_nonmiss_allowed} values per group are candidates for filtering. For the G-test, features that do not have at least \code{min_nonmiss_allowed} values per group are candidates for filtering. The G-test is a test of independence, used here to test the null hypothesis of independence between the number of missing values across groups.


  # check to see whether a imdanovaFilt has already been run on omicsData #
  if("imdanovaFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    min_nonmiss_anova_prev <- attributes(omicsData)$filters$imdanovaFilt$threshold$min_nonmiss_anova
    min_nonmiss_gtest_prev <- attributes(omicsData)$filters$imdanovaFilt$threshold$min_nonmiss_gtest

    stop(paste("A proteomics filter has already been run on this dataset, using a 'min_nonmiss_anova' of ", min_nonmiss_anova_prev, " and 'min_nonmiss_gtest' of ", min_nonmiss_gtest_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous rmdFilt, so go ahead and run it like normal #


    ## initial checks ##

    # check that at least one of min_nonmiss_anova and min_nonmiss_gtest are present #
    if(is.null(min_nonmiss_anova) & is.null(min_nonmiss_gtest)) stop("At least one of min_nonmiss_anova and min_nonmiss_gtest must be present")
    # check that if they aren't NULL, min_nonmiss_anova and min_nonmiss_gtest are numeric, >=2 and >=3, respectively, and neither are bigger than the minimum group size (group_sizes in an attribute of the filter_object, see below) #
    if(!is.null(min_nonmiss_anova)) {
      # check that min_nonmiss_anova is numeric >= 2 #
      if(!(class(min_nonmiss_anova) %in% c("numeric","integer")) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
      # check that min_nonmiss_anova is an integer #
      if(min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
      # check that min_nonmiss_anova is less than the max number of observations #
      if(min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_anova cannot be greater than the minimum group size")
    }
    if(!is.null(min_nonmiss_gtest)) {
      # check that min_nonmiss_gtest is numeric >= 3 #
      if(!(class(min_nonmiss_gtest) %in% c("numeric","integer")) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
      # check that min_nonmiss_gtest is an integer #
      if(min_nonmiss_gtest %% 1 != 0) stop("min_nonmiss_gtest must be an integer >= 3")
      # check that min_nonmiss_gtest is less than the max number of observations #
      if(min_nonmiss_gtest > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_gtest cannot be greater than the minimum group size")
    }

    ## end of initial checks ##


    group_sizes <- attr(filter_object, "group_sizes")
    nonmiss_per_group <- list(nonmiss_totals = filter_object, group_sizes = group_sizes)

    groupDF <- attributes(omicsData)$group_DF
    e_data <- omicsData$e_data
    edata_cname <- attr(omicsData, "cnames")$edata_cname
    samp_cname <- attr(omicsData, "cnames")$fdata_cname
    emeta_cname <- attributes(omicsData)$cnames$emeta_cname

    if(!is.null(min_nonmiss_anova) & !is.null(min_nonmiss_gtest)){
      filter_method <- "combined"
    }else{
      if(!is.null(min_nonmiss_anova) & is.null(min_nonmiss_gtest)){
        filter_method <- "anova"
      }else{
        if(is.null(min_nonmiss_anova) & !is.null(min_nonmiss_gtest)){
          filter_method <- "gtest"
        }
      }
    }
    if(filter_method=="anova"){
      filter.edata <- anova_filter(nonmiss_per_group=nonmiss_per_group, min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_cname)
    }else{
      if(filter_method=="gtest"){
        filter.edata <- gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=groupDF, e_data=e_data, alpha=NULL, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_cname, samp_id = samp_cname)
      }else{
        if(filter_method=="combined"){
          filter.edata.gtest <- gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=groupDF, e_data=e_data, alpha=NULL, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_cname, samp_id = samp_cname)
          #           min.nonmiss.allowed <- 2
          filter.edata.anova <- anova_filter(nonmiss_per_group=nonmiss_per_group, min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_cname)
          filter.edata <- intersect(filter.edata.anova, filter.edata.gtest)
        }
      }
    }
    
    #checking that filter.edata does not specify all of e_data in omicsData
    if(all(omicsData$e_data[[edata_cname]] %in% filter.edata)){stop("filter.edata specifies all of e_data in omicsData")}
    
    if(length(filter.edata) < 1){ # if-statement added by KS 1/12/2018 to catch the case where nothing is removed
      filter_object_new = list(edata_filt = NULL, emeta_filt = NULL, samples_filt = NULL)
    }else{
      filter_object_new = list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL)
    }
    

    # call the function that does the filter application
    results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    omicsData$e_data <- results_pieces$temp.pep2
    omicsData$f_data <- results_pieces$temp.samp2
    omicsData$e_meta <- results_pieces$temp.meta1
    results <- omicsData

    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }

    # set attributes for which filters were run
    attr(results, "filters")$imdanovaFilt <- list(filter_method = NULL, report_text = "", threshold = c(), filtered = c())
    if(filter_method == "anova"){
      attr(results, "filters")$imdanovaFilt$filter_method <- "anova"
      attr(results, "filters")$imdanovaFilt$report_text <- paste("An ANOVA filter was applied to the data, removing ", edata_cname, "s that do not have at least ", min_nonmiss_anova, " non-missing values per group. A total of ", length(filter.edata), " ", edata_cname, "s were filtered out of the dataset by this filter.", sep="")
      attr(results, "filters")$imdanovaFilt$threshold <- data.frame(min_nonmiss_anova=min_nonmiss_anova, min_nonmiss_gtest=NA)
      attr(results, "filters")$imdanovaFilt$filtered <- filter.edata

      # add nonmissing per group, character vectors of the peptides that can be tested with anova and that can be tested with the gtest (default to NULL)
      attr(results, "imdanova")$nonmiss_per_group <- nonmiss_per_group
      attr(results, "imdanova")$test_with_anova <- setdiff(x=omicsData$e_data[, edata_cname], y=filter.edata)
      attr(results, "imdanova")$test_with_gtest <- NULL
    }else{
      if(filter_method == "gtest"){
        attr(results, "filters")$imdanovaFilt$filter_method <- "gtest"
        attr(results, "filters")$imdanovaFilt$report_text <- paste("An IMD (independence of missing data) filter was applied to the data, removing ", edata_cname, "s that do not have at least ", min_nonmiss_gtest, " non-missing values in at least one of the groups specified in the group_DF attribute of the dataset. A total of ", length(filter.edata), " ", edata_cname, "s were filtered out of the dataset by this filter.", sep="")
        attr(results, "filters")$imdanovaFilt$threshold <- data.frame(min_nonmiss_anova=NA, min_nonmiss_gtest=min_nonmiss_gtest)
        attr(results, "filters")$imdanovaFilt$filtered <- filter.edata

        # add nonmissing per group, character vectors of the peptides that can be tested with anova and that can be tested with the gtest (default to NULL)
        attr(results, "imdanova")$nonmiss_per_group <- nonmiss_per_group
        attr(results, "imdanova")$test_with_anova <- NULL
        attr(results, "imdanova")$test_with_gtest <- setdiff(x=omicsData$e_data[, edata_cname], y=filter.edata)
      }else{
        if(filter_method == "combined"){
          attr(results, "filters")$imdanovaFilt$filter_method <- c("anova", "gtest")
          attr(results, "filters")$imdanovaFilt$report_text <- paste("An ANOVA filter was applied to the data, removing ", edata_cname, "s that do not have at least ", min_nonmiss_anova, " non-missing values per group. Additionally, an IMD (independence of missing data) filter was applied to the data, removing ", edata_cname, "s that do not have at least ", min_nonmiss_gtest, " non-missing values in at least one of the groups specified in the group_DF attribute of the dataset. A total of ", length(filter.edata), " ", edata_cname, "s were filtered out of the dataset by this filter.", sep="")
          attr(results, "filters")$imdanovaFilt$threshold <- data.frame(min_nonmiss_anova=min_nonmiss_anova, min_nonmiss_gtest=min_nonmiss_gtest)
          attr(results, "filters")$imdanovaFilt$filtered <- filter.edata

          # add nonmissing per group, character vectors of the peptides that can be tested with anova and that can be tested with the gtest (default to NULL)
          attr(results, "imdanova")$nonmiss_per_group <- nonmiss_per_group
          attr(results, "imdanova")$test_with_anova <- setdiff(x=omicsData$e_data[, edata_cname], y=filter.edata.anova)
          attr(results, "imdanova")$test_with_gtest <- setdiff(x=omicsData$e_data[, edata_cname], y=filter.edata.gtest)
        }
      }
    }


  }

  return(results)
}




# function for customFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.customFilt <- function(filter_object, omicsData){


  edata_cname <- attributes(omicsData)$cnames$edata_cname
  fdata_cname <- attributes(omicsData)$cnames$fdata_cname
  emeta_cname <- attributes(omicsData)$cnames$emeta_cname

  #if filter_object contains removes
  if(!is.null(filter_object$e_data_remove)||!is.null(filter_object$f_data_remove)||!is.null(filter_object$e_meta_remove)){
    
    filter_object_new = list(edata_filt = filter_object$e_data_remove, emeta_filt = filter_object$e_meta_remove, samples_filt = filter_object$f_data_remove)
    
    # check that edata_filt doesn't specify ALL the items in omicsData #
    if(all(omicsData$e_data[, edata_cname] %in% filter_object_new$edata_filt)){stop("edata_filt specifies all the items in the data")}
    
    # check that samples_filt doesn't specify ALL the items in omicsData #
    if(all(omicsData$f_data[, fdata_cname] %in% filter_object_new$samples_filt)){stop("samples_filt specifies all the items in the data")}
    
    # check that emeta_filt doesn't specify ALL the items in omicsData, emeta_filt is present #
    if(!is.null(omicsData$e_meta[, emeta_cname])){
       if(all(omicsData$e_meta[, emeta_cname] %in% filter_object_new$emeta_filt)){stop("emeta_filt specifies all the items in the data")}
    }
   
    
  }
  
  else{
    filter_object_new = list(edata_keep = filter_object$e_data_keep, emeta_keep = filter_object$e_meta_keep, samples_keep = filter_object$f_data_keep)
    
    # check that edata_keep doesn't specify ALL the items in omicsData #
    if(all(omicsData$e_data[, edata_cname] %in% filter_object_new$edata_keep)){stop("edata_keep specifies all the items in the data")}
    
    # check that samples_keep doesn't specify ALL the items in omicsData #
    if(all(omicsData$f_data[, fdata_cname] %in% filter_object_new$samples_keep)){stop("samples_keep specifies all the items in the data")}
    
    # check that emeta_keep doesn't specify ALL the items in omicsData #
    if(!is.null(omicsData$e_meta[, emeta_cname])){
    if(all(omicsData$e_meta[, emeta_cname] %in% filter_object_new$emeta_keep)){stop("emeta_keep specifies all the items in the data")}
    
    }
  }
  
  # call the function that does the filter application
  results_pieces <- MSomics_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

  # return filtered data object #
  results <- omicsData
  results$e_data <- results_pieces$temp.pep2
  results$f_data <- results_pieces$temp.samp2
  if(!is.null(omicsData$e_meta)){ # if-statement added by Kelly 3/24/2017 #
    results$e_meta <- data.frame(results_pieces$temp.meta1)
    names(results$e_meta)[which(names(omicsData$e_meta) == emeta_cname)] <- emeta_cname
  }else{
   # e_meta is null
    results$e_meta <- NULL
  }
  
  

  # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure #
  if(!is.null(attr(results, "group_DF"))){
    results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
  }else{
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if(!is.null(results$e_meta)){
      # number of unique proteins that map to a peptide in e_data #
      if(!is.null(emeta_cname)){
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }


  # check to see whether a customFilt has already been run on omicsData #
  if("customFilt" %in% names(attributes(omicsData)$filters)){
    # yes, so be sure to keep track of this customFilt in a separate $filters attribute #

    # get number of previous customFilts applied #
    n_prev_filts <- length(grep("customFilt", names(attributes(omicsData)$filters)))
    # will need to append this number + 1 to the end of the $filters$customFilt attribute #


    # set attributes for which filters were run
    attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]] <- list(report_text = "", threshold = c(), filtered = c())

    n_edata_filtered <- nrow(omicsData$e_data) - nrow(results$e_data)
    n_fdata_filtered <- nrow(omicsData$f_data) - nrow(results$f_data)
    if(!is.null(omicsData$e_meta)){
      n_emeta_filtered <- nrow(omicsData$e_meta) - nrow(results$e_meta)
    }else{
      n_emeta_filtered = NA
    }
    if(!is.null(omicsData$e_meta)){
      attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data, ", n_emeta_filtered, " ", emeta_cname, "s from e_meta, and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }else{
      attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }
    attr(results, "filters")[[(paste("customFilt", n_prev_filts + 1, sep=""))]]$filtered <- filter_object

  }else{ # no previous customFilt, so go ahead and name the attribute like normal #

    # set attributes for which filters were run
    attr(results, "filters")$customFilt <- list(report_text = "", threshold = c(), filtered = c())
    n_edata_filtered <- nrow(omicsData$e_data) - nrow(results$e_data)
    n_fdata_filtered <- nrow(omicsData$f_data) - nrow(results$f_data)
    if(!is.null(omicsData$e_meta)){
      n_emeta_filtered <- nrow(omicsData$e_meta) - nrow(results$e_meta)
    }else{
      n_emeta_filtered = NA
    }
    if(!is.null(omicsData$e_meta)){
      attr(results, "filters")$customFilt$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data, ", n_emeta_filtered, " ", emeta_cname, "s from e_meta, and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }else{
      attr(results, "filters")$customFilt$report_text <- paste("A custom filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }
    attr(results, "filters")$customFilt$filtered <- filter_object
  }
  return(results)
}









#' Remove items that need to be filtered out
#'
#' This function removes
#'
#' @param omicsData an object of the class \code{pepData}, \code{proData}, \code{lipidData}, or \code{metabData}, usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or \code{\link{as.metabData}}, respectively.
#' @param filter_object a list created by the functions above
#' @return list
#' @author Kelly Stratton, Lisa Bramer
#'
MSomics_filter_worker <- function(filter_object, omicsData){
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname
  
  # pull group_DF attribute #
  group_DF = attr(omicsData, "group_DF")
  
  #check if filter object contains remove arguments
  if(!is.null(filter_object$edata_filt)||!is.null(filter_object$emeta_filt)||!is.null(filter_object$samples_filt)){
    
    ## check to see if e_meta is provided ##
    # if not provided we only need to worry about e_data and f_data #
    if(attr(omicsData, "meta_info") == FALSE){
      
      ## remove entries in edata ##
      if(!is.null(filter_object$edata_filt) & !is.null(edata_cname)){
        
        temp.pep = omicsData$e_data
        
        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)
        
        if(length(edat_ids) > 0){
          # identify which peptides in e_data match filter list and remove #
          temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.pep1 = temp.pep}
        
      }else{ # no entries in edata need to be removed
        temp.pep1 = omicsData$e_data
      }
      
      ## remove samples ##
      if(!is.null(filter_object$samples_filt) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data
        
        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_filt)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)
        
        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[-which(temp.samp[,samp_cname] %in% filter_object$samples_filt),]
        }else{temp.samp2 = temp.samp}
        
        # identify which samples in e_data match filter list and remove #
        if(length(edat_ids2) > 0){
          temp.pep2 = temp.pep1[, -which(names(temp.pep1) %in% filter_object$samples_filt)]
        }else{temp.pep2 = temp.pep1}
        
      }else{ # no entries in f_data need to be removed
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }
      
      temp.meta2 = NULL
      
      
    }else{ # e_meta is present, so we need to work with it as well
      ## remove entries in edata ##
      if(!is.null(filter_object$edata_filt) & !is.null(edata_cname)){
        
        temp.pep = omicsData$e_data
        
        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)
        
        if(length(edat_ids) > 0){
          # identify which peptides in e_data and e_meta match filter list and remove#
          temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.pep1 = temp.pep}
        
        temp.meta = omicsData$e_meta
        
        # check that at least one of the peptides is present in e_meta #
        emeta_ids = which(temp.meta[,edata_cname] %in% filter_object$edata_filt)
        
        if(length(emeta_ids) > 0){
          temp.meta1 = temp.meta[-which(temp.meta[,edata_cname] %in% filter_object$edata_filt),]
        }else{temp.meta1 = temp.meta}
        
      }else{
        temp.pep1 = omicsData$e_data
        temp.meta1 = omicsData$e_meta
      }
      
      ## remove samples ##
      if(!is.null(filter_object$samples_filt) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data
        
        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_filt)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)
        
        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[-which(temp.samp[,samp_cname] %in% filter_object$samples_filt),]
        }else{temp.samp2 = temp.samp}
        
        # identify which samples in e_data match filter list and remove #
        if(length(edat_ids2) > 0){
          inds = which(names(temp.pep1) %in% filter_object$samples_filt)
          temp.pep2 = temp.pep1[, -inds]
        }else{temp.pep2 = temp.pep1}
        
      }else{
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }
      
      ## remove entries in emeta ##
      if(!is.null(filter_object$emeta_filt) & !is.null(emeta_cname)){
        # identify which proteins in data match filter list and remove from e_meta #
        temp.meta = temp.meta1
        
        # check that at least one of the proteins is in e_meta #
        if(!is.null(ncol(temp.meta))){
          emeta_ids2 = which(as.character(temp.meta[,emeta_cname]) %in% filter_object$emeta_filt)
        }else{
          emeta_ids2 = which(temp.meta %in% filter_object$emeta_filt)
        }
        
        
        if(length(emeta_ids2) > 0){
          if(!is.null(ncol(temp.meta))){
            temp.meta2 = temp.meta[-which(temp.meta[,emeta_cname] %in% filter_object$emeta_filt),]
          }else{
            temp.meta2 = temp.meta[-which(temp.meta %in% filter_object$emeta_filt)]
          }
        }else{temp.meta2 = temp.meta}
      }else{
        temp.meta2 = temp.meta1
      }
      
      
      # check for rogue entries in edata #
      if(!is.null(ncol(temp.meta2))){
        edat_ids2 = which(!(temp.pep2[,edata_cname] %in% temp.meta2[,edata_cname]))
      }else{
        edat_ids2 = which(!(temp.pep2[,edata_cname] %in% temp.meta2))
      }
      
      
      # filter out edata entries which no longer have mappings to emeta entries #
      if(length(edat_ids2) > 0){
        #temp.pep2 = temp.pep2[-which(!(temp.pep2[,edata_cname] %in% temp.meta2[,edata_cname])),]
        temp.pep2 = temp.pep2[-edat_ids2,]
      }
      
    }
    
    output <- list(temp.pep2 = temp.pep2, temp.samp2 = temp.samp2, temp.meta1 = temp.meta2, edata_cname = edata_cname, emeta_cname = emeta_cname, samp_cname = samp_cname)
  }
  
  
  #if filter object contains keep arguments
  else{
    
    ## check to see if e_meta is provided ##
    # if not provided we only need to worry about e_data and f_data #
    if(attr(omicsData, "meta_info") == FALSE){
      
      ## keep entries in edata ##
      if(!is.null(filter_object$edata_keep) & !is.null(edata_cname)){
        
        temp.pep = omicsData$e_data
        
        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_keep)
        
        if(length(edat_ids) > 0){
          # identify which peptides in e_data match filter list and keep #
          temp.pep1 = temp.pep[which(temp.pep[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.pep1 = temp.pep}
        
      }else{ # no entries in edata need to be removed
        temp.pep1 = omicsData$e_data
      }
      
      ## keep samples ##
      if(!is.null(filter_object$samples_keep) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data
        
        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_keep)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_keep)
        
        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[which(temp.samp[,samp_cname] %in% filter_object$samples_keep),]
        }else{temp.samp2 = temp.samp}
        
        # identify which samples in e_data match filter list and keep #
        if(length(edat_ids2) > 0){
          edata_cname_id = which(names(temp.pep1) == edata_cname)
          temp.pep2 = temp.pep1[,c(edata_cname_id,(which(names(temp.pep1) %in% filter_object$samples_keep)))]
        }else{temp.pep2 = temp.pep1}
        
      }else{ # no entries in f_data need to be removed
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }
      
      temp.meta2 = NULL
      
      
      
      
    }else{ # e_meta is present, so we need to work with it as well
      ## keep entries in edata ##
      if(!is.null(filter_object$edata_keep) & !is.null(edata_cname)){
        
        temp.pep = omicsData$e_data
        
        # have to check that at least one of the items is present in the data #
        edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_keep)
        
        if(length(edat_ids) > 0){
          # identify which peptides in e_data and e_meta match filter list and keep#
          temp.pep1 = temp.pep[which(temp.pep[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.pep1 = temp.pep}
        
        temp.meta = omicsData$e_meta
        
        # check that at least one of the peptides is present in e_meta #
        emeta_ids = which(temp.meta[,edata_cname] %in% filter_object$edata_keep)
        
        if(length(emeta_ids) > 0){
          temp.meta1 = temp.meta[which(temp.meta[,edata_cname] %in% filter_object$edata_keep),]
        }else{temp.meta1 = temp.meta}
        
      }else{
        temp.pep1 = omicsData$e_data
        temp.meta1 = omicsData$e_meta
      }
      
      ## keep samples ##
      if(!is.null(filter_object$samples_keep) & !is.null(samp_cname)){
        # identify which samples in f_data match filter list #
        temp.samp = omicsData$f_data
        
        # check that at least one sample is in f_data and e_data #
        fdat_ids = which(temp.samp[,samp_cname] %in% filter_object$samples_keep)
        edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_keep)
        
        if(length(fdat_ids) > 0){
          temp.samp2 = temp.samp[which(temp.samp[,samp_cname] %in% filter_object$samples_keep),]
        }else{temp.samp2 = temp.samp}
        
        # identify which samples in e_data match filter list and keep #
        if(length(edat_ids2) > 0){
          edata_cname_id = which(names(temp.pep1) == edata_cname)
          
          inds = which(names(temp.pep1) %in% filter_object$samples_keep)
          temp.pep2 = temp.pep1[,c(edata_cname_id,inds)]
        }else{temp.pep2 = temp.pep1}
        
      }else{
        temp.samp2 = omicsData$f_data
        temp.pep2 = temp.pep1
      }
      
      ## keep entries in emeta ##
      if(!is.null(filter_object$emeta_keep) & !is.null(emeta_cname)){
        # identify which proteins in data match filter list and keep in e_meta #
        temp.meta = temp.meta1
        temp.meta_not_kept = omicsData$e_meta[-which(omicsData$e_meta[[edata_cname]] %in% temp.meta[[edata_cname]]),]
        
        # check that at least one of the proteins is in e_meta (this is e_meta after e_data_keep has been applied) #
        if(!is.null(ncol(temp.meta))){
          emeta_ids2 = which(as.character(temp.meta[,emeta_cname]) %in% filter_object$emeta_keep)
        }else{
          emeta_ids2 = which(temp.meta %in% filter_object$emeta_keep)
        }
        
        # check that at least one of the proteins is in temp_meta_not_kept #
        if(!is.null(ncol(temp.meta_not_kept))){
          emeta_ids_not_kept = which(as.character(temp.meta_not_kept[,emeta_cname]) %in% filter_object$emeta_keep)
        }else{
          emeta_ids_not_kept = which(temp.meta_not_kept %in% filter_object$emeta_keep)
        }
        
        #if there are e_meta_keep (proteins) in the part of e_meta that we previously kept, then these proteins have already been kept
        if(length(emeta_ids2) > 0){
          if(!is.null(ncol(temp.meta))){
            temp.meta2 = temp.meta
          }else{
            temp.meta2 = temp.meta
          }
        }else{temp.meta2 = temp.meta}
        
        #if there are e_meta_keep(proteins) outside of e_meta that we previously kept we will keep these 
        if(length(emeta_ids_not_kept) > 0){
          if(!is.null(ncol(temp.meta_not_kept))){
            temp.meta3 = temp.meta_not_kept[which(temp.meta_not_kept[,emeta_cname] %in% filter_object$emeta_keep),]
            temp.meta2 = rbind(temp.meta2,temp.meta3)
          }else{
            temp.meta3 = temp.meta_not_kept[which(temp.meta_not_kept %in% filter_object$emeta_keep)]
            temp.meta2 = rbind(temp.meta2, temp.meta3)
          }
        }else{temp.meta2 = temp.meta}
        
      }else{
        temp.meta2 = temp.meta1
      }
      
      
      # check for entries in e_meta[,edata_cname] that are not in e_data[,edata_cname]#
      if(!is.null(ncol(temp.meta2))){
        edat_ids2 = which(!temp.meta2[,edata_cname] %in% (temp.pep2[,edata_cname]))
      }else{
        edat_ids2 = which(!(temp.meta2 %in% temp.pep2[,edata_cname]))
      }
      
      
      # add edata entries which were present in emeta but not edata #
      if(length(edat_ids2) > 0){
        additional_peps<- temp.meta2[[edata_cname]][edat_ids2]
        edata_cname_id = which(names(temp.pep1) == edata_cname)
        
        if(is.null(filter_object$samples_keep)){
          inds = which(names(temp.pep1) %in% temp.samp2[,samp_cname]) 
        }
        else inds = which(names(temp.pep1) %in% filter_object$samples_keep)
        
        temp.pep2 = rbind(temp.pep2, omicsData$e_data[which(omicsData$e_data[[edata_cname]] %in% additional_peps) ,c(edata_cname_id,inds)])
      }
      
      
      
    }
    
    output <- list(temp.pep2 = temp.pep2, temp.samp2 = temp.samp2, temp.meta1 = temp.meta2, edata_cname = edata_cname, emeta_cname = emeta_cname, samp_cname = samp_cname)
    
  }
  
  
  # return the pieces needed to assemble a proData/pepData/lipidData/metabData object
  return(output)
}
