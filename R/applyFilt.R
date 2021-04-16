#' Apply a S3 filter  object to a pmartR S3 object
#'
#' This function takes a filter object of class 'cvFilt', 'rmdFilt', 
#'   'moleculeFilt', 'proteomicsFilt', 'imdanovaFilt', or 'customFilt' and 
#'   applies the filter to a dataset of \code{pepData}, \code{proData}, 
#'   \code{lipidData}, \code{metabData}, or \code{nmrData}.
#'   
#' @param filter_object an object of the class 'cvFilt', 'proteomicsFilt',
#'   'rmdFilt', 'moleculeFilt', 'imdanovaFilt', or 'customFilt' created by
#'   \code{cv_filter}, \code{proteomics_filter}, \code{rmd_filter},
#'   \code{molecule_filter}, or \code{imdanova_filter}, respectively.
#' @param omicsData an object of the class \code{pepData}, \code{proData},
#'   \code{lipidData}, \code{metabData}, or \code{nmrData} usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.metabData}}, or 
#'   \code{\link{as.nmrData}}, respectively.
#' @param ... further arguments
#'
#' @return An object of the class \code{pepData}, \code{proData},
#'   \code{lipidData}, \code{metabData}, or \code{nmrData} with specified
#'   cname_ids,
#'   edata_cnames, and emeta_cnames filtered out of the appropriate datasets.
#'
#' @details Various further arguments can be specified depending on the class of
#'   the \code{filter_object} being applied. For a \code{filter_object} of type
#'   'moleculeFilt': \tabular{ll}{ \code{min_num} \tab an integer value
#'   specifying the minimum number of times each biomolecule must be observed
#'   across all samples in order to retain the biomolecule. Default value is 2.
#'   \cr } For a \code{filter_object} of type 'cvFilt': \tabular{ll}{
#'   \code{cv_threshold} \tab an integer value specifying the maximum
#'   coefficient of variation (CV) threshold for the biomolecules. Biomolecules
#'   with CV greater than this threshold will be filtered. Default is 150. \cr }
#'   For a \code{filter_object} of type 'rmdFilt': \tabular{ll}{
#'   \code{pvalue_threshold} \tab numeric value between 0 and 1, specifying the
#'   p-value, below which samples will be removed from the dataset. Default is
#'   0.001. \cr \code{min_num_biomolecules} \tab numeric value greater than 10
#'   (preferably greater than 50) that specifies the minimum number of
#'   biomolecules that must be present in order to create an rmdFilt object.
#'   Using values less than 50 is not advised. \cr } For a \code{filter_object}
#'   of type 'proteomicsFilt': \tabular{ll}{ \code{min_num_peps} \tab an
#'   optional integer value between 1 and the maximum number of peptides that
#'   map to a protein in omicsData. The value specifies the minimum number of
#'   peptides that must map to a protein. Any protein with less than
#'   \code{min_num_peps} mapping to it will be returned as a protein that should
#'   be filtered. Default value is NULL. \cr \code{degen_peps} \tab logical
#'   indicator of whether to filter out degenerate peptides (TRUE) or not
#'   (FALSE). Default value is FALSE.\cr } For a \code{filter_object} of type
#'   'imdanovaFilt': \tabular{ll}{ \code{min_nonmiss_anova} \tab integer value
#'   specifying the minimum number of non-missing feature values allowed per
#'   group for \code{anova_filter}. Default value is 2. \cr
#'   \code{min_nonmiss_gtest} \tab integer value specifying the minimum number
#'   of non-missing feature values allowed per group for \code{gtest_filter}.
#'   Default value is 3.\cr } There are no further arguments for a
#'   \code{filter_object} of type 'customFilt'.
#'   
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' to_filter <- molecule_filter(omicsData = pep_object)
#' pep_object2 <- applyFilt(filter_object = to_filter,
#'                          omicsData = pep_object, min_num = 2)
#' print(str(attributes(pep_object2)$filters))
#' pep_object2 <- group_designation(pep_object2, main_effects = "Condition")
#' to_filter2 <- imdanova_filter(omicsData = pep_object2)
#' pep_object3 <- applyFilt(filter_object = to_filter2,
#'                          omicsData = pep_object2,
#'                          min_nonmiss_anova = 3)
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
#' 
applyFilt <- function (filter_object, omicsData, ...) {
  
  # check that omicsData is of pmartR S3 class#
  if (!inherits(omicsData, c("pepData", "proData", "lipidData", "metabData",
                             "nmrData"))) {
    
    # Throw down an error for blatant abuse of the pmartR standards.
    stop (paste("omicsData must be of class 'pepData', 'proData', 'lipidData',",
                "'metabData' or 'nmrData'",
                sep = " "))
    
  }
  
  # check that filter_object is of an appropriate class#
  if (!inherits(filter_object, c("cvFilt", "proteomicsFilt", "moleculeFilt",
                                 "rmdFilt", "imdanovaFilt", "customFilt"))) {
    
    # Throw an error for unholy argument specification.
    stop (paste("filter_object must be of  'cvFilt', 'proteomicsFilt',",
                "'moleculeFilt', 'rmdFilt', 'imdanovaFilt', or 'customFilt.",
                sep = " "))
    
  }
  
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname
  
  # generate warnings if data type is not present but user asks to filter #
  if (!is.null(filter_object$filter_edata) & is.null(edata_cname)) {
    
    warning (paste("e_data identifier column specified in filter_object is not",
                   "present in omicsData$e_data. Specified e_data features",
                   "cannot be filtered from data.",
                   sep = " "))
    
  }
  
  if (!is.null(filter_object$filter_emeta) & is.null(emeta_cname)) {
    
    warning (paste("e_meta identifier column specified in filter_object is not",
                   "present in omicsData$e_meta. Specified e_meta features",
                   "cannot be filtered from data.",
                   sep = " "))
    
  }
  
  if (!is.null(filter_object$filter_samples) & is.null(samp_cname)) {
    
    warning (paste("Sample identifier column specified in filter_object is not",
                   "present in omicsData. Specified samples cannot be filtered",
                   "from data.",
                   sep = " "))
    
  }
  
  
  UseMethod ("applyFilt")
  
}

# function for moleculeFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.moleculeFilt <- function(filter_object, omicsData, min_num=2){
  
  # Perform some initial checks on the input arguments -------------------------
  
  # Check if a molecule filter has already been applied.
  if ('moleculeFilt' %in% get_filter_type(omicsData)) {
    
    # Gently tell the user they have already applied a molecule filter.
    warning ('A molecule filter has already been applied to this data set.')
    
  }
  
  # check that min_num is of length 1 #
  if (length(min_num) != 1) {
    
    # Warn the user of their treachery with an error.
    stop ("min_num must be of length 1")
    
  }
  
  # The min_num argument must be an integer and >= 1.
  if (min_num %% 1 != 0 || min_num < 1) {
    
    # Denounce the users heinous actions with an error.
    stop ("min_num must be an integer greater than or equal to 1")
    
  }
  
  # check that min_num is less than the total number of samples #
  if (min_num >= attr(filter_object, "num_samps")) {
    
    # Throw an error for the users audacity to besmirch the arguments in such a
    # foul manner.
    stop ("min_num cannot be greater than or equal to the number of samples")
    
  }
  
  # Prepare the information needed to filter the data --------------------------
  
  # Extract the column number containing the identifiers.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Sniff out the indices that fall below the threshold.
  inds <- which(filter_object$Num_Observations < min_num)
  
  # Compute the length of the inds vector and specify filter.edata accordingly.
  if (length(inds) < 1) {
    
    # Set filter.edata to NULL because no rows in omicsData$e_data will be
    # filtered out.
    filter.edata <- NULL
    
  } else {
    
    # Fish out the identifiers that will be filtered.
    filter.edata <- omicsData$e_data[, id_col][inds]
    
  }
  
  # Create a list that is used in the pmartR_filter_worker function.
  filter_object_new <- list(e_data_remove = filter.edata,
                            e_meta_remove = NULL,
                            f_data_remove = NULL)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # Check if group_DF attribute is present.
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
      
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = min_num,
    filtered = filter.edata,
    method = NA
  )
  
  # RETURN THE FILTERED OMICSDATA OBJECT! YAY!!!
  return(omicsData)
  
}

# function for cvFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.cvFilt <- function (filter_object, omicsData, cv_threshold = 150) {
  
  # Perform some initial checks on the input arguments -------------------------
  
  # Check if a CV filter has already been applied.
  if ('cvFilt' %in% get_filter_type(omicsData)) {
    
    # Slap the users wrist with a warning.
    warning ('A CV filter has already been applied to this data set.')
    
  }
  
  # Check that cv_threshold is a scalar.
  if (length(cv_threshold) != 1) {
    
    # Stop the treacherous snake with an error!
    stop ("cv_threshold must be a scalar.")
    
  }
  
  # check that cv_threshold is greater than 0 #
  if (cv_threshold <= 0) stop ("cv_threshold must be greater than 0.")
  
  # Check that cv_threshold is less than the largest CV value.
  if (cv_threshold > max(na.omit(filter_object$CV))) {
    
    # Throw an error because the user does not understand basic number line
    # protocol.
    stop ("cv_threshold cannot be greater than the largest CV value.")
    
  }
  
  # Prepare the information needed to filter the data --------------------------
  
  # Extract the column number containing the identifiers.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Sniff out the indices that fall below the threshold.
  inds <- which(filter_object$CV > cv_threshold)
  
  # Compute the length of the inds vector and specify filter.edata accordingly.
  if (length(inds) < 1) {
    
    # Set filter.edata to NULL because no rows in omicsData$e_data will be
    # filtered out.
    filter.edata <- NULL
    
  } else {
    
    # Fish out the identifiers that will be filtered.
    filter.edata <- omicsData$e_data[, id_col][inds]
    
  }
  
  # Create a list that is used in the pmartR_filter_worker function.
  filter_object_new <- list(e_data_remove = filter.edata,
                            e_meta_remove = NULL,
                            f_data_remove = NULL)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # Check if group_DF attribute is present.
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = cv_threshold,
    filtered = filter.edata,
    method = NA
  )
  
  # We did it! Kudos to us!!
  return (omicsData)
  
}

# function for rmdFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.rmdFilt <- function (filter_object, omicsData,
                               pvalue_threshold = 0.0001,
                               min_num_biomolecules = 50) {
  
  # Perform some initial checks on the input arguments -------------------------
  
  # Check if an RMD filter has already been applied.
  if ('rmdFilt' %in% get_filter_type(omicsData)) {
    
    # Slap the users wrist with a warning.
    warning ('An RMD filter has already been applied to this data set.')
    
  }
  
  # check that pvalue_threshold is between 0 and 1 #
  if (pvalue_threshold < 0 || pvalue_threshold > 1) {
    
    # Figuratively smack the user with an error for not knowing what a p-value
    # is (or what values it can take).
    stop ("pvalue_threshold must be between 0 and 1.")
    
  }
  
  # Check that min_num_biomolecules is less than the number of observations.
  if (get_data_info(omicsData)$num_edata < min_num_biomolecules) {
    
    # Stop the user for trying to filter their entire sample. Maybe we should
    # also ask them why they went to all that effort in gathering the data if
    # they are just going to try to throw it all away.
    stop (paste("There are fewer biomolecules in omicsData than",
                "min_num_biomolecules",
                paste0("(", min_num_biomolecules, ")."),
                "See applyFilt for details.",
                sep = " "))
  }
  
  # Prepare the information needed to filter the data --------------------------
  
  # Extract the column number containing the sample IDs.
  id_col <- which(names(omicsData$f_data) == get_fdata_cname(omicsData))
  
  # Sniff out the indices that fall below the threshold.
  inds <- which(filter_object$pvalue < pvalue_threshold)
  
  # Compute the length of the inds vector and specify filter.samp accordingly.
  if (length(inds) < 1) {
    
    # Set filter.samp to NULL because no columns in omicsData$e_data will be
    # filtered out.
    filter.samp <- NULL
    
    # Check if the number of indices to filter is equal to the total number of
    # samples.
  } else if (length(inds) == dim(filter_object)[1]) {
    
    # Throw an error for trying to remove every sample.
    stop ("With the current p-value threshold all samples will be filtered.")
    
  } else {
    
    # Fish out the sample IDs that will be filtered.
    filter.samp <- as.character(filter_object[inds, id_col])
    
  }
  
  # Create a list that is used in the pmartR_filter_worker function.
  filter_object_new <- list(e_data_remove = NULL,
                            e_meta_remove = NULL,
                            f_data_remove = filter.samp)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # Check if group_DF attribute is present.
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
      
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = pvalue_threshold,
    filtered = filter.samp,
    method = NA
  )
  
  # Return the filtered data! High fives all around!!!
  return(omicsData)
  
}

# function for proteomicsFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
#' @export
applyFilt.proteomicsFilt <- function (filter_object,
                                      omicsData,
                                      min_num_peps = NULL,
                                      degen_peps = FALSE) {
  
  # Perform initial checks on the input arguments ------------------------------
  
  # Check if a proteomics filter has already been applied.
  if ('proteomicsFilt' %in% get_filter_type(omicsData)) {
    
    # Slap the users wrist with a warning.
    warning ('A proteomics filter has already been applied to this data set.')
    
  }
  
  # Make sure at least some of the peptides or proteins will be filtered. This
  # means min_num_peps must not be null or degen_peps must not be FALSE.
  if (is.null(min_num_peps) && !degen_peps) {
    
    # Warn the user that no filtering will actually occur.
    stop (paste("No peptides or proteins will be filtered. Either change",
                "min_num_peps to an integer > 1 or change degen_peps to TRUE.",
                sep = " "))
    
  }
  
  # error checks for min_num_peps, if not NULL #
  if (!is.null(min_num_peps)) {
    
    # check that min_num_peps is of length 1 #
    if (length(min_num_peps) != 1) {
      
      # Warn the user of their treachery with an error.
      stop ("min_num_peps must be of length 1")
      
    }
    
    # The min_num_peps argument must be an integer and > 1.
    if (min_num_peps %% 1 != 0 || min_num_peps <= 1) {
      
      # Denounce the users heinous actions with an error.
      stop ("min_num_peps must be an integer greater than 1")
      
    }
    
    # Check that min_num_peps is less than the total number of peptides.
    if (min_num_peps >= nrow(omicsData$e_data)) {
      
      # Throw an error for the users inability to perform basic counting
      # operations.
      stop (paste("min_num_peps cannot be greater than or equal to the number",
                  "of peptides",
                  sep = " "))
      
    }
    
  }
  
  # check that degen_peps is logical #
  if (!inherits(degen_peps, "logical")) {
    
    # Warn the illogical user that their inputs are also illogical.
    stop ("degen_peps must be either TRUE or FALSE")
    
  }
  
  # degenerate peptides portion ------------------------------------------------
  
  # Extract the name of the column that contains peptide IDs.
  pep_id = attr(omicsData, "cnames")$edata_cname
  
  # Extract the column name containing the protein IDs.
  pro_id = attr(omicsData, "cnames")$emeta_cname
  
  # Check if degen_peps is TRUE.
  if (degen_peps) {
    
    count_bypep <- filter_object$counts_by_pep
    
    # pull peptides with more than one row in e_meta.  
    degen_peptides <- as.character(
      data.frame(count_bypep[which(count_bypep$n > 1),])[, pep_id]
    )
    
    ## identify any proteins that now will not have peptides mapping to them ##
    ## find rows in e_meta that correspond to peptides to be filtered ##
    pepfilt.ids = which(omicsData$e_meta[, pep_id] %in% degen_peptides)
    
    # Fish out the indices for the proteins that will be filtered.
    add_prots <- as.character(setdiff(omicsData$e_meta[pepfilt.ids, pro_id],
                                      omicsData$e_meta[-pepfilt.ids, pro_id]))
    
    if (length(add_prots) == 0) {
      
      add_prots = NULL
      
    }
    
    if (length(degen_peptides) == 0) {
      
      degen_peptides = NULL
      
    }
    
    pepe <- list(e_data_remove = degen_peptides,
                 e_meta_remove = add_prots)
    
  } else {
    
    pepe <- list(e_data_remove = NULL,
                 e_meta_remove = NULL)
    
  }
  
  
  # protein filter portion -----------------------------------------------------
  
  # Check if min_num_peps is not NULL.
  if (!is.null(min_num_peps)) {
    
    # Find all rows with proteins that map to fewer peptides than min_num_peps.
    pro_idx <- which(filter_object$counts_by_pro$n < min_num_peps)
    
    # pull proteins with less than min_num_peps.
    pro_filt <- as.character(filter_object$counts_by_pro[pro_idx, pro_id])
    
    # Find rows in e_meta that correspond to proteins to be filtered.
    protfilt.ids = which(omicsData$e_meta[,pro_id] %in% pro_filt)
    
    # Extricate the indices for the peptides that will be filtered.
    pep_filt <- as.character(setdiff(omicsData$e_meta[protfilt.ids, pep_id],
                                     omicsData$e_meta[-protfilt.ids, pep_id]))
    
    if (length(pep_filt) == 0) {
      
      pep_filt = NULL
      
    }
    
    if (length(pro_filt) == 0) {
      
      pro_filt = NULL
      
    }
    
    pepe2 <- list(e_meta_remove = pro_filt,
                  e_data_remove = pep_filt)
    
  } else {
    
    pepe2 <- list(e_data_remove = NULL,
                  e_meta_remove = NULL)
    
  }
  
  # Consolidate pepe and pepe2 to pass to the pmartR_filter_worker function.
  filter_object_new <- list(e_meta_remove = unique(c(pepe$e_meta_remove, 
                                                     pepe2$e_meta_remove)),
                            e_data_remove = unique(c(pepe$e_data_remove, 
                                                     pepe2$e_data_remove)))
  
  # checking that filter_object_new does not specify all of e_data(peps) or all
  # of e_meta(protiens) in omicsData
  if (all(omicsData$e_meta[[pro_id]] %in% filter_object_new$e_meta_remove)) {
    
    stop ("filter_object specifies all proteins in e_meta")
    
  }
  
  if (all(omicsData$e_data[[pep_id]] %in% filter_object_new$e_data_remove)) {
    
    stop ("filter_object specifies all peps in e_data")
    
  }
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # Check if group_DF attribute is present.
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = data.frame(min_num_peps = min_num_peps,
                           degen_peps = as.character(degen_peps)),
    filtered = filter_object_new,
    method = NA
  )
  
  # Return the filtered data! Look at us go!!
  return(omicsData)
  
}

# function for imdanovaFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.imdanovaFilt <- function (filter_object,
                                    omicsData,
                                    min_nonmiss_anova = NULL, 
                                    min_nonmiss_gtest = NULL, 
                                    remove_singleton_groups = TRUE) {
  
  # #' @details If filter_method="combined" is specified, then both the
  # \code{anova_filter} and \code{gtest_filter} are applied to the data, and the
  # intersection of features from the two filters is the set returned. For
  # ANOVA, features that do not have at least \code{min_nonmiss_allowed} values
  # per group are candidates for filtering. For the G-test, features that do not
  # have at least \code{min_nonmiss_allowed} values per group are candidates for
  # filtering. The G-test is a test of independence, used here to test the null
  # hypothesis of independence between the number of missing values across
  # groups.
  
  
  # Perform initial checks on the input arguments ------------------------------
  
  # Check if an imdanova filter has already been applied.
  if ('imdanovaFilt' %in% get_filter_type(omicsData)) {
    
    # Slap the users wrist with a warning.
    warning ('An imdanova filter has already been applied to this data set.')
    
  }
  
  # Check the number of columns in filter_object. If there are only two columns
  # then there is only one non-singleton group and an imdanova filter cannot be
  # applied.
  if (dim(filter_object)[2] == 2) {
    
    # Throw an error for too few samples to run an imdanova filter.
    stop (paste("An IMD-ANOVA filter cannot be used because there is only one",
                "non-singleton group.",
                sep = " "))
    
  }
  
  # verify remove_singleton_groups is logical T/F #
  if (class(remove_singleton_groups) != "logical") {
    
    # Stop the illogical user in their tracks for using illogical inputs.
    stop ("remove_singletong_groups must be TRUE or FALSE")
    
  }
  
  # check that at least one of min_nonmiss_anova and min_nonmiss_gtest are
  # present #
  if (is.null(min_nonmiss_anova) & is.null(min_nonmiss_gtest)) {
    
    # Throw an error for not providing any filter criteria.
    stop (paste("At least one of min_nonmiss_anova and min_nonmiss_gtest must",
                "be present",
                sep = " "))
    
  }
  
  #### The following objects will be used in the remainder of the checks. ####
  
  # Fish out the group_DF attribute.
  groupDF <- attributes(omicsData)$group_DF
  
  # Extricate the nonsingleton group names.
  nonsingleton_groups <- attr(groupDF, 'nonsingleton_groups')
  
  # Extract the group sizes attribute.
  group_sizes <- attr(filter_object, "group_sizes")
  
  # Find the minimum group size.
  min_grp_size <- min(group_sizes$n_group[which(group_sizes$Group %in%
                                                  nonsingleton_groups)])
  
  #### End object creation for preliminary checking purposes. ####
  
  # Check if min_nonmiss_anova is present.
  if (!is.null(min_nonmiss_anova)) {
    
    # check that min_nonmiss_anova is of length 1 #
    if (length(min_nonmiss_anova) != 1) {
      
      # Warn the user of their treachery with an error.
      stop ("min_nonmiss_anova must be of length 1")
      
    }
    
    # The min_nonmiss_anova argument must be an integer and > 1.
    if (min_nonmiss_anova %% 1 != 0 || min_nonmiss_anova < 2) {
      
      # Denounce the users heinous actions with an error.
      stop ("min_nonmiss_anova must be an integer > 1")
      
    }
    
    # Ensure min_nonmiss_anova is not greater than the smallest group size.
    if (min_nonmiss_anova > min_grp_size) {
      
      # Throw an error for an illegal use of integers.
      stop ("min_nonmiss_anova cannot be greater than the minimum group size")
      
    }
    
  }
  
  # Check if min_nonmiss_gtest is present.
  if (!is.null(min_nonmiss_gtest)) {
    
    # check that min_nonmiss_gtest is of length 1 #
    if (length(min_nonmiss_gtest) != 1) {
      
      # Warn the user of their treachery with an error.
      stop ("min_nonmiss_gtest must be of length 1")
      
    }
    
    # The min_nonmiss_gtest argument must be an integer and > 1.
    if (min_nonmiss_gtest %% 1 != 0 || min_nonmiss_gtest < 2) {
      
      # Denounce the users heinous actions with an error.
      stop ("min_nonmiss_gtest must be an integer > 1")
      
    }
    
    # Ensure min_nonmiss_gtest is not greater than the smallest group size.
    if (min_nonmiss_gtest > min_grp_size) {
      
      # Throw an error for an illegal use of integers.
      stop ("min_nonmiss_gtest cannot be greater than the minimum group size")
      
    }
    
  }
  
  # Prepare the data to be filtered --------------------------------------------
  
  # Create a variable that will be used to determine what should be included in
  # the filtered attribute at the end of the function. If there are singleton
  # groups and if they were removed this value will be changed to TRUE.
  singletons <- FALSE
  
  # if remove_singleton_groups is TRUE, filter those out now #
  if (remove_singleton_groups) {
    
    # Check for singleton groups.
    if (any(group_sizes$n_group == 1)) {
      
      # Change the singleton variable to TRUE.
      singletons <- TRUE
      
      # which group(s) #
      singleton_groups <- group_sizes$Group[group_sizes$n_group == 1]
      
      # which sample name(s) #
      samps_to_rm <- as.character(
        groupDF[which(groupDF$Group %in% singleton_groups),
                get_fdata_cname(omicsData)]
      )
      
      # use custom_filter to remove the sample(s) from the omicsData object #
      my_cust_filt <- custom_filter(omicsData,
                                    f_data_remove = samps_to_rm)
      omicsData <- applyFilt(my_cust_filt, omicsData)
      
      # Check the length of the filters attribute.
      if (length(attr(omicsData, "filters")) == 1) {
        
        # Reset the filters attribute to an empty list because no other filters
        # have been applied previously.
        attr(omicsData, "filters") <- list()
        
        # The following will run if the filters attribute is greater than 1. It
        # cannot be 0 because a custom filter was just applied.
      } else {
        
        # Get the total number of filters previously applied. This is the length
        # of the filters attribute.
        n_filtas <- length(attr(omicsData, "filters"))
        
        # Remove the custom filter from the filters attribute.
        attr(omicsData, "filters") <- attr(omicsData, "filters")[-n_filtas]
        
      }
      
      # Fish out the updated group_DF attribute.
      groupDF <- attributes(omicsData)$group_DF
      
      # update the IMD-ANOVA filter attributes.
      attributes(filter_object)$group_sizes <-
        group_sizes[-which(group_sizes$Group %in% singleton_groups), ]
      
      # give a message to user saying which samples have been removed #
      message(paste("You have specified remove_single_groups = TRUE, so we",
                    "have removed the following sample(s) that correspond to",
                    "singleton groups prior to proceeding with the IMD-ANOVA",
                    "filter: ",
                    sep = " "),
              paste(samps_to_rm, collapse = ", "))
      
      # Will run if no singleton groups are present.
    } else {
      
      message(paste("You have specified remove_singleton_groups = TRUE, but", 
                    "there are no singleton groups to remove.",
                    "Proceeding with application of the IMD-ANOVA filter.",
                    sep = " "))
    }
    
    # Will run if remove_singleton_gropus is FALSE.
  } else {
    
    # Check for singleton groups.
    if (any(attributes(filter_object)$group_sizes$n_group == 1)) {
      
      # Give a message that there are singleton groups present but the user has
      # elected to not remove them .
      message(paste("You have specified remove_singleton_groups = FALSE, so",
                    "the sample(s) corresponding to the singleton group(s)",
                    "were not utilized in the IMD-ANOVA filter and will be",
                    "retained in the resulting omicsData object.",
                    sep = " "))
      
    }
    
  }
  
  # Create a list with the appropriate information for the anova_filter and
  # gtest filter functions.
  nonmiss_per_group <- list(nonmiss_totals = filter_object,
                            group_sizes = attributes(filter_object)$group_sizes)
  
  edata_cname <- get_edata_cname(omicsData)
  samp_cname <- get_fdata_cname(omicsData)
  emeta_cname <- get_emeta_cname(omicsData)
  
  # Set the filter method according to the inputs min_nonmiss_anova and
  # min_nonmiss_gtest.
  if (!is.null(min_nonmiss_anova) && !is.null(min_nonmiss_gtest)) {
    
    filter_method <- "combined"
    
  } else if (!is.null(min_nonmiss_anova) && is.null(min_nonmiss_gtest)) {
    
    filter_method <- "anova"
    
  } else if (is.null(min_nonmiss_anova) && !is.null(min_nonmiss_gtest)) {
    
    filter_method <- "gtest"
    
  }
  
  # Apply the filter according to the filter type.
  if (filter_method == "anova") {
    
    filter.edata <- anova_filter(nonmiss_per_group = nonmiss_per_group,
                                 min_nonmiss_anova = min_nonmiss_anova,
                                 cname_id = edata_cname)
    
  } else if (filter_method == "gtest") {
    
    filter.edata <- gtest_filter(nonmiss_per_group = nonmiss_per_group,
                                 omicsData = omicsData,
                                 min_nonmiss_gtest = min_nonmiss_gtest)
    
  } else if(filter_method == "combined") {
    
    # Run the gtest filter.
    filter.edata.gtest <- gtest_filter(nonmiss_per_group = nonmiss_per_group,
                                       omicsData = omicsData,
                                       min_nonmiss_gtest = min_nonmiss_gtest)
    
    # Run the anova filter.
    filter.edata.anova <- anova_filter(nonmiss_per_group = nonmiss_per_group, 
                                       min_nonmiss_anova = min_nonmiss_anova, 
                                       cname_id = edata_cname)
    
    # Only filter samples found in both filters.
    filter.edata <- intersect(filter.edata.anova, filter.edata.gtest)
    
  }
  
  #checking that filter.edata does not specify all of e_data in omicsData
  if (all(omicsData$e_data[[edata_cname]] %in% filter.edata)) {
    
    stop ("All samples will be filtered. Try reducing the thresholds.")
    
  }
  
  # if-statement added by KS 1/12/2018 to catch the case where nothing is
  # removed
  if (length(filter.edata) < 1) { 
    
    # Stop the remainder of the function from running because none of the
    # samples will be filtered.
    stop ("None of the samples will be removed with the current thresholds.")
    
  } else {
    
    # Prepare the objects for filtering.
    filter_object_new = list(e_data_remove = filter.edata,
                             e_meta_remove = NULL,
                             f_data_remove = NULL)
    
  }
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(omicsData = omicsData,
                                         filter_object = filter_object_new)
  
  # return filtered data object #
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # Check if group_DF attribute is present.
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = data.frame(
      min_nonmiss_anova = if (is.null(min_nonmiss_anova)) {
        # Report NA because min_nonmiss_anova is NULL.
        NA
      } else {
        # Report the vaule used for min_nonmiss_anova.
        min_nonmiss_anova
      },
      min_nonmiss_gtest = if (is.null(min_nonmiss_gtest)) {
        # Report NA because min_nonmiss_gtest is NULL.
        NA
      } else {
        # Report the vaule used for min_nonmiss_gtest.
        min_nonmiss_gtest
      }
    ),
    filtered = if (singletons) {
      # Include the samples that were filtered with the custom filter along with
      # the biomolecules that were filtered with the IMD-ANOVA filter.
      list(samples = samps_to_rm,
           biomolecules = filter.edata)
    } else {
      # Only report the biomolecules because no samples were filtered.
      filter.edata
    },
    method = filter_method
  )
  
  # Create an additional attribute to the omicsData object that is specific to
  # the imdanova filter.
  attr(omicsData, "imdanova")$nonmiss_per_group <- nonmiss_per_group
  
  # Add information to the imdanova attribute according to the method used.
  if (filter_method == "anova") {
    
    # Add the biomolecule IDs that will NOT be filtered. These IDs will be used
    # for tests in later functions.
    attr(omicsData, "imdanova")$test_with_anova <- setdiff(
      x = omicsData$e_data[, edata_cname],
      y = filter.edata
    )
    
    attr(omicsData, "imdanova")$test_with_gtest <- NULL
    
  } else if (filter_method == "gtest") {
    
    attr(omicsData, "imdanova")$test_with_anova <- NULL
    
    # Add the biomolecule IDs that will NOT be filtered. These IDs will be used
    # for tests in later functions.
    attr(omicsData, "imdanova")$test_with_gtest <- setdiff(
      x = omicsData$e_data[, edata_cname],
      y = filter.edata
    )
    
  } else if (filter_method == "combined") {
    
    # Add the biomolecule IDs that will NOT be filtered. These IDs will be used
    # for tests in later functions.
    attr(omicsData, "imdanova")$test_with_anova <- setdiff(
      x = omicsData$e_data[, edata_cname],
      y = filter.edata.anova
    )
    
    # Add the biomolecule IDs that will NOT be filtered. These IDs will be used
    # for tests in later functions.
    attr(omicsData, "imdanova")$test_with_gtest <- setdiff(
      x = omicsData$e_data[, edata_cname],
      y = filter.edata.gtest
    )
    
  }
  
  # It was rough but we made it through!
  return(omicsData)
  
}

# function for customFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.customFilt <- function (filter_object, omicsData) {
  
  # Perform initial checks on the input arguments ------------------------------
  
  # Check if a custom filter has already been applied.
  if ('customFilt' %in% get_filter_type(omicsData)) {
    
    # Slap the users wrist with a warning.
    warning ('A custom filter has already been applied to this data set.')
    
  }
  
  # Prepare the data to be filtered --------------------------------------------
  
  # if filter_object contains 'removes' #
  if (!is.null(c(filter_object$e_data_remove,
                 filter_object$f_data_remove,
                 filter_object$e_meta_remove))) {
    
    filter_object_new = list(e_data_remove = filter_object$e_data_remove,
                             e_meta_remove = filter_object$e_meta_remove,
                             f_data_remove = filter_object$f_data_remove)
    
    # filter_object contains 'keeps' #
  } else {
    
    filter_object_new = list(e_data_keep = filter_object$e_data_keep,
                             e_meta_keep = filter_object$e_meta_keep,
                             f_data_keep = filter_object$f_data_keep)
    
  }
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(omicsData = omicsData,
                                         filter_object = filter_object_new)
  
  # return filtered data object #
  omicsData$e_data <- results_pieces$temp.edata
  omicsData$f_data <- results_pieces$temp.fdata
  omicsData$e_meta <- results_pieces$temp.emeta
  
  # if group attribute is present, re-run group_designation in case filtering
  # any items impacted the group structure #
  if (!is.null(attr(omicsData, "group_DF"))) {
    
    # Re-run group_designation in case filtering any items impacted the group
    # structure. The attributes will also be updated in this function.
    omicsData <- group_designation(omicsData = omicsData,
                                   main_effects = attr(get_group_DF(omicsData),
                                                       "main_effects"),
                                   covariates = attr(get_group_DF(omicsData),
                                                     "covariates"),
                                   time_course = attr(get_group_DF(omicsData),
                                                      "time_course"))
    
  } else {
    
    # Extract data_info attribute from omicsData. Some of the elements will be
    # used to update this attribute.
    dInfo <- attr(omicsData, 'data_info')
    
    # Update the data_info attribute.
    attr(omicsData, 'data_info') <- set_data_info(
      e_data = omicsData$e_data,
      edata_cname = get_edata_cname(omicsData),
      data_scale = get_data_scale(omicsData),
      data_types = dInfo$data_types,
      norm_info = dInfo$norm_info,
      is_normalized = dInfo$norm_info$is_normalized
    )
    
    # Update the meta_info attribute.
    attr(omicsData, 'meta_info') <- set_meta_info(
      e_meta = omicsData$e_meta,
      emeta_cname = get_emeta_cname(omicsData)
    )
    
  }
  
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData, 'filters')[[n_filters]] <- set_filter(
    type = class(filter_object)[[1]],
    threshold = NA,
    filtered = filter_object_new,
    method = NA
  )
  
  # Return the filtered data! Woot!!
  return(omicsData)
  
}

#' Remove items that need to be filtered out
#'
#' This function removes rows and columns in e_data, f_data, and e_meta based
#' on either remove or keep criteria.
#'
#' @param omicsData an object of the class \code{pepData}, \code{proData},
#' \code{lipidData}, or \code{metabData}, usually created by
#' \code{\link{as.pepData}}, \code{\link{as.proData}},
#' \code{\link{as.lipidData}}, or \code{\link{as.metabData}}, respectively.
#'    
#' @param filter_object A list of three elements. Each element contains a set of
#'        names to either remove or keep from e_data, f_data, and e_meta.
#' 
#' @return A list with three elements: first is the filtered e_data object,
#'         second is the filtered f_data object, and third is the filtered 
#'         e_meta object.
#' 
#' @author Kelly Stratton, Lisa Bramer
#'
pmartR_filter_worker <- function (filter_object, omicsData) {
  
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname
  
  # Filter with remove arguments ---------------
  
  if (!is.null(c(filter_object$e_data_remove,
                 filter_object$e_meta_remove,
                 filter_object$f_data_remove))) {
    
    # remove any samples from f_data and e_data #
    if (!is.null(filter_object$f_data_remove)) {
      
      # Find the row indices of the sample names that will be removed from the
      # f_data object.
      inds <- which(
        omicsData$f_data[, which(names(omicsData$f_data) == samp_cname)]
        %in% filter_object$f_data_remove
      )
      
      # Remove the rows in f_data that correspond to the samples that should be
      # removed.
      omicsData$f_data <- omicsData$f_data[-inds, ]
      
      # Find the column indices of the sample names that will be removed from
      # the e_data object.
      inds <- which(names(omicsData$e_data) %in% filter_object$f_data_remove)
      
      # Remove the columns in e_data that correspond to the samples that should
      # be removed.
      omicsData$e_data <- omicsData$e_data[, -inds]
      
      
    } 
    
    # remove any edata molecules from e_data and e_meta #
    if (!is.null(filter_object$e_data_remove)) {
      
      # Find the row indices of the biomolecule names that will be removed from
      # the e_data object.
      inds <- which(
        omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)]
        %in% filter_object$e_data_remove
      )
      
      # Remove the rows in e_data that correspond to the biomolecules that
      # should be removed.
      omicsData$e_data <- omicsData$e_data[-inds, ]
      
      # Check if e_meta is present.
      if (!is.null(omicsData$e_meta)) {
        
        # Find the row indices of the biomolecule names that will be removed
        # from the e_meta object.
        inds <- which(
          omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
          %in% filter_object$e_data_remove
        )
        
        # Remove the rows in e_meta corresponding to the biomolecules that
        # should be removed.
        omicsData$e_meta <- omicsData$e_meta[-inds, ]
        
      }
      
    }
    
    # remove any emeta molecules from e_meta and e_data #
    if (!is.null(filter_object$e_meta_remove)) {
      
      # Find the row indices of the mapping variable names that will be removed
      # from the e_meta object.
      inds <- which(
        omicsData$e_meta[, which(names(omicsData$e_meta) == emeta_cname)]
        %in% filter_object$e_meta_remove
      )
      
      # Check if there is at least one mapping variable that will be removed.
      if (length(inds) > 0) {
        
        # Remove the rows in e_meta corresponding the the mapping variables that
        # should be removed.
        omicsData$e_meta <- omicsData$e_meta[-inds, ]
        
      }
      
      # subset to the intersection of the edata_molecules in both e_data and 
      # e_meta, in case more were removed in one than the other #
      
      # Find the names of the biomolecules that are in both e_data and the
      # reduced e_meta object.
      mols <- intersect(
        omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)],
        omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
      )
      
      # Find the row indices in e_data of the biomolecules that will be kept in
      # both e_data and e_meta.
      inds <- which(
        omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)]
        %in% mols
      )
      
      # Only keep the rows in e_data corresponding to the biomolucules that
      # should not be removed.
      omicsData$e_data <- omicsData$e_data[inds, ]
      
      # Find the row indices in e_meta of the biomolecules that will be kept in
      # both e_data and e_meta.
      inds <- which(
        omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
        %in% mols
      )
      
      # Only keep the rows in e_meta corresponding to the biomolucules that
      # should not be removed.
      omicsData$e_meta <- omicsData$e_meta[inds, ]
      
    }
    
    # Filter with keep arguments ---------------
    
  } else {
    
    # keep samples in f_data and e_data #
    if (!is.null(filter_object$f_data_keep)) {
      
      # Find the row indices in f_data of the samples that will be kept.
      inds <- which(
        omicsData$f_data[, which(names(omicsData$f_data) == samp_cname)]
        %in% filter_object$f_data_keep
      )
      
      # Only keep the rows in f_data that correspond to the samples that
      # should be kept.
      omicsData$f_data <- omicsData$f_data[inds, ]
      
      # Find the column indices in e_data corresponding to the ID column and
      # the sample names that should be kept. The first element in the vector
      # is the column index where the ID column is.
      inds <- c(which(names(omicsData$e_data) == edata_cname), 
                which(names(omicsData$e_data) %in% filter_object$f_data_keep))
      
      # Only keep the columns in e_data that correspond to the samples that
      # should be kept.
      omicsData$e_data <- omicsData$e_data[ , inds]
      
    }
    
    # keep edata molecules in e_data and e_meta if it is present #
    if (!is.null(filter_object$e_data_keep)) {
      
      # Find the row indices in e_data that correspond to the biomolecules
      # that should be kept.
      inds <- which(
        omicsData$e_data[ , which(names(omicsData$e_data) == edata_cname)]
        %in% filter_object$e_data_keep
      )
      
      # Only keep the rows in e_data corresponding to the biomolecules that
      # should be kept.
      omicsData$e_data <- omicsData$e_data[inds, ]
      
      # if e_meta is present and we aren't explicitly specifying to keep
      # anything in it, also keep these e_data molecules in e_meta #
      if (!is.null(omicsData$e_meta) && is.null(filter_object$e_meta_keep)) {
        
        # Find the row indices of e_meta from the biomolecule IDs specified in
        # e_data_keep.
        inds <- which(
          omicsData$e_meta[ , which(names(omicsData$e_meta) == edata_cname)]
          %in% filter_object$e_data_keep
        )
        
        # Only keep the rows of e_meta that contain the biomolecule IDs from
        # the e_data_keep input.
        omicsData$e_meta <- omicsData$e_meta[inds, ]
        
      }
      
    }
    
    # keep emeta molecules in e_meta (here, we are explicitly specifying things
    # to keep).
    if (!is.null(filter_object$e_meta_keep)) {
      
      # Find the row indices in e_meta corresponding to the IDs of the mapping
      # variable. These are the rows that will be kept.
      inds <- which(
        omicsData$e_meta[ , which(names(omicsData$e_meta) == emeta_cname)]
        %in% filter_object$e_meta_keep
      )
      
      # Only keep the IDs of the mapping variable specified in e_meta_keep.
      omicsData$e_meta <- omicsData$e_meta[inds, ]
      
      # keep the union of the edata_molecules in both e_data and e_meta, in case 
      # more were kept in one than the other #
      if (is.null(filter_object$e_data_keep)) {
        
        # Use intersection here, since nothing was explicitly specified to
        # keep from e_data, and if we use the union, then e_data doesn't
        # actually get filtered at all.
        # Find the biomolecule names that occur in both e_data and the reduced
        # e_meta data frame.
        mols <- intersect(
          omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)],
          omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
        )
        
        # Find the row indices in e_data that correspond to the biomolecule
        # IDs that occur in both e_data and the reduced e_meta.
        inds <- which(
          omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)]
          %in% mols
        )
        
        # Only keep the rows in e_data that correspond to the biomolecule IDs
        # found in both e_data and the reduced e_meta data frames.
        omicsData$e_data <- omicsData$e_data[inds, ]
        
        # Find the row indices in e_meta that match the biomolecules found in
        # both e_data and the reduced e_meta data frames.
        inds <- which(
          omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
          %in% mols
        )
        
        # Only keep the rows in e_meta corresponding to the biomolecule IDs
        # from the previous line.
        omicsData$e_meta <- omicsData$e_meta[inds, ]
        
        # Keep arguments are specified for both e_data and e_meta.
      } else {
        
        # use union here, since there WERE things explicitly specified to keep
        # from e_data.
        # Find the names of the biomolecule IDs found in both the reduced
        # e_data and reduced e_meta data frames.
        mols <- union(
          omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)],
          omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
        )
        
        # Find the row indices in e_data corresponding to the biomolecule IDs
        # that occur in both e_data and e_meta. Both data frames were
        # previously filtered separately.
        inds <- which(
          omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)]
          %in% mols
        )
        
        # Only keep the rows in e_data that were extracted from the previous
        # line.
        omicsData$e_data <- omicsData$e_data[inds, ]
        
        # Find the row indices in e_meta corresponding to the biomolecule IDs
        # that occur in both e_data and e_meta. Both data frames were
        # previously filtered separately.
        inds <- which(
          omicsData$e_meta[, which(names(omicsData$e_meta) == edata_cname)]
          %in% mols
        )
        
        # Only keep the rows in e_meta that correspond to the indices
        # extracted in the previous line.
        omicsData$e_meta <- omicsData$e_meta[inds, ]
        
      }
      
    }
    
  }
  
  # Return the filtered omicsData pieces!!!
  return (list(temp.edata = omicsData$e_data,
               temp.fdata = omicsData$f_data,
               temp.emeta = omicsData$e_meta))
  
}
