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
#'   \code{\link{as.lipidData}}, \code{\link{as.metabData}}, or \code{\link{as.nmrData}}, respectively.
#' @param ... further arguments
#'
#' @return An object of the class \code{pepData}, \code{proData},
#'   \code{lipidData}, \code{metabData}, or \code{nmrData} with specified cname_ids,
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
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' to_filter <- molecule_filter(omicsData = pep_object)
#' pep_object2 <- applyFilt(filter_object = to_filter, omicsData = pep_object, min_num = 2)
#' print(str(attributes(pep_object2)$filters))
#' pep_object2 <- group_designation(pep_object2, main_effects = "Condition")
#' to_filter2 <- imdanova_filter(omicsData = pep_object2)
#' pep_object3 <- applyFilt(filter_object = to_filter2, omicsData = pep_object2, min_nonmiss_anova = 3)
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

  # check that omicsData is of pmartR S3 class#
  if(!inherits(omicsData, c("pepData", "proData", "lipidData", "metabData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'lipidData', 'metabData' or 'nmrData'")

  # check that filter_object is of an appropriate class#
  if(!inherits(filter_object, c("cvFilt", "proteomicsFilt", "moleculeFilt", "rmdFilt", "imdanovaFilt", "customFilt"))) stop("filter_object must be of  'cvFilt', 'proteomicsFilt', 'moleculeFilt', 'rmdFilt', 'imdanovaFilt', or 'customFilt.")

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
  filter_object_new <- list(edata_filt = filter.edata,
                            emeta_filt = NULL,
                            samples_filt = NULL)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.pep2
  omicsData$f_data <- results_pieces$temp.samp2
  omicsData$e_meta <- results_pieces$temp.meta1
  
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
    attr(omicsData,
         'data_info') <- set_data_info(e_data = omicsData$e_data,
                                       edata_cname = get_edata_cname(omicsData),
                                       data_scale = get_data_scale(omicsData),
                                       data_types = dInfo$data_types,
                                       norm_info = dInfo$norm_info,
                                       is_normalized = dInfo$norm_info$is_normalized)
    
    # Update the meta_info attribute.
    attr(omicsData,
         'meta_info') <- set_meta_info(e_meta = omicsData$e_meta,
                                       emeta_cname = get_emeta_cname(omicsData))
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData,
       'filters')[[n_filters]] <- set_filter(type = class(filter_object)[[1]],
                                             threshold = min_num,
                                             filtered = filter.edata,
                                             method = NA)

  # RETURN THE FILTERED OMICSDATA OBJECT! YAY!!!
  return(omicsData)
  
}

# function for cvFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.cvFilt <- function(filter_object, omicsData, cv_threshold = 150){
  
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
  filter_object_new <- list(edata_filt = filter.edata,
                            emeta_filt = NULL,
                            samples_filt = NULL)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.pep2
  omicsData$f_data <- results_pieces$temp.samp2
  omicsData$e_meta <- results_pieces$temp.meta1
  
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
    attr(omicsData,
         'data_info') <- set_data_info(e_data = omicsData$e_data,
                                       edata_cname = get_edata_cname(omicsData),
                                       data_scale = get_data_scale(omicsData),
                                       data_types = dInfo$data_types,
                                       norm_info = dInfo$norm_info,
                                       is_normalized = dInfo$norm_info$is_normalized)
    
    # Update the meta_info attribute.
    attr(omicsData,
         'meta_info') <- set_meta_info(e_meta = omicsData$e_meta,
                                       emeta_cname = get_emeta_cname(omicsData))
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData,
       'filters')[[n_filters]] <- set_filter(type = class(filter_object)[[1]],
                                             threshold = cv_threshold,
                                             filtered = filter.edata,
                                             method = NA)

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
  filter_object_new <- list(edata_filt = NULL,
                            emeta_filt = NULL,
                            samples_filt = filter.samp)
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.pep2
  omicsData$f_data <- results_pieces$temp.samp2
  omicsData$e_meta <- results_pieces$temp.meta1
  
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
    attr(omicsData,
         'data_info') <- set_data_info(e_data = omicsData$e_data,
                                       edata_cname = get_edata_cname(omicsData),
                                       data_scale = get_data_scale(omicsData),
                                       data_types = dInfo$data_types,
                                       norm_info = dInfo$norm_info,
                                       is_normalized = dInfo$norm_info$is_normalized)
    
    # Update the meta_info attribute.
    attr(omicsData,
         'meta_info') <- set_meta_info(e_meta = omicsData$e_meta,
                                       emeta_cname = get_emeta_cname(omicsData))
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData,
       'filters')[[n_filters]] <- set_filter(type = class(filter_object)[[1]],
                                             threshold = pvalue_threshold,
                                             filtered = filter.samp,
                                             method = NA)

  # Return the filtered data! High fives all around!!!
  return(omicsData)
  
}

# function for proteomicsFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
#' @export
applyFilt.proteomicsFilt <- function(filter_object,
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
    degen_peptides = as.character(data.frame(count_bypep[which(count_bypep$n > 1),])[, pep_id])
    
    ## identify any proteins that now will not have peptides mapping to them ##
    ## find rows in e_meta that correspond to peptides to be filtered ##
    pepfilt.ids = which(omicsData$e_meta[, pep_id] %in% degen_peptides)
    
    # Fish out the indices for the proteins that will be filtered.
    add_prots <- as.character(omicsData$e_meta[pepfilt.ids, pro_id])
    # add_prots = as.character(setdiff(omicsData$e_meta[pepfilt.ids, pro_id],
    #                                  omicsData$e_meta[-pepfilt.ids, pro_id]))
    
    if(length(add_prots)==0){add_prots = NULL}
    if(length(degen_peptides)==0){degen_peptides = NULL}
    
    pepe <- list(edata_filt = degen_peptides, emeta_filt = add_prots)
    
  } else {
    pepe <- list(edata_filt = NULL, emeta_filt = NULL)
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
    pep_filt <- as.character(omicsData$e_meta[protfilt.ids, pep_id])
    # pep_filt = as.character(setdiff(omicsData$e_meta[protfilt.ids, pep_id],
    #                                 omicsData$e_meta[-protfilt.ids, pep_id]))
    
    if(length(pep_filt)==0){pep_filt = NULL}
    if(length(pro_filt)==0){pro_filt = NULL}
    
    pepe2 <- list(emeta_filt = pro_filt, edata_filt = pep_filt)
    
  } else {
    pepe2 <- list(edata_filt = NULL, emeta_filt = NULL)
  }
  
  # Consolidate pepe and pepe2 to pass to the pmartR_filter_worker function.
  filter_object_new <- list(emeta_filt = unique(c(pepe$emeta_filt, 
                                                  pepe2$emeta_filt)),
                            edata_filt = unique(c(pepe$edata_filt, 
                                                  pepe2$edata_filt)))
  
  # checking that filter_object_new does not specify all of e_data(peps) or all
  # of e_meta(protiens) in omicsData
  if(all(omicsData$e_meta[[pro_id]] %in% filter_object_new$emeta_filt)) {
    
    stop("filter_object specifies all proteins in e_meta")
    
  }
  if(all(omicsData$e_data[[pep_id]] %in% filter_object_new$edata_filt)) {
    
    stop("filter_object specifies all peps in e_data")
    
  }
  
  # Filter the data and update the attributes ----------------------------------
  
  # call the function that does the filter application
  results_pieces <- pmartR_filter_worker(filter_object = filter_object_new,
                                         omicsData = omicsData)
  
  # Update the omicsData data frames.
  omicsData$e_data <- results_pieces$temp.pep2
  omicsData$f_data <- results_pieces$temp.samp2
  omicsData$e_meta <- results_pieces$temp.meta1
  
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
    attr(omicsData,
         'data_info') <- set_data_info(e_data = omicsData$e_data,
                                       edata_cname = get_edata_cname(omicsData),
                                       data_scale = get_data_scale(omicsData),
                                       data_types = dInfo$data_types,
                                       norm_info = dInfo$norm_info,
                                       is_normalized = dInfo$norm_info$is_normalized)
    
    # Update the meta_info attribute.
    attr(omicsData,
         'meta_info') <- set_meta_info(e_meta = omicsData$e_meta,
                                       emeta_cname = get_emeta_cname(omicsData))
    
  }
  
  # Determine the number of filters applied previously and add one to it. This
  # will include the current filter object in the next available space of the
  # filters attribute list.
  n_filters <- length(attr(omicsData, 'filters')) + 1
  
  # Update the filters attribute.
  attr(omicsData,
       'filters')[[n_filters]] <- set_filter(type = class(filter_object)[[1]],
                                             threshold = data.frame(min_num_peps = min_num_peps,
                                                                    degen_peps = as.character(degen_peps)),
                                             filtered = filter_object_new,
                                             method = NA)
  
  # Return the filtered data! Look at us go!!
  return(omicsData)


}

# function for imdanovaFilt
#' @export
#' @name applyFilt
#' @rdname applyFilt
applyFilt.imdanovaFilt <- function(filter_object, omicsData, min_nonmiss_anova=NULL, min_nonmiss_gtest=NULL, remove_singleton_groups = TRUE){
  # #' @details If filter_method="combined" is specified, then both the \code{anova_filter} and \code{gtest_filter} are applied to the data, and the intersection of features from the two filters is the set returned. For ANOVA, features that do not have at least \code{min_nonmiss_allowed} values per group are candidates for filtering. For the G-test, features that do not have at least \code{min_nonmiss_allowed} values per group are candidates for filtering. The G-test is a test of independence, used here to test the null hypothesis of independence between the number of missing values across groups.


  # check to see whether a imdanovaFilt has already been run on omicsData #
  if("imdanovaFilt" %in% names(attributes(omicsData)$filters)){
    # get previous threshold #
    min_nonmiss_anova_prev <- attributes(omicsData)$filters$imdanovaFilt$threshold$min_nonmiss_anova
    min_nonmiss_gtest_prev <- attributes(omicsData)$filters$imdanovaFilt$threshold$min_nonmiss_gtest

    stop(paste("An IMD-ANOVA filter has already been run on this dataset, using a 'min_nonmiss_anova' of ", min_nonmiss_anova_prev, " and 'min_nonmiss_gtest' of ", min_nonmiss_gtest_prev, ". See Details for more information about how to choose a threshold before applying the filter.", sep=""))


  }else{ # no previous imdanovaFilt, so go ahead and run it like normal #


    ## initial checks ##

    # verify remove_singleton_groups is logical T/F #
    if(class(remove_singleton_groups) != "logical"){stop("remove_singletong_groups must be of class logical")}
    
    group_sizes <- attr(filter_object, "group_sizes")
    nonmiss_per_group <- list(nonmiss_totals = filter_object, group_sizes = group_sizes)
    groupDF <- attributes(omicsData)$group_DF
    
    # if remove_singleton_groups is TRUE, filter those out now #
    if(remove_singleton_groups == TRUE){
      if(any(group_sizes$n_group == 1)){
        # which group(s) #
        singleton_groups <- group_sizes$Group[group_sizes$n_group == 1]
        
        # which sample name(s) #
        samps_to_rm <- as.character(groupDF[which(groupDF$Group %in% singleton_groups), get_fdata_cname(omicsData)])
        
        # use custom_filter to remove the sample(s) from the omicsData object #
        my_cust_filt <- custom_filter(omicsData, f_data_remove = samps_to_rm)
        omicsData <- applyFilt(my_cust_filt, omicsData)
        
        # update the IMD-ANOVA filter attributes #
        attributes(filter_object)$names <- as.character(attributes(filter_object)$names[-which(attributes(filter_object)$names %in% singleton_groups)])
        # might need to explicitly remove any NAs... #
        # attributes(filter_object)$names <- attributes(filter_object)$names[-which(is.na(attributes(filter_object)$names))]
        
        attributes(filter_object)$group_sizes <- attributes(filter_object)$group_sizes[-which(attributes(filter_object)$group_sizes$Group %in% singleton_groups), ]
        
        # give a message to user saying which samples have been removed #
        message(paste("You have specified remove_single_groups = TRUE, so we have removed the following sample(s) that correspond to singleton groups prior to proceeding with the IMD-ANOVA filter:", 
                      paste(samps_to_rm, sep = ", ")))
        
      }else{
        message("You have specified remove_singleton_groups = TRUE, but there are no singleton groups to remove. 
                Proceeding with application of the IMD-ANOVA filter.")
      }
    }else{
      if(any(group_sizes$n_group == 1)){
        message("You have specified remove_singleton_groups = FALSE, so the sample(s) corresponding to the singleton group(s) were not utilized in the IMD-ANOVA filter and will be retained in the resulting omicsData object.")
      }
    }
    
    # check that at least one of min_nonmiss_anova and min_nonmiss_gtest are present #
    if(is.null(min_nonmiss_anova) & is.null(min_nonmiss_gtest)) stop("At least one of min_nonmiss_anova and min_nonmiss_gtest must be present")
    # check that if they aren't NULL, min_nonmiss_anova and min_nonmiss_gtest are numeric, >=2 and >=3, respectively, and neither are bigger than the minimum group size (group_sizes in an attribute of the filter_object, see below) #
    if(!is.null(min_nonmiss_anova)) {
      # check that min_nonmiss_anova is numeric >= 2 #
      if(!inherits(min_nonmiss_anova, c("numeric","integer")) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
      # check that min_nonmiss_anova is an integer #
      if(min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
      # check that min_nonmiss_gtest is less than the minimum group size #
      nonsingleton_groups <- attributes(attributes(omicsData)$group_DF)$nonsingleton_groups
      if(min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group[which(attributes(filter_object)$group_sizes$Group %in% nonsingleton_groups)])) stop("min_nonmiss_anova cannot be greater than the minimum group size")
    }
    if(!is.null(min_nonmiss_gtest)) {
      # check that min_nonmiss_gtest is numeric >= 3 #
      if(!inherits(min_nonmiss_gtest, c("numeric","integer")) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
      # check that min_nonmiss_gtest is an integer #
      if(min_nonmiss_gtest %% 1 != 0) stop("min_nonmiss_gtest must be an integer >= 3")
      # check that min_nonmiss_gtest is less than the minimum group size #
      nonsingleton_groups <- attributes(attributes(omicsData)$group_DF)$nonsingleton_groups
      if(min_nonmiss_gtest > min(attributes(filter_object)$group_sizes$n_group[which(attributes(filter_object)$group_sizes$Group %in% nonsingleton_groups)])) stop("min_nonmiss_gtest cannot be greater than the minimum group size")
    }
    
    ## end of initial checks ##

    
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
    results_pieces <- pmartR_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

    # return filtered data object #
    omicsData$e_data <- results_pieces$temp.pep2
    omicsData$f_data <- results_pieces$temp.samp2
    omicsData$e_meta <- results_pieces$temp.meta1
    results <- omicsData

    # if group attribute is present (and it MUST be), re-run group_designation in case filtering any items impacted the group structure
    if(!is.null(attr(results, "group_DF"))){
      results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
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

  # if filter_object contains 'removes' #
  if(!is.null(filter_object$e_data_remove)|!is.null(filter_object$f_data_remove)|!is.null(filter_object$e_meta_remove)){
    
    filter_object_new = list(edata_filt = filter_object$e_data_remove, emeta_filt = filter_object$e_meta_remove, samples_filt = filter_object$f_data_remove)
    
    # check that edata_filt doesn't specify ALL the items in omicsData #
    if(all(omicsData$e_data[, edata_cname] %in% filter_object_new$edata_filt)){stop("e_data_remove specifies all the items in the data")}
    
    # check that samples_filt doesn't specify ALL the items in omicsData #
    if(all(omicsData$f_data[, fdata_cname] %in% filter_object_new$samples_filt)){stop("f_data_remove specifies all the items in the data")}
    
    # check that emeta_filt doesn't specify ALL the items in omicsData, emeta_filt is present #
    if(!is.null(omicsData$e_meta)){
       if(all(omicsData$e_meta[, emeta_cname] %in% filter_object_new$emeta_filt)){stop("e_meta_remove specifies all the items in the data")}
    }
    
  }
  
  else{ # filter_object contains 'keeps' #
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
  results_pieces <- pmartR_filter_worker(omicsData = omicsData, filter_object = filter_object_new)

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
    attributes(results)$data_info$prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
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
pmartR_filter_worker <- function(filter_object, omicsData){
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname
  
  # pull group_DF attribute #
  group_DF = attr(omicsData, "group_DF")
  
  # initialize the new omicsData parts #
  temp.edata <- omicsData$e_data
  temp.fdata <- omicsData$f_data
  temp.emeta <- omicsData$e_meta
  
  #check if filter object contains remove arguments
  if(!is.null(filter_object$edata_filt)|!is.null(filter_object$emeta_filt)|!is.null(filter_object$samples_filt)){
    
    # remove any samples from f_data and e_data #
    if(!is.null(filter_object$samples_filt)){
      inds <- which(temp.fdata[, which(names(temp.fdata) == samp_cname)] %in% filter_object$samples_filt)
      temp.fdata <- temp.fdata[-inds, ]
      
      inds <- which(names(temp.edata) %in% filter_object$samples_filt)
      temp.edata <- temp.edata[ ,-inds]
    }
    
    # remove any edata molecules from e_data and e_meta #
    if(!is.null(filter_object$edata_filt)){
      inds <- which(temp.edata[ , which(names(temp.edata) == edata_cname)] %in% filter_object$edata_filt)
      temp.edata <- temp.edata[-inds, ]
      
      # also remove these from e_meta, if it is present #
      if(!is.null(temp.emeta)){
        inds <- which(temp.emeta[ , which(names(temp.emeta) == edata_cname)] %in% filter_object$edata_filt)
        temp.emeta <- temp.emeta[-inds, ]
      }
    }
    
    # remove any emeta molecules from e_meta and e_data #
    if(!is.null(filter_object$emeta_filt)){
      inds <- which(temp.emeta[ , which(names(temp.emeta) == emeta_cname)] %in% filter_object$emeta_filt)
      if(length(inds) > 0){
        temp.emeta <- temp.emeta[-inds, ] 
      }
      
      # subset to the intersection of the edata_molecules in both e_data and 
      # e_meta, in case more were removed in one than the other #
      mols <- intersect(temp.edata[, which(names(temp.edata) == edata_cname)],
                        temp.emeta[, which(names(temp.emeta) == edata_cname)])
      inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
      temp.edata <- temp.edata[inds, ]
      inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
      temp.emeta <- temp.emeta[inds, ]
    }
    
  }else{ # filter object contains keep arguments #
    
    # keep samples in f_data and e_data #
    if(!is.null(filter_object$samples_keep)){
      inds <- which(temp.fdata[, which(names(temp.fdata) == samp_cname)] %in% filter_object$samples_keep)
      temp.fdata <- temp.fdata[inds, ]
      
      inds <- c(which(names(temp.edata) == edata_cname), 
                which(names(temp.edata) %in% filter_object$samples_keep))
      temp.edata <- temp.edata[ , inds]
    }
    
    # keep edata molecules in e_data #
    if(!is.null(filter_object$edata_keep)){
      inds <- which(temp.edata[ , which(names(temp.edata) == edata_cname)] %in% filter_object$edata_keep)
      temp.edata <- temp.edata[inds, ]
      
      # if e_meta is present and we aren't explicitly specifying to keep anything
      # in it, also keep these e_data molecules in e_meta #
      if(!is.null(temp.emeta) & is.null(filter_object$emeta_keep)){
        inds <- which(temp.emeta[ , which(names(temp.emeta) == edata_cname)] %in% filter_object$edata_keep)
        temp.emeta <- temp.emeta[inds, ]
      }
    }
    
    # keep emeta molecules in e_meta (here, we are explicitly specifying things to keep) #
    if(!is.null(filter_object$emeta_keep)){
      inds <- which(temp.emeta[ , which(names(temp.emeta) == emeta_cname)] %in% filter_object$emeta_keep)
      temp.emeta <- temp.emeta[inds, ]
      
      # keep the union of the edata_molecules in both e_data and e_meta, in case 
      # more were kept in one than the other #
      if(is.null(filter_object$edata_keep)){
        # use intersection here, since nothing was explicitly specified to keep
        # from edata, and if we use the union, then edata doesn't actually get filtered at all #
        mols <- intersect(temp.edata[, which(names(temp.edata) == edata_cname)],
                          temp.emeta[, which(names(temp.emeta) == edata_cname)])
        inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
        temp.edata <- temp.edata[inds, ]
        inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
        temp.emeta <- temp.emeta[inds, ]
      }else{
        # use union here, since there WERE things explicitly specified to keep from edata #
        mols <- union(temp.edata[, which(names(temp.edata) == edata_cname)],
                      temp.emeta[, which(names(temp.emeta) == edata_cname)])
        inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
        temp.edata <- temp.edata[inds, ]
        inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
        temp.emeta <- temp.emeta[inds, ]
      }
    }
  }
  
  # return the pieces needed to assemble a proData/pepData/lipidData/metabData object #
  output <- list(temp.pep2 = temp.edata,
                 temp.samp2 = temp.fdata,
                 temp.meta1 = temp.emeta,
                 edata_cname = edata_cname,
                 emeta_cname = emeta_cname,
                 samp_cname = samp_cname)
  
  return(output)
}