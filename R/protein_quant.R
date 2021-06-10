#' protein_quant wrapper function
#'
#' This function takes in a pepData object, method (quantification method, mean,
#' median or rrollup), and the optional argument isoformRes (defaults to NULL).
#' An object of the class 'proData' is returned.
#'
#' @param pepData an omicsData object of the class 'pepData'.
#' @param method is one of four protein quantification methods, 'rollup',
#'   'rrollup', 'qrollup' and 'zrollup'. When 'rollup' is selected, combine_fn
#'   must also be provided and will determine whether pquant_mean or
#'   pquant_median function will be used.
#' @param isoformRes is a list of data frames, the result of applying the
#'   'bpquant' function to original pepData object. Defaults to NULL.
#' @param qrollup_thresh is a numeric value; is the peptide abundance cutoff
#'   value. Is an argument to qrollup function.
#' @param single_pep logical indicating whether or not to remove proteins that
#'   have just a single peptide mapping to them, defaults to FALSE.
#' @param single_observation logical indicating whether or not to remove
#'   peptides that have just a single observation, defaults to FALSE.
#' @param combine_fn can either be 'mean' or 'median'.
#' @param use_parallel logical indicating whether or not to use "doParallel"
#'   loop in applying rollup functions. Defaults to TRUE. Is an argument of
#'   rrollup, qrollup and zrollup functions.
#' @param emeta_cols A character vector indicating additional columns of e_meta
#'   that should be kept after rolling up to the protein level. The default,
#'   NULL, only keeps the column containing the mapping variable along with the
#'   new columns created (peps_per_pro and n_peps_used).
#' 
#' @return an omicsData object of the class 'proData'
#'
#' @details If isoformRes is provided then, a temporary pepData object is formed
#'   using the isoformRes information as the e_meta component and the original
#'   pepData object will be used for e_data and f_data components. The
#'   emeta_cname for the temporary pepData object will be the 'protein_isoform'
#'   column of isoformRes. Then one of the three 'method' functions can be
#'   applied to the temporary pepData object to return a proData object. If
#'   isofromRes is left NULL, then depending on the input for 'method', the
#'   correct 'method' function is applied directly to the input pepData object
#'   and a proData object is returned.
#'
#' @references Webb-Robertson, B.-J. M., Matzke, M. M., Datta, S., Payne, S. H.,
#'   Kang, J., Bramer, L. M., ... Waters, K. M. (2014). \emph{Bayesian
#'   Proteoform Modeling Improves Protein Quantification of Global Proteomic
#'   Measurements}. Molecular & Cellular Proteomics.: MCP, 13(12), 3639-3646.
#'
#' @examples
#' \dontrun{
#' library(pmartR)
#' library(pmartRdata)
#'
#' mypepData <- group_designation(omicsData = pep_object, main_effects = c("Condition"))
#' mypepData = edata_transform(mypepData, "log2")
#'
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
#'
#' imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='bon')
#'
#' isoformRes = bpquant(statRes = imd_anova_res, pepData = mypepData)
#'
#' #case where isoformRes is NULL:
#' results<- protein_quant(pepData = mypepData, method = 'rollup', combine_fn = 'median', isoformRes = NULL)
#'
#' #case where isoformRes is provided:
#' results2 = protein_quant(pepData = pep_object, method = 'rollup', combine_fn = 'mean', isoformRes = isoformRes)
#' }
#'
#' @rdname protein_quant
#' @export
#' 
protein_quant <- function (pepData, method, isoformRes = NULL,
                           qrollup_thresh = NULL, single_pep = FALSE,
                           single_observation = FALSE, combine_fn = "median",
                           use_parallel = TRUE, emeta_cols = NULL) {
  
  # Preflight checks -----------------------------------------------------------
  
  # Make sure the data are on one of the log scales.
  if (get_data_scale(pepData) == "abundance") {
    
    stop (paste("The data must be log transformed. Use edata_transform to",
                "convert to the log scale.",
                sep = " "))
    
  }
  
  if(!inherits(pepData, "pepData")) stop("pepData must be an object of class pepData")
  if(!(method %in% c('rollup', 'rrollup', 'qrollup', 'zrollup'))) stop("method must be one of, rollup, rrollup, qrollup, zrollup")
  if(!(combine_fn %in% c('median', 'mean'))) stop("combine_fn must be either 'mean' or 'median'")
  
  #gives message if single_pep and single_observation are TRUE and method is not zrollup
  if (method != 'zrollup' && (single_pep == TRUE || single_observation == TRUE)) message("single_pep and single_observation will be ignored, as they are only applicable if method is zrollup")
  
  #gives message if qrollup_thresh is not NULL and method is not qrollup
  if(method != 'qrollup' && !is.null(qrollup_thresh)) message("qrollup_thresh argument will be ignored, as it is only applicable if method is qrollup")
  
  # Set the combine_fn input to the appropriate function.
  if(combine_fn == "median"){
    chosen_combine_fn <- combine_fn_median
  }else{chosen_combine_fn <- combine_fn_mean}
  
  # Extract attribute info to be used throughout the function ------------------
  
  # Pull out column names from e_data, f_data, and e_meta.
  edata_cname <- attr(pepData, "cnames")$edata_cname
  fdata_cname <- attr(pepData, "cnames")$fdata_cname
  emeta_cname <- attr(pepData, "cnames")$emeta_cname
  
  # Extricate e_data column name index.
  edata_cname_id<- which(names(pepData$e_data) == edata_cname)
  
  e_meta<- pepData$e_meta
  
  # Prepare attribute info when isoformRes is present.
  if (!is.null(isoformRes)) {
    
    # The following attributes are reset when the as.pepData function is called
    # and will need to be manually updated after creating a new pepData object
    # with the isoformRes_subset attribute.
    filtas <- attr(pepData, "filters")
    groupies <- attr(pepData, "group_DF")
    scales <- attr(pepData, "data_scale_orig")
    inovas <- attr(pepData, "imdanova")
    
    # we will extract 'isoformRes_subset' attribute from isoformRes, which is
    # all the proteins that mapped to a nonzero proteoformID
    isoformRes2 <- attr(isoformRes, "isoformRes_subset")
    
    # Extract more attributes to create the new pepData object.
    data_scale <- get_data_scale(pepData)
    is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized
    
    # Find the peptides that occur in both e_data and isoformRes_subset.
    peptides <- which(
      pepData$e_data[, edata_cname] %in% isoformRes2[, edata_cname]
    )
    
    # Produce a pepData object with the reduced e_data object.
    pepData <- as.pepData(e_data = pepData$e_data[peptides, ],
                          f_data = pepData$f_data,
                          e_meta = isoformRes2,
                          edata_cname = edata_cname,
                          fdata_cname = fdata_cname,
                          emeta_cname = "Protein_Isoform",
                          data_scale = data_scale,
                          is_normalized = is_normalized)
    
    # Update the attributes that are reset in as.pepData.
    attr(pepData, "filters") <- filtas
    attr(pepData, "group_DF") <- groupies
    attr(pepData, "data_info")$data_scale_orig <- scales
    attr(pepData, "imdanova") <- inovas
    
  }
  
  # Quantitate the heck out of the peptides ------------------------------------
  
  if (method == 'rollup') {
    
    results <- pquant(pepData = pepData,
                      combine_fn = chosen_combine_fn)
    
  }
  
  if (method == 'rrollup') {
    results <- rrollup(pepData,
                       combine_fn = chosen_combine_fn,
                       parallel = use_parallel)
  }
  
  if (method == 'qrollup') {
    
    if (is.null(qrollup_thresh)) {
      
      stop ("qrollup_thresh parameter value must be specified") 
      
    }
    
    results <- qrollup(pepData,
                       qrollup_thresh = qrollup_thresh,
                       combine_fn = chosen_combine_fn,
                       parallel = use_parallel)
  }
  
  if (method == 'zrollup') {
    
    # Create both a proteomics and molecule filter object. These filter
    # objects will be used based on the input to single_pep and
    # single_observation arguments.
    proteomicsfilt <- proteomics_filter(pepData)
    moleculefilt <- molecule_filter(pepData)
    
    # Check for single pepes and single observations.
    if (single_pep == FALSE && single_observation == FALSE) {
      
      # Throw errors for unruly data.
      if(any(moleculefilt$Num_Observations == 1) && any(proteomicsfilt$counts_by_pro == 1)) stop("There are peptides with less than 2 observations and proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides/proteins and then run zrollup, set both 'single_observation' and 'single_pep' arguments to TRUE")
      if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
      if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE")
      
    }
    
    if (single_pep == TRUE && single_observation == FALSE) {
      
      # since single_pep is TRUE we will remove proteins with single peptide
      # mapped to them
      pepData = applyFilt(proteomicsfilt, pepData, min_num_peps = 2)
      
      if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
      
    }
    
    if(single_pep == FALSE && single_observation == TRUE){
      
      # since single_observation is TRUE we will remove peptides with single
      # observation
      pepData = applyFilt(moleculefilt, pepData)
      
      if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE.")
      
    }
    
    if (single_pep == TRUE && single_observation == TRUE) {
      
      pepData0 = applyFilt(moleculefilt, pepData, min_num = 2)
      pepData = applyFilt(proteomicsfilt, pepData0, min_num_peps = 2)
      
    }
    
    results <- zrollup(pepData,
                       combine_fn = chosen_combine_fn,
                       parallel = use_parallel)
  }
  
  # Update e_meta after quantitation -------------------------------------------
  
  if (is.null(isoformRes)) {
    
    # Update e_meta with peptide counts.
    results$e_meta <- results$e_meta %>%
      dplyr::group_by(!!rlang::sym(emeta_cname)) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      # Only qrollup will cause n_peps_used != peps_per_protein when isoformRes
      # is NULL.
      dplyr::mutate(
        n_peps_used = if ("n_peps_used" %in% colnames(results$e_meta)) {
          n_peps_used
        } else {
          peps_per_pro
        }
      ) %>%
      # Move n_pep_used to the end. This line will only make a change to the
      # order of the columns when qrollup is selected.
      dplyr::relocate(n_peps_used, .after = dplyr::last_col()) %>%
      dplyr::select(dplyr::any_of(c(emeta_cname, "peps_per_pro",
                                    "n_peps_used", emeta_cols))) %>%
      dplyr::distinct(dplyr::all_of(.)) %>%
      data.frame()
    
  } else {
    
    # if isoform is specified, proteins used in rollup will be unique peptides/protein found in isoformRes2
    peps_used <- isoformRes2 %>%
      dplyr::group_by(!!rlang::sym(emeta_cname)) %>%
      dplyr::mutate(n_peps_used = dplyr::n()) %>%
      dplyr::select(!!rlang::sym(emeta_cname), n_peps_used) %>%
      dplyr::slice(1)
    
    # store total number of peptides mapping to each protein (different from above)
    peps_per_protein = e_meta %>%
      dplyr::group_by(!!rlang::sym(emeta_cname)) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::select(!!rlang::sym(emeta_cname), peps_per_pro) %>%
      dplyr::slice(1)
    
    # join the two count columns to the output e_meta
    results$e_meta <- results$e_meta %>%
      dplyr::left_join(peps_per_protein, by = emeta_cname) %>%
      dplyr::left_join(peps_used, by = emeta_cname)
    
  }
  
  # Update proData attributes --------------------------------------------------
  
  # Update the original data scale for the proData object. This needs to be
  # manually updated because the original data scale for the proData object will
  # be set to the current data scale of the pepData object when the as.proData
  # function was called. If these two scales are different this is the only way
  # to set the original data scale for the proData object to the original data
  # scale of the pepData object.
  attr(results, "data_info")$data_scale_orig <- get_data_scale_orig(pepData)
  
  # Update the group_DF attribute (if it exists). This attribute will be "reset"
  # when the as.proData function is called in the rollup functions. It will need
  # to be manually updated to reflect anything done to the peptide data before
  # protein_quant.
  attr(results, "group_DF") <- attr(pepData, "group_DF")
  
  return(results)
  
}
