#' Protein Quantification
#'
#' This function takes in a pepData object, method (quantification method, mean,
#' median or rrollup), and the optional argument isoformRes (defaults to NULL).
#' An object of the class 'proData' is returned.
#'
#' @param pepData an omicsData object of the class 'pepData'
#' @param method character string specifying one of four protein quantification methods, 'rollup',
#'   'rrollup', 'qrollup' and 'zrollup'
#' @param isoformRes list of data frames, the result of applying the
#'   'bpquant' function to original pepData object. Defaults to NULL.
#' @param qrollup_thresh numeric value; is the peptide abundance cutoff
#'   value. Is an argument to qrollup function.
#' @param single_pep logical indicating whether or not to remove proteins that
#'   have just a single peptide mapping to them, defaults to FALSE.
#' @param single_observation logical indicating whether or not to remove
#'   peptides that have just a single observation, defaults to FALSE.
#' @param combine_fn character string specifying either be 'mean' or 'median'
#' @param use_parallel logical indicating whether or not to use "doParallel"
#'   loop in applying rollup functions. Defaults to TRUE. Is an argument of
#'   rrollup, qrollup and zrollup functions.
#' @param emeta_cols character vector indicating additional columns of e_meta
#'   that should be kept after rolling up to the protein level. The default,
#'   NULL, only keeps the column containing the mapping variable along with the
#'   new columns created (peps_per_pro and n_peps_used).
#'
#' @return omicsData object of the class 'proData'
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
#' library(pmartRdata)
#'
#' mypepData <- group_designation(omicsData = pep_object, main_effects = c("Phenotype"))
#' mypepData = edata_transform(omicsData = mypepData, "log2")
#'
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
#'
#' imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust_a ='bon', pval_adjust_g ='bon')
#'
#' isoformRes = bpquant(statRes = imd_anova_res, pepData = mypepData)
#'
#' #case where isoformRes is NULL:
#' results<- protein_quant(pepData = mypepData, method = 'rollup', combine_fn = 'median', isoformRes = NULL)
#'
#' #case where isoformRes is provided:
#' # results2 = protein_quant(pepData = mypepData, method = 'rollup', combine_fn = 'mean', isoformRes = isoformRes)
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

  # Check if isoformRes is actually and isoformRes.
  if (!is.null(isoformRes) && class(isoformRes) != "isoformRes") {
    stop ("The input for isoformRes must be of class 'isoformRes'.")
  }

  # Why must you be so illogical?! WWHHHHYYYYYY??
  if (!is.logical(single_pep)) stop ("sinlge_pep must be either TRUE or FALSE.")
  if (!is.logical(single_observation)) {
    stop ("sinlge_observation must be either TRUE or FALSE.")
  }

  # Make sure emeta_cols is a character vector.
  if (!is.null(emeta_cols) && !is.character(emeta_cols)) {

    stop ("emeta_cols must be a character vector.")

  }

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

  # Nab the check.names attribute. This needs to be used anytime another
  # omicsData object is created. There are cases where the names will not match
  # between the old and new omicsData objects which leads to errors.
  fijate <- get_check_names(pepData)

  # Grab more attributes that will be used at some point somewhere.
  data_scale <- get_data_scale(pepData)
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Prepare attribute info when isoformRes is present.
  if (!is.null(isoformRes)) {

    # Keep a copy of the original e_meta data frame. This will be used to
    # compute the number of peptides used per protein at the end of the
    # function.
    e_meta <- pepData$e_meta

    # The following attributes are reset when the as.pepData function is called
    # and will need to be manually updated after creating a new pepData object
    # with the isoformRes_subset attribute.
    filtas <- attr(pepData, "filters")
    groupies <- attr(pepData, "group_DF")
    scales <- get_data_scale_orig(pepData)
    inovas <- attr(pepData, "imdanova")

    # we will extract 'isoformRes_subset' attribute from isoformRes, which is
    # all the proteins that mapped to a nonzero proteoformID
    isoformRes2 <- attr(isoformRes, "isoformRes_subset")

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
                          is_normalized = is_normalized,
                          check.names = fijate)

    # Update the attributes that are reset in as.pepData.
    attr(pepData, "filters") <- filtas
    attr(pepData, "group_DF") <- groupies
    attr(pepData, "data_info")$data_scale_orig <- scales
    attr(pepData, "imdanova") <- inovas

    # Grab the new e_meta because it has the columns we need to create the
    # proData object at the end of this function. Only keep the protein ID and
    # Protein_Isoform columns. More will be added later if they pass the new
    # pre_flight check that is causing all sorts of trouble.
    e_meta_iso <- pepData$e_meta %>%
      dplyr::select(
        dplyr::any_of(c(emeta_cname, "Protein_Isoform"))
      )

    # Update emeta_cname because it is different when using the isoform crap.
    emeta_cname_iso <- "Protein_Isoform"

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

    # Update pepData -----------------------------------------------------------

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

  # Check if emeta_cols is NULL. If it is everything can proceed as usual. If it
  # is not NULL we have to check if the number of unique rows in e_meta is the
  # same as the number of rows in e_data. If they aren't the pre_flight function
  # will throw an error. To avoid this we will set emeta_cols to NULL. This will
  # only keep the protein ID, n_peps_used, and peps_per_pro columns.
  if (!is.null(emeta_cols)) {

    # Nab number of rows in e_data to compare with number of rows in e_meta.
    n_row_edata <- nrow(results$e_data)

    # Grab the number of rows in e_meta depending on the rollup method used.
    n_row_emeta <- if (is.null(results$e_meta)) {
      # Use either the original or isoform e_meta depending on the input.
      if (is.null(isoformRes)) {
        nrow(unique(pepData$e_meta))
      } else {
        nrow(e_meta_iso)
      }
    } else {
      nrow(unique(results$e_meta))
    }

    # Change emeta_cols to NULL if the number of e_data and unique e_meta rows
    # do not match.
    if (n_row_edata != n_row_emeta) {

      emeta_cols <- NULL

    }

  }

  # The rollup functions now return a list containing e_data and e_meta. The
  # e_meta element of this list is NULL except for the case when qrollup is
  # used. Check if e_meta is NULL and do stuff accordingly.
  if (is.null(results$e_meta)) {

    # The isoformRes portion of the code changes the pepData object. This causes
    # trouble when isoformRes is present. If isoformRes is not null use the
    # emeta object created previously instead of pepData$e_meta.
    if (!is.null(isoformRes)) {

      results$e_meta <- e_meta_iso

    } else {

      results$e_meta <- pepData$e_meta

    }

  }

  # Check if isoformRes is NULL. results$e_meta will be updated differently
  # depending on whether isoformRes is present.
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
      # Keep the mapping variable and new columns created plus any columns
      # specified by the user.
      dplyr::select(dplyr::any_of(c(emeta_cname, "peps_per_pro",
                                    "n_peps_used", emeta_cols))) %>%
      # Only keep distinct combinations of the columns that are kept.
      dplyr::distinct(dplyr::all_of(.)) %>%
      data.frame()

    # The following runs when isoformRes is present. In this case n_peps_used
    # will be calculated based on protein isoform instead of protein (which
    # includes all isoforms).
  } else {

    # Count the number of peptides per isoform. If qrollup was used this was
    # already done in the qrollup function otherwise the peptides per isoform
    # will be counted from the isoformRes_subset data frame.
    peps_used <- if ("n_peps_used" %in% colnames(results$e_meta)) {

      # Extract counts from e_meta because peptides were counted in qrollup.
      results$e_meta %>%
        dplyr::select(!!rlang::sym(emeta_cname), n_peps_used) %>%
        # Only keep unique combinations of emeta_cname and n_peps_used
        dplyr::distinct(dplyr::all_of(.))

    } else {

      # Count peptides per isoform from the bpquant output.
      isoformRes2 %>%
        dplyr::group_by(Protein_Isoform) %>%
        dplyr::mutate(n_peps_used = dplyr::n()) %>%
        # Only keep the first row of each group.
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!rlang::sym(emeta_cname), n_peps_used)

    }

    # store total number of peptides mapping to each protein (different from
    # above)
    peps_per_protein = e_meta %>%
      dplyr::group_by(!!rlang::sym(emeta_cname)) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      # Keep the mapping variable, pep_per_pro, and any columns specified by the
      # user.
      dplyr::select(dplyr::any_of(c(emeta_cname, "peps_per_pro", emeta_cols)))

    # Check if n_peps_used is a column in e_meta. This will only be the case if
    # qrollup was used.
    if ("n_peps_used" %in% colnames(results$e_meta)) {

      # Remove the n_peps_used column.
      results$e_meta <- results$e_meta %>%
        dplyr::select(-n_peps_used)

    }

    # join the two count columns to the output e_meta
    results$e_meta <- results$e_meta %>%
      dplyr::left_join(peps_per_protein, by = emeta_cname) %>%
      dplyr::left_join(peps_used, by = emeta_cname) %>%
      # Move any columns specified by the user after n_peps_used.
      dplyr::relocate(dplyr::any_of(emeta_cols), .after = n_peps_used) %>%
      # Only keep distinct combinations of the columns that are kept.
      dplyr::distinct(dplyr::all_of(.)) %>%
      data.frame()

  }

  # Create a proData object ----------------------------------------------------

  # The isoform crap is making this really difficult. When creating the
  # data frames above emeta_cname needs to be the protein ID column. However,
  # when creating the proData object the column name changes to
  # "Protein_isoform".
  if (!is.null(isoformRes)) {

    emeta_cname <- emeta_cname_iso

  }

  prodata <- as.proData(e_data = results$e_data,
                        f_data = pepData$f_data,
                        e_meta = results$e_meta,
                        edata_cname = emeta_cname,
                        fdata_cname = fdata_cname,
                        emeta_cname = emeta_cname,
                        data_scale = get_data_scale(pepData),
                        is_normalized = is_normalized,
                        check.names = fijate)

  # Update proData attributes --------------------------------------------------

  # Update the original data scale for the proData object. This needs to be
  # manually updated because the original data scale for the proData object will
  # be set to the current data scale of the pepData object when the as.proData
  # function was called. If these two scales are different this is the only way
  # to set the original data scale for the proData object to the original data
  # scale of the pepData object.
  attr(prodata, "data_info")$data_scale_orig <- get_data_scale_orig(pepData)

  # Update the group_DF attribute (if it exists). This attribute will be "reset"
  # when the as.proData function is called in the rollup functions. It will need
  # to be manually updated to reflect anything done to the peptide data before
  # protein_quant.
  attr(prodata, "group_DF") <- attr(pepData, "group_DF")

  # Update the pro_quant_info attribute to reflect which rollup method was used.
  attr(prodata, "pro_quant_info")$method <- method

  return (prodata)

}

#' Protein Quantitation using Mean or Median Peptide Abundances
#'
#' This function takes in a pepData object and returns a proData object
#'
#' @param pepData omicsData object of class 'pepData'
#'
#' @param combine_fn A character string that can either be 'mean' or 'median'.
#'
#' @return An omicsData object of class 'proData'
#'
pquant <- function (pepData,
                    combine_fn) {

  # Nab the check names attribute because data.frame conventions are making my
  # life miserable. In fact, if it weren't for check.names and COVID my life
  # would be pretty awesome right now.
  check_names <- get_check_names(pepData)

  # check that pepData is of appropraite class #
  if (!inherits(pepData, "pepData")) {

    stop ("pepData is not an object of the appropriate class")

  }

  # check that a protein mapping is provided #
  if (is.null(pepData$e_meta)) {

    stop (paste("A mapping to proteins must be provided in order to use the",
                "protein_filter function.",
                sep = " "))

  }

  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname

  # Quantitate the heck out of the peptides!
  res <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
               y = pepData$e_data,
               by = pep_id,
               all.x = FALSE,
               all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    dplyr::group_by(!!rlang::sym(pro_id)) %>%
    dplyr::mutate(dplyr::across(.cols = -dplyr::any_of(pro_id),
                                .fns = combine_fn)) %>%
    dplyr::distinct() %>%
    data.frame(check.names = check_names)


  # Extricate attribute info for creating the proData object.
  # data_scale <- attr(pepData, "data_info")$data_scale
  # is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Add check here from Lisa's email -------------------------------------------

  # Create a proData object with the quantitated proteins.
  # prodata <- as.proData(e_data = res,
  #                       f_data = pepData$f_data,
  #                       e_meta = dplyr::select(pepData$e_meta,
  #                                              -rlang::sym(pep_id)),
  #                       edata_cname = pro_id,
  #                       fdata_cname = samp_id,
  #                       emeta_cname = pro_id,
  #                       data_scale = data_scale,
  #                       is_normalized = is_normalized,
  #                       check.names = check_names)

  return (
    list(e_data = res,
         e_meta = NULL)
  )

}

#' Applies rrollup function
#'
#' This function applies the rrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#'
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#'
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying rrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the rrollup method, peptides are scaled based on a reference
#'   peptide and protein abundance is set as the mean of these scaled peptides.
#'
#' @references Matzke, M. M., Brown, J. N., Gritsenko, M. A., Metz, T. O.,
#'   Pounds, J. G., Rodland, K. D., ... Webb-Robertson, B.-J. (2013). \emph{A
#'   comparative analysis of computational approaches to relative protein
#'   quantification using peptide peak intensities in label-free LC-MS
#'   proteomics experiments}. Proteomics, 13(0), 493-503.
#'
#'   Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins, J. N.,
#'   Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical tool for
#'   quantitative analysis of -omics data}. Bioinformatics (Oxford, England),
#'   24(13), 1556-1558.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = rrollup(pepData = pep_object)
#' }
#'
#' @rdname rrollup
#'
rrollup <- function (pepData, combine_fn, parallel = TRUE) {

  # Nab the check names attribute because data.frame conventions are making my
  # life miserable. In fact, if it weren't for check.names and COVID my life
  # would be pretty awesome right now.
  check_names <- get_check_names(pepData)

  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) {

    stop("pepData is not an object of the appropriate class")

  }

  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }

  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname

  # Combine e_data and e_meta (just the peptide and protein ID columns) into one
  # data frame by peptide ID. This is a right join with e_meta being the data
  # frame on the right.
  temp <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
                y = pepData$e_data,
                by = pep_id,
                all.x = FALSE,
                all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    data.frame(check.names = check_names)

  #pull protein column from temp and apply unique function
  unique_proteins <- unique(temp[[pro_id]])

  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  r <- foreach::foreach(i=1:length(unique_proteins))%dopar%{

    row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
    current_subset<- temp[row_ind,]
    current_subset<- current_subset[,-which(names(temp) == pro_id)]

    #### Perform R_Rollup ####
    ## store number of peptides ##
    num_peps = nrow(current_subset)

    # Create a matrix with one row and the same number of columns as the number
    # of samples. This will be used to create a data frame with only one row
    # after using the apply function columnwise.
    res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))

    ## if only 1 peptide, set the protein value to the peptide ##
    if(num_peps==1){
      protein_val = unlist(current_subset)
    }else{
      ## Step 1: Select Reference Peptide -- peptide with least amount of
      ## missing data ##
      na.cnt = apply(is.na(current_subset),1,sum)
      least.na = which(na.cnt == min(na.cnt))

      ## If tied, select one with highest median abundance##
      if(length(least.na)>1){
        mds = apply(current_subset,1,median,na.rm=T)[least.na]
        least.na = least.na[which(mds==max(mds))]
      }
      prot_val = unlist(current_subset[least.na,])

      ## Step 2: Ratio all peptides to the reference.  Since the data is on the
      ## log scale, this is the difference ##
      scaling_factor = apply(rep(as.numeric(prot_val),
                                 each = nrow(current_subset)) - current_subset,
                             1,
                             median,
                             na.rm=T)

      ## Step 3: Use the median of the ratio as a scaling factor for each
      ## peptide ##
      x_scaled = current_subset + rep(scaling_factor, ncol(current_subset))

      ## Step 4: Set Abundance as Median Peptide Abundance ##
      protein_val = apply(x_scaled, 2, combine_fn)

    }

    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] = protein_val

    # Convert the single row matrix to a single row data frame and rename the
    # columns to the original column names.
    res <- data.frame(res)
    names(res) <- names(current_subset)

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case res will be assigned to the ith element of r.
    res

  }

  # Combine the protein abundances (or is it abundanci?).
  final_result <- data.frame(unique_proteins,
                             data.table::rbindlist(r),
                             check.names = check_names)
  names(final_result)[1] <- pro_id

  # Extricate attribute info for creating the proData object.
  # data_scale <- attr(pepData, "data_info")$data_scale
  # is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Add check here from Lisa's email -------------------------------------------

  # Create a proData object with the quantitated proteins.
  # prodata <- as.proData(e_data = final_result,
  #                       f_data = pepData$f_data,
  #                       e_meta = dplyr::select(pepData$e_meta,
  #                                              -rlang::sym(pep_id)),
  #                       edata_cname = pro_id,
  #                       fdata_cname = samp_id,
  #                       emeta_cname = pro_id,
  #                       data_scale = data_scale,
  #                       is_normalized = is_normalized,
  #                       check.names = check_names)

  return (
    list(e_data = final_result,
         e_meta = NULL)
  )

}

#' Applies qrollup function
#'
#' This function applies the qrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#' @param qrollup_thresh numeric value between 0 and 1 inclusive. Peptides above
#'   this threshold are used to roll up to the protein level
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying qrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the qrollup method, peptides are selected according to a user
#'   selected abundance cutoff value (qrollup_thresh), and protein abundance is
#'   set as the mean of these selected peptides.
#'
#' @references Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins,
#'   J. N., Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical
#'   tool for quantitative analysis of -omics data}. Bioinformatics (Oxford,
#'   England), 24(13), 1556-1558.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = qrollup(pepData = pep_object, qrollup_thresh = 2)
#' }
#'
#' @rdname qrollup
#'
qrollup <- function (pepData, qrollup_thresh,
                     combine_fn, parallel = TRUE) {

  # Nab the check names attribute because data.frame conventions are making my
  # life miserable. In fact, if it weren't for check.names and COVID my life
  # would be pretty awesome right now.
  check_names <- get_check_names(pepData)

  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")

  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  # check that a qrollup_thresh is numeric and between 0 and 1#
  if(!is.numeric(qrollup_thresh) || qrollup_thresh > 1 || qrollup_thresh < 0){
    stop("qrollup_thresh must be Numeric and between 0 and 1")
  }

  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname

  # Combine e_data and e_meta (just the peptide and protein ID columns) into one
  # data frame by peptide ID. This is a right join with e_meta being the data
  # frame on the right.
  temp <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
                y = pepData$e_data,
                by = pep_id,
                all.x = FALSE,
                all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    data.frame(check.names = check_names)

  #pull protein column from temp and apply unique function
  unique_proteins<- unique(temp[[pro_id]])

  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  r<-foreach::foreach(i=1:length(unique_proteins))%dopar%{

    row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
    current_subset<- temp[row_ind,]
    current_subset<- current_subset[,-which(names(temp) == pro_id)]

    #### Perform Q_Rollup ####
    ## store number of peptides ##
    num_peps = nrow(current_subset)

    res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
    ## if only 1 peptide, set the protein value to the peptide ##
    if(num_peps==1){
      protein_val = unlist(current_subset)
      peps_used <- 1
    }else{
      ## Step 1: Subset peptides whose abundance is >= to qrollup_thresh ##
      means = apply(current_subset,1,mean,na.rm=T)
      quantil = quantile(means, probs = qrollup_thresh, na.rm = T)

      new_subset = current_subset[which(means >= quantil), ]
      peps_used <- nrow(new_subset)

      #after step 1 if only 1 peptide, set the protein value to the peptide
      if(nrow(new_subset) == 1){
        protein_val = unlist(new_subset)
      }else{
        ## Step 2: Set protein abundance as the mean/median of peptide abundances
        protein_val = apply(new_subset, 2, combine_fn)
      }
    }

    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] <- protein_val

    # Convert the single row matrix to a single row data frame and rename the
    # columns to the original column names.
    res <- data.frame(res)
    names(res) <- names(current_subset)

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case the list containing res and peps_used will be assigned to the ith
    # element of r.
    list(res, peps_used)

  }

  # Combine the protein abundances (or is it abundanci?).
  final_result <- data.frame(unique_proteins,
                             data.table::rbindlist(
                               lapply(r, function(x) x[[1]])
                             ),
                             check.names = check_names)
  names(final_result)[1] <- pro_id

  # Combine the peptide counts with their corresponding proteins.
  temp_pepes <- data.frame(final_result[, 1],
                           n_peps_used = sapply(r, function(x) x[[2]]),
                           check.names = check_names)
  names(temp_pepes)[1] <- pro_id

  # Combine the peptide counts with pepData$e_meta by protein. This is done to
  # preserve information in the rows of e_meta. For example, there could be
  # fewer unique proteins than there are rows of e_meta. If we do not combine
  # the peptide counts this way then only one row per protein will be kept and
  # the additional information will be lost.
  temp_emeta <- dplyr::left_join(x = pepData$e_meta,
                                 y = temp_pepes,
                                 by = pro_id) %>%
    # Remove peptide id column.
    dplyr::select(-rlang::sym(pep_id))

  # Extricate attribute info for creating the proData object.
  # data_scale <- attr(pepData, "data_info")$data_scale
  # is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Add check here from Lisa's email -------------------------------------------

  # Create a proData object with the quantitated proteins.
  # prodata <- as.proData(e_data = final_result,
  #                       f_data = pepData$f_data,
  #                       e_meta = temp_emeta,
  #                       edata_cname = pro_id,
  #                       fdata_cname = samp_id,
  #                       emeta_cname = pro_id,
  #                       data_scale = data_scale,
  #                       is_normalized = is_normalized,
  #                       check.names = check_names)

  return (
    list(e_data = final_result,
         e_meta = temp_emeta)
  )

}

#' Applies zrollup function
#'
#' This function applies the zrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying zrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the zrollup method, peptides are scaled as, pep_scaled = (pep -
#'   median)/sd, and protein abundance is set as the mean of these scaled
#'   peptides.
#'
#' @references Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins,
#'   J. N., Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical
#'   tool for quantitative analysis of -omics data}. Bioinformatics (Oxford,
#'   England), 24(13), 1556-1558.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = zrollup(pepData = pep_object)
#' }
#'
#' @rdname zrollup
#'
zrollup <- function (pepData, combine_fn, parallel = TRUE) {

  # Preliminary checks ---------------------------------------------------------

  # Nab the check names attribute because data.frame conventions are making my
  # life miserable. In fact, if it weren't for check.names and COVID my life
  # would be pretty awesome right now.
  check_names <- get_check_names(pepData)

  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")

  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }

  # Grab useful stuff ----------------------------------------------------------

  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname

  # Combine e_data and e_meta (just the peptide and protein ID columns) into one
  # data frame by peptide ID. This is a right join with e_meta being the data
  # frame on the right.
  temp <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
                y = pepData$e_data,
                by = pep_id,
                all.x = FALSE,
                all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    data.frame(check.names = check_names)

  #pull protein column from temp and apply unique function
  unique_proteins <- unique(temp[[pro_id]])

  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  r <- foreach::foreach(i=1:length(unique_proteins))%dopar%{

    row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
    current_subset<- temp[row_ind,]
    current_subset<- current_subset[,-which(names(temp) == pro_id)]

    #### Perform Z_Rollup ####
    ## store number of peptides ##
    num_peps = nrow(current_subset)

    res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
    ## Step 1: Compute mean and sd of peptides ##
    mds = apply(current_subset, 1, median, na.rm = T)
    sds = apply(current_subset, 1, sd, na.rm = T)

    ## Step 2: Scale peptide data as pep_scaled = (pep - median)/sd
    medians_mat = matrix(mds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)
    standiv_mat = matrix(sds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)

    peptides_scaled = apply((current_subset - medians_mat)/standiv_mat,
                            2, combine_fn)

    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] = peptides_scaled

    # Convert the single row matrix to a single row data frame and rename the
    # columns to the original column names.
    res <- data.frame(res, check.names = check_names)
    names(res) <- names(current_subset)

    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case res will be assigned to the ith element of r.
    res

  }

  # Combine the protein abundances (or is it abundanci?).
  final_result <- data.frame(unique_proteins,
                             data.table::rbindlist(r),
                             check.names = check_names)
  names(final_result)[1] <- pro_id

  # Extricate attribute info for creating the proData object.
  # data_scale <- attr(pepData, "data_info")$data_scale
  # is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Add check here from Lisa's email -------------------------------------------

  # Create a proData object with the quantitated proteins.
  # prodata <- as.proData(e_data = final_result,
  #                       f_data = pepData$f_data,
  #                       e_meta = dplyr::select(pepData$e_meta,
  #                                              -rlang::sym(pep_id)),
  #                       edata_cname = pro_id,
  #                       fdata_cname = samp_id,
  #                       emeta_cname = pro_id,
  #                       data_scale = data_scale,
  #                       is_normalized = is_normalized,
  #                       check.names = check_names)

  return (
    list(e_data = final_result,
         e_meta = NULL)
  )

}

# Functions for combining data in rollup methods -------------------------------

combine_fn_mean <- function (x) {

  if (all(is.na(x)))
    mean(x) else
      mean(x, na.rm = T)

}

combine_fn_median <- function (x) median(x, na.rm = T)
