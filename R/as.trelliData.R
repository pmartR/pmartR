#' @name .is_edata
#'
#' @title Test if a file is an edata file
#'
#' @param edata Must be a dataframe. Required.
#'
#' @return A boolean where TRUE means the file is an acceptable edata file.
#' 
.is_edata <- function(edata) {
  # If edata is NULL, return FALSE
  if (is.null(edata) || is.data.frame(edata) == FALSE) {
    message("edata must be a data.frame.")
    return(FALSE)
  }

  # edata must have at least 2 samples
  if (ncol(edata) < 3) {
    message("edata must have at least 3 columns, a 'descriptor' column and at least 2 samples.")
    return(FALSE)
  }

  ## All columns with the exception of one column MUST be numeric

  # Get counts of the number of numeric columns
  LogicCounts <- lapply(1:ncol(edata), function(col) {
    is.numeric(edata[, col])
  }) %>%
    unlist() %>%
    table()

  # If more than one column is not numeric, return error
  if ("FALSE" %in% names(LogicCounts) && LogicCounts[["FALSE"]] > 1) {
    message("All columns in edata must be numeric, with exception of a 'descriptor' column.")
    return(FALSE)
  }

  # Otherwise, return True
  return(TRUE)
}

#' @name as.trelliData.edata
#'
#' @title Generate an object from edata to pass to trelliscope building
#'   functions
#'
#' @description The only acceptable input file type is a single edata file.
#'   Transformation and normalization must be specified. Isobaric protein or NMR
#'   data does not need to be normalized.
#'
#' @param e_data a \eqn{p * (n + 1)} data.frame of expression data, where
#'   \eqn{p} is the number of biomolecules observed and \eqn{n} is the number of
#'   samples (an additional biomolecule identifier/name column should also be
#'   present anywhere in the data.frame). Each row corresponds to data for one
#'   biomolecule. One column specifying a unique identifier for each biomolecule
#'   (row) must be present. We do not recommend passing data that requires
#'   reference normalization (isobaric, nmr, etc.)
#' @param edata_cname character string specifying the name of the column
#'   containing the biomolecule identifiers. It should be the only non-numeric
#'   colummn in edata.
#' @param omics_type A string specifying the data type. Acceptable options are
#'   "pepData", "isobaricpepData", "proData", "metabData", "lipidData", 
#'   "nmrData", or "seqData". 
#' @param data_scale_original A character string indicating original scale
#'   of the data. Valid values are: 'log2', 'log', 'log10', or 'abundance'.
#'   Default is abundance. This parameter is ignored if the data is "seqData". 
#' @param data_scale A character string indicating the scale to transform the
#'   data to. Valid values are: 'log2', 'log', 'log10', or 'abundance'. If the
#'   value is the same as data_scale_original, then transformation is not
#'   applied. Default is log2. This parameter is ignored if the data is "seqData". 
#' @param normalization_fun A character string indicating the pmartR
#'   normalization function to use on the data, if is_normalized is FALSE.
#'   Acceptable choices are 'global', 'loess', and 'quantile'. This parameter is
#'   ignored if the data is "seqData". 
#' @param normalization_params A vector or list where the normalization
#'   parameters are the names, and the parameter values are the list values. For
#'   example, an acceptable entry for 'normalize_global' would be
#'   list("subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE,
#'   "backtransform" = TRUE). This parameter is ignored if the data is "seqData". 
#' @param is_normalized A logical indicator of whether the data is already
#'   normalized (and will therefore skip the normalization step). This parameter is
#'   ignored if the data is "seqData". 
#' @param force_normalization A logical indicator to force normalization that is
#'   not required for both isobaric protein and NMR data. This parameter is ignored
#'   if the data is "seqData." 
#'
#' @return An object of class 'trelliData' containing the raw data.
#'  To be passed to trelliscope building functions.  
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' 
#' ###########################
#' ## MS/NMR OMICS EXAMPLES ##
#' ###########################
#' 
#' # Simple MS/NMR Example 
#' trelliData1 <- as.trelliData.edata(e_data = pep_edata,
#'                                    edata_cname = "Peptide",
#'                                    omics_type = "pepData")
#'  
#' ######################
#' ## RNA-SEQ EXAMPLES ##  
#' ######################
#'                                     
#' # RNA-seq Example
#' trelliData_seq1 <- as.trelliData.edata(e_data = rnaseq_edata, 
#'                                       edata_cname = "Transcript",
#'                                       omics_type = "seqData")
#'                                   
#'           
#' @author David Degnan, Daniel Claborne, Lisa Bramer
#'
#' @export
as.trelliData.edata <- function(e_data,
                                edata_cname,
                                omics_type,
                                data_scale_original = "abundance",
                                data_scale = "log2",
                                normalization_fun = "global",
                                normalization_params = list(
                                  "subset_fn" = "all", "norm_fn" = "median",
                                  "apply_norm" = TRUE, "backtransform" = TRUE
                                ),
                                is_normalized = FALSE,
                                force_normalization = FALSE) {
  edata <- e_data

  # Initial checks -------------------------------------------------------------

  # Run the .is_edata check to confirm the file is an acceptable edata file
  if (!.is_edata(edata)) {
    stop("")
  }
  
  # Skip most checks if the data is seqData
  if (omics_type != "seqData") {
    
    # Check that the object type is in the acceptable class
    if (omics_type %in% c("pepData", "isobaricpepData", "proData", "metabData", "lipidData", "nmrData") == FALSE) {
      stop(paste(omics_type, "is not an acceptable omics_type."))
    }
    
    # Check that both data_scale_original and data scale are an acceptable option.
    log_transforms <- c("abundance", "log2", "log", "log10")
    if (data_scale_original %in% log_transforms == FALSE) {
      stop(paste(data_scale_original, "is not an acceptable data scale."))
    }
    if (data_scale %in% log_transforms == FALSE) {
      stop(paste(data_scale, "is not an acceptable data scale."))
    }
    
    # Check the normalization function. Isobaric is not included since it requires f_data.
    if (normalization_fun %in% c("global", "loess", "quantile") == FALSE) {
      stop(paste(normalization_fun, "is not an acceptable normalization function type."))
    }
    
    # Normalization parameters should have apply_norm in it
    if (!is.null(normalization_params) && normalization_fun == "global") {
      if ("apply_norm" %in% names(normalization_params) == FALSE || (normalization_params$apply_norm == FALSE)) {
        stop("apply_norm must be TRUE to apply normalization parameters.")
      }
    } 
    
  } else {
    message(paste("Notice: seqData will be log count per million (LCPM) transformed for visualization purposes.",
                  "No other transformation method or normalization will be applied to the data.",
                  "Therefore, the statistics results may not perfectly match visualized patterns of expressed data.",
                  "We also suggest filtering out biomolecules before applying visualizations, as seqData tends to be quite large."))
  }

  # Build an omics data object--------------------------------------------------

  # Generate a f data frame
  fdata <- data.frame("Sample" = colnames(edata)[colnames(edata) != edata_cname], "Condition" = NA)
  fdata_cname <- "Sample"

  # Diverge building paths depending on input data type 
  if (omics_type != "seqData") {
    
    # Build omics object - this also checks that edata_cname is in edata 
    omicsData <- eval(parse(text = paste0(
      "as.", omics_type,
      "(e_data = edata, edata_cname = edata_cname, f_data = fdata, fdata_cname = fdata_cname,
      data_scale = data_scale_original)"
    )))
    
    # Transform if appropriate
    if (data_scale_original != data_scale) {
      omicsData <- edata_transform(omicsData, data_scale)
    }
    
    # If the data is already normalized, skip this step
    if (is_normalized == FALSE) {
      # Normalization is not required for both isobaric protein and NMR data, but can
      # be forced with .force_normalization
      if (omics_type %in% c("isobaricpepData", "nmrData") == FALSE | force_normalization) {
        # Get the normalization function
        norm_fun <- switch(normalization_fun,
                           "global" = normalize_global,
                           "global_basic" = normalize_global_basic,
                           "loess" = normalize_loess,
                           "quantile" = normalize_quantile
        )
        
        # Add omics data to normalization parameters
        normalization_params[["omicsData"]] <- omicsData
        
        # Apply normalization with its parameters
        omicsData <- do.call(norm_fun, normalization_params)
      }
    } else {
      attr(omicsData, "data_info")$norm_info$is_normalized <- TRUE
    }
    
    # Finally, generate the trelliData object-------------------------------------
    
    # Put the edata into the trelliData omics slot
    trelliData <- list(
      trelliData = omicsData$e_data %>%
        tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
        dplyr::rename(Sample = name, Abundance = value),
      omicsData = omicsData,
      statRes = NULL
    )
    
  } else {
    
    # Build omics object - this also checks that edata_cname is in edata 
    omicsData <- eval(parse(text = paste0(
      "as.", omics_type,
      "(e_data = edata, edata_cname = edata_cname, f_data = fdata, fdata_cname = fdata_cname,
      data_scale = 'counts')"
    )))
    
    # Conduct the lcpm transformation
    biomolecules <- omicsData$e_data[[edata_cname]]
    temp_data <- omicsData$e_data %>% dplyr::select(-dplyr::all_of(edata_cname))
    samp_sum <- apply(temp_data, 2, sum, na.rm = TRUE) + 1
    div_sum <- sweep((temp_data + .5), 2, samp_sum, `/`)
    lcpm <- log2(div_sum * 10^6)
    lcpm <- lcpm %>% dplyr::mutate(!!edata_cname := biomolecules)
    
    # Pivot the counts data.frame longer
    counts_pivoted <- omicsData$e_data %>%
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Sample = name, Count = value)
    
    # Pivot the lcpm data.frame longer
    lcpm_pivoted <- lcpm %>%
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Sample = name, LCPM = value)
    
    # Join datasets
    pivoted <- dplyr::left_join(counts_pivoted, lcpm_pivoted, by = c(edata_cname, "Sample"))
    
    # Finally, generate the trelliData object-------------------------------------
    
    # Put the edata into the trelliData omics slot
    trelliData <- list(
      trelliData = pivoted,
      omicsData = omicsData,
      statRes = NULL
    )
    
  }

  # Save Panel By information and set class."panel_by_options" list the potential
  # inputs for the panel_by function. "panel_by_omics"/"panel_by_stat" will hold
  # the column name of the trelliData/trelliData.stat that the data has
  # been grouped by. And "panel_by" tracks whether the panel_by function has been
  # applied or not.
  attr(trelliData, "edata_col") <- edata_cname
  attr(trelliData, "fdata_col") <- "Sample"
  attr(trelliData, "emeta_col") <- NULL
  attr(trelliData, "panel_by_options") <- c(edata_cname, fdata_cname)
  attr(trelliData, "panel_by") <- FALSE
  attr(trelliData, "panel_by_col") <- NA
  class(trelliData) <- c("trelliData", "trelliData.edata")
  
  # Add a special label for seqData
  if (omics_type == "seqData") {
    class(trelliData) <- c(class(trelliData), "trelliData.seqData")
  }

  return(trelliData)
}

#' @name as.trelliData
#'
#' @title Generate an object from omicsData and/or statRes objects to pass to
#'   trelliscope building functions
#'
#' @description Either an omicData and/or a statRes object are accepted.
#'   omicData must be transformed and normalized, unless the data is isobaric
#'   protein or NMR data. If group_designation() has been run on the omicData
#'   object to add "main_effects", the resulting plots will include groups. The
#'   main effects group_designation and e_meta columns are merged to the e_data
#'   in long format to create the trelliData dataframe, and e_meta is
#'   merged to statRes in long format to create trelliData.stat dataframe.
#'
#' @param omicsData an object of the class 'pepData', 'isobaricpepData',
#'   proData', 'metabData', 'lipidData', or 'nmrData', created by
#'   \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param statRes statRes an object of the class 'statRes', created by
#'   \code{\link{imd_anova}}
#'
#' @return An object of class 'trelliData' containing the raw data and optionally, statRes.  
#'  To be passed to trelliscope building functions.
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' library(pmartRdata)
#' 
#' ###########################
#' ## MS/NMR OMICS EXAMPLES ##
#' ###########################
#' 
#' # Transform the data
#' omicsData <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' 
#' # Group the data by condition
#' omicsData <- group_designation(omicsData = omicsData, main_effects = c("Phenotype"))
#' 
#' # Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = omicsData)
#' omicsData <- applyFilt(filter_object = imdanova_Filt, omicsData = omicsData,
#'                        min_nonmiss_anova = 2)
#'
#' # Normalize my pepData
#' omicsData <- normalize_global(omicsData, "subset_fn" = "all", "norm_fn" = "median",
#'                              "apply_norm" = TRUE, "backtransform" = TRUE)
#'
#' # Implement the IMD ANOVA method and compute all pairwise comparisons 
#' # (i.e. leave the `comparisons` argument NULL)
#' statRes <- imd_anova(omicsData = omicsData, test_method = 'combined')
#'
#' # Generate the trelliData object
#' trelliData2 <- as.trelliData(omicsData = omicsData)
#' trelliData3 <- as.trelliData(statRes = statRes)
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' ######################
#' ## RNA-SEQ EXAMPLES ##  
#' ######################
#' 
#' # Group data by condition
#' omicsData_seq <- group_designation(omicsData = rnaseq_object, main_effects = c("Virus"))
#' 
#' # Filter low transcript counts
#' omicsData_seq <- applyFilt(filter_object = total_count_filter(omicsData = omicsData_seq), 
#'  omicsData = omicsData_seq, min_count = 15)
#' 
#' # Select a normalization and statistics method (options are 'edgeR', 'DESeq2', and 'voom').
#' # See ?difexp_seq for more details
#' statRes_seq <- diffexp_seq(omicsData = omicsData_seq, method = "voom")
#' 
#' # Generate the trelliData object
#' trelliData_seq2 <- as.trelliData(omicsData = omicsData_seq)
#' trelliData_seq3 <- as.trelliData(statRes = statRes_seq)
#' trelliData_seq4 <- as.trelliData(omicsData = omicsData_seq, statRes = statRes_seq)
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
as.trelliData <- function(omicsData = NULL, statRes = NULL) {
  
  require_normalization <- TRUE

  # Initial checks--------------------------------------------------------------

  # Either omicsData or statRes must be included
  if (is.null(omicsData) & is.null(statRes)) {
    stop("At least 1 omicsData or 1 statRes object must be provided.")
  }

  # If omicsData is provided...
  if (!is.null(omicsData)) {
    
    # ...it must be an omics data object
    if (any(class(omicsData) %in% c("pepData", "isobaricpepData", "proData", "metabData", "lipidData", "nmrData", "seqData")) == FALSE) {
      stop(paste(class(omicsData), "is not a supported omicsData class."))
    }

    # ...it must be log transformed if it's not NMR, isobaric, or seqData
    if (any(class(omicsData) %in% c("nmrData", "isobaricpepData", "seqData")) == FALSE &
      get_data_scale(omicsData) %in% c("log2", "log", "log10") == FALSE) {
      stop("omicsData must be log transformed.")
    }

    # ...it must be normalized if it's not NMR, isobaric, or seqData
    if (require_normalization & any(class(omicsData) %in% c("isobaricpepData", "nmrData", "seqData")) == FALSE) {
      if (!get_data_norm(omicsData)) {
        stop("omicsData must be normalized.")
      }
    }

    # ...and if group_designation is set, the first main effect must not have any singletons
    if (!is.null(attributes(omicsData)$group_DF)) {
      group_counts <- attributes(omicsData)$group_DF$Group %>%
        table(dnn = "Group") %>%
        data.frame()
      if (1 %in% group_counts$Freq) {
        warning(paste(
          "singleton groups found in group_designation:",
          paste0(group_counts[group_counts$Freq == 1, "Group"], collapse = ", ")
        ))
      }
    }
  }

  # If statRes is provided
  if (!is.null(statRes)) {
    # ...it must be a statRes object
    if ("statRes" %in% class(statRes) == FALSE) {
      stop("statRes must be an object of the statRes class. See ?imd_anova.")
    }
  }

  # If both omicsData and staRes are provided...
  if (!is.null(omicsData) & !is.null(statRes)) {
    # ...they should be from the same dataset.

    # Get edata cname
    edata_cname <- pmartR::get_edata_cname(omicsData)

    # Confirm that all statRes biomolecules are in omicsData.
    if (is.null(statRes[[edata_cname]])) {
      stop("omicsData and statRes are from different datasets.")
    }
  }
  
  # Generate the trelliData object----------------------------------------------

  # Create placeholders for trelliData objects
  pre_trelliData <- NULL
  fdata_cname <- NULL
  emeta_cname <- NULL

  # Format omicsData if applicable
  if (!is.null(omicsData)) {
    
    # Get edata_cname, edata, and fdata_cname
    edata_cname <- pmartR::get_edata_cname(omicsData)
    edata <- omicsData$e_data
    fdata_cname <- pmartR::get_fdata_cname(omicsData)
    fdata <- omicsData$f_data

    # Generate telliData.omics object
    pre_trelliData <- edata %>%
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Abundance = value, !!fdata_cname := name)

    # Add group_designation if it exists
    if (!is.null(attributes(omicsData)$group_DF)) {
      pre_trelliData <- dplyr::left_join(pre_trelliData, attributes(omicsData)$group_DF, by = fdata_cname)
    }

    # Add emeta columns if emeta exists
    if (!is.null(omicsData$e_meta)) {
      # Pull emeta
      emeta <- omicsData$e_meta

      # Add emeta columns
      pre_trelliData <- dplyr::left_join(pre_trelliData, emeta, by = edata_cname)
    }
    
    # Clean up if seqData
    if (inherits(omicsData, "seqData")) {
      
      # Rename abundance as count
      pre_trelliData <- pre_trelliData %>% dplyr::rename(Count = Abundance)
      
      # Generate log counts per million (lcpm)
      biomolecules <- omicsData$e_data[[edata_cname]]
      temp_data <- omicsData$e_data %>% dplyr::select(-dplyr::all_of(edata_cname))
      samp_sum <- apply(temp_data, 2, sum, na.rm = TRUE) + 1
      div_sum <- sweep((temp_data + .5), 2, samp_sum, `/`)
      lcpm <- log2(div_sum * 10^6)
      lcpm <- lcpm %>% dplyr::mutate(!!edata_cname := biomolecules)
      
      # Pivot the lcpm data.frame longer
      lcpm_pivoted <- lcpm %>%
        tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
        dplyr::rename(!!fdata_cname := name, LCPM = value)
      
      # Add LCPM to the dataframe 
      pre_trelliData <- dplyr::left_join(pre_trelliData, lcpm_pivoted, by = c(edata_cname, fdata_cname))
      
    }
    
    
  } else {
    omicsData <- NULL
  }

  # Format statRes if applicable
  if (!is.null(statRes)) {
    
    # Get edata cname
    edata_cname <- pmartR::get_edata_cname(statRes)

    # Get column names of all fold changes, as well as p-values
    pvalue_anova_cols <- colnames(statRes)[grepl("P_value_A", colnames(statRes))]
    pvalue_gtest_cols <- colnames(statRes)[grepl("P_value_G", colnames(statRes))]
    fold_change_cols <- colnames(statRes)[grepl("Fold_change", colnames(statRes))]
    
    # Change class to prevent dplyr issues
    class(statRes) <- "data.frame"
    
    # trelliData is simply the statRes object if that is all that was provided 
    if (is.null(omicsData)) {pre_trelliData <- statRes} else {
      pre_trelliData <- dplyr::left_join(pre_trelliData, statRes, by = edata_cname)
    }
    
    # Fix statRes
    class(statRes) <- c("statRes", "data.frame")
    
  } else {
    statRes <- NULL
  }
  
  # Generate a trelliData object
  trelliData <- list(
    trelliData = pre_trelliData,
    omicsData = omicsData,
    statRes = statRes
  )
  
  # Save edata column name
  attr(trelliData, "edata_col") <- edata_cname

  # Add fdata_cname and emeta column names as attributes
  if (!is.null(fdata_cname)) {
    attr(trelliData, "fdata_col") <- fdata_cname
  } 

  if (!is.null(omicsData$e_meta)) {
    attr(trelliData, "emeta_col") <- colnames(omicsData$e_meta)[colnames(omicsData$e_meta) != edata_cname]
  }

  # Save Panel By information and set class."panel_by_options" list the potential
  # inputs for the panel_by function. "panel_by" will hold
  # the column name of the trelliData that the data has
  # been grouped by. And "panel_by" tracks whether the panel_by function has been
  # applied or not.
  group_options <- colnames(pre_trelliData)[grepl("P_value_|Fold_change_|Flag_|Count_|Mean_", colnames(pre_trelliData)) == FALSE]
  group_nonoptions <- c("Abundance", "Comparison", "Group", "Count", "LCPM")
  group_options <- group_options[group_options %in% group_nonoptions == FALSE]
  attr(trelliData, "panel_by_options") <- group_options
  attr(trelliData, "panel_by") <- FALSE
  attr(trelliData, "panel_by_col") <- NA
  class(trelliData) <- c("trelliData")
  
  if (inherits(omicsData, "seqData")) {
    
    class(trelliData) <- c(class(trelliData), "trelliData.seqData")
    
    if (!is.null(omicsData)) {
      message(paste("Notice: seqData will be log count per million (LCPM) transformed for visualization purposes.",
                    "No other transformation method or normalization will be applied to the data.",
                    "Therefore, the statistics results may not perfectly match visualized patterns of expressed data.",
                    "We also suggest filtering out biomolecules before applying visualizations, as seqData tends to be quite large."))
    } else {
      message(paste("Notice: seqData tends to be quite large, so we suggest filtering out biomolecules before applying visualizations."))  
    }
    
  }

  return(trelliData)
}

#' @name trelli_panel_by
#'
#' @title Set the "panel_by" variable for a trelliData object
#'
#' @description Allows for grouping omics or stats data for downstream plotting
#'     and cognostic functions
#'
#' @param trelliData A trelliscope data object made by as.trelliData or as.trelliData.edata. Required.
#' @param panel The name of a column in trelliData to panel the data by. Required.
#' 
#'
#' @return A trelliData object with attributes "panel_by_omics" or "panel_by_stat" to determine 
#' which columns to divide the data by.
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' library(pmartRdata)
#' 
#' trelliData1 <- as.trelliData.edata(e_data = pep_edata,
#'                                    edata_cname = "Peptide",
#'                                    omics_type = "pepData")
#' # Transform the data
#' omicsData <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' 
#' # Group the data by condition
#' omicsData <- group_designation(omicsData = omicsData, main_effects = c("Phenotype"))
#'
#' # Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = omicsData)
#' omicsData <- applyFilt(filter_object = imdanova_Filt, omicsData = omicsData,
#'                        min_nonmiss_anova = 2)
#'
#' # Normalize my pepData
#' omicsData <- normalize_global(omicsData, "subset_fn" = "all", "norm_fn" = "median",
#'                              "apply_norm" = TRUE, "backtransform" = TRUE)
#'
#' # Implement the IMD ANOVA method and compute all pairwise comparisons 
#' # (i.e. leave the `comparisons` argument NULL)
#' statRes <- imd_anova(omicsData = omicsData, test_method = 'combined')
#'
#' # Generate the trelliData object
#' trelliData2 <- as.trelliData(omicsData = omicsData)
#' trelliData3 <- as.trelliData(statRes = statRes)
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#'
#' ## "panel_by" with an edata file. 
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide")
#' trelli_panel_by(trelliData = trelliData1, panel = "Sample")
#' 
#' ## "panel_by" with trelliData containing omicsData. 
#' ## Generate trelliData2 using the example code for as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide")
#' trelli_panel_by(trelliData = trelliData2, panel = "RazorProtein")
#' 
#' ## "panel_by" with trelliData containing statRes. 
#' ## Generate trelliData3 using the example code for as.trelliData
#' trelli_panel_by(trelliData = trelliData3, panel = "Peptide")
#' 
#' ## "panel_by" with trelliData containing both omicsData and statRes. 
#' ## Generate trelliData4 using the example code for as.trelliData
#' trelli_panel_by(trelliData = trelliData4, panel = "Peptide")
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein")
#' trelli_panel_by(trelliData = trelliData4, panel = "SampleID")
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_panel_by <- function(trelliData, panel) {
  # Run initial checks----------------------------------------------------------

  # Confirm that trelliData is a trelliData object
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData from as.trelliData or as.trelliData.edata.")
  }

  # Confirm that panel is an acceptable option
  if (panel %in% attr(trelliData, "panel_by_options") == FALSE) {
    stop(paste(
      "panel is not an acceptable option. The following can be selected:",
      paste(attr(trelliData, "panel_by_options"), collapse = ", ")
    ))
  }
  
  # Confirm that panel_by is false
  if (attr(trelliData, "panel_by")) {
    stop(paste("trelliData has already been paneled by", attr(trelliData, "panel_by_col")))
  }
  
  # Give feedback based on the number of panels---------------------------------
  
  # Count the number of panels
  panel_count <- trelliData$trelliData[[panel]] %>% unique() %>% length()
  
  # If the number of panels is only a few, a trelliscope display doesn't make much sense
  if (panel_count <= 6) {
    message(paste0("Paneling by '", panel, "' results in only ", panel_count,
            " panels, which is a small number for trelliscope."))
  }
  
  # If the number of panels is many, encourage filtering 
  if (panel_count >= 10000) {
    message(paste0("Paneling by '", panel, "' results in ", panel_count,
                   " panels. Consider filtering the number of panels down",
                   " before building the trelliscope display if a faster build time",
                   " is desired."))
  }

  # Group and nest samples------------------------------------------------------
  
  # If the data is only statRes, there's no need to nest 
  if (!is.null(trelliData$omicsData)) {
    
    # Identify fold change columns
    foldchanges <- NULL 
    if (any(grepl("Fold_change_", colnames(trelliData$trelliData)))) {
      foldchanges <- colnames(trelliData$trelliData)[grepl("Fold_change_", colnames(trelliData$trelliData))]
    }
    
    # Identify p value columns
    pvalues <- NULL 
    if (any(grepl("P_value_", colnames(trelliData$trelliData)))) {
      pvalues <- colnames(trelliData$trelliData)[grepl("P_value_", colnames(trelliData$trelliData))]
    }
    
    # Identify expression columns (different depending on whether the data is MS/NMR or seq)
    if (inherits(trelliData, "trelliData.seqData")) {
      expression_cols <-  c("LCPM", "Count")
    } else {
      expression_cols <- c("Abundance")
    }
    
    # Identify which columns are needed, which columns should be tossed, and which
    # columns could stay based on whether the data is paneled by the edata_col, 
    # fdata_col, or emeta_col 
    if (panel == attr(trelliData, "edata_col")) {
  
      # In this case, only fold changes need to be in the nested columns, if applicable.
      needed_cols <- unique(c(attr(trelliData, "edata_col"), 
                              attr(trelliData, "fdata_col"), 
                              panel, "Group", expression_cols, 
                              foldchanges))
      
    } else if (panel == attr(trelliData, "fdata_col")) {
      
      # Only edata_col, fdata_col, and abundance are needed when grouping by sample
      needed_cols <- unique(c(attr(trelliData, "edata_col"), 
                              attr(trelliData, "fdata_col"), 
                              expression_cols))
      
      # Toss unnecessary columns
      trelliData$trelliData <- trelliData$trelliData %>% dplyr::select_at(needed_cols)
      
    } else if (panel %in% attr(trelliData, "emeta_col")) {
      
      # Count, Mean, and Flag columns are not needed 
      trelliData$trelliData <- trelliData$trelliData[!grepl("Count_|Mean_|Flag_", colnames(trelliData$trelliData))]
      
      # In this case, fold changes and pvalues should go into the nested column, if applicable
      needed_cols <- unique(c(attr(trelliData, "edata_col"), 
                              attr(trelliData, "fdata_col"), 
                              panel, "Group", expression_cols, 
                              foldchanges, pvalues))
      
    }
    
    # Nest the data
    nested <- trelliData$trelliData %>% 
      dplyr::select(dplyr::any_of(needed_cols)) %>%
      dplyr::group_by_at(panel) 
    
    # Add data that does not need to be nested
    if (all(colnames(trelliData$trelliData) %in% needed_cols) == FALSE) {
      
      browser()
      
      # Pull columns to append
      cols_2_append <- c(panel, colnames(trelliData$trelliData)[colnames(trelliData$trelliData) %in% needed_cols == FALSE])
      
      # See if data can be kept numeric. If there's any duplicates, then no. 
      as_numeric <- TRUE
      
      # Collapse by panel
      to_append <- trelliData$trelliData[,cols_2_append] %>%
        dplyr::group_by_at(panel) %>%
        dplyr::summarize_all(function(x) {
          
          # Grab unique entries
          unique_entries <- unique(x)
          
          # Return as is if numerics are acceptable
          if (length(unique_entries) == 1 & as_numeric) {
            return(unique_entries)
          } else if (length(unique_entries) == 1) {
            return(as.character(unique_entries))
          } else {
            as_numeric <- FALSE
            return(paste0(unique_entries, collapse = ";"))
          }
        })
      
      # Append collapsed data.frame
      nested <- dplyr::left_join(nested, to_append, by = panel)
    
    }
    
    # Add nested data
    trelliData$trelliData <- nested
    
  }
  
  


  # Export results--------------------------------------------------------------
  attr(trelliData, "panel_by") <- TRUE
  attr(trelliData, "panel_by_col") <- panel
  return(trelliData)
  
}
