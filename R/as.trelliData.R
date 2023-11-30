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
#'   "pepData", "isobaricpepData", "proData", "metabData", "lipidData", or
#'   "nmrData".
#' @param data_scale_original A character string indicating original scale
#'   of the data. Valid values are: 'log2', 'log', 'log10', or 'abundance'.
#'   Default is abundance.
#' @param data_scale A character string indicating the scale to transform the
#'   data to. Valid values are: 'log2', 'log', 'log10', or 'abundance'. If the
#'   value is the same as data_scale_original, then transformation is not
#'   applied. Default is log2.
#' @param normalization_fun A character string indicating the pmartR
#'   normalization function to use on the data, if is_normalized is FALSE.
#'   Acceptable choices are 'global', 'loess', and 'quantile'.
#' @param normalization_params A vector or list where the normalization
#'   parameters are the names, and the parameter values are the list values. For
#'   example, an acceptable entry for 'normalize_global' would be
#'   list("subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE,
#'   "backtransform" = TRUE).
#' @param is_normalized A logical indicator of whether the data is already
#'   normalized (and will therefore skip the normalization step)
#' @param force_normalization A logical indicator to force normalization that is
#'   not required for both isobaric protein and NMR data
#'
#' @return An object of class 'trelliData' containing the raw data.
#'  To be passed to trelliscope building functions.  
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' 
#' # Simple example 
#' trelliData1 <- as.trelliData.edata(e_data = pep_edata,
#'                                    edata_cname = "Peptide",
#'                                    omics_type = "pepData")
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

  # Build an omics data object--------------------------------------------------

  # Generate a f data frame
  fdata <- data.frame("Sample" = colnames(edata)[colnames(edata) != edata_cname], "Condition" = NA)
  fdata_cname <- "Sample"

  # Build the omics object
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
    trelliData.omics = omicsData$e_data %>%
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Sample = name, Abundance = value),
    trelliData.stat = NULL,
    omicsData = omicsData,
    statRes = NULL
  )

  # Save Panel By information and set class."panel_by_options" list the potential
  # inputs for the panel_by function. "panel_by_omics"/"panel_by_stat" will hold
  # the column name of the trelliData.omics/trelliData.stat that the data has
  # been grouped by. And "panel_by" tracks whether the panel_by function has been
  # applied or not.
  attr(trelliData, "fdata_col") <- "Sample"
  attr(trelliData, "emeta_col") <- NULL
  attr(trelliData, "panel_by_options") <- c(edata_cname, fdata_cname)
  attr(trelliData, "panel_by_omics") <- NA
  attr(trelliData, "panel_by_stat") <- NA
  attr(trelliData, "panel_by") <- FALSE
  class(trelliData) <- c("trelliData", "trelliData.edata")

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
#'   in long format to create the trelliData.omics dataframe, and e_meta is
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
    if (any(class(omicsData) %in% c("pepData", "isobaricpepData", "proData", "metabData", "lipidData", "nmrData")) == FALSE) {
      stop(paste(class(omicsData), "is not a supported omicsData class."))
    }

    # ...it must be log transformed if it's not NMR or isobaric
    if (any(class(omicsData) %in% c("nmrData", "isobaricpepData")) == FALSE &
      get_data_scale(omicsData) %in% c("log2", "log", "log10") == FALSE) {
      stop("omicsData must be log transformed.")
    }

    # ...it must be normalized if it's not NMR nor isobaric
    if (require_normalization & any(class(omicsData) %in% c("isobaricpepData", "nmrData")) == FALSE) {
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
  trelliData.omics <- NULL
  trelliData.stat <- NULL
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
    trelliData.omics <- edata %>%
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Abundance = value, !!fdata_cname := name)

    # Add group_designation if it exists
    if (!is.null(attributes(omicsData)$group_DF)) {
      trelliData.omics <- dplyr::left_join(trelliData.omics, attributes(omicsData)$group_DF, by = fdata_cname)
    }

    # Add emeta columns if emeta exists
    if (!is.null(omicsData$e_meta)) {
      # Pull emeta
      emeta <- omicsData$e_meta

      # Add emeta columns
      trelliData.omics <- dplyr::left_join(trelliData.omics, emeta, by = edata_cname)
    }
  } else {
    omicsData <- NULL
  }

  # Format statRes if applicable
  if (!is.null(statRes)) {
    # Get edata cname
    edata_cname <- pmartR::get_edata_cname(statRes)

    # Get column names of all fold changes, as well as p-values
    pvalue_cols <- colnames(statRes)[grepl("P_value", colnames(statRes))]
    fold_change_cols <- colnames(statRes)[grepl("Fold_change", colnames(statRes))]
    
    # Change class to prevent dplyr issues
    class(statRes) <- "data.frame"
    
    # Pivot longer so that the first column is the edata_cname, extract comparison,
    # panel_by comparison, nest dataframes, and then extract the p_value and fold_change
    # for each group
    trelliData.stat <- statRes %>%
      dplyr::select(dplyr::all_of(c(edata_cname, pvalue_cols, fold_change_cols))) %>%
      tidyr::pivot_longer(dplyr::all_of(c(pvalue_cols, fold_change_cols))) %>%
      dplyr::mutate(
        Comparison = gsub("P_value_A_|P_value_G_|Fold_change_", "", name),
        Metric = lapply(name, function(x) {
          if (grepl("P_value_A", x)) {
            return("p_value_anova")
          } else if (grepl("P_value_G", x)) {
            return("p_value_gtest")
          } else {
            return("fold_change")
          }
        }) %>% unlist()
      ) %>%
      dplyr::select(-name) %>%
      dplyr::group_by(dplyr::across(c(Comparison, !!dplyr::sym(edata_cname)))) %>%
      dplyr::summarise(
        "p_value_anova" = ifelse(length(value[which(Metric == "p_value_anova")]) == 0, NA, value[which(Metric == "p_value_anova")]),
        "p_value_gtest" = ifelse(length(value[which(Metric == "p_value_gtest")]) == 0, NA, value[which(Metric == "p_value_gtest")]),
        "fold_change" = value[which(Metric == "fold_change")]
      ) %>%
      dplyr::relocate(!!dplyr::sym(edata_cname))

    # Add emeta columns if emeta exists
    if (!is.null(omicsData$e_meta)) {
      # Add emeta columns
      trelliData.stat <- dplyr::left_join(trelliData.stat, emeta, by = emeta_cname)
    }
    
    # Fix statRes
    class(statRes) <- c("statRes", "data.frame")
    
  } else {
    statRes <- NULL
  }
  
  # Generate a trelliData object
  trelliData <- list(
    trelliData.omics = trelliData.omics,
    trelliData.stat = trelliData.stat,
    omicsData = omicsData,
    statRes = statRes
  )

  # Add fdata_cname and emeta column names as attributes
  if (!is.null(fdata_cname)) {
    attr(trelliData, "fdata_col") <- fdata_cname
  }

  if (!is.null(omicsData$e_meta)) {
    attr(trelliData, "emeta_col") <- colnames(omicsData$e_meta)[colnames(omicsData$e_meta) != edata_cname]
  } else {
    attr(trelliData, "emeta_col") <- NULL
  }

  # Save Panel By information and set class."panel_by_options" list the potential
  # inputs for the panel_by function. "panel_by_omics"/"panel_by_stat" will hold
  # the column name of the trelliData.omics/trelliData.stat that the data has
  # been grouped by. And "panel_by" tracks whether the panel_by function has been
  # applied or not.
  group_options <- c(colnames(trelliData.omics), colnames(trelliData.stat)) %>% unique()
  group_nonoptions <- c("Abundance", "Comparison", "p_value_anova", "p_value_gtest", "fold_change", "Group")
  group_options <- group_options[group_options %in% group_nonoptions == FALSE]
  attr(trelliData, "panel_by_options") <- group_options
  attr(trelliData, "panel_by_omics") <- NA
  attr(trelliData, "panel_by_stat") <- NA
  attr(trelliData, "panel_by") <- FALSE
  class(trelliData) <- c("trelliData")

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
    if (is.na(attr(trelliData, "panel_by_omics"))) {
      paneled_by <- attr(trelliData, "panel_by_stat")
    } else {
      paneled_by <- attr(trelliData, "panel_by_omics")
    }
    stop(paste("trelliData has already been paneled by", paneled_by))
  }

  # Determine which dataframes this grouping variable applies to----------------

  # Test if grouping applies to omicsData
  if (!is.null(trelliData$trelliData.omics) &&
    panel %in% colnames(trelliData$trelliData.omics)) {
    apply_to_omics <- TRUE
  } else {
    apply_to_omics <- FALSE
  }

  # Test if grouping applies to statRes
  if (!is.null(trelliData$trelliData.stat) &&
    panel %in% colnames(trelliData$trelliData.stat)) {
    apply_to_stat <- TRUE
  } else {
    apply_to_stat <- FALSE
  }

  # Test panel_by option--------------------------------------------------------

  # If there are instances of groups with less than 3 data points, it is not a
  # good grouping option.
  .test_grouping <- function(trelliData_subclass, trelliData_subclass_name) {
    # Get the smallest group size
    smallest_group_size <- trelliData_subclass %>%
      dplyr::group_by_at(panel) %>%
      dplyr::summarize(N = dplyr::n()) %>%
      dplyr::select(N) %>%
      min()

    # If the smallest group size is less than 3, then give warning
    if (smallest_group_size < 3) {
      warning(paste0(
        "Grouping by ", panel,
        " results in panels with less than 3 ",
        "data points in ", trelliData_subclass_name, "."
      ))
    }
  }

  # Onlt test grouping if it applies to omics
  if (apply_to_omics) {
    .test_grouping(
      trelliData$trelliData.omics,
      "trelliData.omics"
    )
  }

  # Group and nest samples------------------------------------------------------

  .group_samples <- function(trelliData_subclass) {
    trelliData_subclass %>%
      dplyr::group_by_at(panel) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::rename(Nested_DF = data)
  }

  # Group omicsData
  if (apply_to_omics) {
    trelliData$trelliData.omics <- .group_samples(trelliData$trelliData.omics)
    attr(trelliData, "panel_by_omics") <- panel
  }

  # Group statRes
  if (apply_to_stat) {
    trelliData$trelliData.stat <- .group_samples(trelliData$trelliData.stat)
    attr(trelliData, "panel_by_stat") <- panel
  }

  # Export results--------------------------------------------------------------
  attr(trelliData, "panel_by") <- TRUE
  return(trelliData)
}

#' @name trelli_pvalue_filter
#' 
#' @title Filter a paneled trelliData object by a p-value 
#' 
#' @description This use-case-specific function allows users to filter down their plots to a 
#'    specified p-value IF statistics data has been included. This function is mostly
#'    relevant to the MODE application. 
#'    
#' @param trelliData A trelliData object with statistics results (statRes). Required.
#' @param p_value_test A string to indicate which p_values to plot. Acceptable
#'    entries are "anova" or "gtest". Default is "anova". Unlike the
#'    plotting functions, here p_value_test cannot be null. Required.
#' @param p_value_thresh A value between 0 and 1 to indicate the p-value threshold 
#'    at which to keep plots. Default is 0.05. Required.
#' @param comparison The specific comparison to filter significant values to. Can
#'    be null. See attr(statRes, "comparisons") for the available options. Optional.
#' 
#' @return A paneled trelliData object with only plots corresponding to significant
#' p-values from a statistical test.
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' library(pmartRdata)
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
#' trelliData3 <- as.trelliData(statRes = statRes)
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' # Filter a trelliData object with only statistics results, while not caring about a comparison
#' trelli_pvalue_filter(trelliData3, p_value_test = "anova", p_value_thresh = 0.1)
#' 
#' # Filter a trelliData object with only statistics results, while caring about a specific comparison
#' trelli_pvalue_filter(
#'  trelliData3, p_value_test = "anova", p_value_thresh = 0.1, comparison = "Phenotype3_vs_Phenotype2")
#' 
#' # Filter both a omicsData and statRes object, while not caring about a specific comparison
#' trelli_pvalue_filter(trelliData4, p_value_test = "anova", p_value_thresh = 0.001)
#' 
#' # Filter both a omicsData and statRes object, while caring about a specific comparison
#' trelli_pvalue_filter(
#'  trelliData4, p_value_test = "gtest", p_value_thresh = 0.25, 
#'  comparison = "Phenotype3_vs_Phenotype2"
#' )
#' 
#' 
#' }
#'
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_pvalue_filter <- function(trelliData, 
                                 p_value_test = "anova", 
                                 p_value_thresh = 0.05,
                                 comparison = NULL) {
  
  # Run initial checks----------------------------------------------------------
  
  # trelliData object must be of the trelliData class
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData.")
  }
  
  # trelliData must contain statistics 
  if (is.null(trelliData$statRes)) {
    stop("trelliData must contain a statRes object.")
  }
  
  # Ensure the p_value_test is anova, gtest, or combined
  if (p_value_test %in% c("anova", "gtest", "combined") == FALSE) {
    stop("p_value_test must be anova, gtest, or combined.")
  }
  
  # Ensure that p_value threshold is a single numeric
  if (!is.numeric(p_value_thresh)) {
    stop("p_value_threshold must be a number.")
  }
  p_value_thresh <- abs(p_value_thresh)
  
  # If comparison is not NULL...
  if (!is.null(comparison)) {
    # Ensure that comparison is an acceptable input
    if (!is.character(comparison) | length(comparison) > 1) {
      stop("comparison must be a string of length 1.")
    }
    
    # Ensure that comparison is from the comparisons lists
    Comparisons <- attr(trelliData$statRes, "comparisons")
    if (comparison %in% Comparisons == FALSE) {
      stop(paste0(comparison, "is not an acceptable comparison"))
    }
  }
  
  # Filter by p_value-----------------------------------------------------------
  
  # Pull the biomolecule name 
  biomolecule <- attr(trelliData$statRes, "cnames")$edata_cname
  
  # Filter down to only a comparison if the user specifies one
  if (!is.null(comparison)) {
    trelliData$trelliData.stat <- trelliData$trelliData.stat[trelliData$trelliData.stat$Comparison == comparison,]
  }
  
  # Filter statRes
  if (p_value_test == "anova") {
    trelliData$trelliData.stat <- trelliData$trelliData.stat[trelliData$trelliData.stat$p_value_anova <= p_value_thresh,]
  } else if (p_value_test == "gtest") {
    trelliData$trelliData.stat <- trelliData$trelliData.stat[trelliData$trelliData.stat$p_value_gtest <= p_value_thresh,] 
  } 
  
  # Filter omics down 
  if (!is.null(trelliData$trelliData.omics)) {
    biomolecule_subset <-  trelliData$trelliData.stat[,biomolecule] %>% unlist()
    trelliData$trelliData.omics <- trelliData$trelliData.omics[trelliData$trelliData.omics[[biomolecule]] %in% biomolecule_subset,]
  }
  
  return(trelliData)
  
}
