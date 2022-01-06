#' @name .is_edata
#' 
#' @title Test if a file is an edata file
#' 
#' @param edata Must be a dataframe. Required. 
#' 
#' @return A boolean where TRUE means the file is an acceptable edata file. 
#' @examples 
#' \dontrun{
#' 
#' .is_edata(pmartRdata::lipid_edata)
#' 
#' }
#' @export
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
  LogicCounts <- lapply(1:ncol(edata), function(col) {is.numeric(edata[,col])}) %>% 
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
#' @title Generate an object from edata to pass to trelliscope building functions 
#' 
#' @description The only acceptable input file type is a single edata file. Transformation
#'    and normalization must be specified. Isobaric protein or NMR data does not
#'    need to be normalized. 
#'  
#' @param e_data a \eqn{p * (n + 1)} data.frame of expression data, where
#'   \eqn{p} is the number of peptides observed and \eqn{n} is the number of
#'   samples (an additional peptide identifier/name column should also be
#'   present anywhere in the data.frame). Each row corresponds to data for each
#'   peptide. One column specifying a unique identifier for each peptide (row)
#'   must be present.
#' @param edata_cname character string specifying the name of the column
#'   containing the protein identifiers (or other mapping variable) in
#'   \code{e_meta} (if applicable). Defaults to NULL. If \code{e_meta} is NULL,
#'   then either do not specify \code{emeta_cname} or specify it as NULL. If
#'   \code{e_meta} is NULL, then specify \code{emeta_cname} as NULL.
#' @param omics_type A string of the data type. Acceptable options are "pepData",
#'   "isobaricpepData", "proData", "metabData", "lipidData", or "nmrData".
#' @param data_scale_original A character string indicating original scale 
#'   of the data. Valid values are: 'log2', 'log', 'log10', or 'abundance'.
#'   Default is abundance. 
#' @param data_scale A character string indicating the scale to transform 
#'   the data to. Valid values are: 'log2', 'log', 'log10', or 'abundance'. If 
#'   the value is the same as data_scale_original, then transformation is not applied. 
#'   Default is log2.
#' @param normalization_fun A character string indicating the pmartR normalization
#'   function to use on the data, if is_normalized is FALSE. Acceptable choices are 
#'   'global', 'loess', and 'quantile'.
#' @param normalization_params A vector or list where the normalization parameters
#'   are the names, and the parameter values are the list values. For example, an 
#'   acceptable entry for 'normalize_global' would be list("subset_fn" = "all", "norm_fn" = "median", 
#'   "apply_norm" = TRUE, "backtransform" = TRUE).
#'  
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("lipid_edata_pos")
#' 
#' # Simple example 
#' trelliData <- as.trelliData.edata(e_data = lipid_edata_pos,
#'                                   edata_cname = "LipidCommonName",
#'                                   omics_type = "lipidData")
#'                                   
#' }   
#'           
#' @author David Degnan, Daniel Claborne
#' 
#' @export
as.trelliData.edata <- function(e_data, 
                                edata_cname,
                                omics_type,
                                data_scale_original = "abundance",
                                data_scale = "log2", 
                                normalization_fun = "global", 
                                normalization_params = list("subset_fn" = "all", "norm_fn" = "median", 
                                                            "apply_norm" = TRUE, "backtransform" = TRUE),
                                ...
) {
  
  .as.trelliData.edata(e_data, edata_cname, omics_type, data_scale_original, 
                       data_scale, normalization_fun, normalization_params, ...)
  
}

.as.trelliData.edata <- function(edata, edata_cname, omics_type, data_scale_original,
                                 data_scale, normalization_fun, normalization_params,
                                 is_normalized = FALSE, force_normalization = FALSE) {
  
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
  if (data_scale_original %in% log_transforms == FALSE) {stop(paste(data_scale_original, "is not an acceptable data scale."))}
  if (data_scale %in% log_transforms == FALSE) {stop(paste(data_scale, "is not an acceptable data scale."))}
  
  # Check the normalization function. Isobaric is not included since it requires f_data.
  if (normalization_fun %in% c("global", "loess", "quantile") == FALSE) {
    stop(paste(normalization_fun, "is not an acceptable normalization function type."))
  }
  
  # Build an omics data object--------------------------------------------------
  
  # Generate a f data frame 
  fdata <- data.frame("Sample" = colnames(edata), "Condition" = NA)
  fdata_cname <- "Sample"
  
  # Build the omics object 
  omicsData <- eval(parse(text = paste0("as.", omics_type, 
                                        "(e_data = edata, edata_cname = edata_cname, f_data = fdata, fdata_cname = fdata_cname,
      data_scale = data_scale_original)")))
  
  # Transform if appropriate
  if (data_scale_original != data_scale) {
    omicsData <- edata_transform(omicsData, data_scale)
  }
  
  # If the data is already normalized, skip this step 
  if (is_normalized == FALSE ) {
    
    # Normalization is not required for both isobaric protein and NMR data, but can
    # be forced with .force_normalization
    if (omics_type %in% c("isobaricpepData", "nmrData") == FALSE | force_normalization) {
      
      # Get the normalization function
      norm_fun <- switch(normalization_fun, "global" = normalize_global, 
                         "global_basic" = normalize_global_basic, "loess" = normalize_loess,
                         "quantile" = normalize_quantile)
      
      # Add omics data to normalization parameters
      normalization_params[["omicsData"]] <- omicsData
      
      # Apply normalization with its parameters
      omicsData <- do.call(norm_fun, normalization_params)
      
    }
    
  } else {
    attr(omicsData, "data_info")$norm_info$is_normalized <- TRUE
  }
  
  # Finally, generate the trelliData object-----------------------------------
  
  # Put the edata into the trelliData omics slot 
  trelliData <- list(
    trelliData.omics = omicsData$e_data %>% 
      tidyr::pivot_longer(colnames(edata)[colnames(edata) != edata_cname]) %>%
      dplyr::rename(Sample = name, Abundance = value),
    trelliData.stat = NULL,
    omicsData = omicsData,
    statRes = NULL
  )
  
  # Save group by information and set class."group_by_options" list the potential
  # inputs for the group_by function. "group_by_omics"/"group_by_stat" will hold 
  # the column name of the trelliData.omics/trelliData.stat that the data has 
  # been grouped by. And "group_by" tracks whether the group_by function has been
  # applied or not. 
  attr(trelliData, "group_by_options") <- c(edata_cname, fdata_cname)
  attr(trelliData, "group_by_omics") <- NA
  attr(trelliData, "group_by_stat") <- NA
  attr(trelliData, "group_by") <- FALSE
  class(trelliData) <- c("trelliData", "trelliData.edata")
  
  return(trelliData)
  
}

#' @name as.trelliData
#' 
#' @title Generate an object from omicsData and/or statRes objects to pass to trelliscope building functions
#' 
#' @description Either an omicsData and/or a statRes object is accepted. omicsData must
#'    be transformed and normalized, unless the data is isobaric protein or NMR data.
#'    Adding a "main_effects" group_desgination() will result in plots with groups. 
#'    The main effects group_designation and e_meta columns are merged to the e_data in long format to
#'    create the trelliData.omics dataframe, and e_meta is merged to statRes in long
#'    format to create trelliData.stat dataframe. 
#'    
#' @param omicsData an object of the class 'pepData', 'isobaricpepData', proData', 
#'    'metabData', 'lipidData', 'nmrData', created by as.pepData, as.isobaricpepData.
#'    as.proData, as.metabData, as.lipidData, as.nmrData, respectively. 
#' @param statRes statRes an object of the class 'statRes', created by \code{\link{summarize}}
#' 
#' @examples
#' \dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' # Generate an example e_meta file for lipid data 
#' lipid_emeta <- data.frame("LipidCommonName" = lipid_edata_pos$LipidCommonName, 
#'   "LipidFamily" = lipid_edata_pos$LipidCommonName %>% as.character() %>% 
#'     strsplit("(", fixed = T) %>% lapply(function(el) {el[1]}) %>% unlist())
#'     
#' # Extend fdata to have two infection groups 
#' lipid_fdata_pos2 <- data.frame("Sample_Name" = lipid_fdata_pos$Sample_Name, 
#'   "Condition" = c(rep("Mock", 3), rep("Infection_A", 4), rep("Infection_B", 4)),
#'   "Weight" = runif(11))
#' 
#' # Build lipid data object
#' lipids <- as.lipidData(e_data = lipid_edata_pos, f_data = lipid_fdata_pos2, e_meta = lipid_emeta,
#'   edata_cname = "LipidCommonName", fdata_cname = "Sample_Name", 
#'   emeta_cname = "LipidCommonName")
#' 
#' # Transform the data
#' omicsData <- edata_transform(omicsData = lipids, data_scale = "log2")
#' 
#' # Group the data by condition
#' omicsData <- group_designation(omicsData = omicsData, main_effects = c("Condition"))
#' 
#' # Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = omicsData)
#' omicsData <- applyFilt(filter_object = imdanova_Filt, omicsData = omicsData, min_nonmiss_anova=2)
#' 
#' # Normalize my pepData
#' omicsData <- normalize_global(omicsData, "subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE, "backtransform" = TRUE)
#' 
#' # Implement the IMD ANOVA method and compute all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
#' statRes <- imd_anova(omicsData = omicsData, test_method = 'anova')
#' 
#' # Generate the trelliData object 
#' trelliData2 <- as.trelliData(omicsData = omicsData)
#' trelliData3 <- as.trelliData(statRes = statRes)   
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' }
#' 
#' @author David Degnan
#' 
#' @export
as.trelliData <- function(omicsData = NULL, statRes = NULL, ...) {
  
  .as.trelliData(omicsData, statRes, require_normalization = TRUE)
  
}

.as.trelliData <- function(omicsData, statRes, require_normalization) {
  
  # Initial checks--------------------------------------------------------------
  
  # Either omicsData or statRes must be included
  if (is.null(omicsData) & is.null(statRes)) {
    stop("At least 1 omicsData or 1 statRes object must be provided.")
  }
  
  # If omicsData is provided...
  if (!is.null(omicsData)) {
    
    # ...it must be an omics data object
    if (class(omicsData) %in% c("pepData", "isobaricpepData", "proData", "metabData", "lipidData", "nmrData") == FALSE) {
      stop(paste(class(omicsData), "is not a supported omicsData class."))
    }
    
    # ...it must be log transformed if it's not NMR or isobaric
    if (class(omicsData) %in% c("nmrData", "isobaricpepData") == FALSE & 
        get_data_scale(omicsData) %in% c("log2", "log", "log10") == FALSE) {
      stop("omicsData must be log transformed.")
    }
    
    # ...it must be normalized if it's not NMR nor isobaric 
    if (require_normalization & class(omicsData) %in% c("isobaricpepData", "nmrData") == FALSE) {
      if (!get_data_norm(omicsData)) {
        stop("omicsData must be normalized.")
      }
    }
    
    # ...and if group_designation is set, the first main effect must not have any singletons
    if (!is.null(attributes(omicsData)$group_DF)) {
      group_counts <- attributes(omicsData)$group_DF$Group %>% table(dnn = "Group") %>% data.frame()
      if (1 %in% group_counts$Freq) {
        warning(paste("singleton groups found in group_designation:", 
                      paste0(group_counts[group_counts$Freq == 1, "Group"], collapse = ", ")))
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
    if (!all(statRes[[edata_cname]] %in% omicsData$e_data[[edata_cname]])) {
      stop("omicsData and statRes are from different datasets.")
    }
    
  }
  
  # Generate the trelliData object----------------------------------------------
  
  # Create placeholders for trelliData objects
  trelliData.omics <- NULL
  trelliData.stat <- NULL
  
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
      trelliData.omics <- merge(trelliData.omics, attributes(omicsData)$group_DF, by = fdata_cname)
    }
    
    # Add emeta columns if emeta exists
    if (!is.null(omicsData$e_meta)) {
      
      # Pull emeta and emeta_colname
      emeta_cname <- pmartR::get_emeta_cname(omicsData)
      emeta <- omicsData$e_meta
      
      # Add emeta columns
      trelliData.omics <- merge(trelliData.omics, emeta, by = emeta_cname)
      
    }
    
  } else {omicsData <- NULL}
  
  # Format statRes if applicable
  if (!is.null(statRes)) {
    
    # Get edata cname
    edata_cname <- pmartR::get_edata_cname(statRes)
    
    # Get column names of all fold changes, as well as p-values
    pvalue_cols <- colnames(statRes)[grepl("P_value", colnames(statRes))]
    fold_change_cols <- colnames(statRes)[grepl("Fold_change", colnames(statRes))]
    
    # Pivot longer so that the first column is the edata_cname, extract comparison,
    # group_by comparison, nest dataframes, and then extract the p_value and fold_change
    # for each group
    trelliData.stat <- statRes %>%
      dplyr::select(c(edata_cname, pvalue_cols, fold_change_cols))  %>%
      tidyr::pivot_longer(c(pvalue_cols, fold_change_cols)) %>%
      dplyr::mutate(
        Comparison = gsub("P_value_A_|Fold_change_", "", name),
        name = lapply(1:length(name), function(el) {
          gsub(paste0("_", Comparison[el]), "", name[el])
        }) %>% unlist()
      ) %>%
      dplyr::group_by(dplyr::across(c(Comparison, !!rlang::sym(edata_cname)))) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        "p_value" = purrr::map(data, function(x) {unlist(x[x$name == "P_value_A", "value"])}) %>% unlist(),
        "fold_change" = purrr::map(data, function(x) {unlist(x[x$name == "Fold_change", "value"])}) %>% unlist()
      ) %>%
      dplyr::select(-data)
    
    # Add emeta columns if emeta exists
    if (!is.null(omicsData$e_meta)) {
      
      # Add emeta columns
      trelliData.stat <- merge(trelliData.stat, emeta, by = emeta_cname)
      
    }
    
  } else {statRes <- NULL}
  
  # Generate a trelliData object
  trelliData <- list(
    trelliData.omics = trelliData.omics,
    trelliData.stat = trelliData.stat,
    omicsData = omicsData,
    statRes = statRes
  )
  
  # Save group by information and set class."group_by_options" list the potential
  # inputs for the group_by function. "group_by_omics"/"group_by_stat" will hold 
  # the column name of the trelliData.omics/trelliData.stat that the data has 
  # been grouped by. And "group_by" tracks whether the group_by function has been
  # applied or not. 
  group_options <- c(colnames(trelliData.omics), colnames(trelliData.stat)) %>% unique()
  group_nonoptions <- c("Abundance", "Comparison", "P_value", "Fold_change", "Group")
  group_options <- group_options[group_options %in% group_nonoptions == FALSE]
  attr(trelliData, "group_by_options") <- group_options
  attr(trelliData, "group_by_omics") <- NA
  attr(trelliData, "group_by_stat") <- NA
  attr(trelliData, "group_by") <- FALSE
  class(trelliData) <- c("trelliData")
  
  return(trelliData)
  
}

#' @name trelli_group_by
#' 
#' @title Set the "group_by" variable for a trelliData object 
#'    
#' @description Allows for grouping omics or stats data for downstream plotting 
#'     and cognostic functions
#'     
#' @param trelliData A trelliscope data object made by as.trelliData or as.trelliData.edata. Required.
#' @param group The name of a column in trelliData to group the data by. Required. 
#' 
#' @examples
#' \dontrun{
#' 
#' library(pmartRdata)
#' library(pmartR)
#' 
#' ## "group_by" with an edata file. Generate with example code in as.trelliData.edata
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName")
#' trelli_group_by(trelliData = trelliData, group = "Sample")
#' 
#' ## "group_by" with trelliData containing omicsData. Generate with example code in as.trelliData
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName")
#' trelli_group_by(trelliData = trelliData2, group = "LipidFamily")
#' 
#' ## "group_by" with trelliData containing statRes. Generate with example code in as.trelliData
#' trelli_group_by(trelliData = trelliData3, group = "LipidCommonName")
#' 
#' ## "group_by" with trelliData containing both omicsData and statRes. Generate with example code in as.trelliData
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily")
#' trelli_group_by(trelliData = trelliData4, group = "LipidCommonName")
#' trelli_group_by(trelliData = trelliData4, group = "Sample_Name")
#' 
#' }  
#' 
#' @author David Degnan
#' 
#' @export
trelli_group_by <- function(trelliData, group) {
  
  # Run initial checks----------------------------------------------------------
  
  # Confirm that trelliData is a trelliData object
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData from as.trelliData or as.trelliData.edata.")
  }
  
  # Confirm that group is an acceptable option
  if (group %in% attr(trelliData, "group_by_options") == FALSE) {
    stop(paste("groups is not an acceptable option. The following can be selected:", 
               paste(attr(trelliData, "group_by_options"), collapse = ", ")))
  }
  
  # Confirm that group_by is false
  if (attr(trelliData, "group_by")) {
    
    # If group_by has already been conducted with this group, write a message and return trelliData
    if (attr(trelliData, "group_by") == group) {
      stop(paste("trelliData has already been grouped by", group))
    } 
    
  }
  
  # Determine which dataframes this grouping variable applies to----------------
  
  # Test if grouping applies to omicsData 
  if (!is.null(trelliData$trelliData.omics) && group %in% colnames(trelliData$trelliData.omics)) {
    apply_to_omics <- TRUE} else {apply_to_omics <- FALSE}
  
  # Test if grouping applies to statRes 
  if (!is.null(trelliData$trelliData.stat) && group %in% colnames(trelliData$trelliData.stat)) {
    apply_to_stat <- TRUE} else {apply_to_stat <- FALSE}
  
  # Test group_by option--------------------------------------------------------
  
  # If there are instances of groups with less than 3 data points, it is not a good grouping option. 
  .test_grouping <- function(trelliData_subclass, trelliData_subclass_name) {
    
    # Get the smallest group size
    smallest_group_size <- trelliData_subclass %>% 
      dplyr::group_by_at(group) %>% 
      dplyr::summarize(N = dplyr::n()) %>%
      dplyr::select(N) %>%
      min()
    
    # If the smallest group size is less than 3, then don't proceed 
    if (smallest_group_size < 3) {
      warning(paste0("Grouping by ", group, " results in groups with less than 3 ",
                     "data points in ", trelliData_subclass_name, "."))
    }
    
  }
  
  # Onlt test grouping if it applies to omics
  if (apply_to_omics) {.test_grouping(trelliData$trelliData.omics, "trelliData.omics")}
  
  # Group and nest samples------------------------------------------------------
  
  .group_samples <- function(trelliData_subclass) {
    
    trelliData_subclass %>%
      dplyr::group_by_at(group) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::rename(Nested_DF = data)
  }
  
  # Group omicsData
  if (apply_to_omics) {
    trelliData$trelliData.omics <- .group_samples(trelliData$trelliData.omics)
    attr(trelliData, "group_by_omics") <- group
  }
  
  # Group statRes
  if (apply_to_stat) {
    trelliData$trelliData.stat <- .group_samples(trelliData$trelliData.stat)
    attr(trelliData, "group_by_stat") <- group
  }
  
  # Export results--------------------------------------------------------------
  attr(trelliData, "group_by") <- TRUE
  return(trelliData)
  
}

#' @name trelli_abundance_boxplot
#' 
#' @title Boxplot trelliscope building function for abundance data   
#' 
#' @description Specify a boxplot design and cognostics for the abundance boxplot trelliscope.
#'    Each boxplot will have its own groups as specified by the first main effect in group_designation.
#' 
#' @param trelliData A trelliscope data object made by as.trelliData or as.trelliData.edata,
#'    and grouped by trelli_group_by. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    n, mean, median, sd, and skew for abundance. If statRes data is included, 
#'    p_value and fold_change cognostics can be added. If no cognostics are desired,
#'    set to NULL. 
#' @param ggplot_params An optional string of ggplot parameters to the backend ggplot
#'    function. For example "ylab('Group') + xlab('Abundance') + ylim(c(2,20))".
#'    Default is NULL. 
#' @param interactive A logical argument indicating whether the plots should be interactive
#'    or not. Interactive plots are ggplots piped to ggplotly (for now). Default is FALSE.  
#' @param path The base directory of the trelliscope application. Default is Downloads. 
#' @param name The name of the display. Default is Trelliscope.
#' @param test_mode A logical to return a smaller trelliscope to confirm plot and design.
#'    Default is FALSE.
#' @param test_example The index number of the plot to return for test_mode. Default is 1. 
#' 
#'    
#' @examples
#' \dontrun{
#' 
#' ## Build the abundance boxplot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName") %>% 
#'    trelli_abundance_boxplot(test_mode = T, test_example = 1:10)
#' trelli_group_by(trelliData = trelliData, group = "Sample") %>% trelli_abundance_boxplot()
#' 
#' ## Build the abundance boxplot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'    trelli_abundance_boxplot(test_mode = T, test_example = 1:10)
#' trelli_group_by(trelliData = trelliData2, group = "LipidFamily") %>% trelli_abundance_boxplot()
#'     
#' ## Build the abundance boxplot with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_group_by(trelliData = trelliData4, group = "LipidCommonName") %>%
#'    trelli_abundance_boxplot(test_mode = T, test_example = 1:10)
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily") %>% trelli_abundance_boxplot()
#'    
#' 
#' }
#' 
#' @author David Degnan
#' 
#' @export
trelli_abundance_boxplot <- function(trelliData,
                           cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
                           ggplot_params = NULL,
                           interactive = FALSE,
                           path = "~/Downloads/Trelliscope",
                           name = "Trelliscope",
                           test_mode = FALSE,
                           test_example = 1,
                           ...) {
  
  # Run initial checks----------------------------------------------------------
  
  # trelliData object must be of the trelliData class
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData from as.trelliData or as.trelliData.edata.")
  }
  
  # Check that trelliData has been passed to the "trelli_group_by" function
  if (!attr(trelliData, "group_by")) {
    stop("trelliData must be grouped with trelli_group_by.")
  }
  
  # Assert that trelliData has omicsData
  if (is.null(trelliData$trelliData.omics)) {
    stop("trelliData must have omicsData for this abundance boxplot function.")
  }
  
  # If cognostics are not NULL...
  if (!is.null(cognostics)) {
    
    # ...assert that cognostics are acceptable options 
    cog_options <- c("n", "mean", "median", "sd", "skew", "p_value", "fold_change")
    if (all(cognostics %in% cog_options) == FALSE) {
      stop(paste("Unacceptable cognostic option included. Acceptable options are: ", 
                 paste(cog_options, collapse = ", ")))
    }
    
    # ...remove stat specific options 
    if (any(c("p_value", "fold_change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[-match(c("p_value", "fold_change"), cognostics, nomatch = 0)]
    }
    
  }
  
  # If no cognostics, set to NULL
  if (length(cognostics) == 0) {cognostics <- NULL}
  
  # If ggplot_params is not null...
  if (!is.null(ggplot_params)) {
    
    # ...and ggplot_params is not a character or is a vector of greater than 1, 
    if (!is.character(ggplot_params) | length(ggplot_params) > 1) {
      stop("ggplot_params must be a single character string or NULL.")
    }
    
  }
  
  # If interactive is not TRUE/FALSE, inform the user
  if (!is.logical(interactive)) {
    stop("interactive must be a true or false.")
  }
  if (is.na(interactive)) {interactive <- FALSE}
  
  # test_mode must be a TRUE/FALSE
  if (!is.logical(test_mode)) {
    stop("test_mode must be a true or false")
  }
  if (is.na(test_mode)) {test_mode <- FALSE}
  
  # Ensure that test_example is an integer
  if (!is.numeric(test_example) | 0 %in% test_example) {
    "test_example should be non-zero integers."
  }
  test_example <- unique(abs(round(test_example)))
  
  # Ensure that test_example is in the range of possibilities 
  if (max(test_example) > nrow(trelliData$trelliData.omics)) {
    stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.omics)))
  }
  
  # Make boxplot function-------------------------------------------------------
  
  # First, generate the boxplot function
  box_plot_fun <- function(DF, title) {
    
    # Add a blank group if no group designation was given
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {DF$Group <- "x"} 
    
    # Build plot 
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = Abundance)) + 
      ggplot2::geom_boxplot() + ggplot2::geom_point() + ggplot2::theme_bw() + 
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) + 
      ggplot2::ylab(paste(attr(trelliData$omicsData, "data_info")$data_scale, "Abundance"))
    
    # Remove x axis if no groups
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      boxplot <- boxplot + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                          axis.ticks.x = ggplot2::element_blank(), 
                                          axis.text.x = ggplot2::element_blank())
    }
    
    # Add additional parameters
    if (!is.null(ggplot_params)) {
      boxplot <- boxplot + unquote(ggplot_params)
    }
    
    return(boxplot)
  }
  
  # Create cognostic function---------------------------------------------------
  
  # Second, create function to return cognostics
  box_cog_fun <- function(DF, group) {
    
    # Set basic cognostics for ungrouped data or in case when data is not split by fdata_cname
    cog <- list(
     "n" = dplyr::tibble(`Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Biomolecule Count")),
     "mean" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = T), 4), desc = "Mean Abundance")), 
     "median" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = T), 4), desc = "Median Abundance")), 
     "sd" = dplyr::tibble(`Standard Deviation Abundance` = trelliscopejs::cog(round(sd(DF$Abundance, na.rm = T), 4), desc = "Abundance Standard Deviation")), 
     "skew" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = T), 4), desc= "Abundance Skewness"))
    )
    
    # Start list of cogs
    cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% tibble::tibble()
    
    # Get fdata cname and group_by selection
    fdata_cname <- pmartR::get_fdata_cname(trelliData$omicsData)
    group_by_choice <- attr(trelliData, "group_by_omics")
    
    # Additional group cognostics can be added only if group_designation was set and
    # trelli_group_by is not the fdata_cname
    if (!is.null(attributes(trelliData$omicsData)$group_DF) & fdata_cname != group_by_choice) {
      
      # A quick cognostic function 
      quick_cog <- function(name, value) {
        dplyr::tibble(!!rlang::sym(name) := trelliscopejs::cog(value, desc = name))
      }
      
      # Create a list to convert from short name to long
      name_converter <- list("n" = "Count", "mean" = "Mean Abundance", 
        "median" = "Median Abundance", "sd" = "Standard Deviation Abundance", 
        "skew" = "Skew Abundance", "p_value" = "P Value", "fold_change" = "Fold Change")
      
      # Since the number of groups is unknown, first group_by the Groups,
      # then calculate all summary statistics, pivot to long format,
      # subset down to requested statistics, switch name to a more specific name, 
      # combine group and name, and generate the cognostic tibble
      cogs_to_add <- DF %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(
          "n" = sum(!is.na(Abundance)), 
          "mean" = round(mean(Abundance, na.rm = T), 4),
          "median" = round(median(Abundance, na.rm = T), 4),
          "sd" = round(sd(Abundance, na.rm = T), 4),
          "skew" = round(e1071::skewness(Abundance, na.rm = T), 4)
        ) %>%
        tidyr::pivot_longer(c(n, mean, median, sd, skew)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Group, lapply(name, function(x) {name_converter[[x]]}) %>% unlist())
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Group) 
      
      # Add new cognostics 
      cog_to_trelli <- cbind(cog_to_trelli, do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
        quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
      })) %>% tibble::tibble()) %>% tibble::tibble()
        
    }
    
    # Add cognostics that only apply when there is stats data and it is the 
    # same column as omicsData
    if (!is.null(trelliData$trelliData.stat) && attr(trelliData, "group_by_omics") == attr(trelliData, "group_by_stat")) {
      
      # Get edata cname
      edata_cname <- pmartR::get_edata_cname(trelliData$statRes)
      
      # Subset down the dataframe down to group, unnest the dataframe, 
      # pivot_longer to comparison, subset columns to requested statistics, 
      # switch name to a more specific name
      cogs_to_add <- trelliData$trelliData.stat %>%
        dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == group) %>%
        dplyr::select(Nested_DF) %>%
        tidyr::unnest(cols = c(Nested_DF)) %>%
        dplyr::select(c(Comparison, p_value, fold_change)) %>%
        tidyr::pivot_longer(c(p_value, fold_change)) %>%
        dplyr::mutate(
          name = paste(Comparison, lapply(name, function(x) {name_converter[[x]]}) %>% unlist()),
          value = round(value, 4)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Comparison)
      
      # Add new cognostics 
      cog_to_trelli <- cbind(cog_to_trelli, do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
        quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
      })) %>% tibble::tibble()) %>% tibble::tibble()
      
    }

    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # If test_mode is on, then just build the required panels
  if (test_mode) {toBuild <- trelliData$trelliData.omics[test_example,]} 
  else {toBuild <- trelliData$trelliData.omics}
  
  # Build trelliscope without cognostics if none are provided. Otherwise, build with cognostics.
  if (is.null(cognostics)) {
    
    toBuild %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), box_plot_fun)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = T, ...)
    
  } else {
    
    toBuild %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), box_plot_fun),
        cog = trelliscopejs::map2_cog(Nested_DF, as.character(unlist(toBuild[,1])), box_cog_fun)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = T, ...)

  }

}