#' Summarizes potential plotting options for a trelliData object
#'
#' @param object An object from the as.trelliData.edata or as.trelliData
#'   functions
#' @param ... further arguments passed to or from other methods.
#'
#' @return A data.frame containing panel plot options for this trelliData object.
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#'
#' library(dplyr)
#' library(pmartRdata)
#' 
#' trelliData <- as.trelliData.edata(e_data = pep_edata,
#'                                    edata_cname = "Peptide",
#'                                    omics_type = "pepData")
#' 
#' # Use an edata example. Build with as.trelliData.edata.
#' summary(trelliData)
#' summary(trelliData %>% trelli_panel_by("Peptide"))
#' summary(trelliData %>% trelli_panel_by("Sample"))
#'
#' }
#'
#' @export
#' @rdname summary-trelliData
#' @name summary-trelliData
summary.trelliData <- function(object, ...) {
  trelliData <- object

  #######################################
  ## EXTRACT ATTRIBUTES AND PROPERTIES ##
  #######################################

  # First, we must know if this data has been grouped at all
  panel_by <- attr(trelliData, "panel_by")
  panel_by_col <- attr(trelliData, "panel_by_col")

  # Second, let's determine if there's omicsData or statRes or both in this
  # object
  omics <- !is.null(trelliData$omicsData) 
  stat <- !is.null(trelliData$statRes)

  # Get edata_cname
  edata_cname <- attr(trelliData, "edata_col")

  # Extract fdata_cname
  fdata_cname <- attr(trelliData, "fdata_col")
  fdata_cname_missing <- is.null(fdata_cname)

  # Extract emeta colnames
  emeta_cols <- attr(trelliData, "emeta_col")
  emeta_cols_missing <- is.null(emeta_cols)

  #####################
  ## BUILD DATAFRAME ##
  #####################
  
  # Create a base data.frame which can be filtered
  All_Options <- data.table::data.table(
    `Panel By Choice` = c(
      rep("e_data cname", 4),
      rep("f_data cname", 2),
      rep("e_meta column", 5)
    ),
    Plot = c(
      "abundance boxplot", "abundance histogram", "missingness bar",
      "fold change bar", "abundance boxplot", "missingness bar",
      "abundance heatmap", "missingness bar",
      "fold change boxplot", "fold change heatmap", "fold change volcano"
    ),
    `Data Type` = c(
      "omics", "omics", NA, "stat", "omics", NA, "omics",
      NA, "stat", "stat", "stat"
    )
  )
  
  # Update plot names when data is seqData
  if (inherits(trelliData, "trelliData.seqData")) {
    
    All_Options <- All_Options %>%
      dplyr::mutate(
        Plot = gsub("abundance", "rnaseq", .data$Plot),
        Plot = gsub("missingness", "rnaseq nonzero", .data$Plot)
      )
  }
  
  #######################
  ## FILTER DATA.FRAME ##
  #######################
  
  # Filter by "Data Type"
  if (omics == FALSE) {
    All_Options <- All_Options %>% 
      dplyr::filter(`Data Type` != "omics" | is.na(`Data Type`))
  }
  if (stat == FALSE) {
    All_Options <- All_Options %>% 
      dplyr::filter(`Data Type` != "stat" | is.na(`Data Type`))
  }
  
  # Filter by Panel choice 
  if (!panel_by) {
    if (fdata_cname_missing) {
      All_Options <- All_Options %>% 
        dplyr::filter(`Panel By Choice` != "f_data cname")
    }
    if (emeta_cols_missing) {
      All_Options <- All_Options %>% 
        dplyr::filter(`Panel By Choice` != "e_meta column")
    }
  } else {
    
    if (panel_by_col == edata_cname) {
      All_Options <- All_Options %>% dplyr::filter(`Panel By Choice` == "e_data cname")
    }
    if (!is.null(fdata_cname) && panel_by_col == fdata_cname) {
      All_Options <- All_Options %>% dplyr::filter(`Panel By Choice` == "f_data cname")
    }
    if (!is.null(emeta_cols) && panel_by_col %in% emeta_cols) {
      All_Options <- All_Options %>% dplyr::filter(`Panel By Choice` == "e_meta column")
    }
    
  }

  ######################
  ## CALCULATE COUNTS ##
  ######################

  # Remove data type and add a holder for count
  All_Options <- All_Options %>%
    dplyr::mutate(`Number of Plots` = 0)
  
  # Replace names and get counts TODO: apply styler
  if ("e_data cname" %in% All_Options$`Panel By Choice`) {
    bio_count <- trelliData$trelliData[[edata_cname]] %>% unique() %>% length()
    All_Options[All_Options$`Panel By Choice` == "e_data cname", "Number of Plots"] <- bio_count %>% as.character()
    All_Options[All_Options$`Panel By Choice` == "e_data cname", "Panel By Choice"] <- edata_cname
  }

  if ("f_data cname" %in% All_Options$`Panel By Choice`) {
    sample_count <- trelliData$trelliData[[fdata_cname]] %>% unique() %>% length()
    All_Options[All_Options$`Panel By Choice` == "f_data cname", "Number of Plots"] <- sample_count %>% as.character()
    All_Options[All_Options$`Panel By Choice` == "f_data cname", "Panel By Choice"] <- fdata_cname
  }

  if ("e_meta column" %in% All_Options$`Panel By Choice`) {

    # Add counts
    if (panel_by && attr(trelliData, "panel_by_col") %in% emeta_cols) {
      
      emeta_var <- attr(trelliData, "panel_by_col")
      emeta_count <- trelliData$trelliData[[emeta_var]] %>% unique() %>% length()
      All_Options[All_Options$`Panel By Choice` == "e_meta column", "Number of Plots"] <- emeta_count %>% as.character()
      All_Options[All_Options$`Panel By Choice` == "e_meta column", "Panel By Choice"] <- emeta_var
    
    } else {
      
      # Get counts per e_meta variable
      emeta_counts <- lapply(emeta_cols, function(name) {
        trelliData$trelliData[[name]] %>% unique() %>% length()
      }) %>% 
        unlist() %>% 
        paste(collapse = ", ")
      
      
      All_Options[
        All_Options$`Panel By Choice` == "e_meta column", "Number of Plots"
      ] <- emeta_counts
      All_Options[
        All_Options$`Panel By Choice` == "e_meta column", "Panel By Choice"
      ] <- paste(emeta_cols, collapse = ", ")
    }
  }
  
  # Remove data type column
  All_Options <- All_Options %>% dplyr::select(-`Data Type`)
  
  # If paneled, filter Panel By Choice
  if (panel_by) {
    All_Options <- All_Options %>% dplyr::filter(`Panel By Choice` == panel_by_col)
  }

  return(All_Options)
}
