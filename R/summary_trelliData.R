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

  # Second, let's determine if there's omicsData or statRes or both in this
  # object
  omics <- !is.null(trelliData$trelliData.omics)
  stat <- !is.null(trelliData$trelliData.stat)

  # Third, let's determine if this object is a trelliData.edata
  edata_only <- inherits(trelliData, "trelliData.edata")

  # Get edata_cname
  if (omics) {
    edata_cname <- get_edata_cname(trelliData$omicsData)
  } else {
    edata_cname <- get_edata_cname(trelliData$statRes)
  }

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
      rep("e_meta column", 6)
    ),
    Plot = c(
      "abundance boxplot", "abundance histogram", "missingness bar",
      "fold change bar", "abundance boxplot", "missingness bar",
      "abundance boxplot", "abundance heatmap", "missingness bar",
      "fold change boxplot", "fold change heatmap", "fold change volcano"
    ),
    `Data Type` = c(
      "omics", "omics", NA, "stat", "omics", NA, "omics", "omics",
      NA, "stat", "stat", "stat"
    )
  )
  
  # Update plot names when data is seqData
  if (inherits(trelliData, "trelliData.seqData")) {
    
    All_Options <- All_Options %>%
      dplyr::mutate(
        Plot = gsub("abundance", "rnaseq", Plot),
        Plot = gsub("missingness", "rnaseq nonzero", Plot)
      )
    
  }

  #################################
  ## SUBSET AND RETURN DATAFRAME ##
  #################################

  # Filter by "Data Type"
  if (omics == FALSE) {
    All_Options <- All_Options %>% 
      dplyr::filter(`Data Type` != "omics" | is.na(`Data Type`))
  }
  if (stat == FALSE) {
    All_Options <- All_Options %>% 
      dplyr::filter(`Data Type` != "stat" | is.na(`Data Type`))
  }

  # If there is no grouping information, we should suggest potential plots.
  if (panel_by == FALSE) {
    # Filter by "Panel By" choices
    if (fdata_cname_missing) {
      All_Options <- All_Options %>% 
        dplyr::filter(`Panel By Choice` != "f_data cname")
    }
    if (emeta_cols_missing) {
      All_Options <- All_Options %>% 
        dplyr::filter(`Panel By Choice` != "e_meta column")
    }

    # Remove data type and add a holder for count
    All_Options <- All_Options %>%
      dplyr::select(-`Data Type`) %>%
      dplyr::mutate(`Number of Plots` = 0)

    # Replace names and get counts TODO: apply styler
    if ("e_data cname" %in% All_Options$`Panel By Choice`) {
      bio_var <- ifelse(omics, attr(trelliData$omicsData, "cnames")$edata_cname, attr(trelliData$statRes, "cnames")$edata_cname)
      bio_count <- ifelse(omics, trelliData$trelliData.omics[[bio_var]] %>% unique() %>% length(), trelliData$trelliData.stat[[bio_var]] %>% unique() %>% length())
      All_Options[All_Options$`Panel By Choice` == "e_data cname", "Number of Plots"] <- bio_count %>% as.character()
      All_Options[All_Options$`Panel By Choice` == "e_data cname", "Panel By Choice"] <- edata_cname
    }

    if ("f_data cname" %in% All_Options$`Panel By Choice`) {
      sample_var <- ifelse(omics, attr(trelliData$omicsData, "cnames")$fdata_cname, attr(trelliData$statRes, "cnames")$fdata_cname)
      sample_count <- ifelse(omics, trelliData$trelliData.omics[[sample_var]] %>% unique() %>% length(), trelliData$trelliData.stat[[sample_var]] %>% unique() %>% length())
      All_Options[All_Options$`Panel By Choice` == "f_data cname", "Number of Plots"] <- sample_count %>% as.character()
      All_Options[All_Options$`Panel By Choice` == "f_data cname", "Panel By Choice"] <- fdata_cname
    }

    if ("e_meta column" %in% All_Options$`Panel By Choice`) {
      # Get counts per e_meta variable
      emeta_counts <- lapply(emeta_cols, function(name) {
        trelliData$trelliData.omics[[name]] %>% unique() %>% length()
      }) %>% 
        unlist() %>% 
        paste(collapse = ", ")

      # Add counts
      All_Options[
        All_Options$`Panel By Choice` == "e_meta column", "Number of Plots"
      ] <- emeta_counts
      All_Options[
        All_Options$`Panel By Choice` == "e_meta column", "Panel By Choice"
      ] <- paste(emeta_cols, collapse = ", ")
    }
  } else {
    # Determine what the data has been grouped by
    Grouped <- ifelse(
      omics,
      attr(trelliData, "panel_by_omics"),
      attr(trelliData, "panel_by_stat")
    )

    # Determine if group variable is edata_cname, fdata_cname, or an emeta
    # column, and get a count
    theMatch <- match(Grouped, c(edata_cname, fdata_cname, emeta_cols))
    if (theMatch == 1) {
      panel_by_choice <- "e_data cname"
      count <- ifelse(
        omics, 
        nrow(trelliData$trelliData.omics), 
        nrow(trelliData$trelliData.stat)
      )
    } else if (theMatch == 2) {
      panel_by_choice <- "f_data cname"
      count <- ifelse(
        omics,
        nrow(trelliData$omicsData$f_data),
        attr(trelliData$statRes, "group_DF") %>% nrow()
      )
    } else {
      panel_by_choice <- "e_meta column"
      count <- trelliData$omicsData$e_meta[[Grouped]] %>%
        unique() %>%
        length()
    }

    # Finally, subset down the table, remove data type, and add a count
    All_Options <- All_Options %>%
      dplyr::filter(`Panel By Choice` == panel_by_choice) %>%
      dplyr::select(-`Data Type`) %>%
      dplyr::mutate(`Number of Plots` = count %>% as.numeric()) %>%
      dplyr::mutate(`Panel By Choice` = Grouped)
  }

  return(All_Options)
}
