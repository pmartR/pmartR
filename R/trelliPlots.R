#' Performs initial checks for trelliData objects
#'
#' This function runs necessary checks for pmartR trelliscope plotting
#' functions. It cleans any parameters (rounding numerics to integers, etc.),
#' and returns them.
#'
#' @param trelliData trelliData object the user passed to a plotting
#'   function
#' @param trelliCheck Check if the object type is supposed to be "omics",
#'   "statRes" or put a vector of both
#' @param cognostics A vector of the user provided cognstics
#' @param acceptable_cognostics The acceptable cognostics for this plot
#' @param ggplot_params The vector of user provided ggplots
#' @param interactive The user provided logical for whether the plot should be
#'   interactive
#' @param build_trelliscope The user provided value for building trelliscope
#' @param test_mode The user provided logical for whether a smaller trelliscope
#'   should be returned
#' @param test_example The user provided vector of plot indices
#' @param single_plot The user provided logical for whether a single plot should
#'   be returned
#' @param seqDataCheck Whether seqData is permitted for this plot. "no" means that 
#'   seqData cannot be used at all, "permissible" means that seqData can be used, 
#'   and "required" means that seqData is required for the plotting function.
#' @param seqText Text that should appear when seqDataCheck is violated. 
#' @param p_value_skip Whether to skip specific p_value checks. 
#' @param p_value_thresh The user provided threshold for plotting significant
#'   p-values.
#' 
#' @return No return value, validates a trelliData object before passing it to builder functions.
#' 
trelli_precheck <- function(trelliData, 
                            trelliCheck,
                            cognostics,
                            acceptable_cognostics,
                            ggplot_params,
                            interactive,
                            build_trelliscope, 
                            test_mode,
                            test_example,
                            single_plot,
                            seqDataCheck,
                            seqText = NULL,
                            p_value_skip = FALSE,
                            p_value_thresh = NULL) {
  
  ####################
  ## SEQDATA CHECKS ##
  ####################
  
  seq_test <- inherits(trelliData, "trelliData.seqData")
  if (seqDataCheck == "no" & seq_test) {stop(paste("seqData is not permitted for this plotting function.", seqText))}
  if (seqDataCheck == "required" & !seq_test) {stop(paste("seqData is required for this plotting function.", seqText))}
  
  #########################
  ## TEST EXAMPLE CHECKS ##
  #########################

  # Assert that test mode is a true or false
  if (!is.logical(test_mode) & !is.na(test_mode)) {
    stop("test_mode must be a TRUE or FALSE")
  }

  # Ensure that test_example is an integer
  if (test_mode) {
    if (!is.numeric(test_example) | 0 %in% test_example) {
      stop("test_example should be a non-zero integer.")
    }
    test_example <- unique(abs(round(test_example)))
  }

  #######################
  ## trelliData checks ##
  #######################

  # trelliData object must be of the trelliData class
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData.")
  }

  # Check that trelliData has been passed to the "trelli_panel_by" function.
  if (!attr(trelliData, "panel_by")) {
    stop("trelliData must be paneled with trelli_panel_by.")
  }

  # Check that omics data exists
  if ("omics" %in% trelliCheck) {
    # Assert that trelliData has omicsData
    if (is.null(trelliData$omicsData)) {
      stop("trelliData must have omicsData for this plotting function.")
    }

    # Ensure that test_example is in the range of possibilities
    if (test_mode) {
      if (max(test_example) > nrow(trelliData$trelliData)) {
        stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData)))
      }
    }
  }

  # Check that statRes data exists
  if ("stat" %in% trelliCheck) {
    # Assert that trelliData has statRes
    if (is.null(trelliData$statRes)) {
      stop("trelliData must have statRes for this plotting function.")
    }

    # Ensure that test_example is in the range of possibilities
    if (test_mode) {
      if (max(test_example) > nrow(trelliData$statRes)) {
        stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.stat)))
      }
    }
  }

  ######################
  ## COGNOSTIC CHECKS ##
  ######################
  
  # Add emeta columns in acceptable cognostics 
  if (!is.null(attr(trelliData, "emeta_col"))) {
    acceptable_cognostics <- c(acceptable_cognostics, attr(trelliData, "emeta_col"))
  }

  # If cognostics are not NULL...
  if (!is.null(cognostics)) {
    # ...assert that cognostics are acceptable options
    if (all(cognostics %in% acceptable_cognostics) == FALSE) {
      stop(paste(
        "Unacceptable cognostic option included. Acceptable options are: ",
        paste(acceptable_cognostics, collapse = ", ")
      ))
    }
  }

  #########################
  ## GGPLOT PARAMS CHECK ##
  #########################

  # If ggplot_params is not null...
  if (!is.null(ggplot_params)) {
    # ...and ggplot_params is not a character or is a vector of greater than 1,
    if (!is.character(ggplot_params)) {
      stop("ggplot_params must be a string, vector of strings, or NULL.")
    }
  }

  ####################
  ## LOGICAL CHECKS ##
  ####################

  # If interactive is not TRUE/FALSE, inform the user
  if (!is.logical(interactive) & !is.na(interactive)) {
    stop("interactive must be a TRUE or FALSE.")
  }
  
  # If build_trelliscope is not TRUE/FALSE, inform the user
  if (!is.logical(build_trelliscope) & !is.na(build_trelliscope)) {
    stop("build_trelliscope must be a TRUE or FALSE.")
  }

  # single_plot must be a TRUE/FALSE
  if (!is.logical(single_plot) & !is.na(single_plot)) {
    stop("single_plot must be a TRUE or FALSE.")
  }

  #############################
  ## P-VALUE SPECIFIC CHECKS ##
  #############################

  # Check p_value test
  if ("stat" %in% trelliCheck & p_value_skip != TRUE) {
      
    # Check p_value threshold
    if (!is.numeric(p_value_thresh)) {
      stop("p_value_thresh must be a numeric.")
    }
    if (p_value_thresh < 0 | p_value_thresh > 1) {
      stop("p_value_thresh must be between 0 and 1.")
    }
  }
  
}

# This function builds a trelliscope display from a trelliscope data.frame using lazy panel
trelli_builder_lazy <- function(toBuild, path, name, ...) {
  
  trelliscope::as_trelliscope_df(
    df = toBuild,
    name = name,
    path = path,
    force_plot = TRUE,
    ...
  ) %>%
    trelliscope::view_trelliscope()
  
}

#' @name trelli_abundance_boxplot
#'
#' @title Boxplot trelliscope building function for abundance data 
#'
#' @description Specify a boxplot design and cognostics for the abundance
#'   boxplot trelliscope. Each boxplot will have its own groups as specified by
#'   the first main effect in group_designation. Use "trelli_rnaseq_boxplot"
#'   for RNA-Seq data. 
#'
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by trelli_panel_by. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are "count", "mean abundance", "median abundance", and "cv abundance". 
#'   If data are paneled by a biomolecule, the count will be "sample count".
#'   If data are paneled by a sample or a biomolecule class, the count will be "biomolecule count". 
#'   If statRes data is included, "anova p-value" and "fold change" data per comparisons
#'   may be added. If grouping information is included, only "sample count" and 
#'   "mean abundance" will be calculated, along with "anova p-value" and "fold change"
#'   if specified. "anova p-value" will not be included if paneling a trelliscope
#'   display by a biomolecule class. Default is "sample count" and "mean abundance".
#'   Any additional cognostics from e_meta variables may be added if their column
#'   name is specified, and the data has been paneled by the edata_cname or by 
#'   an e_meta variable.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "ylim(c(2,20))").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param include_points Add points as a geom_jitter. Default is TRUE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of boxplots that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#'
#' # Build the abundance boxplot with an edata file where each panel is a biomolecule. 
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'
#' # Build the abundance boxplot wher each panel is a sample.
#' # Include all applicable cognostics. Remove points. 
#' trelli_panel_by(trelliData = trelliData1, panel = "Sample") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, 
#'                             include_points = FALSE,
#'                             cognostics = c("count", 
#'                                            "mean abundance", 
#'                                            "median abundance", 
#'                                            "cv abundance"),
#'                              path = tempdir()
#'                            )
#'
#' # Build the abundance boxplot with an omicsData object.
#' # Let the panels be biomolecules. Here, grouping information is included.
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'
#' # Build the abundance boxplot with an omicsData object. 
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'
#' # Build the abundance boxplot with an omicsData and statRes object.
#' # Panel by a biomolecule, and add statistics data to the cognostics
#' trelli_panel_by(trelliData = trelliData4, panel = "Peptide") %>%
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir(),
#'                             cognostics = c("mean abundance", "anova p-value", "fold change"))
#'
#' # Other options include modifying the ggplot  
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir(),
#'      ggplot_params = c("ylab('')", "ylim(c(20,30))"))
#'
#' # Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData4, panel = "Peptide") %>% 
#'     trelli_abundance_boxplot(
#'      interactive = TRUE, test_mode = TRUE, test_example = 1:10, path = tempdir())
#' 
#' \dontshow{closeAllConnections()}
#' }
#' }
#' 
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_abundance_boxplot <- function(trelliData,
                                     cognostics = c("count", "mean abundance"),
                                     ggplot_params = NULL,
                                     interactive = FALSE,
                                     include_points = TRUE,
                                     path = tempdir(),
                                     name = "Trelliscope",
                                     build_trelliscope = TRUE,
                                     test_mode = FALSE,
                                     test_example = 1,
                                     single_plot = FALSE,
                                     ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("count", "mean abundance", "median abundance", "cv abundance", "anova p-value", "fold change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_boxplot instead.",
                  p_value_thresh = NULL)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Paneling by e_meta col is not allowed
  if (panel %in% attr(trelliData, "emeta_col")) {
    stop("Paneling by an e_meta column is not allowed for trelli_abundance_boxplot. Try `trelli_abundance_heatmap`")
  }
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$statRes)) {
    if (any(c("anova p-value", "fold change") %in% cognostics) & is.null(trelliData$statRes)) {
      cognostics <- cognostics[-match(c("anova p-value", "fold change"), cognostics, nomatch = 0)]
      message(paste("'anova p-value' and/or 'fold change' were listed as cognostics, but not provided in the trelliData object.",
                    "Did you forget to include a statRes object?")
             )
    }    
  }
  
  # Remove median and cv as cognostics if data is grouped
  if (!inherits(trelliData, "trelliData.edata") & panel == attr(trelliData, "edata_col")) {
    if (any(c("median abundance", "cv abundance") %in% cognostics)) {
      cognostics <- cognostics[-match(c("median abundance", "cv abundance"), cognostics, nomatch = 0)]
      message("'median abundance' and 'cv abundance' are not permitted when groups have been specified.")
    }
  }
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Make sure include_points is a true or false
  if (!is.logical(include_points) & !is.na(include_points)) {
    stop("include_points must be a TRUE or FALSE.")
  }
  
  if (any(c("anova p-value", "fold change") %in% cognostics) & panel != get_edata_cname(trelliData$omicsData)) {
    message(paste("Please panel by", get_edata_cname(trelliData$omicsData), "to get 'anova p-value' and 'fold change' as cognostics in the trelliscope display."))
    cognostics <- cognostics[cognostics %in% c("anova p-value", "fold change") == FALSE]
  }

  # Start summary toBuild data.frame
  toBuild <- trelliData$trelliData
  
  # First, add any missing metas------------------------------------------------
  
  if (!single_plot) {
    
    # Add metas without groups
    if ("Group" %in% colnames(toBuild) == FALSE) {
      
      toBuild <- toBuild %>% 
        dplyr::group_by_at(panel) %>%
        dplyr::summarize(
          count = sum(!is.na(Abundance)),
          `mean abundance` = mean(Abundance, na.rm = T),
          `median abundance` = median(Abundance, na.rm = T),
          `cv abundance` = sd(Abundance, na.rm = T) / `mean abundance` * 100
        )
      
      # Now, select only the requested metas
      toBuild <- toBuild %>%
        dplyr::select(dplyr::all_of(c(panel, cognostics)))
      
    } else {
      
      # Make a group variable name just to make group_by_at work properly
      theGroup <- "Group"
      
      # Add metas per group
      toBuild <- toBuild %>% 
        dplyr::group_by_at(c(panel, theGroup)) %>%
        dplyr::summarise(
          count = sum(!is.na(Abundance)),
          `mean abundance` = mean(Abundance, na.rm = T), 
        ) %>%
        tidyr::pivot_wider(id_cols = panel, names_from = Group, values_from = c(count, `mean abundance`), names_sep = " ")
      
      # Now remove unwanted cognostics
      if ("mean abundance" %in% cognostics == FALSE) {
        toBuild <- toBuild[,grepl("mean abundance", colnames(toBuild)) == FALSE]
      }
      if (any(grepl("count", cognostics)) == FALSE) {
        toBuild <- toBuild[,grepl("count", colnames(toBuild)) == FALSE]
      }
      
      # Add emeta columns
      if (!is.null(attr(trelliData, "emeta_col"))) {
        
        # Pull emeta uniqued columns that should have been prepped in the pivot_longer section
        emeta <- trelliData$trelliData[,c(panel, attr(trelliData, "emeta_col"))] %>% unique()
        
        # Add emeta cognostics
        toBuild <- dplyr::left_join(toBuild, emeta, by = panel)
        
      }
      
      # Now, select only the requested cognostics 
      toSelect <- lapply(colnames(toBuild), function(column) {
        lapply(c(panel, cognostics), function(choice) {grepl(pattern = choice, x = column)}) %>%
          unlist() %>%
          any()
      }) %>% unlist()
      toBuild <- toBuild[,toSelect]
      
    }
    
    # Add p-values and fold changes 
    if ("anova p-value" %in% cognostics) {
      anova_cols <- colnames(trelliData$trelliData)[grepl("P_value_A", colnames(trelliData$trelliData))]
      anova_data <- trelliData$trelliData[,c(panel, anova_cols)] %>% unique()
      toBuild <- dplyr::left_join(toBuild, anova_data, by = panel)
    }
    if ("fold change" %in% cognostics) {
      fc_cols <- colnames(trelliData$trelliData)[grepl("Fold_change", colnames(trelliData$trelliData))]
      fc_data <- trelliData$trelliData[,c(panel, fc_cols)] %>% unique()
      toBuild <- dplyr::left_join(toBuild, fc_data, by = panel)
    }
    
    # Filter down if test mode
    if (test_mode) {
      toBuild <- toBuild[test_example,]
    }
    
  }
  
  # Second, make boxplot function-----------------------------------------------

  box_plot_fun <- function(Panel) {
    
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}})
    
    # Build plot
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = Abundance)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylab(paste(attr(trelliData$omicsData, "data_info")$data_scale, "Abundance"))

    # Add include_points
    if (include_points) {
      boxplot <- boxplot + ggplot2::geom_jitter(height = 0, width = 0.25)
    }

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        boxplot <- boxplot + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      boxplot <- boxplot %>% plotly::ggplotly()
    }

    return(boxplot)
  }
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(box_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Build trelliscope display---------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
    
  } else {
    
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
    
  }
}

#' @name trelli_abundance_histogram
#'
#' @title Histogram trelliscope building function for abundance data
#'
#' @description Specify a plot design and cognostics for the abundance histogram
#'   trelliscope. Main_effects grouping are ignored. Data must be grouped by
#'   edata_cname. For RNA-Seq data, use "trelli_rnaseq_histogram". 
#'
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by edata_cname in trelli_panel_by.
#'   Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are "count", "mean abundance", "median abundance", "cv abundance", 
#'   and "skew abundance". All are included by default. Any additional cognostics 
#'   from e_meta variables may be added if their column name is specified, and the 
#'   data has been paneled by the edata_cname or by an e_meta variable.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "ylim(c(1,2))").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of histograms that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' # Build the abundance histogram with an edata file. 
#' # Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' 
#' # Build the abundance histogram with an omicsData object. 
#' # Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'     
#' # Build the abundance histogram with an omicsData and statRes object. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData4, panel = "Peptide") %>%
#'    trelli_abundance_histogram(
#'      test_mode = TRUE, test_example = 1:10, cognostics = "sample count", path = tempdir())
#'    
#' # Users can modify the plotting function with ggplot parameters and interactivity, 
#' # and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "xlab('Abundance')"), interactive = TRUE,
#'      cognostics = c("mean abundance", "median abundance"), path = tempdir())  
#'  
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_abundance_histogram <- function(trelliData,
                                       cognostics = c("sample count", "mean abundance", "median abundance", "cv abundance", "skew abundance"),
                                       ggplot_params = NULL,
                                       interactive = FALSE,
                                       path = tempdir(),
                                       name = "Trelliscope",
                                       build_trelliscope = TRUE,
                                       test_mode = FALSE,
                                       test_example = 1,
                                       single_plot = FALSE,
                                       ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("sample count", "mean abundance", "median abundance", "cv abundance", "skew abundance"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_histogram instead.",
                  p_value_thresh = NULL)

  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
  if (edata_cname != panel) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Start summary toBuild data.frame
  toBuild <- trelliData$trelliData
  
  # First, add any missing cognostics-------------------------------------------
  toBuild <- toBuild %>% 
    dplyr::group_by_at(panel) %>%
    dplyr::summarize(
      `sample count` = sum(!is.na(Abundance)),
      `mean abundance` = mean(Abundance, na.rm = T),
      `median abundance` = median(Abundance, na.rm = T),
      `cv abundance` = sd(Abundance, na.rm = T) / `mean abundance` * 100,
      `skew abundance` = e1071::skewness(Abundance, na.rm = T)
    )
  
  # Add emeta columns
  if (!is.null(attr(trelliData, "emeta_col"))) {
    
    # Pull emeta uniqued columns that should have been prepped in the pivot_longer section
    emeta <- trelliData$trelliData[,c(panel, attr(trelliData, "emeta_col"))] %>% unique()
    
    # Add emeta cognostics
    toBuild <- dplyr::left_join(toBuild, emeta, by = panel)
    
  }
  
  # Now, select only the requested cognostics 
  toSelect <- lapply(colnames(toBuild), function(column) {
    lapply(c(panel, cognostics), function(choice) {grepl(pattern = choice, x = column)}) %>%
      unlist() %>%
      any()
  }) %>% unlist()
  toBuild <- toBuild[,toSelect]
  
  # Filter down if test mode
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }

  # Make histogram function-----------------------------------------------------

  hist_plot_fun <- function(Panel) {
    
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}})
    
    # Remove NAs
    DF <- DF[!is.na(DF$Abundance), ]

    # Build plot
    histogram <- ggplot2::ggplot(DF, ggplot2::aes(x = Abundance)) +
      ggplot2::geom_histogram(bins = 10, fill = "steelblue", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::xlab(paste(attr(trelliData$omicsData, "data_info")$data_scale, "Abundance")) +
      ggplot2::ylab("Frequency")

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        histogram <- histogram + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      histogram <- histogram %>% plotly::ggplotly()
    }

    return(histogram)
  }

  # Build trelliscope display---------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(hist_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
    
  } else {
    
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
    
  }
}

#' @name trelli_abundance_heatmap
#'
#' @title Heatmap trelliscope building function for abundance data
#'
#' @description Specify a plot design and cognostics for the abundance heatmap
#'   trelliscope. Data must be grouped by an e_meta column. Main_effects order
#'   the y-variables. All statRes data is ignored. For RNA-Seq data, use "trelli_rnaseq_heatmap".
#'   Note that z-scores are plotted per biomolecule for ease of comparison.
#'
#' @param trelliData A trelliscope data object made by as.trelliData, and
#'   grouped by an emeta variable. Required.
#' @param cognostics A vector of cognostic options. Defaults are "sample count", 
#'   "mean abundance" and "biomolecule count". "sample count" and "mean abundance"
#'   are reported per group, and "biomolecule count" is the total number of biomolecules
#'   in the biomolecule class (e_meta column).Any additional cognostics from e_meta 
#'   variables may be added if their column name is specified, and the data has 
#'   been paneled by the edata_cname  or by an e_meta variable.
#' @param ggplot_params An optional vector of strings ofggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   tempdir.
#' @param name The name of the display. Default is Trelliscope
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of heatmaps that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' # Build the abundance heatmap with an omicsData object with emeta variables. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData2, panel = "RazorProtein") %>%
#'    trelli_abundance_heatmap(test_mode = TRUE, test_example = 1:3, path = tempdir())
#'    
#' # Users can modify the plotting function with ggplot parameters and interactivity, 
#' # and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>% 
#'    trelli_abundance_heatmap(
#'      test_mode = TRUE, test_example = 1:5, 
#'      ggplot_params = c("ylab('')", "xlab('')"), 
#'      interactive = TRUE, cognostics = c("biomolecule count"),
#'      path = tempdir()
#'    )
#'
#' \dontshow{closeAllConnections()}
#' }
#' }
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_abundance_heatmap <- function(trelliData,
                                     cognostics = c("sample count", "mean abundance", "biomolecule count"),
                                     ggplot_params = NULL,
                                     interactive = FALSE,
                                     path = tempdir(),
                                     name = "Trelliscope",
                                     build_trelliscope = TRUE,
                                     test_mode = FALSE,
                                     test_example = 1,
                                     single_plot = FALSE,
                                     ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("sample count", "mean abundance", "biomolecule count"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_heatmap instead.",
                  p_value_thresh = NULL)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is grouped by an e_meta variable
  if (panel %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }

  # If no group designation, set cognostics to NULL.
  if (is.null(attributes(trelliData$omicsData)$group_DF)) {
    cognostics <- NULL
  }
  
  # Get the column names of edata and fdata
  edata_cname <- get_edata_cname(trelliData$omicsData)
  fdata_cname <- get_fdata_cname(trelliData$omicsData)
  
  # Start builder dataframe
  toBuild <- trelliData$trelliData
  
  # First, add any missing cognostics-------------------------------------------
  
  if (!single_plot) {
  
    # Make a group variable name just to make group_by_at work properly
    theGroup <- "Group"
    
    # Add cognostics per group
    toBuild <- toBuild %>% 
      dplyr::group_by_at(c(panel, theGroup)) %>%
      dplyr::summarise(
        count = sum(!is.na(Abundance)),
        `mean abundance` = mean(Abundance, na.rm = T)
      ) %>%
      tidyr::pivot_wider(id_cols = panel, names_from = Group, values_from = c(count, `mean abundance`), names_sep = " ")
    
    # Remove unwanted cognostics 
    if ("sample count" %in% cognostics == FALSE) {
      toBuild <- toBuild[,!grepl("count", colnames(toBuild))]
    }
    if ("mean abundance" %in% cognostics == FALSE) {
      toBuild <- toBuild[,!grepl("mean abundance", colnames(toBuild))]
    }
    
    # Add biomolecule count if requested 
    if ("biomolecule count" %in% cognostics) {
      bio_counts <- trelliData$trelliData %>% 
        dplyr::select(panel, edata_cname) %>%
        unique() %>%
        dplyr::group_by_at(panel) %>%
        dplyr::summarise(`biomolecule count` = dplyr::n())
      toBuild <- dplyr::left_join(toBuild, bio_counts, by = panel)
    }
    
    # Now, select only the requested cognostics 
    toSelect <- lapply(colnames(toBuild), function(column) {
      lapply(c(panel, cognostics), function(choice) {grepl(pattern = choice, x = column)}) %>%
        unlist() %>%
        any()
    }) %>% unlist()
    toBuild <- toBuild[,toSelect]
    
    # Filter down if test mode
    if (test_mode) {
      toBuild <- toBuild[test_example,]
    }
  }

  # Make heatmap function-------------------------------------------------------

  hm_plot_fun <- function(Panel) {
    
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}})

    # If group designation was set, then convert Group to a factor variable
    if (!is.null(attributes(trelliData$omicsData)$group_DF)) {
      theLevels <- attr(trelliData$omicsData, "group_DF") %>%
        dplyr::arrange(Group) %>%
        dplyr::select(dplyr::all_of(fdata_cname)) %>%
        unlist()
      DF[, fdata_cname] <- factor(unlist(DF[, fdata_cname]), levels = theLevels)
    } else {
      DF[, fdata_cname] <- as.factor(unlist(DF[, fdata_cname]))
    }
    
    # Calculate z score
    DF <- DF %>%
      dplyr::group_by_at(edata_cname) %>%
      dplyr::mutate(`Z-Score` = (Abundance - mean(Abundance, na.rm = T)) / sd(Abundance, na.rm = T))
    
    # Build plot: this should be edata_cname
    hm <- ggplot2::ggplot(DF, ggplot2::aes(x = as.factor(.data[[edata_cname]]), y = .data[[fdata_cname]], fill = `Z-Score`)) +
      ggplot2::geom_tile() +
      ggplot2::theme_bw() +
      ggplot2::ylab("Sample") +
      ggplot2::xlab("Biomolecule") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "white") 

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        hm <- hm + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      hm <- hm %>% plotly::ggplotly()
    }

    return(hm)
  }

  # Build trelliscope display---------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(hm_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
    
  } else {
    
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
    
  }
}

#' @name trelli_missingness_bar
#'
#' @title Bar chart trelliscope building function for missing data
#'
#' @description Specify a plot design and cognostics for the missing barchart
#'    trelliscope. Missingness is displayed per panel_by variable. Main_effects
#'    data is used to split samples when applicable. For RNA-Seq data, use 
#'    "trelli rnaseq nonzero bar". 
#'
#' @param trelliData A trelliscope data object made by as.trelliData.edata or
#'    as.trelliData. Required.
#' @param cognostics A vector of cognostic options for each plot. Defaults are "total count",
#'    "observed count", and "observed proportion". If grouping
#'    data is included, all cognostics will be reported per group. If the 
#'    trelliData is paneled by a biomolecule, the counts and proportion we be 
#'    samples. If paneled by a sample or biomolecule class, the counts and proportions
#'    will be biomolecules. If statRes data is included, "g-test p-value" may 
#'    be included. 
#' @param proportion A logical to determine whether plots should display counts
#'    or proportions. Default is TRUE.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'    the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'    Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'    interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'    now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'    a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'    and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'    is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'    trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of missingness bar charts that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' # Build the missingness bar plot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData1, panel = "Peptide") %>% 
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' trelli_panel_by(trelliData = trelliData1, panel = "Sample") %>% 
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10, 
#'    cognostics = "observed proportion", path = tempdir())
#' 
#' # Build the missingness bar plot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' 
#' # Build the missingness bar plot with a statRes object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData3, panel = "Peptide") %>%
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10, path = tempdir(),
#'                          cognostics = c("observed proportion", "g-test p-value"))
#' 
#' # Build the missingness bar plot with an omicsData and statRes object. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>%
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10, path = tempdir()) 
#' 
#' # Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'    trelli_missingness_bar(
#'      test_mode = TRUE, test_example = 1:5, interactive = TRUE, path = tempdir())
#'    
#' # Or visualize only count data 
#' trelli_panel_by(trelliData = trelliData2, panel = "Peptide") %>% 
#'    trelli_missingness_bar(
#'      test_mode = TRUE, test_example = 1:5, 
#'      cognostics = "observed count", proportion = FALSE,
#'      path = tempdir()
#'    )
#'
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_missingness_bar <- function(trelliData,
                                   cognostics = c("total count", "observed count", "observed proportion"),
                                   proportion = TRUE,
                                   ggplot_params = NULL,
                                   interactive = FALSE,
                                   path = tempdir(),
                                   name = "Trelliscope",
                                   build_trelliscope = TRUE,
                                   test_mode = FALSE,
                                   test_example = 1,
                                   single_plot = FALSE,
                                   ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("either"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("total count", "observed count", "observed proportion", "g-test p-value", attr(trelliData, "emeta_col")),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_nonzero_bar instead.",
                  p_value_thresh = NULL)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Check that proportion is a non NA logical
  if (!is.logical(proportion) | is.na(proportion)) {
    stop("proportion must be a TRUE or FALSE.")
  }
  
  # Test whether statRes data is included to use g-test p-value
  if ("g-test p-value" %in% cognostics) {
    if (is.null(trelliData$statRes)) {
      cognostics <- cognostics[cognostics != "g-test p-value"]
      message("'g-test p-value' can only be included if stats data (statRes) is included")
    } else if (panel != attr(trelliData, "edata_col")) {
      cognostics <- cognostics[cognostics != "g-test p-value"]
      message("'g-test p-value' can only be included if the data has been paneled by the biomolecule column 'edata_cname'")
    }
  }
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Determine if the trelliData is paneled by the edata column
  paneled_by_edata <- panel == attr(trelliData, "edata_col")
  
  # Set stats mode 
  stats_mode <- FALSE
  if (is.null(trelliData$omicsData)) {
    stats_mode <- TRUE
  }

  # Start builder dataframe
  toBuild <- trelliData$trelliData
  theGroupName <- "Group"
  
  # First, add any missing cognostics-------------------------------------------
  
  # Extract observed counts 
  if (!is.null(trelliData$omicsData)) {
    
    if ("Group" %in% colnames(trelliData$trelliData) == FALSE) {
      toBuild <- toBuild %>% 
        dplyr::group_by_at(panel) %>%
        dplyr::summarize(
          `total count` = dplyr::n(),
          `observed count` = sum(!is.na(Abundance)),
          `observed proportion` = round(`observed count` / `total count`, 4)
        )
    } else {
      toBuild <- toBuild %>%
        dplyr::group_by_at(c(panel, theGroupName)) %>%
        dplyr::summarize(
          `total count` = dplyr::n(),
          `observed count` = sum(!is.na(Abundance)),
          `observed proportion` = round(`observed count` / `total count`, 4)
        ) %>% 
        tidyr::pivot_wider(id_cols = panel, names_from = theGroupName, names_sep = " ",
                           values_from = c("total count", "observed count", "observed proportion"))
    }

  } else {
    
    # Get count columns 
    cols2rename <- colnames(toBuild)[grepl("Count_", colnames(toBuild))]
    
    # Make totals dataframe
    totalsDF <- attr(trelliData$statRes, "group_DF")$Group %>% 
      table(dnn = "Group") %>% 
      data.frame() %>% 
      dplyr::rename(`total count` = Freq)
    
    # Extract panel and observed counts. Format with group information. Add
    # totals and observed proportion. 
    toBuild <- toBuild %>%
      dplyr::select_at(c(panel, cols2rename)) %>%
      tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
      dplyr::mutate(name = gsub("Count_", "", name, fixed = T)) %>%
      dplyr::rename(Group = name, `observed count` = value) %>%
      dplyr::left_join(totalsDF, by = "Group") %>%
      dplyr::mutate(`observed proportion` = round(`observed count` / `total count`, 4)) %>%
      tidyr::pivot_wider(id_cols = panel, names_from = theGroupName, names_sep = " ",
                         values_from = c("total count", "observed count", "observed proportion"))
  
  }
  
  # Add p-values if possible
  if ("g-test p-value" %in% cognostics) {
    
    # Determine if g-tests were conducted
    gcols_pos <- grepl("P_value_G", colnames(trelliData$trelliData))
    
    # If no g-tests, skip 
    if (!any(gcols_pos)) {
      message("No g-test p-value columns found. Skipping this cognostic.")
      cognostics <- cognostics[!(cognostics == "p-value g-test")]
    } else {
      gcols <- colnames(trelliData$trelliData)[gcols_pos]
      toBuild <- dplyr::left_join(toBuild, 
        trelliData$trelliData %>% 
          dplyr::select_at(c(panel, gcols)), 
        by = panel
      )
    }
  }
  
  # Add emeta columns
  if (!is.null(trelliData) && panel != attr(trelliData, "fdata_col") && !is.null(attr(trelliData, "emeta_col"))) {
    
    # Pull emeta uniqued columns that should have been prepped in the pivot_longer section
    emeta <- trelliData$trelliData[,c(panel, attr(trelliData, "emeta_col"))] %>% unique()
    
    # Add emeta cognostics
    toBuild <- dplyr::left_join(toBuild, emeta, by = panel) %>% 
      unique()
    
  }
  
  # Now, select only the requested cognostics 
  toSelect <- lapply(colnames(toBuild), function(column) {
    lapply(c(panel, cognostics), function(choice) {grepl(pattern = choice, x = column)}) %>%
      unlist() %>%
      any()
  }) %>% unlist()
  toBuild <- toBuild[,toSelect]
  
  # Filter down if test mode --> must be here to get the right selection of panels
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }
  
  # Second, write the plotting function-----------------------------------------

  # First, generate the boxplot function
  missing_bar_plot_fun <- function(Panel) {
    
    # WARNING: trelliData$trelliData must be used every time. There is no escape. 
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}}) 
    
    # Set group information
    if ("Group" %in% colnames(DF) == FALSE) {DF$Group <- "X"}
    
    if (is.null(trelliData$omicsData)) {
      
      # Due to the way the new trelliscope functions run the lazy_plot function,
      # we have to re-do how we calculated the columns for the summary toBuild file,
      # but we can use the same variables defined globally.
      MissPlotDF <- DF %>%
        dplyr::select_at(c(panel, cols2rename)) %>%
        tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
        dplyr::mutate(name = gsub("Count_", "", name, fixed = T)) %>%
        dplyr::rename(Group = name, `observed count` = value) %>%
        dplyr::left_join(totalsDF, by = "Group") %>%
        dplyr::mutate(
          `observed proportion` = round(`observed count` / `total count`, 4),
          `missing count` = `total count` - `observed count`,
          `missing proportion` = round(`missing count` / `total count`, 4)
        ) %>%
        tidyr::pivot_longer(cols = c(3:7))
      
    } else {
      
      # Split out group names from summary toBuild table, pivot_wider to add missing
      # columns, an then pivot longer 
      MissPlotDF <- DF %>% 
        dplyr::group_by_at(c(panel, theGroupName)) %>%
        tidyr::nest() %>%
        dplyr::mutate(
          `total count` = purrr::map_int(data, nrow),
          `observed count` = purrr::map_int(data, function(x) {sum(!is.na(x$Abundance))}),
          `missing count` = `total count` - `observed count`,
          `observed proportion` = round(`observed count` / `total count`, 4),
          `missing proportion` = round(`missing count` / `total count`, 4)
        ) %>%
        dplyr::select(-data) %>%
        tidyr::pivot_longer(cols = c(3:7))
      
    }
    
    # Subset based on count or proportion
    if (proportion) {
      MissPlotDF <- MissPlotDF %>%
        dplyr::filter(name %in% c("observed proportion", "missing proportion")) %>%
        dplyr::mutate(name = 
          factor(ifelse(name == "observed proportion", "Present", "Absent"), levels = c("Absent", "Present"))
        )
      ylab <- "Proportion"
    } else {
      MissPlotDF <- MissPlotDF %>%
        dplyr::filter(name %in% c("observed count", "missing count")) %>%
        dplyr::mutate(name = 
          factor(ifelse(name == "observed count", "Present", "Absent"), levels = c("Absent", "Present"))
        )
      ylab <- "Count"
    }

    # Build plot
    missing_bar <- ggplot2::ggplot(MissPlotDF, ggplot2::aes(x = Group, y = value, fill = name)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::ylab(ylab) +
      ggplot2::scale_fill_manual(values = c("Absent" = "black", "Present" = "steelblue")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank()
      )

    # Remove x axis if no groups
    if ("Group" %in% colnames(trelliData$trelliData) == FALSE & stats_mode == FALSE) {
      missing_bar <- missing_bar + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      )
    }

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        missing_bar <- missing_bar + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      missing_bar <- missing_bar %>% plotly::ggplotly()
    }

    return(missing_bar)
  }


  # Build trelliscope display---------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(missing_bar_plot_fun)) %>%
    dplyr::select(-Panel)

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
    
  } else {
    
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
    
  }
  
}

determine_significance <- function(DF, p_value_thresh) {
  
  # Dermine significance
  DF$Significance <- DF$p_value <= p_value_thresh
  
  # Filter out NA values 
  DF <- DF[!(is.nan(DF$Significance) | is.nan(DF$fold_change)), ]
  
  # Return NULL if no rows left
  if (nrow(DF) == 0) {
    return(NULL)
  }

  # Add variables
  GreaterThan <- paste0(">", p_value_thresh)
  LessThan <- paste0("<", p_value_thresh)
  DF$Significance <- ifelse(DF$Significance, LessThan, GreaterThan)

  # Add attributes for plotting functions
  attr(DF, "GreaterThan") <- GreaterThan
  attr(DF, "LessThan") <- LessThan

  return(DF)
}

#' @name trelli_foldchange_bar
#'
#' @title Bar chart trelliscope building function for fold_change
#'
#' @description Specify a plot design and cognostics for the fold_change
#'   barchart trelliscope. Fold change must be grouped by edata_cname.
#'
#' @param trelliData A trelliscope data object with statRes results. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   and the defaults are "fold change" and "p-value". If the omics data is MS/NMR,
#'   p-value will be the results from the ANOVA test. If the omics data is seqData,
#'   the p-value will be the results from the function "diffexp_seq".
#' @param p_value_thresh A value between 0 and 1 to indicate significant
#'   biomolecules for p_value_test. Default is 0.05.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of fold_change bar plots that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' # Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData3, panel = "Peptide") %>% 
#'   trelli_foldchange_bar(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_foldchange_bar <- function(trelliData,
                                  cognostics = c("fold change", "p-value"),
                                  p_value_thresh = 0.05,
                                  ggplot_params = NULL,
                                  interactive = FALSE,
                                  path = tempdir(),
                                  name = "Trelliscope",
                                  build_trelliscope = TRUE,
                                  test_mode = FALSE,
                                  test_example = 1,
                                  single_plot = FALSE,
                                  ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "stat",
                  cognostics = cognostics,
                  acceptable_cognostics = c("fold change", "p-value"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is edata_cname
  edata_cname <- attr(trelliData, "edata_col")
  if (panel != edata_cname) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Start builder dataframe
  toBuild <- trelliData$trelliData
  
  # First, generate the cognostics----------------------------------------------
  
  # Extract fold change and p-value ANOVA (omics) or p-value columns
  needed_cols <- colnames(toBuild)[grepl("P_value|Fold_change", colnames(toBuild))]
  if (any(grepl("P_value_G_", needed_cols))) {
    needed_cols <- needed_cols[!grepl("P_value_G", needed_cols)]
  }
  
  # Pull required columns 
  toBuild <- toBuild %>% dplyr::select_at(c(panel, needed_cols))
  
  # Add emeta columns: TO FIX
  if (!is.null(trelliData) && panel != attr(trelliData, "fdata_col") && !is.null(attr(trelliData, "emeta_col"))) {
    
    # Pull emeta uniqued columns that should have been prepped in the pivot_longer section
    emeta <- trelliData$trelliData[,c(panel, attr(trelliData, "emeta_col"))] %>% unique()
    
    # Add emeta cognostics
    toBuild <- dplyr::left_join(toBuild, emeta, by = panel) %>% 
      unique()
    
  }
  
  # Filter down if test mode --> must be here to get the right selection of panels
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }
  
  # Make foldchange bar function------------------------------------------------

  fc_bar_plot_fun <- function(Panel) {
    
    # Pull data.frame, and use needed cols from earlier. Then pull the comparisons,
    # fold changes, and p-values. Make it flexible for seqData as well. 
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}}) %>%
      dplyr::select_at(c(panel, needed_cols)) %>%
      unique() %>%
      dplyr::mutate_at(colnames(.)[2:ncol(.)], as.numeric) %>%
      tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
      dplyr::mutate(
        Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
        Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
      ) %>% 
      dplyr::select(-name) %>%
      tidyr::pivot_wider(id_cols = c(panel, Comparison), names_from = Type, values_from = value)
    
    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh)
      if (is.null(DF)) {return(NULL)}
      
      # Make bar plot 
      bar <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison, color = Significance)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_color_manual(values = structure(c("black", NA), .Names = c(attr(DF, "LessThan"), attr(DF, "GreaterThan"))), na.value = NA)
    } else {
      # Make bar plot
      bar <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
        ggplot2::geom_bar(stat = "identity")
    }

    # Extend bar plot
    bar <- bar + ggplot2::theme_bw() + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      ) + ggplot2::ylab("Fold Change") + ggplot2::guides(fill = "none")

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        bar <- bar + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      bar <- bar %>% plotly::ggplotly()
    }

    return(bar)
  }


  # Build trelliscope function--------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(fc_bar_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
    
  } else {
    
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
  }
}

#' @name trelli_foldchange_boxplot
#'
#' @title Boxplot trelliscope building function for fold_changes
#'
#' @description Specify a plot design and cognostics for the fold_change boxplot
#'   trelliscope. Fold change must be grouped by an emeta column, which means
#'   both an omicsData object and statRes are required to make this plot.
#'
#' @param trelliData A trelliscope data object with omicsData and statRes
#'   results. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries are 
#'   "biomolecule count", "proportion significant", "mean fold change",
#'   and "sd fold change". Default is "biomolecule count". 
#' @param p_value_thresh A value between 0 and 1 to indicate significant
#'   biomolecules for the anova (MS/NMR) or diffexp_seq (RNA-seq) test. Default is 0.05.
#' @param include_points Add points. Default is TRUE.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'
#' @return No return value, builds a trelliscope display of fold_change boxplots that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{ 
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' # Build fold_change box plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>% 
#'   trelli_foldchange_boxplot(test_mode = TRUE, 
#'                             test_example = 1:10,
#'                             cognostics = c("biomolecule count", 
#'                                            "proportion significant", 
#'                                            "mean fold change",
#'                                            "sd fold change"),
#'                             path = tempdir()
#'                            )
#'                            
#'                            
#' #####################
#' ## RNA-SEQ EXAMPLE ##                            
#' #####################
#'
#' # Build fold_change box plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") %>% 
#'   trelli_foldchange_boxplot(test_mode = TRUE, 
#'                             test_example = c(16823, 16890, 17680, 17976, 17981, 19281),
#'                             cognostics = c("biomolecule count", 
#'                                            "proportion significant", 
#'                                            "mean fold change",
#'                                            "sd fold change"),
#'                             path = tempdir()
#'                            )
#' 
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_foldchange_boxplot <- function(trelliData,
                                      cognostics = "biomolecule count",
                                      p_value_thresh = 0.05,
                                      include_points = TRUE,
                                      ggplot_params = NULL,
                                      interactive = FALSE,
                                      path = tempdir(),
                                      name = "Trelliscope",
                                      build_trelliscope = TRUE,
                                      test_mode = FALSE,
                                      test_example = 1,
                                      single_plot = FALSE,
                                      ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("omics", "stat"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("biomolecule count", 
                                            "proportion significant", 
                                            "mean fold change",
                                            "sd fold change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is an emeta column
  if (!(panel %in% attr(trelliData, "emeta_col"))) {
    stop("trelliData must be paneled_by an e_meta column.")
  }

  # Make sure include_points is a true or false
  if (!is.logical(include_points) & !is.na(include_points)) {
    stop("include_points must be a TRUE or FALSE.")
  }
  
  # Start builder dataframe
  toBuild <- trelliData$trelliData
  
  # Make cognostic function-----------------------------------------------------
  
  # Extract fold change and p-value ANOVA (omics) or p-value columns (seqData)
  needed_cols <- colnames(toBuild)[grepl("P_value|Fold_change", colnames(toBuild))]
  if (any(grepl("P_value_G_", needed_cols))) {
    needed_cols <- needed_cols[!grepl("P_value_G", needed_cols)]
  }
  
  # Pull required columns 
  toBuild <- toBuild %>% dplyr::select_at(c(panel, needed_cols))
  toBuild <- unique(toBuild)
  theComparison <- "Comparison"
  
  # Pivot longer to extract comparisons. Then split into p-value and fold change columns
  preBuild <- toBuild %>%
    tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
    dplyr::mutate(
      Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
      Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
    ) %>%
    dplyr::select(-name)
  
  # Pull fold change information and calculate metas 
  foldBuild <- preBuild %>%
    dplyr::filter(Type == "fold_change") %>%
    dplyr::group_by_at(c(panel, theComparison)) %>%
    dplyr::summarize(
      `biomolecule count` = sum(!is.na(value)), 
      `mean fold change` = mean(value, na.rm = T),
      `sd fold change` = sd(value, na.rm = T)
    ) 
  
  # Pull p-value information and calculate metas 
  pvalBuild <- preBuild %>%
    dplyr::filter(Type == "p_value") %>%
    dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>%
    dplyr::group_by_at(c(panel, theComparison)) %>%
    dplyr::summarize(SignificantCount = sum(value <= p_value_thresh))
  
  # Merge data.frames, calculate proportion significant, pivot wider by comparison again
  toBuild <- dplyr::left_join(foldBuild, pvalBuild, by = c(panel, theComparison)) %>%
    dplyr::mutate(`proportion significant` = SignificantCount / `biomolecule count`) %>%
    dplyr::select(-SignificantCount) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = panel, names_from = theComparison, names_sep = " ",
                       values_from = c("biomolecule count", "mean fold change", "sd fold change", "proportion significant"))
  
  # Add additiona e_meta columns
  emeta <- trelliData$trelliData[,unique(c(panel, attr(trelliData, "emeta_col")))] %>% unique()
  
  # Add emeta cognostics
  toBuild <- dplyr::left_join(toBuild, emeta, by = panel) %>% unique()
  
  # Filter down if test mode --> must be here to get the right selection of panels
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }
  
  # Make foldchange boxplot function--------------------------------------------

  fc_box_plot_fun <- function(Panel) {
    
    # Pull data.frame, and use needed cols from earlier. Then pull the comparisons,
    # fold changes, and p-values. Make it flexible for seqData as well. 
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}}) %>%
      dplyr::select_at(c(panel, needed_cols)) %>%
      unique() %>%
      dplyr::mutate_at(colnames(.)[2:ncol(.)], as.numeric) %>%
      tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
      dplyr::mutate(
        Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
        Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
      ) %>% 
      dplyr::select(-name) %>%
      tidyr::pivot_wider(id_cols = c(panel, theComparison), names_from = Type, 
                         values_from = value, values_fn = list) %>%
      tidyr::unnest(cols = c(fold_change, p_value))
    
    if (p_value_thresh != 0) {
    
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh)
      if (is.null(DF)) {
        return(NULL)
      }
    
      # Make boxplot
      boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
        ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5), 
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
        ) + ggplot2::ylab("Fold Change") + 
        ggplot2::guides(fill = "none")
      
      # Add include_points
      if (include_points) {
        boxplot <- boxplot + ggplot2::geom_jitter(ggplot2::aes(shape = Significance), height = 0, width = 0.25) +
          ggplot2::scale_shape_manual(values = structure(c(17, 16), .Names = c(attr(DF, "LessThan"), attr(DF, "GreaterThan"))))
      }
      
    } else {
      
      # Make boxplot
      boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
        ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5), 
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
        ) + ggplot2::ylab("Fold Change") + 
        ggplot2::guides(fill = "none")
      
      # Add include_points
      if (include_points) {
          boxplot <- boxplot + ggplot2::geom_jitter(height = 0, width = 0.25)
      }
    }

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        boxplot <- boxplot + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      boxplot <- boxplot %>% plotly::ggplotly()
    }

    return(boxplot)
  }

  # Build the trelliscope-------------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(fc_box_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
  } else {
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
  }
}

#' @name trelli_foldchange_heatmap
#' 
#' @title Heatmap trelliscope building function for fold_change
#'
#' @description Specify a plot design and cognostics for the fold_change heatmap
#'   trelliscope. Fold change must be grouped by an emeta column, which means
#'   both an omicsData object and statRes are required to make this plot.
#'
#' @param trelliData A trelliscope data object with omicsData and statRes
#'   results. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries are 
#'   "biomolecule count", "proportion significant", "mean fold change",
#'   and "sd fold change". Default is "biomolecule count". 
#' @param p_value_thresh A value between 0 and 1 to indicate significant biomolecules 
#'   for the anova (MS/NMR) or diffexp_seq (RNA-seq) test. Default is 0.05.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#' 
#' @return No return value, builds a trelliscope display of fold-change heatmaps that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' 
#' ##########################
#' ## MS/NMR OMICS EXAMPLE ##
#' ##########################
#' 
#' # Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>% 
#'   trelli_foldchange_heatmap(test_mode = TRUE, 
#'                             test_example = 1:10,
#'                             path = tempdir())
#' 
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_foldchange_heatmap <- function(trelliData,
                                      cognostics = "biomolecule count",
                                      p_value_thresh = 0.05,
                                      ggplot_params = NULL,
                                      interactive = FALSE,
                                      path = tempdir(),
                                      name = "Trelliscope",
                                      build_trelliscope = TRUE,
                                      test_mode = FALSE,
                                      test_example = 1,
                                      single_plot = FALSE,
                                      ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("omics", "stat"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("biomolecule count", 
                                            "proportion significant", 
                                            "mean fold change",
                                            "sd fold change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is an emeta column
  if (!(panel %in% attr(trelliData, "emeta_col"))) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Start builder dataframe
  toBuild <- trelliData$trelliData
  
  # First, build the metas------------------------------------------------------
  
  # Extract fold change and p-value ANOVA (omics) or p-value columns (seqData)
  needed_cols <- colnames(toBuild)[grepl("P_value|Fold_change", colnames(toBuild))]
  if (any(grepl("P_value_G_", needed_cols))) {
    needed_cols <- needed_cols[!grepl("P_value_G", needed_cols)]
  }
  
  # Pull required columns 
  toBuild <- toBuild %>% dplyr::select_at(c(panel, needed_cols))
  toBuild <- unique(toBuild)
  theComparison <- "Comparison"
  
  # Pivot longer to extract comparisons. Then split into p-value and fold change columns
  preBuild <- toBuild %>%
    tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
    dplyr::mutate(
      Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
      Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
    ) %>%
    dplyr::select(-name)
  
  # Pull fold change information and calculate metas 
  foldBuild <- preBuild %>%
    dplyr::filter(Type == "fold_change") %>%
    dplyr::group_by_at(c(panel, theComparison)) %>%
    dplyr::summarize(
      `biomolecule count` = sum(!is.na(value)), 
      `mean fold change` = mean(value, na.rm = T),
      `sd fold change` = sd(value, na.rm = T)
    ) 
  
  # Pull p-value information and calculate metas 
  pvalBuild <- preBuild %>%
    dplyr::filter(Type == "p_value") %>%
    dplyr::mutate(value = ifelse(is.na(value), 1, value)) %>%
    dplyr::group_by_at(c(panel, theComparison)) %>%
    dplyr::summarize(SignificantCount = sum(value <= p_value_thresh))
  
  # Merge data.frames, calculate proportion significant, pivot wider by comparison again
  toBuild <- dplyr::left_join(foldBuild, pvalBuild, by = c(panel, theComparison)) %>%
    dplyr::mutate(`proportion significant` = SignificantCount / `biomolecule count`) %>%
    dplyr::select(-SignificantCount) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = panel, names_from = theComparison, names_sep = " ",
                       values_from = c("biomolecule count", "mean fold change", "sd fold change", "proportion significant"))
  
  # Add additiona e_meta columns
  emeta <- trelliData$trelliData[,unique(c(panel, attr(trelliData, "emeta_col")))] %>% unique()
  
  # Add emeta cognostics
  toBuild <- dplyr::left_join(toBuild, emeta, by = panel) %>% unique()
  
  # Filter down if test mode --> must be here to get the right selection of panels
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }
  
  # Make foldchange boxplot function--------------------------------------------
  
  fc_hm_plot_fun <- function(Panel) {
    
    # Get edata cname
    edata_cname <- attr(trelliData, "edata_col")
    
    # Pull data.frame, and use needed cols from earlier. Then pull the comparisons,
    # fold changes, and p-values. Make it flexible for seqData as well. 
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}}) %>%
      dplyr::select_at(c(panel, edata_cname, needed_cols)) %>%
      unique() %>%
      dplyr::mutate_at(colnames(.)[3:ncol(.)], as.numeric) %>%
      tidyr::pivot_longer(cols = c(3:ncol(.))) %>%
      dplyr::mutate(
        Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
        Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
      ) %>% 
      dplyr::select(-name) %>%
      tidyr::pivot_wider(id_cols = c(panel, edata_cname, theComparison), names_from = Type, 
                         values_from = value, values_fn = list) %>%
      tidyr::unnest(cols = c(fold_change, p_value))
    
    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh)
      if (is.null(DF)) {
        return(NULL)
      }
      
      # Make heatmap with significance
      hm <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = as.factor(.data[[edata_cname]]), fill = fold_change)) +
        ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ylab("Biomolecule") + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggplot2::labs(fill = "Fold Change") + 
        ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
        ggplot2::geom_point(ggplot2::aes(color = Significance)) + 
        ggplot2::scale_color_manual(values = structure(c("black", NA), .Names = c(attr(DF, "LessThan"), attr(DF, "GreaterThan"))), na.value = NA)
      
    } else {
    
      # Make heatmap
      hm <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = as.factor(.data[[edata_cname]]), fill = fold_change)) +
        ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ylab("Biomolecule") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggplot2::labs(fill = "Fold Change") +
        ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "white")
      
    }

    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        hm <- hm + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      hm <- hm %>% plotly::ggplotly()
    }
    
    return(hm)
  }

  # Build the trelliscope-------------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(fc_hm_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
  } else {
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
  }

}

#' @name trelli_foldchange_volcano
#' 
#' @title Volcano trelliscope building function for fold_change   
#' 
#' @description Specify a plot design and cognostics for the fold_change volcano
#'   trelliscope. Fold change must be grouped by an emeta column, which means
#'   both an omicsData object and statRes are required to make this plot.
#'
#' @param trelliData A trelliscope data object with omicsData and statRes
#'   results. Required.
#' @param comparison The specific comparison to visualize in the fold_change
#'   volcano. See attr(statRes, "comparisons") for the available options.
#'   If all comparisons are desired, the word "all" can be used, which is the
#'   default. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries are 
#'   "biomolecule count", "proportion significant", "proportion significant up",
#'   and "proportion significant down". Default is "biomolecule count". 
#' @param p_value_thresh A value between 0 and 1 to indicate significant
#'   biomolecules for p_value_test. Default is 0.05.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param path The base directory of the trelliscope application. Default is
#'   a temporary directory.
#' @param name The name of the display. Default is Trelliscope.
#' @param build_trelliscope If TRUE, a trelliscope display will be built. Otherwise,
#'   a dataframe will be returned. Default is TRUE. 
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'   and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'   is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'   trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#' 
#' @return No return value, builds a trelliscope display of fold-change volcano plots that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' if (interactive()) {
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
#' trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>% 
#'   trelli_foldchange_volcano(comparison = "all", test_mode = TRUE, test_example = 1:10,
#'                             cognostics = c("biomolecule count", "proportion significant"),
#'                             path = tempdir())
#' 
#' \dontshow{closeAllConnections()}
#' }
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_foldchange_volcano <- function(trelliData,
                                      comparison = "all", 
                                      cognostics = "biomolecule count",
                                      p_value_thresh = 0.05,
                                      ggplot_params = NULL,
                                      interactive = FALSE,
                                      path = tempdir(),
                                      name = "Trelliscope",
                                      build_trelliscope = TRUE,
                                      test_mode = FALSE,
                                      test_example = 1,
                                      single_plot = FALSE,
                                      ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  p_value_test <- trelli_precheck(trelliData = trelliData, 
                                  trelliCheck = c("omics", "stat"),
                                  cognostics = cognostics,
                                  acceptable_cognostics = c("biomolecule count",
                                                            "proportion significant",
                                                            "proportion significant up",
                                                            "proportion significant down"),
                                  ggplot_params = ggplot_params,
                                  interactive = interactive,
                                  build_trelliscope = build_trelliscope,
                                  test_mode = test_mode, 
                                  test_example = test_example,
                                  seqDataCheck = "permissible",
                                  single_plot = single_plot,
                                  p_value_thresh = p_value_thresh)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Ensure that comparison is an acceptable input
  if (!is.character(comparison)) {
    stop("comparison must be a string.")
  }
  
  # Ensure that comparison is from the comparisons lists
  Comparisons <- attr(trelliData$statRes, "comparisons")
  if ("all" %in% comparison) {comparison <- Comparisons}
  
  if (any(comparison %in% Comparisons) == FALSE) {
    stop(paste0(comparison, "is not an acceptable comparison"))
  }
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is an emeta column
  if (!(panel %in% attr(trelliData, "emeta_col"))) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Start builder dataframe
  toBuild <- trelliData$trelliData
  
  # First, build the metas------------------------------------------------------
  
  # Extract fold change and p-value ANOVA (omics) or p-value columns (seqData)
  needed_cols <- colnames(toBuild)[grepl("P_value|Fold_change", colnames(toBuild))]
  if (any(grepl("P_value_G_", needed_cols))) {
    needed_cols <- needed_cols[!grepl("P_value_G", needed_cols)]
  }
  
  # Pull required columns 
  toBuild <- toBuild %>% dplyr::select_at(c(panel, needed_cols))
  toBuild <- unique(toBuild)
  theComparison <- "Comparison"
  
  # Pivot longer to extract comparisons. Then split into p-value and fold change columns
  toBuild <- toBuild %>%
    tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
    dplyr::mutate(
      Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
      Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
    ) %>%
    dplyr::select(-name) %>%
    dplyr::filter(Comparison == comparison) %>%
    dplyr::select(-Comparison) %>%
    tidyr::pivot_wider(id_cols = panel, names_from = Type, values_from = value, values_fn = list) %>%
    tidyr::unnest(cols = c(fold_change, p_value)) %>%
    dplyr::mutate(fold_change = ifelse(p_value <= p_value_thresh, 
                                ifelse(fold_change >= 0, "up", "down"), "NS")) %>%
    dplyr::filter(!is.na(fold_change)) %>%
    dplyr::group_by_at(panel) %>%
    dplyr::summarize(
      `biomolecule count` = dplyr::n(),
      `proportion significant` = sum(fold_change %in% c("up", "down")),
      `proportion significant up` = sum(fold_change == "up"),
      `proportion significant down` = sum(fold_change == "down")
    )
    
  # Add additiona e_meta columns
  emeta <- trelliData$trelliData[,unique(c(panel, attr(trelliData, "emeta_col")))] %>% unique()
  
  # Add emeta cognostics
  toBuild <- dplyr::left_join(toBuild, emeta, by = panel) %>% unique()
  
  # Filter down if test mode --> must be here to get the right selection of panels
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }

  # Make foldchange volcano function--------------------------------------------
  
  fc_volcano_plot_fun <- function(Panel) {

    # Get edata cname
    edata_cname <- attr(trelliData, "edata_col")
    
    # Pull data.frame, and use needed cols from earlier. Then pull the comparisons,
    # fold changes, and p-values. Make it flexible for seqData as well. 
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}}) %>%
      dplyr::select_at(c(panel, edata_cname, needed_cols)) %>%
      unique() %>%
      dplyr::mutate_at(colnames(.)[3:ncol(.)], as.numeric) %>%
      tidyr::pivot_longer(cols = c(3:ncol(.))) %>%
      dplyr::mutate(
        Comparison = gsub("Fold_change_|P_value_|P_value_A_", "", name),
        Type = ifelse(grepl("Fold_change", name), "fold_change", "p_value")
      ) %>% 
      dplyr::filter(Comparison == comparison) %>%
      dplyr::select(-name) %>%
      tidyr::pivot_wider(id_cols = c(panel, edata_cname, theComparison), names_from = Type, 
                         values_from = value, values_fn = list) %>%
      tidyr::unnest(cols = c(fold_change, p_value))
    
    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh) 
      if (is.null(DF)) {return(NULL)}
      
      GreaterThan <- attr(DF, "GreaterThan")
      LessThan <- attr(DF, "LessThan")
      DF <- DF %>% dplyr::filter(!is.na(fold_change))
      
      # Indicate which comparisons should be highlighted
      DF$Significance <- lapply(1:nrow(DF), function(row) {
        if (DF$Significance[row] == LessThan) {
          ifelse(DF$fold_change[row] > 0, 
                 paste(DF$Significance[row], "& High"), 
                 paste(DF$Significance[row], "& Low")
          )
        } else {GreaterThan}
      }) %>% unlist()
      
      # Indicate what symbols depict low, high, and no significance
      LowSig <- paste(LessThan, "& Low")
      HighSig <- paste(LessThan, "& High")
      NoSig <- GreaterThan
      
      # Make volcano plot 
      volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value), color = Significance)) +
        ggplot2::geom_point() + ggplot2::theme_bw() + 
        ggplot2::scale_color_manual(values = structure(c("blue", "red", "black"), 
                                                       .Names = c(LowSig, HighSig, NoSig))) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value") 

    } else {
        volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value))) +
          ggplot2::geom_point() + ggplot2::theme_bw() + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value")
    }
      
    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        volcano <- volcano + eval(parse(text = paste0("ggplot2::", param)))
      }
    }

    # If interactive, pipe to ggplotly
    if (interactive) {
      volcano <- volcano %>% plotly::ggplotly()
    }
    
    return(volcano)
  }

  # Build the trelliscope-------------------------------------------------------

  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(fc_volcano_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- toBuild[test_example[1], "plots"]
    return(singleData$plots)
  } else {
    # If build_trelliscope is true, then build the display. Otherwise, return 
    if (build_trelliscope) {
      trelli_builder_lazy(toBuild, path, name, ...)
    } else {
      return(toBuild)
    }
  }
  
}
