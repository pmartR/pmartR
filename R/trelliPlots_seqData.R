#' @name trelli_rnaseq_boxplot
#'
#' @title Boxplot trelliscope building function for RNA-seq data 
#'
#' @description Specify a boxplot design and cognostics for the RNA-Seq
#'   boxplot trelliscope. Each boxplot will have its own groups as specified by
#'   the first main effect in group_designation. Use "trelli_abundance_boxplot"
#'   for MS/NMR-based omics. 
#'
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by trelli_panel_by. Must be built using 
#'   seqData. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are "count", "mean lcpm", "median lcpm", and "cv lcpm". 
#'   If data are paneled by a biomolecule, the count will be "sample count".
#'   If data are paneled by a sample or a biomolecule class, the count will be "biomolecule count". 
#'   If statRes data is included, "p-value" and "fold change" data per comparisons
#'   may be added. If grouping information is included, only "sample count" and 
#'   "mean lcpm" will be calculated, along with "p-value" and "fold change"
#'   if specified. "p-value" will not be included if paneling a trelliscope
#'   display by a biomolecule class. Default is "sample count" and "mean lcpm".
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "ylim(c(2,20))").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param build_trelliscope The user provided value for building trelliscope
#' @param include_points Add points as a geom_jitter. Default is TRUE.
#' @param path The base directory of the trelliscope application. Default is
#'   Downloads.
#' @param name The name of the display. Default is Trelliscope.
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
#' \dontrun{
#' library(pmartRdata)
#' 
#' trelliData_seq1 <- as.trelliData.edata(e_data = rnaseq_edata,
#'                                       edata_cname = "Transcript",
#'                                       omics_type = "seqData")
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
#' 
#' ## Generate trelliData objects using the as.trelliData.edata example code.
#' 
#' # Build the RNA-seq boxplot with an edata file where each panel is a biomolecule. 
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'    
#' # Build the RNA-seq boxplot where each panel is a sample.
#' # Include all applicable cognostics. Remove points. 
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Sample") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, 
#'                             include_points = FALSE,
#'                             cognostics = c("count", 
#'                                            "mean lcpm", 
#'                                            "median lcpm", 
#'                                            "cv lcpm"),
#'                             path = tempdir()
#'                            )
#' 
#' # Build the RNA-seq boxplot with an omicsData object.
#' # Let the panels be biomolecules. Here, grouping information is included.
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'    
#' # Build the RNA-seq boxplot with an omicsData object. The panel is a biomolecule class,
#' # which is proteins in this case.
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Gene") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'     
#' # Build the RNA-seq boxplot with an omicsData and statRes object.
#' # Panel by a biomolecule, and add statistics data to the cognostics
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Transcript") %>%
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10,
#'      cognostics = c("mean lcpm", "p-value", "fold change"), path = tempdir())
#'  
#' # Other options include modifying the ggplot  
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>%
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10,
#'      ggplot_params = c("ylab('')", "xlab('')"), path = tempdir())
#' 
#' # Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") %>%
#'     trelli_rnaseq_boxplot(interactive = TRUE, test_mode = TRUE, 
#'      test_example = 1:10, path = tempdir())
#' 
#' \dontshow{closeAllConnections()}
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_rnaseq_boxplot <- function(trelliData,
                                  cognostics = c("count", "mean lcpm"),
                                  ggplot_params = NULL,
                                  interactive = FALSE,
                                  build_trelliscope = TRUE, 
                                  include_points = TRUE,
                                  path = tempdir(),
                                  name = "Trelliscope",
                                  test_mode = FALSE,
                                  test_example = 1,
                                  single_plot = FALSE,
                                  ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("count", "mean lcpm", "median lcpm", "cv lcpm", "p-value", "fold change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "required",
                  seqText = "Use trelli_abundance_boxplot instead.",
                  p_value_thresh = NULL)
  
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("p-value", "fold change") %in% cognostics) & is.null(trelliData$statRes)) {
      cognostics <- cognostics[-match(c("p-value", "fold change"), cognostics, nomatch = 0)]
      message(paste("'p-value' and/or 'fold change' were listed as cognostics, but not provided in the trelliData object.",
                    "Did you forget to include a statRes object?")
      )
    }    
  }
  
  # Remove median and cv as cognostics if data is grouped
  if (!inherits(trelliData, "trelliData.edata")) {
    if (any(c("median lcpm", "cv lcpm") %in% cognostics)) {
      cognostics <- cognostics[-match(c("median lcpm", "cv lcpm"), cognostics, nomatch = 0)]
      message("'median lcpm' and 'cv lcpm' are not permitted when groups have been specified.")
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
  
  if (any(c("p-value", "fold change") %in% cognostics) & panel != get_edata_cname(trelliData$omicsData)) {
    message(paste("Please panel by", get_edata_cname(trelliData$omicsData), "to get 'p-value' and 'fold change' as cognostics in the trelliscope display."))
  }

  # Start summary toBuild data.frame
  toBuild <- trelliData$trelliData
    
  # First, generate the cognostics----------------------------------------------
  
  # Add cognostics without groups
  if ("Group" %in% colnames(toBuild) == FALSE) {
    
    toBuild <- toBuild %>% 
      dplyr::group_by_at(panel) %>%
      dplyr::summarize(
        count = sum(Count != 0),
        `mean lcpm` = mean(LCPM, na.rm = T),
        `median lcpm` = median(LCPM, na.rm = T),
        `cv lcpm` = sd(LCPM, na.rm = T) / `mean lcpm` * 100
      )

    # Now, select only the requested cognostics 
    toBuild <- toBuild %>%
      dplyr::select(dplyr::all_of(c(panel, cognostics)))
    
  } else {
    
    # Make a group variable name just to make group_by_at work properly
    theGroup <- "Group"
    
    # Add cognostics per group
    toBuild <- toBuild %>% 
      dplyr::group_by_at(c(panel, theGroup)) %>%
      dplyr::summarise(
        count = sum(Count != 0),
        `mean lcpm` = mean(LCPM, na.rm = T),
      ) %>%
      tidyr::pivot_wider(id_cols = panel, names_from = Group, values_from = c(count, `mean lcpm`), names_sep = " ")
    
    # Now remove unwanted cognostics
    if ("mean lcpm" %in% cognostics == FALSE) {
      toBuild <- toBuild[,grepl("mean lcpm", colnames(toBuild)) == FALSE]
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
    
  }
  
  # Add p-values and fold changes 
  if ("p-value" %in% cognostics) {
    anova_cols <- colnames(trelliData$trelliData)[grepl("P_value_", colnames(trelliData$trelliData))]
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
  
  # Make boxplot function-------------------------------------------------------
  
  # First, generate the boxplot function
  box_plot_fun <- function(Panel) {
    
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}})
    
    # Add a blank group if no group designation was given
    if (!("Group" %in% colnames(DF))) {
      DF$Group <- "x"
    }
    
    # Build plot
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = LCPM)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylab("Log Counts per Million (LCPM)") 
    
    # Add include_points
    if (include_points) {
      boxplot <- boxplot + ggplot2::geom_jitter(height = 0, width = 0.25)
    }
    
    # Remove x axis if no groups
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      boxplot <- boxplot + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      )
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

#' @name trelli_rnaseq_histogram
#'
#' @title Histogram trelliscope building function for RNA-Seq data
#'
#' @description Specify a plot design and cognostics for the abundance histogram
#'   trelliscope. Main_effects grouping are ignored. Data must be grouped by
#'   edata_cname. For MS/NMR data, use "trelli_abundance_histogram". 
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by edata_cname in trelli_panel_by.
#'   Must be built using seqData. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are "sample count", "mean lcpm", "median lcpm", "cv lcpm", 
#'   and "skew lcpm". All are included by default. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "ylim(c(1,2))").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
#' @param build_trelliscope The user provided value for building trelliscope
#' @param path The base directory of the trelliscope application. Default is
#'   Downloads.
#' @param name The name of the display. Default is Trelliscope.
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
#' \dontrun{
#' library(pmartRdata)
#' 
#' trelliData_seq1 <- as.trelliData.edata(e_data = rnaseq_edata, 
#'                                       edata_cname = "Transcript",
#'                                       omics_type = "seqData")
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
#' 
#' # Build the RNA-seq histogram with an edata file. 
#' # Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'    trelli_rnaseq_histogram(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' 
#' # Build the RNA-seq histogram with an omicsData object. 
#' # Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'    trelli_rnaseq_histogram(test_mode = TRUE, test_example = 1:10, path = tempdir())
#'     
#' # Build the RNA-seq histogram with an omicsData and statRes object. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Transcript") %>%
#'    trelli_rnaseq_histogram(test_mode = TRUE, test_example = 1:10, 
#'      cognostics = "sample count", path = tempdir())
#'    
#' # Users can modify the plotting function with ggplot parameters and interactivity, 
#' # and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'    trelli_rnaseq_histogram(test_mode = TRUE, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "xlab('')"), interactive = TRUE,
#'      cognostics = c("mean lcpm", "median lcpm"), path = tempdir())  
#' 
#' \dontshow{closeAllConnections()}
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_rnaseq_histogram <- function(trelliData,
                                    cognostics = c("sample count", "mean lcpm", "median lcpm", "cv lcpm", "skew lcpm"),
                                    ggplot_params = NULL,
                                    interactive = FALSE,
                                    build_trelliscope = TRUE,
                                    path = tempdir(),
                                    name = "Trelliscope",
                                    test_mode = FALSE,
                                    test_example = 1,
                                    single_plot = FALSE,
                                    ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("sample count", "mean lcpm", "median lcpm", "cv lcpm", "skew lcpm"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "required",
                  seqText = "Use trelli_abundance_histogram instead.",
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
  
  # Convert sample count to sample nonzero count
  if ("sample count" %in% cognostics) {
    cognostics[grepl("sample count", cognostics, fixed = T)] <- "sample nonzero count"
  }
  
  # Start summary toBuild data.frame
  toBuild <- trelliData$trelliData
  
  # First, add any missing cognostics-------------------------------------------
  
  toBuild <- toBuild %>% 
    dplyr::group_by_at(panel) %>%
    dplyr::summarize(
      `sample nonzero count` = sum(Count != 0),
      `mean lcpm` = mean(LCPM, na.rm = T),
      `median lcpm` = median(LCPM, na.rm = T),
      `cv lcpm` = sd(LCPM, na.rm = T) / `mean lcpm` * 100,
      `skew lcpm` = e1071::skewness(LCPM, na.rm = T)
    )
  
  # Add emeta columns
  if (!is.null(attr(trelliData, "emeta_col"))) {
    
    # Pull emeta uniqued columns that should have been prepped in the pivot_longer section
    emeta <- trelliData$trelliData[,c(panel, attr(trelliData, "emeta_col"))] %>% unique()
    
    # Add emeta cognostics
    toBuild <- dplyr::left_join(toBuild, emeta, by = panel)
    
  }
  
  # Now, select only the requested cognostics 
  toBuild <- toBuild %>%
    dplyr::select(dplyr::all_of(c(panel, cognostics)))
  
  # Filter down if test mode
  if (test_mode) {
    toBuild <- toBuild[test_example,]
  }
  
  # Make histogram function-----------------------------------------------------
  
  hist_plot_fun <- function(Panel) {
    
    DF <- dplyr::filter(trelliData$trelliData, Panel == {{Panel}})
    
    # Build plot
    histogram <- ggplot2::ggplot(DF, ggplot2::aes(x = LCPM)) +
      ggplot2::geom_histogram(bins = 10, fill = "steelblue", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::xlab("Log Counts Per Million") +
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

#' @name trelli_rnaseq_heatmap
#'
#' @title Heatmap trelliscope building function for RNA-seq data
#'
#' @description Specify a plot design and cognostics for the RNA-seq heatmap
#'   trelliscope. Data must be grouped by an e_meta column. Main_effects order
#'   the y-variables. All statRes data is ignored. 
#'   For MS/NMR data, use "trelli_abundance_heatmap".
#'
#' @param trelliData A trelliscope data object made by as.trelliData, and
#'   grouped by an emeta variable. Must be built using seqData. Required.
#' @param cognostics A vector of cognostic options. Defaults are "sample count", 
#'   "mean lcpm" and "biomolecule count". "sample count" and "mean lcpm"
#'   are reported per group, and "biomolecule count" is the total number of biomolecules
#'   in the biomolecule class (e_meta column).
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly. Default is FALSE.
#' @param build_trelliscope The user provided value for building trelliscope
#' @param path The base directory of the trelliscope application. Default is
#'   Downloads.
#' @param name The name of the display. Default is Trelliscope
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
#' \dontrun{
#' library(pmartRdata)
#' 
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
#' trelliData_seq4 <- as.trelliData(omicsData = omicsData_seq, statRes = statRes_seq)
#' 
#' # Build the RNA-seq heatmap with an omicsData object with emeta variables. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Gene") %>% 
#'    trelli_rnaseq_heatmap(test_mode = TRUE, test_example = c(1532, 1905, 6134), path = tempdir())
#'    
#' # Users can modify the plotting function with ggplot parameters and interactivity, 
#' # and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") %>% 
#'    trelli_rnaseq_heatmap(test_mode = TRUE, test_example = c(1532, 1905, 6134), 
#'      ggplot_params = c("ylab('')", "xlab('')"), 
#'      interactive = TRUE, cognostics = c("biomolecule count"), path = tempdir())  
#' 
#' \dontshow{closeAllConnections()}
#' }
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_rnaseq_heatmap <- function(trelliData,
                                  cognostics = c("sample count", "mean lcpm", "biomolecule count"),
                                  ggplot_params = NULL,
                                  interactive = FALSE,
                                  build_trelliscope = TRUE,
                                  path = tempdir(),
                                  name = "Trelliscope",
                                  test_mode = FALSE,
                                  test_example = 1,
                                  single_plot = FALSE,
                                  ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "omics",
                  cognostics = cognostics,
                  acceptable_cognostics = c("sample count", "mean lcpm", "biomolecule count"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "required",
                  seqText = "Use trelli_abundance_heatmap instead.",
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
  
  # Start summary toBuild data.frame
  toBuild <- trelliData$trelliData
  
  # First, add any missing cognostics-------------------------------------------
  
  # Make a group variable name just to make group_by_at work properly
  theGroup <- "Group"
  
  # Add cognostics per group
  toBuild <- toBuild %>% 
    dplyr::group_by_at(c(panel, theGroup)) %>%
    dplyr::summarise(
      count = sum(Count != 0),
      `mean lcpm` = mean(LCPM, na.rm = T)
    ) %>%
    tidyr::pivot_wider(id_cols = panel, names_from = Group, values_from = c(count, `mean lcpm`), names_sep = " ")
  
  # Remove unwanted cognostics 
  if ("sample count" %in% cognostics == FALSE) {
    toBuild <- toBuild[,!grepl("count", colnames(toBuild))]
  }
  if ("mean lcpm" %in% cognostics == FALSE) {
    toBuild <- toBuild[,!grepl("mean lcpm", colnames(toBuild))]
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
  
  # Filter down if test mode
  if (test_mode) {
    toBuild <- toBuild[test_example,]
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
      dplyr::mutate(`Z-Score` = (LCPM - mean(LCPM, na.rm = T)) / sd(LCPM, na.rm = T))
    
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

#' @name trelli_rnaseq_nonzero_bar
#'
#' @title Bar chart trelliscope building function for Non-Zero counts in RNA-seq data
#'
#' @description Specify a plot design and cognostics for the Non-Zero barchart
#'    trelliscope. Non-Zeroes are displayed per panel_by variable. Main_effects
#'    data is used to split samples when applicable. For MS/NMR data, use 
#'    "trelli missingness bar". 
#'
#' @param trelliData A trelliscope data object made by as.trelliData.edata or
#'    as.trelliData. Must be built using seqData. Required.
#' @param cognostics A vector of cognostic options for each plot. Defaults are "total count",
#'    "non-zero count", and "non-zero proportion". If grouping
#'    data is included, all cognostics will be reported per group. If the 
#'    trelliData is paneled by a biomolecule, the counts and proportion we be 
#'    samples. If paneled by a sample or biomolecule class, the counts and proportions
#'    will be biomolecules.
#' @param proportion A logical to determine whether plots should display counts
#'    or proportions. Default is TRUE.
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'    the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'    Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'    interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'    now). Default is FALSE.
#' @param build_trelliscope The user provided value for building trelliscope
#' @param path The base directory of the trelliscope application. Default is
#'    Downloads.
#' @param name The name of the display. Default is Trelliscope.
#' @param test_mode A logical to return a smaller trelliscope to confirm plot
#'    and design. Default is FALSE.
#' @param test_example A vector of plot indices to return for test_mode. Default
#'    is 1.
#' @param single_plot A TRUE/FALSE to indicate whether 1 plot (not a
#'    trelliscope) should be returned. Default is FALSE.
#' @param ... Additional arguments to be passed on to the trelli builder
#'   
#' @return No return value, builds a trelliscope display of bar charts that is stored in `path`
#' 
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \dontrun{
#' library(pmartRdata)
#' 
#' trelliData_seq1 <- as.trelliData.edata(e_data = rnaseq_edata, 
#'                                       edata_cname = "Transcript",
#'                                       omics_type = "seqData")
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
#' 
#' # Build the non-zero bar plot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'   trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Sample") %>% 
#'   trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:10, 
#'    cognostics = "non-zero proportion", path = tempdir())
#' 
#' # Build the non-zero bar plot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'   trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:10, path = tempdir())
#' 
#' # Build the non-zero bar plot with a statRes object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData_seq3, panel = "Transcript") %>%
#'   trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:10,
#'                          cognostics = c("non-zero proportion"), path = tempdir())
#' 
#' # Build the non-zero bar plot with an omicsData and statRes object. 
#' # Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") %>%
#'   trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:10, path = tempdir()) 
#' 
#' # Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'    trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:5, 
#'      interactive = TRUE, path = tempdir())
#'    
#' # Or visualize only count data 
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'    trelli_rnaseq_nonzero_bar(test_mode = TRUE, test_example = 1:5, 
#'      cognostics = "non-zero count", proportion = FALSE, path = tempdir())
#' 
#' \dontshow{closeAllConnections()}   
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_rnaseq_nonzero_bar <- function(trelliData,
                                   cognostics = c("total count", "non-zero count", "non-zero proportion"),
                                   proportion = TRUE,
                                   ggplot_params = NULL,
                                   interactive = FALSE,
                                   build_trelliscope = TRUE,
                                   path = tempdir(),
                                   name = "Trelliscope",
                                   test_mode = FALSE,
                                   test_example = 1,
                                   single_plot = FALSE,
                                   ...) {
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("either"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("total count", "non-zero count", "non-zero proportion"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  build_trelliscope = build_trelliscope,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "required",
                  seqText = "Use trelli_missingness_bar instead.",
                  p_value_thresh = NULL)
  
  # Extract panel column
  panel <- attr(trelliData, "panel_by_col")
  
  # Check that proportion is a non NA logical
  if (!is.logical(proportion) | is.na(proportion)) {
    stop("proportion must be a TRUE or FALSE.")
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
  
  # Extract non-zero counts 
  if (!is.null(trelliData$omicsData)) {
    
    if ("Group" %in% colnames(trelliData$trelliData) == FALSE) {
      toBuild <- toBuild %>% 
        dplyr::group_by_at(panel) %>%
        dplyr::summarize(
          `total count` = dplyr::n(),
          `non-zero count` = sum(Count != 0),
          `non-zero proportion` = round(`non-zero count` / `total count`, 4)
        )
    } else {
      toBuild <- toBuild %>%
        dplyr::group_by_at(c(panel, theGroupName)) %>%
        dplyr::summarize(
          `total count` = dplyr::n(),
          `non-zero count` = sum(Count != 0),
          `non-zero proportion` = round(`non-zero count` / `total count`, 4)
        ) %>% 
        tidyr::pivot_wider(id_cols = panel, names_from = theGroupName, names_sep = " ",
                           values_from = c("total count", "non-zero count", "non-zero proportion"))
    }
    
  } else {
    
    # Get count columns 
    cols2rename <- colnames(toBuild)[grepl("Count_", colnames(toBuild))]
    
    # Make totals dataframe
    totalsDF <- attr(trelliData$statRes, "group_DF")$Group %>% 
      table(dnn = "Group") %>% 
      data.frame() %>% 
      dplyr::rename(`total count` = Freq)
    
    # Extract panel and non-zero counts. Format with group information. Add
    # totals and non-zero proportion. 
    toBuild <- toBuild %>%
      dplyr::select_at(c(panel, cols2rename)) %>%
      tidyr::pivot_longer(cols = c(2:ncol(.))) %>%
      dplyr::mutate(name = gsub("Count_", "", name, fixed = T)) %>%
      dplyr::rename(Group = name, `non-zero count` = value) %>%
      dplyr::left_join(totalsDF, by = "Group") %>%
      dplyr::mutate(`non-zero proportion` = round(`non-zero count` / `total count`, 4)) %>%
      tidyr::pivot_wider(id_cols = panel, names_from = theGroupName, names_sep = " ",
                         values_from = c("total count", "non-zero count", "non-zero proportion"))
    
  }
  
  # Add emeta columns
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
  
  # Second, write the plotting function-----------------------------------------
  
  # First, generate the bar plot function
  zero_bar_plot_fun <- function(Panel) {
    
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
        dplyr::rename(Group = name, `non-zero count` = value) %>%
        dplyr::left_join(totalsDF, by = "Group") %>%
        dplyr::mutate(
          `non-zero proportion` = round(`non-zero count` / `total count`, 4),
          `zero count` = `total count` - `non-zero count`,
          `zero proportion` = round(`zero count` / `total count`, 4)
        ) %>%
        tidyr::pivot_longer(cols = c(3:7))
      
    } else {
      
      # Split out group names from summary toBuild table, pivot_wider to add zero
      # columns, an then pivot longer 
      MissPlotDF <- DF %>% 
        dplyr::group_by_at(c(panel, theGroupName)) %>%
        tidyr::nest() %>%
        dplyr::mutate(
          `total count` = purrr::map_int(data, nrow),
          `non-zero count` = purrr::map_int(data, function(x) {sum(x$Count != 0)}),
          `zero count` = `total count` - `non-zero count`,
          `non-zero proportion` = round(`non-zero count` / `total count`, 4),
          `zero proportion` = round(`zero count` / `total count`, 4)
        ) %>%
        dplyr::select(-data) %>%
        tidyr::pivot_longer(cols = c(3:7))
      
    }
    
    # Subset based on count or proportion
    if (proportion) {
      MissPlotDF <- MissPlotDF %>%
        dplyr::filter(name %in% c("non-zero proportion", "zero proportion")) %>%
        dplyr::mutate(
          name = factor(ifelse(name == "non-zero proportion", "Non-Zero", "Zero"), levels = c("Zero", "Non-Zero"))
        )
      ylab <- "Proportion"
    } else {
      MissPlotDF <- MissPlotDF %>%
        dplyr::filter(name %in% c("non-zero count", "zero count")) %>%
        dplyr::mutate( 
          name = factor(ifelse(name == "non-zero count", "Non-Zero", "Zero"), levels = c("Zero", "Non-Zero"))
        )
      ylab <- "Count"
    }
    
    # Build plot
    zero_bar <- ggplot2::ggplot(MissPlotDF, ggplot2::aes(x = Group, y = value, fill = name)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::ylab(ylab) +
      ggplot2::scale_fill_manual(values = c("Non-Zero" = "steelblue", "Zero" = "black")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank()
      )
    
    # Remove x axis if no groups
    if ("Group" %in% colnames(trelliData$trelliData) == FALSE & stats_mode == FALSE) {
      zero_bar <- zero_bar + ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      )
    }
    
    # Add additional parameters
    if (!is.null(ggplot_params)) {
      for (param in ggplot_params) {
        zero_bar <- zero_bar + eval(parse(text = paste0("ggplot2::", param)))
      }
    }
    
    # If interactive, pipe to ggplotly
    if (interactive) {
      zero_bar <- zero_bar %>% plotly::ggplotly()
    }
    
    return(zero_bar)
  }
  
  
  # Build trelliscope display---------------------------------------------------
  
  # Add a panel column for plotting
  trelliData$trelliData$Panel <- trelliData$trelliData[[panel]]
  toBuild$Panel <- toBuild[[panel]]
  
  # Add plots and remove that panel column
  toBuild <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plots = trelliscope::panel_lazy(zero_bar_plot_fun)) %>%
    dplyr::select(-Panel)
  
  # Now, select only the requested cognostics 
  selected <- colnames(toBuild)[lapply(colnames(toBuild), function(x) {grepl("total count|non-zero count|non-zero proportion", x)}) %>% unlist()]
  selected <- c(selected, panel, "plots")
  
  toBuild <- toBuild %>%
    dplyr::select(dplyr::all_of(selected))
  
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