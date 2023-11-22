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
#'   as.trelliData.edata, and grouped by trelli_panel_by. Required.
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
#' @examples
#' \dontrun{
#' 
#' ## Generate trelliData objects using the as.trelliData.edata example code.
#' 
#' # Build the RNA-seq boxplot with an edata file where each panel is a biomolecule. 
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10)
#'    
#' # Build the RNA-seq boxplot where each panel is a sample.
#' # Include all applicable cognostics. Remove points. 
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Sample") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, 
#'                             include_points = FALSE,
#'                             cognostics = c("count", 
#'                                            "mean lcpm", 
#'                                            "median lcpm", 
#'                                            "cv lcpm")
#'                            )
#' 
#' # Build the RNA-seq boxplot with an omicsData object.
#' # Let the panels be biomolecules. Here, grouping information is included.
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Transcript") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10)
#'    
#' # Build the RNA-seq boxplot with an omicsData object. The panel is a biomolecule class,
#' # which is proteins in this case.
#' trelli_panel_by(trelliData = trelliData_seq2, panel = "Gene") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10)
#'     
#' # Build the RNA-seq boxplot with an omicsData and statRes object.
#' # Panel by a biomolecule, and add statistics data to the cognostics
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Transcript") %>%
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10,
#'                             cognostics = c("mean lcpm", "p-value", "fold change"))
#'  
#' # Other options include modifying the ggplot  
#' trelli_panel_by(trelliData = trelliData_seq1, panel = "Transcript") %>% 
#'    trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "xlab('')"))
#' 
#' # Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData_seq4, panel = "Gene") %>% 
#'     trelli_rnaseq_boxplot(interactive = TRUE, test_mode = TRUE, test_example = 1:10)
#' 
#' }
#'
#' @author David Degnan, Lisa Bramer
#'
#' @export
trelli_rnaseq_boxplot <- function(trelliData,
                                     cognostics = c("count", "mean lcpm"),
                                     ggplot_params = NULL,
                                     interactive = FALSE,
                                     include_points = TRUE,
                                     path = .getDownloadsFolder(),
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
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "required",
                  seqText = "Use trelli_abundance_boxplot instead.",
                  p_value_thresh = NULL)
  
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("p-value", "fold change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
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
  
  # Determine if count is biomolecule count or sample count 
  if ("count" %in% cognostics) {
    
    if (attr(trelliData, "panel_by_omics") == get_edata_cname(trelliData$omicsData)) {
      cognostics[which(cognostics == "count")] <- "sample nonzero count"
    } else {
      cognostics[which(cognostics == "count")] <- "biomolecule nonzero count"
    }
    
  }
  
  if (any(c("p-value", "fold change") %in% cognostics) & attr(trelliData, "panel_by_omics") != get_edata_cname(trelliData$omicsData)) {
    message(paste("Please panel by", get_edata_cname(trelliData$omicsData), "to get 'p-values' and 'fold change' as cognostics in the trelliscope display."))
  }
  
  # Make boxplot function-------------------------------------------------------
  
  # First, generate the boxplot function
  box_plot_fun <- function(DF, title) {
    # Add a blank group if no group designation was given
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      DF$Group <- "x"
    }
    
    # Build plot
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = LCPM)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title) +
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
  
  # Create cognostic function---------------------------------------------------
  
  # Second, create function to return cognostics
  box_cog_fun <- function(DF, biomolecule) {
    # Set basic cognostics for ungrouped data or in case when data is not split by fdata_cname
    cog <- list(
      "sample nonzero count" = dplyr::tibble(`Sample Non-Zero Count` = trelliscopejs::cog(sum(DF$Count != 0), desc = "Sample Non-Zero Count")),
      "biomolecule nonzero count" = dplyr::tibble(`Biomolecule Non-Zero Count` = trelliscopejs::cog(sum(DF$Count != 0), desc = "Biomolecule Non-Zero Count")),
      "mean lcpm" = dplyr::tibble(`Mean LCPM` = trelliscopejs::cog(round(mean(DF$LCPM, na.rm = TRUE), 4), desc = "Mean LCPM")), 
      "median lcpm" = dplyr::tibble(`Median LCPM` = trelliscopejs::cog(round(median(DF$LCPM, na.rm = TRUE), 4), desc = "Median LCPM")), 
      "cv lcpm" = dplyr::tibble(`CV LCPM` = trelliscopejs::cog(round((sd(DF$LCPM, na.rm = TRUE) / mean(DF$LCPM, na.rm = T)) * 100, 4), desc = "CV LCPM"))
    )
    
    # If cognostics are any of the cog, then add them 
    if (any(cognostics %in% c("sample nonzero count", "biomolecule nonzero count", "mean lcpm", "median lcpm", "cv lcpm"))) {
      cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% dplyr::tibble()
    } else {
      cog_to_trelli <- NULL
    }
    
    # Get fdata cname and panel_by selection
    fdata_cname <- pmartR::get_fdata_cname(trelliData$omicsData)
    panel_by_choice <- attr(trelliData, "panel_by_omics")
    
    # Additional group cognostics can be added only if group_designation was set and
    # trelli_panel_by is not the fdata_cname
    if (!is.null(attributes(trelliData$omicsData)$group_DF) & fdata_cname != panel_by_choice & !is.null(cog_to_trelli)) {
      # Since the number of groups is unknown, first panel_by the Groups,
      # then calculate all summary statistics, pivot to long format,
      # subset down to requested statistics, switch name to a more specific name,
      # combine group and name, and generate the cognostic tibble
      cogs_to_add <- DF %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(
          "sample nonzero count" = sum(Count != 0), 
          "biomolecule nonzero count" = sum(Count != 0), 
          "mean lcpm" = round(mean(LCPM, na.rm = TRUE), 4),
          "median lcpm" = round(median(LCPM, na.rm = TRUE), 4),
          "cv lcpm" = round((sd(LCPM, na.rm = T) / mean(LCPM, na.rm =T)) * 100, 4),
        ) %>%
        tidyr::pivot_longer(c(`sample nonzero count`, `biomolecule nonzero count`, `mean lcpm`, `median lcpm`, `cv lcpm`)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(name = paste(Group, name)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Group) 
      
      # Add new cognostics 
      cog_to_trelli <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
        quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
      })) %>% dplyr::tibble() 
      
    }
    
    # Add cognostics that only apply when stats data is grouped by edata_cname
    edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
    
    if (!is.null(trelliData$trelliData.stat) && !is.na(attr(trelliData, "panel_by_stat")) && edata_cname == attr(trelliData, "panel_by_stat")) {
      
      # Downselect to only stats 
      stat_cogs <- cognostics[cognostics %in% c("fold change", "p-value")]
      
      if (length(stat_cogs) != 0) {
        # Subset down the dataframe down to group, unnest the dataframe,
        # pivot_longer to comparison, subset columns to requested statistics,
        # switch name to a more specific name
        cogs_to_add <- trelliData$trelliData.stat %>%
          dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == biomolecule) %>%
          dplyr::select(Nested_DF) %>%
          tidyr::unnest(cols = c(Nested_DF)) %>%
          dplyr::rename(
            `p-value` = p_value,
            `fold change` = fold_change
          ) %>%
          dplyr::select(c(Comparison, stat_cogs)) %>%
          tidyr::pivot_longer(stat_cogs) %>%
          dplyr::mutate(
            name = paste(Comparison, name),
            value = round(value, 4)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-Comparison)
        
        # Generate new cognostics
        new_cogs <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
          quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
        })) %>% dplyr::tibble()
        
        # Add new cognostics, removing when it is NULL 
        cog_to_trelli <- cbind(cog_to_trelli, new_cogs) %>% dplyr::tibble()
        
      }
    }
    
    return(cog_to_trelli)
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- trelliData$trelliData.omics[test_example[1], ]
    return(box_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
  } else {
    # If test_mode is on, then just build the required panels
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example, ]
    } else {
      toBuild <- trelliData$trelliData.omics
    }
    
    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = box_plot_fun,
                   cogFUN = box_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...)
    
  }
}