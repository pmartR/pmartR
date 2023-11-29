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
    if (is.null(trelliData$trelliData.omics)) {
      stop("trelliData must have omicsData for this plotting function.")
    }

    # Ensure that test_example is in the range of possibilities
    if (test_mode) {
      if (max(test_example) > nrow(trelliData$trelliData.omics)) {
        stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.omics)))
      }
    }
  }

  # Check that statRes data exists
  if ("stat" %in% trelliCheck) {
    # Assert that trelliData has statRes
    if (is.null(trelliData$trelliData.stat)) {
      stop("trelliData must have statRes for this plotting function.")
    }

    # Ensure that test_example is in the range of possibilities
    if (test_mode) {
      if (max(test_example) > nrow(trelliData$trelliData.stat)) {
        stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.stat)))
      }
    }
  }

  ######################
  ## COGNOSTIC CHECKS ##
  ######################

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

# A quick cognostic function
quick_cog <- function(name, value) {
  dplyr::tibble(!!dplyr::sym(name) := trelliscopejs::cog(value, desc = name))
}

# This function builds all trelliscopes.
trelli_builder <- function(toBuild, cognostics, plotFUN, cogFUN, path, name, remove_nestedDF, ...) {
  
  # Remove any blank names 
  if ("" %in% unlist(toBuild[,1])) {
    message(paste("Removing", length(sum("" %in% toBuild[, 1])), "blank biomolecule names."))
    toBuild <- toBuild[toBuild[, 1] != "", ]
  }

  if (nrow(toBuild) == 0) {
    stop("No data to build trelliscope with.")
  }
  
  # Plots will always be included 
  preLaunch <- toBuild %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[, 1])), plotFUN)
    ) 
  
  # Cognostics are optional
  if (!is.null(cognostics)) {
    
    preLaunch <- preLaunch %>%
      dplyr::mutate(
        cog = trelliscopejs::map2_cog(Nested_DF, as.character(unlist(toBuild[, 1])), cogFUN)
      )
    
  }
  
  # Some dataframes are only one row, which ends up as cognostics. Removing the nestedDF
  # is recommended
  if (remove_nestedDF) {preLaunch <- preLaunch %>% dplyr::select(-Nested_DF)}
  
  # Finally, build the diplay  
  preLaunch %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = TRUE, ...) 
    
}

# Get downloads folder
.getDownloadsFolder <- function(.test_mode = FALSE) {
  if (Sys.info()['sysname'] == "Windows" | .test_mode) {
    folder <- dirname("~")
    folder <- file.path(folder, "Downloads")
    return(folder)
  } else {
    folder <- path.expand("~")
    folder <- file.path(folder, "Downloads")
    folder <- paste0(folder, .Platform$file.sep)
    return(folder)
  }
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
#' # Build the abundance boxplot with an omicsData object. The panel is a biomolecule class,
#' # which is proteins in this case.
#' trelli_panel_by(trelliData = trelliData2, panel = "RazorProtein") %>% 
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
#' trelli_panel_by(trelliData = trelliData4, panel = "RazorProtein") %>% 
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
                  acceptable_cognostics = c("count", "mean abundance", "median abundance", "cv abundance", "anova p-value", "fold change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_boxplot instead.",
                  p_value_thresh = NULL)
  
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("anova p-value", "fold change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[-match(c("anova p-value", "fold change"), cognostics, nomatch = 0)]
      message(paste("'anova p-value' and/or 'fold change' were listed as cognostics, but not provided in the trelliData object.",
                    "Did you forget to include a statRes object?")
             )
    }    
  }
  
  # Remove median and cv as cognostics if data is grouped
  if (!inherits(trelliData, "trelliData.edata")) {
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
  
  # Determine if count is biomolecule count or sample count 
  if ("count" %in% cognostics) {
    
    if (attr(trelliData, "panel_by_omics") == get_edata_cname(trelliData$omicsData)) {
      cognostics[which(cognostics == "count")] <- "sample count"
    } else {
      cognostics[which(cognostics == "count")] <- "biomolecule count"
    }
    
  }
  
  if (any(c("anova p-value", "fold change") %in% cognostics) & attr(trelliData, "panel_by_omics") != get_edata_cname(trelliData$omicsData)) {
    message(paste("Please panel by", get_edata_cname(trelliData$omicsData), "to get 'anova p-value' and 'fold change' as cognostics in the trelliscope display."))
  }
  
  # Make boxplot function-------------------------------------------------------

  # First, generate the boxplot function
  box_plot_fun <- function(DF, title) {
    # Add a blank group if no group designation was given
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      DF$Group <- "x"
    }

    # Build plot
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = Abundance)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::ylab(paste(attr(trelliData$omicsData, "data_info")$data_scale, "Abundance"))

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
      "sample count" = dplyr::tibble(`Sample Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Sample Count")),
      "biomolecule count" = dplyr::tibble(`Biomolecule Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Biomolecule Count")),
      "mean abundance" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = TRUE), 4), desc = "Mean Abundance")), 
      "median abundance" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = TRUE), 4), desc = "Median Abundance")), 
      "cv abundance" = dplyr::tibble(`CV Abundance` = trelliscopejs::cog(round((sd(DF$Abundance, na.rm = TRUE) / mean(DF$Abundance, na.rm = T)) * 100, 4), desc = "CV Abundance"))
    )
    
    # If cognostics are any of the cog, then add them 
    if (any(cognostics %in% c("sample count", "biomolecule count", "mean abundance", "median abundance", "cv abundance"))) {
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
          "sample count" = sum(!is.na(Abundance)), 
          "biomolecule count" = sum(!is.na(Abundance)), 
          "mean abundance" = round(mean(Abundance, na.rm = TRUE), 4),
          "median abundance" = round(median(Abundance, na.rm = TRUE), 4),
          "cv abundance" = round((sd(Abundance, na.rm = T) / mean(Abundance, na.rm =T)) * 100, 4),
        ) %>%
        tidyr::pivot_longer(c(`sample count`, `biomolecule count`, `mean abundance`, `median abundance`, `cv abundance`)) %>%
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
      stat_cogs <- cognostics[cognostics %in% c("fold change", "anova p-value")]
      
      if (length(stat_cogs) != 0) {
        # Subset down the dataframe down to group, unnest the dataframe,
        # pivot_longer to comparison, subset columns to requested statistics,
        # switch name to a more specific name
        cogs_to_add <- trelliData$trelliData.stat %>%
          dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == biomolecule) %>%
          dplyr::select(Nested_DF) %>%
          tidyr::unnest(cols = c(Nested_DF)) %>%
          dplyr::rename(
            `anova p-value` = p_value_anova,
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
#'   are "sample count", "mean abundance", "median abundance", "cv abundance", 
#'   and "skew abundance". All are included by default. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "ylim(c(1,2))").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
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
                  acceptable_cognostics = c("sample count", "mean abundance", "median abundance", "cv abundance", "skew abundance"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_histogram instead.",
                  p_value_thresh = NULL)

  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
  if (edata_cname != attr(trelliData, "panel_by_omics")) {
    stop("trelliData must be grouped by edata_cname.")
  }

  # Make histogram function-----------------------------------------------------

  # First, generate the histogram function
  hist_plot_fun <- function(DF, title) {
    # Remove NAs
    DF <- DF[!is.na(DF$Abundance), ]

    # Build plot
    histogram <- ggplot2::ggplot(DF, ggplot2::aes(x = Abundance)) +
      ggplot2::geom_histogram(bins = 10, fill = "steelblue", color = "black") +
      ggplot2::ggtitle(title) +
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
  
  # Create cognostic function---------------------------------------------------

  # Second, create function to return cognostics
  hist_cog_fun <- function(DF, biomolecule) {
    # Set basic cognostics for ungrouped data or in case when data is not split by fdata_cname
    cog <- list(
      "sample count" = dplyr::tibble(`Sample Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Sample Count")),
      "mean abundance" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = TRUE), 4), desc = "Mean Abundance")), 
      "median abundance" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = TRUE), 4), desc = "Median Abundance")), 
      "cv abundance" = dplyr::tibble(`CV Abundance` = trelliscopejs::cog(round((sd(DF$Abundance, na.rm = TRUE) / mean(DF$Abundance, na.rm = T)) * 100, 4), desc = "CV Abundance")), 
      "skew abundance" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = TRUE), 4), desc= "Skew Abundance"))
    )
    
    # If cognostics are any of the cog, then add them 
    cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% dplyr::tibble()
    
    return(cog_to_trelli)
  }

  # Build trelliscope display---------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- trelliData$trelliData.omics[test_example[1], ]
    return(hist_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
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
                   plotFUN = hist_plot_fun,
                   cogFUN = hist_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...)
    
  }
}

#' @name trelli_abundance_heatmap
#'
#' @title Heatmap trelliscope building function for abundance data
#'
#' @description Specify a plot design and cognostics for the abundance heatmap
#'   trelliscope. Data must be grouped by an e_meta column. Main_effects order
#'   the y-variables. All statRes data is ignored. For RNA-Seq data, use "trelli_rnaseq_heatmap".
#'
#' @param trelliData A trelliscope data object made by as.trelliData, and
#'   grouped by an emeta variable. Required.
#' @param cognostics A vector of cognostic options. Defaults are "sample count", 
#'   "mean abundance" and "biomolecule count". "sample count" and "mean abundance"
#'   are reported per group, and "biomolecule count" is the total number of biomolecules
#'   in the biomolecule class (e_meta column).
#' @param ggplot_params An optional vector of strings of ggplot parameters to
#'   the backend ggplot function. For example, c("ylab('')", "xlab('')").
#'   Default is NULL.
#' @param interactive A logical argument indicating whether the plots should be
#'   interactive or not. Interactive plots are ggplots piped to ggplotly (for
#'   now). Default is FALSE.
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
                  acceptable_cognostics = c("sample count", "mean abundance", "biomolecule count"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_heatmap instead.",
                  p_value_thresh = NULL)
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is grouped by an e_meta variable
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }

  # If no group designation, set cognostics to NULL.
  if (is.null(attributes(trelliData$omicsData)$group_DF)) {
    cognostics <- NULL
  }
  
  # Get the edata variable name
  edata_cname <- get_edata_cname(trelliData$omicsData)
  
  # Make heatmap function-------------------------------------------------------

  # First, generate the heatmap function
  hm_plot_fun <- function(DF, title) {
    # Get fdata_cname
    fdata_cname <- get_fdata_cname(trelliData$omicsData)

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
    
    # Build plot: this should be edata_cname
    hm <- ggplot2::ggplot(DF, ggplot2::aes(x = as.factor(.data[[edata_cname]]), y = .data[[fdata_cname]], fill = Abundance)) +
      ggplot2::geom_tile() +
      ggplot2::theme_bw() +
      ggplot2::ylab("Sample") +
      ggplot2::xlab("Biomolecule") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
      ) +
      ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
      ggplot2::ggtitle(title)

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

  # Create cognostic function---------------------------------------------------

  hm_cog_fun <- function(DF, emeta_var) {
    # Subset down the dataframe down to group, unnest the dataframe,
    # pivot_longer to comparison, subset columns to requested statistics,
    # switch name to a more specific name
    cogs_to_add <- DF %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        "sample count" = sum(!is.na(Abundance)), 
        "mean abundance" = round(mean(Abundance, na.rm = TRUE), 4),
      ) %>%
      tidyr::pivot_longer(c(`sample count`, `mean abundance`)) %>%
      dplyr::filter(name %in% cognostics) %>%
      dplyr::mutate(name = paste(Group, name)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Group) 
    
    if ("biomolecule count" %in% cognostics) {
      cogs_to_add <- rbind(
        cogs_to_add, 
        data.frame(
          name = "Biomolecule Count",
          value = DF[[edata_cname]] %>% unique() %>% length()
        )
      )
      
    }
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
      quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
    })) %>% dplyr::tibble()

    return(cog_to_trelli)
  }

  # Build trelliscope display---------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- trelliData$trelliData.omics[test_example[1], ]
    return(hm_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
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
                   plotFUN = hm_plot_fun,
                   cogFUN = hm_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...) 
    
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
                                   path = .getDownloadsFolder(),
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
                  acceptable_cognostics = c("total count", "observed count", "observed proportion", "g-test p-value"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  seqDataCheck = "no",
                  seqText = "Use trelli_rnaseq_nonzero_bar instead.",
                  p_value_thresh = NULL)
  
  # Check that proportion is a non NA logical
  if (!is.logical(proportion) | is.na(proportion)) {
    stop("proportion must be a TRUE or FALSE.")
  }
  
  # Test whether statRes data is included to use g-test p-value
  if ("g-test p-value" %in% cognostics) {
    if (is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[cognostics != "g-test p-value"]
      message("'g-test p-value' can only be included if stats data (statRes) is included")
    } else if (is.na(attr(trelliData, "panel_by_stat")) || get_edata_cname(trelliData$statRes) != attr(trelliData, "panel_by_stat")) {
      cognostics <- cognostics[cognostics != "g-test p-value"]
      message("'g-test p-value' can only be included if the data has been paneled by the biomolecule column 'edata_cname'")
    }
  }
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Determine if the trelliData is paneled by the edata column
  if (is.na(attr(trelliData, "panel_by_omics"))) {
    paneled_by_edata <- attr(trelliData, "panel_by_stat") == get_edata_cname(trelliData$statRes)
  } else {
    paneled_by_edata <- attr(trelliData, "panel_by_omics") == get_edata_cname(trelliData$omicsData)
  }
  
  # Generate a function to make missingness dataframes--------------------------
  get_missing_DF <- function(DF) {
    # If there is no omics data, use statRes
    if (is.null(trelliData$omicsData)) {
      # Add to total counts
      Missingness <- data.table::data.table(
        Group = gsub("Count_", "", colnames(DF)),
        `Present Count` = unlist(DF)
      ) %>%
        merge(totalCounts, by = "Group") %>%
        dplyr::mutate(
          `Absent Count` = Total - `Present Count`,
          `Absent Proportion` = round(`Absent Count` / Total, 4),
          `Present Proportion` = round(`Present Count` / Total, 4)
        ) %>%
        dplyr::select(-Total) %>%
        dplyr::mutate(Group, `Absent Count`, `Present Count`, `Absent Proportion`, `Present Proportion`)
    } else {
      # Add a blank group if no group designation was given
      if (is.null(attributes(trelliData$omicsData)$group_DF)) {
        DF$Group <- "x"
      }

      # Create missingness data.frame
      Missingness <- DF %>%
        dplyr::group_by(Group) %>%
        dplyr::summarise(
          `Absent Count` = sum(is.na(Abundance)),
          `Present Count` = sum(!is.na(Abundance)),
          `Absent Proportion` = round(`Absent Count` / sum(c(`Absent Count`, `Present Count`)), 4),
          `Present Proportion` = round(`Present Count` / sum(c(`Absent Count`, `Present Count`)), 4)
        )
    }

    return(Missingness)
  }

  # Make missingness bar function-----------------------------------------------

  # First, generate the boxplot function
  missing_bar_plot_fun <- function(DF, title) {
    # Get missingness dataframe
    Missingness <- get_missing_DF(DF)

    # Subset based on count or proportion
    if (proportion) {
      MissPlotDF <- Missingness %>%
        dplyr::select(c(Group, `Absent Proportion`, `Present Proportion`)) %>%
        dplyr::rename(Absent = `Absent Proportion`, Present = `Present Proportion`) %>%
        tidyr::pivot_longer(c(Absent, Present))
      ylab <- "Proportion"
    } else {
      MissPlotDF <- Missingness %>%
        dplyr::select(c(Group, `Absent Count`, `Present Count`)) %>%
        dplyr::rename(Absent = `Absent Count`, Present = `Present Count`) %>%
        tidyr::pivot_longer(c(Absent, Present))
      ylab <- "Count"
    }

    # Build plot
    missing_bar <- ggplot2::ggplot(MissPlotDF, ggplot2::aes(x = Group, y = value, fill = name)) +
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(title) +
      ggplot2::ylab(ylab) +
      ggplot2::scale_fill_manual(values = c("black", "steelblue")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank()
      )

    # Remove x axis if no groups
    if (is.null(attributes(trelliData$omicsData)$group_DF) & stats_mode == FALSE) {
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

  # Create the cognostic function-----------------------------------------------

  # Next, generate the cognostic function
  missing_bar_cog_fun <- function(DF, GroupingVar) {
    # Get missingness dataframe
    Missingness <- get_missing_DF(DF)

    # Build cognostics
    Miss_Cog <- Missingness %>%
      dplyr::mutate(`total count` = `Present Count` + `Absent Count`) %>%
      dplyr::select(-c(`Absent Count`, `Absent Proportion`)) %>%
      dplyr::rename(`observed count` = `Present Count`, 
                    `observed proportion` = `Present Proportion`) %>%
      tidyr::pivot_longer(c(`observed count`, `observed proportion`, `total count`)) %>%
      dplyr::filter(name %in% cognostics)
    
    # Expand the name depending on whether the counts are samples or biomolecules
    if (paneled_by_edata & nrow(Miss_Cog) > 0) {
      Miss_Cog <- Miss_Cog %>% 
        dplyr::mutate(
          name = lapply(name, function(x) {
            splitNames <- strsplit(x, " ") %>% unlist()
            return(paste(splitNames[1], "sample", splitNames[2]))
          }) %>% unlist()
        )
    } else if (nrow(Miss_Cog) > 0) {
      Miss_Cog <- Miss_Cog %>% 
        dplyr::mutate(
          name = lapply(name, function(x) {
            splitNames <- strsplit(x, " ") %>% unlist()
            return(paste(splitNames[1], "biomolecule", splitNames[2]))
          }) %>% unlist()
        )
    }
    
    # Add grouping data if there's more than one group 
    if (length(unique(Miss_Cog$Group)) > 1) {
      Miss_Cog <- Miss_Cog %>% dplyr::mutate(name = paste(Group, name))
    }

    # Remove group column
    Miss_Cog <- Miss_Cog %>% dplyr::select(c(name, value))
    
    # Add statistics if applicable
    if ("g-test p-value" %in% cognostics) {
      
      edata_cname <- get_edata_cname(trelliData$statRes)
      
      if (!is.na(attr(trelliData, "panel_by_stat")) && edata_cname == attr(trelliData, "panel_by_stat")) {
        
        Miss_Cog <- rbind(Miss_Cog, 
                          trelliData$trelliData.stat %>%
                            dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == GroupingVar) %>%
                            dplyr::select(Nested_DF) %>%
                            tidyr::unnest(cols = c(Nested_DF)) %>%
                            dplyr::mutate(Comparison = paste(Comparison, "g-test p-value")) %>%
                            dplyr::select(Comparison, p_value_gtest) %>%
                            dplyr::rename(name = Comparison, value = p_value_gtest) %>%
                            dplyr::mutate(value = round(value, 4))
        )
        
      }
      
    }
    
    # Return NULL if there's no cognostics 
    if (nrow(Miss_Cog) == 0) {
      return(NULL)
    }
    
    # Generate cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(Miss_Cog), function(row) {
      quick_cog(name = Miss_Cog$name[row], value = Miss_Cog$value[row])
    })) %>% dplyr::tibble()

    return(cog_to_trelli)
  }

  # Build trelliscope display---------------------------------------------------

  # If test_mode is on, then just build the required panels. If the data is statRes, we
  # will need to restructure the data a bit.
  if (!is.null(trelliData$trelliData.omics)) {
    stats_mode <- FALSE
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example, ]
    } else {
      toBuild <- trelliData$trelliData.omics
    }
  } else {
    stats_mode <- TRUE

    # Get the edata column name
    edata_cname <- get_edata_cname(trelliData$statRes)

    # Get the columns with counts
    count_cols <- colnames(trelliData$statRes)[grepl("Count", colnames(trelliData$statRes))]

    # Build toBuild dataframe
    toBuild <- trelliData$statRes %>%
      dplyr::select(c(edata_cname, count_cols)) %>%
      dplyr::group_by(!!dplyr::sym(edata_cname)) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::rename(Nested_DF = data) 
    
    # Save total counts 
    totalCounts <- attr(trelliData$statRes, "group_DF")$Group %>% 
      table(dnn = "Group") %>% 
      data.frame() %>% 
      dplyr::rename(Total = Freq)

    if (test_mode) {
      toBuild <- toBuild[test_example, ]
    }
  }

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- toBuild[test_example[1], ]
    return(missing_bar_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
  } else {
    # Pass parameters to trelli_builder function
    if (!is.null(trelliData$omicsData)) {
      
      trelli_builder(toBuild = toBuild,
                     cognostics = cognostics, 
                     plotFUN = missing_bar_plot_fun,
                     cogFUN = missing_bar_cog_fun,
                     path = path,
                     name = name,
                     remove_nestedDF = FALSE,
                     ...)
      
      
    } else {
      
      trelli_builder(toBuild = toBuild,
                     cognostics = cognostics, 
                     plotFUN = missing_bar_plot_fun,
                     cogFUN = missing_bar_cog_fun,
                     path = path,
                     name = name,
                     remove_nestedDF = TRUE,
                     ...)
      
    }

    
  }
}

determine_significance <- function(DF, p_value_thresh, is_seq) {
  
  # Subset by significance
  if (is_seq) {
    DF$Significance <- DF$p_value <= p_value_thresh
  } else {
    DF$Significance <- DF$p_value_anova <= p_value_thresh 
  }
  
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
#'   p-value will be the results from the ANOVA test. If the omics data is sedData,
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
                                  path = .getDownloadsFolder(),
                                  name = "Trelliscope",
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
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$statRes)
  if (is.na(attr(trelliData, "panel_by_stat")) || edata_cname != attr(trelliData, "panel_by_stat")) {
    stop("trelliData must be grouped by edata_cname.")
  }

  # Make foldchange bar function------------------------------------------------

  fc_bar_plot_fun <- function(DF, title) {
    
    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh, is_seq = inherits(trelliData, "trelliData.seqData"))
      if (is.null(DF)) {
        return(NULL)
      }
      
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
    bar <- bar + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
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

  # Create the cognostic function-----------------------------------------------

  fc_bar_cog_fun <- function(DF, Biomolecule) {
    
    # Extend cognostics if fold change is in it
    if ("fold change" %in% cognostics) {
      DF <- DF %>% dplyr::rename(`fold change` = fold_change)
    }
    
    # Extend cognostics if p_value is in it
    if ("p-value" %in% cognostics) {
      if (inherits(trelliData, "trelliData.seqData")) {
        DF <- DF %>% dplyr::rename(`p-value` = p_value)
      } else {
        DF <- DF %>% dplyr::rename(`p-value` = p_value_anova)
      }
    }
    
    # Prepare DF for quick_cog function
    PreCog <- DF %>%
      dplyr::select(-p_value_gtest) %>%
      dplyr::select(c(cognostics, Comparison)) %>%
      tidyr::pivot_longer(cognostics) %>%
      dplyr::mutate(Comparison = paste(Comparison, name)) %>%
      dplyr::select(-name)

    # Make quick cognostics
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(PreCog), function(row) {
      quick_cog(PreCog$Comparison[row], round(PreCog$value[row], 4))
    })) %>% dplyr::tibble()
    
    return(cog_to_trelli)
  }

  # Build trelliscope function--------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- trelliData$trelliData.stat[test_example[1], ]
    return(fc_bar_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
  }

  # Subset down to test example if applicable
  if (test_mode) {
    toBuild <- trelliData$trelliData.stat[test_example, ]
  } else {
    toBuild <- trelliData$trelliData.stat
  }

  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = fc_bar_plot_fun,
                 cogFUN = fc_bar_cog_fun,
                 path = path,
                 name = name,
                 remove_nestedDF = FALSE,
                 ...)
  
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
                                      path = .getDownloadsFolder(),
                                      name = "Trelliscope",
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
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is an emeta column
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }

  # Make sure include_points is a true or false
  if (!is.logical(include_points) & !is.na(include_points)) {
    stop("include_points must be a TRUE or FALSE.")
  }
  # Make foldchange boxplot function--------------------------------------------

  fc_box_plot_fun <- function(DF, title) {
    
    if (p_value_thresh != 0) {
    
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh, is_seq = inherits(trelliData, "trelliData.seqData"))
      if (is.null(DF)) {
        return(NULL)
      }
    
      # Make boxplot
      boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
        ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5), 
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
        ) + ggplot2::ylab("Fold Change") + ggplot2::ggtitle(title) +
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
        ) + ggplot2::ylab("Fold Change") + ggplot2::ggtitle(title) +
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

  # Make cognostic function-----------------------------------------------------

  fc_box_cog_fun <- function(DF, Group) {
    
    # Calculate biomolecule count for seqData
    if (inherits(trelliData, "trelliData.seqData")) {
      
      cog_to_trelli <- DF %>%
        dplyr::group_by(Comparison) %>%
        dplyr::summarise(
          "biomolecule count" =  sum(!is.nan(fold_change)), 
          "proportion significant" = round(sum(p_value[!is.na(p_value)] <= p_value_thresh) / `biomolecule count`, 4),
          "mean fold change" = round(mean(fold_change, na.rm = TRUE), 4),
          "sd fold change" = round(sd(fold_change, na.rm = TRUE), 4)
        ) %>%
        tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `mean fold change`, `sd fold change`)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Comparison, name)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Comparison) 
      
    } else {
      
      # Calculate stats and subset to selected choices
      cog_to_trelli <- DF %>%
        dplyr::group_by(Comparison) %>%
        dplyr::summarise(
          "biomolecule count" = sum(!is.nan(fold_change)), 
          "proportion significant" = round(sum(p_value_anova[!is.na(p_value_anova)] <= p_value_thresh) / `biomolecule count`, 4),
          "mean fold change" = round(mean(fold_change, na.rm = TRUE), 4),
          "sd fold change" = round(sd(fold_change, na.rm = TRUE), 4)
        ) %>%
        tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `mean fold change`, `sd fold change`)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Comparison, name)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Comparison) 
      
    }
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cog_to_trelli), function(row) {
      quick_cog(cog_to_trelli$name[row], cog_to_trelli$value[row])
    })) %>% dplyr::tibble()
    
    return(cog_to_trelli)
  }

  # Build the trelliscope-------------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    singleData <- trelliData$trelliData.stat[test_example[1], ]
    return(fc_box_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
  } else {
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- trelliData$trelliData.stat[test_example, ]
    } else {
      toBuild <- trelliData$trelliData.stat
    }

    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = fc_box_plot_fun,
                   cogFUN = fc_box_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...)
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
                                      path = .getDownloadsFolder(),
                                      name = "Trelliscope",
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
                  test_mode = test_mode, 
                  test_example = test_example,
                  seqDataCheck = "permissible",
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh)
  
  # Round test example to integer 
  if (test_mode) {
    test_example <- unique(abs(round(test_example)))
  }
  
  # Check that group data is an emeta column
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Make foldchange boxplot function--------------------------------------------
  
  fc_hm_plot_fun <- function(DF, title) {
    
    # Get edata cname
    edata_cname <- get_edata_cname(trelliData$statRes)
    
    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh, is_seq = inherits(trelliData, "trelliData.seqData"))
      if (is.null(DF)) {
        return(NULL)
      }
      
      # Make heatmap with significance
      hm <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = as.factor(.data[[edata_cname]]), fill = fold_change)) +
        ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ylab("Biomolecule") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggplot2::labs(fill = "Fold Change") + ggplot2::ggtitle(title) + 
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
        ggplot2::scale_fill_gradient(low = "blue", high = "red", na.value = "white") +
        ggplot2::ggtitle(title)
      
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

  # Make cognostic function-----------------------------------------------------
  
  fc_hm_cog_fun <- function(DF, Group) {
    
    # Calculate biomolecule count for seqData
    if (inherits(trelliData, "trelliData.seqData")) {
      
      cog_to_trelli <- DF %>%
        dplyr::group_by(Comparison) %>%
        dplyr::summarise(
          "biomolecule count" =  sum(!is.nan(fold_change)), 
          "proportion significant" = round(sum(p_value[!is.na(p_value)] <= p_value_thresh) / `biomolecule count`, 4),
          "mean fold change" = round(mean(fold_change, na.rm = TRUE), 4),
          "sd fold change" = round(sd(fold_change, na.rm = TRUE), 4)
        ) %>%
        tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `mean fold change`, `sd fold change`)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Comparison, name)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Comparison) 
      
    } else {
      
      # Calculate stats and subset to selected choices
      cog_to_trelli <- DF %>%
        dplyr::group_by(Comparison) %>%
        dplyr::summarise(
          "biomolecule count" = sum(!is.nan(fold_change)), 
          "proportion significant" = round(sum(p_value_anova[!is.na(p_value_anova)] <= p_value_thresh) / `biomolecule count`, 4),
          "mean fold change" = round(mean(fold_change, na.rm = TRUE), 4),
          "sd fold change" = round(sd(fold_change, na.rm = TRUE), 4)
        ) %>%
        tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `mean fold change`, `sd fold change`)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Comparison, name)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Comparison) 
      
    }
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cog_to_trelli), function(row) {
      quick_cog(cog_to_trelli$name[row], cog_to_trelli$value[row])
    })) %>% dplyr::tibble()
    
    return(cog_to_trelli)
  }

  # Build the trelliscope-------------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.stat[test_example[1], ]
    return(fc_hm_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
    
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- trelliData$trelliData.stat[test_example, ]
    } else {
      toBuild <- trelliData$trelliData.stat
    }

    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = fc_hm_plot_fun,
                   cogFUN = fc_hm_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...)
    
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
                                      path = .getDownloadsFolder(),
                                      name = "Trelliscope",
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
                                  test_mode = test_mode, 
                                  test_example = test_example,
                                  seqDataCheck = "permissible",
                                  single_plot = single_plot,
                                  p_value_thresh = p_value_thresh)
  
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
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Make foldchange volcano function--------------------------------------------
  
  fc_volcano_plot_fun <- function(DF, title) {

    if (p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh, is_seq = inherits(trelliData, "trelliData.seqData"))
      if (is.null(DF)) {
        return(NULL)
      }
      
      # Indicate which comparisons should be highlighted
      DF$Significance <- lapply(1:nrow(DF), function(row) {
        if (DF$Significance[row] == attr(DF, "LessThan")) {
          ifelse(DF$fold_change[row] > 0, 
                 paste(DF$Significance[row], "& High"), 
                 paste(DF$Significance[row],"& Low")
          )
        } else {attr(DF, "GreaterThan")}
      }) %>% unlist()
      
      # Indicate what symbols depict low, high, and no significance
      LowSig <- paste(attr(DF, "LessThan"), "& Low")
      HighSig <- paste(attr(DF, "LessThan"), "& High")
      NoSig <- attr(DF, "GreaterThan")
      
      # Make volcano plot 
      if (inherits(trelliData, "trelliData.seqData")) {
        volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value), color = Significance)) +
          ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::scale_color_manual(values = structure(c("blue", "red", "black"), 
                                                         .Names = c(LowSig, HighSig, NoSig))) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value") 
      } else {
        volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value_anova), color = Significance)) +
          ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::scale_color_manual(values = structure(c("blue", "red", "black"), 
                                                         .Names = c(LowSig, HighSig, NoSig))) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value") 
      }
      
    } else {
      
      if (inherits(trelliData, "trelliData.seqData")) {
        volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value))) +
          ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value")
      } else {
        volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value_anova))) +
          ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value")
      }
      
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

  # Make cognostic function-----------------------------------------------------
  
  fc_volcano_cog_fun <- function(DF, Group) {
    
    if (p_value_thresh != 0) {
      
      if (inherits(trelliData, "trelliData.seqData")) {
        
        # Make cognostics
        cog_to_trelli <- DF %>%
          dplyr::summarise(
            "biomolecule count" = sum(!is.nan(fold_change)), 
            "proportion significant" = round(sum(p_value[!is.na(p_value)] <= p_value_thresh) / `biomolecule count`, 4),
            "proportion significant up" = round(sum(p_value[!is.na(p_value)] <= p_value_thresh & fold_change[!is.na(fold_change)] > 0) / `biomolecule count`, 4),
            "proportion significant down" = round(sum(p_value[!is.na(p_value)] <= p_value_thresh & fold_change[!is.na(fold_change)] < 0) / `biomolecule count`, 4)
          ) %>%
          tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `proportion significant up`, `proportion significant down`)) %>%
          dplyr::filter(name %in% cognostics)
        
      } else {
        
        # Make cognostics
        cog_to_trelli <- DF %>%
          dplyr::summarise(
            "biomolecule count" = sum(!is.nan(fold_change)), 
            "proportion significant" = round(sum(p_value_anova[!is.na(p_value_anova)] <= p_value_thresh) / `biomolecule count`, 4),
            "proportion significant up" = round(sum(p_value_anova[!is.na(p_value_anova)] <= p_value_thresh & fold_change[!is.na(fold_change)] > 0) / `biomolecule count`, 4),
            "proportion significant down" = round(sum(p_value_anova[!is.na(p_value_anova)] <= p_value_thresh & fold_change[!is.na(fold_change)] < 0) / `biomolecule count`, 4)
          ) %>%
          tidyr::pivot_longer(c(`biomolecule count`, `proportion significant`, `proportion significant up`, `proportion significant down`)) %>%
          dplyr::filter(name %in% cognostics)
        
      }
   
      # Convert to trelliscope cogs
      cog_to_trelli <- do.call(cbind, lapply(1:nrow(cog_to_trelli), function(row) {
        quick_cog(cog_to_trelli$name[row], cog_to_trelli$value[row])
      })) %>% dplyr::tibble()
      
    } else {
      return(NULL)
    }
    
    return(cog_to_trelli)
  }

  # Build the trelliscope-------------------------------------------------------

  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    if (length(comparison) > 1) {stop("single_plot will only work if 1 comparison has been selected.")}
    
    singleData <- trelliData$trelliData.stat[test_example[1], ]
    return(fc_volcano_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
    
    # Set grouping variable names
    emeta_var <- attr(trelliData, "panel_by_omics")
    theComp <- "Comparison"
    
    # Create specific build variable with comparison information (different than default panel by)
    toBuild <- trelliData$trelliData.stat %>% 
      tidyr::unnest(cols = Nested_DF) %>%
      dplyr::ungroup() %>%
      dplyr::group_by_at(c(emeta_var, theComp)) %>%
      tidyr::nest() %>%
      dplyr::filter(Comparison %in% comparison) %>%
      dplyr::rename(Nested_DF = data)
    
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- toBuild[test_example, ]
    } 
    
    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = fc_volcano_plot_fun,
                   cogFUN = fc_volcano_cog_fun,
                   path = path,
                   name = name,
                   remove_nestedDF = FALSE,
                   ...)
    
  }
  
}
