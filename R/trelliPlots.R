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
#' @param p_value_skip Whether to skip specific p_value checks. Placeholder for
#'   potential future functions.
#' @param p_value_thresh The user provided threshold for plotting significant
#'   p-values.
#' @param p_value_test The user provided provided logical for whether p-value
#'   testing should occur.
trelli_precheck <- function(trelliData, 
                            trelliCheck,
                            cognostics, 
                            acceptable_cognostics,
                            ggplot_params,
                            interactive, 
                            test_mode,
                            test_example,
                            single_plot,
                            p_value_skip = FALSE,
                            p_value_thresh = NULL,
                            p_value_test = NULL) {
  
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
      stop(paste("Unacceptable cognostic option included. Acceptable options are: ", 
                 paste(acceptable_cognostics, collapse = ", ")))
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
    
    # Check p-value test method 
    if (length(p_value_test) > 1 | !is.logical(p_value_test)) {
      stop("p_value_test must be a single TRUE or FALSE.")
    }
    if (is.na(p_value_test)) {p_value_test <- FALSE}
    
    # Check that the data has the requested p_value 
    if (p_value_test & any(grepl("P_value_A", colnames(trelliData$statRes))) == FALSE) {
      message(paste(
        "No imd-anova stats were detected in the statRes object which is",
        "required for p_value_test to be TRUE. Seting P_value_test to FALSE."
      ))
      p_value_test <- FALSE
    } 
    
    # Check p_value threshold
    if (!is.numeric(p_value_thresh)) {
      stop("p_value_thresh must be a numeric.")
    }
    if (p_value_thresh < 0 | p_value_thresh > 1) {
      stop("p_value_thresh must be between 0 and 1.")
    }
    
  }
  
  return(p_value_test)
  
}

# A quick cognostic function 
quick_cog <- function(name, value) {
  dplyr::tibble(!!rlang::sym(name) := trelliscopejs::cog(value, desc = name))
}

# Create a list to convert from short name to long
name_converter_abundance <- list("n" = "Count", "mean" = "Mean Abundance", 
                       "median" = "Median Abundance", "sd" = "Standard Deviation Abundance", 
                       "skew" = "Skew Abundance", "p_value_anova" = "Anova P Value",
                       "p_value_gtest" = "G-Test P Value", "fold_change" = "Fold Change")
name_converter_foldchange <- list("n" = "Count", "mean" = "Mean Fold Change", 
                                  "median" = "Median Fold Change", "sd" = "Standard Deviation Fold Change")

# This function builds all trelliscopes.
trelli_builder <- function(toBuild, cognostics, plotFUN, cogFUN, path, name, ...) {
  
  # Remove any blank names 
  if ("" %in% unlist(toBuild[,1])) {
    message(paste("Removing", length(sum("" %in% toBuild[,1])), "blank biomolecule names."))
    toBuild <- toBuild[toBuild[,1] != "",]
  }
  
  if (nrow(toBuild) == 0) {
    stop("No data to build trelliscope with.")
  }
  
  # Build trelliscope without cognostics if none are provided. Otherwise, build with cognostics.
  if (is.null(cognostics)) {
    
    toBuild %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), plotFUN)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = TRUE, ...) 
    
  } else {
    
    toBuild <- toBuild %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), plotFUN),
        cog = trelliscopejs::map2_cog(Nested_DF, as.character(unlist(toBuild[,1])), cogFUN)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = TRUE, ...) 
    
  }
}

# Get downloads folder
.getDownloadsFolder <- function(.test_mode = FALSE) {
  if (Sys.info()['sysname'] == "Windows" | .test_mode) {
    folder <- dirname("~")
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
#'   the first main effect in group_designation.
#'
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by trelli_panel_by. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are n, mean, median, sd, and skew for abundance. If statRes data is
#'   included, p_value and fold_change cognostics can be added. If no cognostics
#'   are desired, set to NULL.
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
#' 
#' @examples
#' \dontrun{
#' 
#' ## Build the abundance boxplot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData, panel = "Lipid") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
#' trelli_panel_by(trelliData = trelliData, panel = "Sample") %>% trelli_abundance_boxplot()
#' 
#' ## Build the abundance boxplot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Lipid") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
#' trelli_panel_by(trelliData = trelliData2, panel = "LipidFamily") %>% trelli_abundance_boxplot()
#'     
#' ## Build the abundance boxplot with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData4, panel = "Lipid") %>%
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% trelli_abundance_boxplot()
#'    
#' ## Other options include modifying the ggplot  
#' trelli_panel_by(trelliData = trelliData, panel = "Lipid") %>% 
#'    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "ylim(c(2,20))"))
#' 
#' ## Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% trelli_abundance_boxplot(interactive = TRUE)
#' 
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_abundance_boxplot <- function(trelliData,
                                     cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
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
                  acceptable_cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = NULL,
                  p_value_test = NULL)
  
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("p_value", "fold_change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[-match(c("p_value", "fold_change"), cognostics, nomatch = 0)]
    }    
  }
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Make sure include_points is a true or false
  if (!is.logical(include_points) & !is.na(include_points)) {
    stop("include_points must be a TRUE or FALSE.")
  }

  # Make boxplot function-------------------------------------------------------
  
  # First, generate the boxplot function
  box_plot_fun <- function(DF, title) {
    
    # Add a blank group if no group designation was given
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {DF$Group <- "x"} 
    
    # Build plot 
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Group, fill = Group, y = Abundance)) + 
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::theme_bw() + 
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) + 
      ggplot2::ylab(paste(attr(trelliData$omicsData, "data_info")$data_scale, "Abundance")) 
    
    # Add include_points 
    if (include_points) {
      boxplot <- boxplot + ggplot2::geom_jitter(height = 0, width = 0.25)
    } 
    
    # Remove x axis if no groups
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      boxplot <- boxplot + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                          axis.ticks.x = ggplot2::element_blank(), 
                                          axis.text.x = ggplot2::element_blank())
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
      "n" = dplyr::tibble(`Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Biomolecule Count")),
      "mean" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = TRUE), 4), desc = "Mean Abundance")), 
      "median" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = TRUE), 4), desc = "Median Abundance")), 
      "sd" = dplyr::tibble(`Standard Deviation Abundance` = trelliscopejs::cog(round(sd(DF$Abundance, na.rm = TRUE), 4), desc = "Abundance Standard Deviation")), 
      "skew" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = TRUE), 4), desc= "Abundance Skewness"))
    )
    
    # If cognostics are any of the cog, then add them 
    if (any(cognostics %in% c("n", "mean", "median", "sd", "skew"))) {
      cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% tibble::tibble()
    } else {cog_to_trelli <- NULL}
    
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
          "n" = sum(!is.na(Abundance)), 
          "mean" = round(mean(Abundance, na.rm = TRUE), 4),
          "median" = round(median(Abundance, na.rm = TRUE), 4),
          "sd" = round(sd(Abundance, na.rm = TRUE), 4),
          "skew" = round(e1071::skewness(Abundance, na.rm = TRUE), 4)
        ) %>%
        tidyr::pivot_longer(c(n, mean, median, sd, skew)) %>%
        dplyr::filter(name %in% cognostics) %>%
        dplyr::mutate(
          name = paste(Group, lapply(name, function(x) {name_converter_abundance[[x]]}) %>% unlist())
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Group) 
      
      # Add new cognostics 
      cog_to_trelli <- cbind(cog_to_trelli, do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
        quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
      })) %>% tibble::tibble()) %>% tibble::tibble()
      
    }
    
    # Add cognostics that only apply when stats data is grouped by edata_cname
    edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
    
    if (!is.null(trelliData$trelliData.stat) && !is.na(attr(trelliData, "panel_by_stat")) && edata_cname == attr(trelliData, "panel_by_stat")) {
      
      # Downselect to only stats 
      stat_cogs <- cognostics[cognostics %in% c("fold_change", "p_value")]
      
      if (length(stat_cogs) != 0) {
        
        # Subset down the dataframe down to group, unnest the dataframe, 
        # pivot_longer to comparison, subset columns to requested statistics, 
        # switch name to a more specific name
        
        # Update stat cogs to accept the new p-value groups
        if ("p_value" %in% stat_cogs) {
          theNames <- trelliData$trelliData.stat$Nested_DF[[1]] %>% colnames()
          p_value_cols <- theNames[grepl("p_value", theNames)]
          stat_cogs <- stat_cogs[stat_cogs != "p_value"]
          stat_cogs <- c(stat_cogs, p_value_cols)
        }
        
        cogs_to_add <- trelliData$trelliData.stat %>%
          dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == biomolecule) %>%
          dplyr::select(Nested_DF) %>%
          tidyr::unnest(cols = c(Nested_DF)) %>%
          dplyr::select(c(Comparison, stat_cogs)) %>%
          tidyr::pivot_longer(stat_cogs) %>%
          dplyr::mutate(
            name = paste(Comparison, lapply(name, function(x) {name_converter_abundance[[x]]}) %>% unlist()),
            value = round(value, 4)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-Comparison)
        
        # Generate new cognostics 
        new_cogs <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
          quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
        })) %>% tibble::tibble()

        
        # Add new cognostics, removing when it is NULL 
        cog_to_trelli <- cbind(cog_to_trelli, new_cogs) %>% tibble::tibble()
        
      }
      
    }
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.omics[test_example[1],]
    return(box_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
  
    # If test_mode is on, then just build the required panels
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example,]
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
                   ...)
    
  }
  
}

#' @name trelli_abundance_histogram
#' 
#' @title Histogram trelliscope building function for abundance data
#'
#' @description Specify a plot design and cognostics for the abundance histogram
#'   trelliscope. Main_effects grouping are ignored. Data must be grouped by
#'   edata_cname.
#'
#' @param trelliData A trelliscope data object made by as.trelliData or
#'   as.trelliData.edata, and grouped by edata_cname in trelli_panel_by.
#'   Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are n, mean, median, sd, and skew. p_value and fold_change can be added if
#'   statRes is included.
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
#' 
#' @examples
#' \dontrun{
#' 
#' ## Build the abundance histogram with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData, panel = "Lipid") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
#' 
#' ## Build the abundance histogram with an omicsData object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Lipid") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
#'     
#' ## Build the abundance histogram with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData4, panel = "Lipid") %>%
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
#'    
#' ## Users can modify the plotting function with ggplot parameters and interactivity, 
#' ## and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData, panel = "Lipid") %>% 
#'    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "xlab('Abundance')"), interactive = TRUE,
#'      cognostics = c("mean", "median"))  
#'    
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_abundance_histogram <- function(trelliData,
                                       cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
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
                  acceptable_cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = NULL,
                  p_value_test = NULL)

  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("p_value", "fold_change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[-match(c("p_value", "fold_change"), cognostics, nomatch = 0)]
    }    
  }
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
  if (edata_cname != attr(trelliData, "panel_by_omics")) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Make histogram function-----------------------------------------------------
  
  # First, generate the histogram function
  hist_plot_fun <- function(DF, title) {
    
    # Remove NAs
    DF <- DF[!is.na(DF$Abundance),]
    
    # Build plot 
    histogram <- ggplot2::ggplot(DF, ggplot2::aes(x = Abundance)) + 
      ggplot2::geom_histogram(bins = 10, fill = "steelblue", color = "black") + ggplot2::ggtitle(title) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
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
      "n" = dplyr::tibble(`Count` = trelliscopejs::cog(sum(!is.na(DF$Abundance)), desc = "Biomolecule Count")),
      "mean" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = TRUE), 4), desc = "Mean Abundance")), 
      "median" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = TRUE), 4), desc = "Median Abundance")), 
      "sd" = dplyr::tibble(`Standard Deviation Abundance` = trelliscopejs::cog(round(sd(DF$Abundance, na.rm = TRUE), 4), desc = "Abundance Standard Deviation")), 
      "skew" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = TRUE), 4), desc= "Abundance Skewness"))
    )
    
    # If cognostics are any of the cog, then add them 
    cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% tibble::tibble()
    
    # Add statistics if applicable 
    if (!is.null(trelliData$trelliData.stat)) {
      
      # Downselect to only stats 
      stat_cogs <- cognostics[cognostics %in% c("fold_change", "p_value")]
      
      if (length(stat_cogs) != 0) {
        
        # Update stat cogs to accept the new p-value groups
        if ("p_value" %in% stat_cogs) {
          theNames <- trelliData$trelliData.stat$Nested_DF[[1]] %>% colnames()
          p_value_cols <- theNames[grepl("p_value", theNames)]
          stat_cogs <- stat_cogs[stat_cogs != "p_value"]
          stat_cogs <- c(stat_cogs, p_value_cols)
        }
        
        # Subset down the dataframe down to group, unnest the dataframe, 
        # pivot_longer to comparison, subset columns to requested statistics, 
        # switch name to a more specific name
        cogs_to_add <- trelliData$trelliData.stat %>%
          dplyr::filter(trelliData$trelliData.stat[[edata_cname]] == biomolecule) %>%
          dplyr::select(Nested_DF) %>%
          tidyr::unnest(cols = c(Nested_DF)) %>%
          dplyr::select(c(Comparison, stat_cogs)) %>%
          tidyr::pivot_longer(stat_cogs) %>%
          dplyr::mutate(
            name = paste(Comparison, lapply(name, function(x) {name_converter_abundance[[x]]}) %>% unlist()),
            value = round(value, 4)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-Comparison)
        
        # Generate new cognostics 
        new_cogs <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
          quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
        })) %>% tibble::tibble()
        
        # Add new cognostics 
        cog_to_trelli <- cbind(cog_to_trelli, new_cogs) %>% tibble::tibble()
      
        
      }
      
    }
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.omics[test_example[1],]
    return(hist_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
  
    # If test_mode is on, then just build the required panels
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example,]
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
                   ...)
    
  }
  
}

#' @name trelli_abundance_heatmap
#' 
#' @title Heatmap trelliscope building function for abundance data
#'
#' @description Specify a plot design and cognostics for the abundance heatmap
#'   trelliscope. Data must be grouped by an emeta column. Main_effects order
#'   the y-variables. All statRes data is ignored.
#'
#' @param trelliData A trelliscope data object made by as.trelliData, and
#'   grouped by an emeta variable. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are n, mean, median, sd, and skew per "main_effects" group designation.
#'   Otherwise, no cognostics are returned.
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
#' 
#' @examples
#' \dontrun{
#' 
#' ## Build the abundance heatmap with an omicsData object with emeta variables. Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData2, panel = "LipidFamily") %>% 
#'    trelli_abundance_heatmap(test_mode = TRUE, test_example = 1:3)
#'    
#' ## Users can modify the plotting function with ggplot parameters and interactivity, 
#' ## and can also select certain cognostics.     
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
#'    trelli_abundance_heatmap(test_mode = TRUE, test_example = 1:5, 
#'      ggplot_params = c("ylab('')", "xlab('')"), interactive = TRUE, cognostics = c("mean", "median"))  
#'    
#' }
#' 
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_abundance_heatmap <- function(trelliData,
                                     cognostics = c("n", "mean", "median", "sd", "skew"),
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
                  acceptable_cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = NULL,
                  p_value_test = NULL)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is grouped by an e_meta variable
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # If no group designation, set cognostics to NULL. 
  if (is.null(attributes(trelliData$omicsData)$group_DF)) {
    congnostics <- NULL
  }
  
  # Make heatmap function-------------------------------------------------------
  
  # First, generate the heatmap function
  hm_plot_fun <- function(DF, title) {
    
    # Get fdata_cname
    fdata_cname <- get_fdata_cname(trelliData$omicsData)
    
    # If group designation was set, then convert Group to a factor variable 
    if (!is.null(attributes(trelliData$omicsData)$group_DF)) {
      theLevels <- attr(trelliData$omicsData, "group_DF") %>% dplyr::arrange(Group) %>% dplyr::select(dplyr::all_of(fdata_cname)) %>% unlist()
      DF[,fdata_cname] <- factor(unlist(DF[,fdata_cname]), levels = theLevels)
    } else {
      DF[,fdata_cname] <- as.factor(unlist(DF[,fdata_cname]))
    }
    
    # Get edata and fdata cname
    edata_cname <- get_edata_cname(trelliData$omicsData)
    
    # Build plot: this should be edata_cname
    hm <- ggplot2::ggplot(DF, ggplot2::aes(x = as.factor(.data[[edata_cname]]), y = .data[[fdata_cname]], fill = Abundance)) +
      ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ylab("Sample") + ggplot2::xlab("Biomolecule") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + 
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
        "n" = sum(!is.na(Abundance)), 
        "mean" = round(mean(Abundance, na.rm = TRUE), 4),
        "median" = round(median(Abundance, na.rm = TRUE), 4),
        "sd" = round(sd(Abundance, na.rm = TRUE), 4),
        "skew" = round(e1071::skewness(Abundance, na.rm = TRUE), 4)
      ) %>%
      tidyr::pivot_longer(c(n, mean, median, sd, skew)) %>%
      dplyr::filter(name %in% cognostics) %>%
      dplyr::mutate(
        name = paste(Group, lapply(name, function(x) {name_converter_abundance[[x]]}) %>% unlist())
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Group) 
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
      quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
    })) %>% tibble::tibble() 
    
    return(cog_to_trelli)
    
  } 
  
  # Build trelliscope display---------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.omics[test_example[1],]
    return(hm_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
    
    # If test_mode is on, then just build the required panels
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example,]
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
                   ...) 
    
  }
  
}

#' @name trelli_missingness_bar
#' 
#' @title Bar chart trelliscope building function for missing data   
#' 
#' @description Specify a plot design and cognostics for the missing barchart
#'   trelliscope. Missingness is displayed per panel_by variable. Main_effects
#'   data is used to split samples when applicable.
#'
#' @param trelliData A trelliscope data object made by as.trelliData.edata or
#'   as.trelliData. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are is n and proportion.
#' @param proportion A logical to determine whether plots should display counts
#'   or proportions. Default is TRUE.
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
#'   
#' @examples
#' \dontrun{
#' 
#' ## Build the missingness bar plot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_panel_by(trelliData = trelliData, panel = "Lipid") %>% 
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
#' trelli_panel_by(trelliData = trelliData, panel = "Sample") %>% trelli_missingness_bar()
#' 
#' ## Build the missingness bar plot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData2, panel = "Lipid") %>% 
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
#' 
#' ## Build the missingness bar plot with a statRes object. Generate trelliData in as.trelliData
#' trelli_panel_by(trelliData = trelliData3, panel = "Lipid") %>%
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
#' 
#' ## Build the missingness bar plot with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_panel_by(trelliData = trelliData4, panel = "Lipid") %>%
#'   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10) 
#' 
#' ## Or making the plot interactive 
#' trelli_panel_by(trelliData = trelliData2, panel = "Lipid") %>% 
#'    trelli_missingness_bar(test_mode = TRUE, test_example = 1:5, interactive = TRUE)
#'    
#' ## Or visualize only count data 
#' trelli_panel_by(trelliData = trelliData2, panel = "Lipid") %>% 
#'    trelli_missingness_bar(test_mode = TRUE, test_example = 1:5, cognostics = "n", proportion = FALSE)
#'    
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_missingness_bar <- function(trelliData,
                                   cognostics = c("n", "proportion"),
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
                  acceptable_cognostics = c("n", "proportion"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = NULL,
                  p_value_test = NULL)
  
  # Check that proportion is a non NA logical
  if (!is.logical(proportion) | is.na(proportion)) {
    stop("proportion must be a TRUE or FALSE.")
  }
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
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
          `Absent Proportion` = `Absent Count` / Total,
          `Present Proportion` = `Present Count` / Total
        ) %>%
        dplyr::select(-Total) %>%
        dplyr::mutate(Group, `Absent Count`, `Present Count`, `Absent Proportion`, `Present Proportion`)
      
    } else {
      
      # Add a blank group if no group designation was given
      if (is.null(attributes(trelliData$omicsData)$group_DF)) {DF$Group <- "x"} 
      
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
      ggplot2::geom_bar(stat = "identity", position = "stack", color = "black") + ggplot2::theme_bw() + 
      ggplot2::ggtitle(title) + ggplot2::ylab(ylab) +
      ggplot2::scale_fill_manual(values = c("black", "steelblue")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.title = ggplot2::element_blank()) 
    
    # Remove x axis if no groups
    if (is.null(attributes(trelliData$omicsData)$group_DF) & stats_mode == FALSE) {
      missing_bar <- missing_bar + ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                                                  axis.ticks.x = ggplot2::element_blank(), 
                                                  axis.text.x = ggplot2::element_blank())
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
      tidyr::pivot_longer(c(`Absent Count`, `Present Count`, `Absent Proportion`, `Present Proportion`))
    
    # Add grouping data if there's more than one group 
    if (length(unique(Miss_Cog$Group)) > 1) {
      Miss_Cog <- Miss_Cog %>% dplyr::mutate(name = paste(Group, name))
    }
      
    # Remove group column
    Miss_Cog <- Miss_Cog %>% dplyr::select(c(name, value))
    
    # Subset down to selected cognostics
    if (length(cognostics) == 1) {
      if (cognostics == "n") {
        Miss_Cog <- Miss_Cog %>% dplyr::filter(grepl("Count", name))
      } else if (cognostics == "proportion") {
        Miss_Cog <- Miss_Cog %>% dplyr::filter(grepl("Proportion", name))
      }
    }
    
    # Generate cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(Miss_Cog), function(row) {
      quick_cog(name = Miss_Cog$name[row], value = Miss_Cog$value[row])
    })) %>% tibble::tibble()
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # If test_mode is on, then just build the required panels. If the data is statRes, we 
  # will need to restructure the data a bit. 
  if (!is.null(trelliData$trelliData.omics)) {
    stats_mode <- FALSE
    if (test_mode) {
      toBuild <- trelliData$trelliData.omics[test_example,]
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
      dplyr::group_by(!!rlang::sym(edata_cname)) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::rename(Nested_DF = data)
    
    # Save total counts 
    totalCounts <- attr(trelliData$statRes, "group_DF")$Group %>% 
      table(dnn = "Group") %>% 
      data.frame() %>% 
      dplyr::rename(Total = Freq)
    
    if (test_mode) {toBuild <- toBuild[test_example,]} 
  }
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- toBuild[test_example[1],]
    return(missing_bar_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
    
    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = missing_bar_plot_fun,
                   cogFUN = missing_bar_cog_fun,
                   path = path,
                   name = name,
                   ...)
    
  }
  
}

determine_significance <- function(DF, p_value_thresh) {
  
  # Indicate significant values 
  DF$Significance <- DF$p_value_anova <= p_value_thresh
  
  # Filter out NA values 
  DF <- DF[!(is.nan(DF$Significance) | is.nan(DF$fold_change)),]
  
  # Return NULL if no rows left
  if (nrow(DF) == 0) {return(NULL)}
  
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
#'   are are fold_change and p_value.
#' @param p_value_test A logical to indicate whether specific significant
#'   biomolecules are to be incidated with a black outline, if an imd-anova was
#'   run. Default is FALSE.
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
#'   
#' @examples
#' \dontrun{
#' 
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData3, panel = "Lipid") %>% 
#'   trelli_foldchange_bar(test_mode = TRUE, test_example = 1:10, p_value_test = TRUE)
#'   
#' ## Or make the plot interactive  
#' trelli_panel_by(trelliData = trelliData4, panel = "Lipid") %>% 
#'   trelli_foldchange_bar(test_mode = TRUE, test_example = 1:10, p_value_test = TRUE, interactive = TRUE) 
#'    
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_bar <- function(trelliData,
                                  cognostics = c("fold_change", "p_value"),
                                  p_value_test = FALSE,
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
                  trelliCheck = "stat",
                  cognostics = cognostics,
                  acceptable_cognostics = c("p_value", "fold_change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh,
                  p_value_test = p_value_test)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$statRes)
  if (is.na(attr(trelliData, "panel_by_stat")) || edata_cname != attr(trelliData, "panel_by_stat")) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Make foldchange bar function------------------------------------------------
  
  fc_bar_plot_fun <- function(DF, title) {
    
    if (p_value_test && p_value_thresh != 0) {
      
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
    bar <- bar + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5), 
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
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
    
    # Extend cognostics if p_value is in it
    if ("p_value" %in% cognostics) {
      theNames <- trelliData$trelliData.stat$Nested_DF[[1]] %>% colnames()
      p_value_cols <- theNames[grepl("p_value", theNames)]
      cognostics <- cognostics[cognostics != "p_value"]
      cognostics <- c(cognostics, p_value_cols)
    }
    
    # Prepare DF for quick_cog function
    PreCog <- DF %>%
      tidyr::pivot_longer(cognostics) %>%
      dplyr::mutate(Comparison = paste(Comparison, name)) %>%
      dplyr::select(-name)
    
    # Make quick cognostics
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(PreCog), function(row) {
      quick_cog(gsub("_", " ", PreCog$Comparison[row]), round(PreCog$value[row], 4))
    })) %>% tibble::tibble()
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope function--------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.stat[test_example[1],]
    return(fc_bar_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  }
  
  # Subset down to test example if applicable
  if (test_mode) {
    toBuild <- trelliData$trelliData.stat[test_example,]
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
                 ...)
  
}

#' @name trelli_foldchange_boxplot
#' 
#' @title Boxplot trelliscope building function for fold_change   
#' 
#' @description Specify a plot design and cognostics for the fold_change boxplot
#'   trelliscope. Fold change must be grouped by an emeta column, which means
#'   both an omicsData object and statRes are required to make this plot.
#'
#' @param trelliData A trelliscope data object with omicsData and statRes
#'   results. Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are are n, mean, median, and sd.
#' @param p_value_test A logical to indicate whether specific significant
#'   biomolecules are to be indicated with a changed symbol if an imd-anova was
#'   run. Default is FALSE.
#' @param p_value_thresh A value between 0 and 1 to indicate significant
#'   biomolecules for p_value_test. Default is 0.05.
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
#'   
#' @examples
#' \dontrun{ 
#' 
#' ## Build fold_change box plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
#'   trelli_foldchange_boxplot(p_value_test = TRUE)
#'
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_boxplot <- function(trelliData,
                                      cognostics = c("n", "median", "mean", "sd"),
                                      p_value_test = FALSE,
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
  p_value_test <- trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("omics", "stat"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("n", "median", "mean", "sd"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh,
                  p_value_test = p_value_test)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
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
    
    # Get significant values
    DF <- determine_significance(DF, p_value_thresh)
    if (is.null(DF)) {return(NULL)}
  
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
      if (p_value_test && p_value_thresh != 0) {
        boxplot <- boxplot + ggplot2::geom_jitter(ggplot2::aes(shape = Significance), height = 0, width = 0.25) +
          ggplot2::scale_shape_manual(values = structure(c(17, 16), .Names = c(attr(DF, "LessThan"), attr(DF, "GreaterThan"))))
      } else {
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
    
    # Calculate stats and subset to selected choices 
    cog_to_trelli <- DF %>%
      dplyr::group_by(Comparison) %>%
      dplyr::summarise(
        "n" = sum(!is.nan(fold_change)), 
        "mean" = round(mean(fold_change, na.rm = TRUE), 4),
        "median" = round(median(fold_change, na.rm = TRUE), 4),
        "sd" = round(sd(fold_change, na.rm = TRUE), 4)
      ) %>%
     tidyr::pivot_longer(c(n, mean, median, sd)) %>%
     dplyr::filter(name %in% cognostics) %>%
     dplyr::mutate(
      name = paste(Comparison, lapply(name, function(x) {name_converter_foldchange[[x]]}) %>% unlist())
     ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Comparison) 
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cog_to_trelli), function(row) {
      quick_cog(gsub("_", " ", cog_to_trelli$name[row]), cog_to_trelli$value[row])
    })) %>% tibble::tibble()
    
    return(cog_to_trelli)
    
  }
  
  # Build the trelliscope-------------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.stat[test_example[1],]
    return(fc_box_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
  
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- trelliData$trelliData.stat[test_example,]
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
#'   Required.
#' @param cognostics A vector of cognostic options for each plot. Valid entry is
#'   n.
#' @param p_value_test A logical to indicate whether specific significant
#'   biomolecules are to be indicated with a changed color if an imd-anova was
#'   run. Default is FALSE.
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
#'   
#' @examples
#' \dontrun{ 
#' 
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
#'   trelli_foldchange_volcano(comparison = "Mock_vs_Infection_A")
#'
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_volcano <- function(trelliData,
                                      comparison, 
                                      cognostics = c("n"),
                                      p_value_test = TRUE,
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
                  acceptable_cognostics = c("n"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_thresh = p_value_thresh,
                  p_value_test = p_value_test)
  
  # Ensure that comparison is an acceptable input
  if (!is.character(comparison) | length(comparison) > 1) {
    stop("comparison must be a string of length 1.")
  }
  
  # Ensure that comparison is from the comparisons lists
  Comparisons <- attr(trelliData$statRes, "comparisons")
  if (comparison %in% Comparisons == FALSE) {
    stop(paste0(comparison, "is not an acceptable comparison"))
  }
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is an emeta column
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Make foldchange volcano function--------------------------------------------
  
  fc_volcano_plot_fun <- function(DF, title) {
    
    # Subset comparison
    DF <- DF[DF$Comparison == comparison,]
    
    if (p_value_test && p_value_thresh != 0) {

      # Get significant values
      DF <- determine_significance(DF, p_value_thresh)
      if (is.null(DF)) {return(NULL)}
      
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
      
      # Make volcano plot Stopped Here
      volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value_anova), color = Significance)) +
        ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
        ggplot2::scale_color_manual(values = structure(c("blue", "red", "black"), 
          .Names = c(LowSig, HighSig, NoSig))) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value")
      
    } else {
      
      volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value_anova))) +
        ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
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
  
  # Make cognostic function-----------------------------------------------------
  
  fc_volcano_cog_fun <- function(DF, Group) {
    
    # Subset comparison
    DF <- DF[DF$Comparison == comparison,]
    
    # Get significant values
    DF <- determine_significance(DF, p_value_thresh)
    if (is.null(DF)) {return(NULL)}
    
    # Indicate which comparisons should be highlighted
    DF$SigCounts <- lapply(1:nrow(DF), function(row) {
      if (DF$Significance[row] == attr(DF, "LessThan")) {
        ifelse(DF$fold_change[row] > 0, "High", "Low")
      } else {return("Not Significant")}
    }) %>% unlist()
    
    # Get count 
    counts <- DF$SigCounts%>%
      table(dnn = "Cog") %>%
      data.table::data.table()
    
    # Add 0's if necessary
    if (nrow(counts) != 3) {
      all_options <- c("High", "Low", "Not Significant")
      missing <- all_options[all_options %in% counts$Cog == FALSE]
      counts <- rbind(counts, data.frame(Cog = missing, N = 0))
    }
    
    # Set order
    counts <- counts[order(counts$Cog),]
    counts$N <- as.numeric(counts$N)
    
    # Convert to trelliscope cogs
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(counts), function(row) {
      quick_cog(paste(counts$Cog[row], "Count"), counts$N[row])
    })) %>% tibble::tibble()
    
    return(cog_to_trelli)
    
  }
  
  # Build the trelliscope-------------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.stat[test_example[1],]
    return(fc_volcano_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
  
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- trelliData$trelliData.stat[test_example,]
    } else {
      toBuild <- trelliData$trelliData.stat
    }
    
    # Pass parameters to trelli_builder function
    trelli_builder(toBuild = toBuild,
                   cognostics = cognostics, 
                   plotFUN = fc_volcano_plot_fun,
                   cogFUN = fc_volcano_cog_fun,
                   path = path,
                   name = name,
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
#' @param cognostics A vector of cognostic options for each plot. Valid entries
#'   are are n, mean, median, and sd.
#' @param p_value_test A logical to indicate whether specific significant
#'   biomolecules are to be indicated with a dot if an imd-anova was run.
#'   Default is FALSE.
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
#' 
#' @examples
#' \dontrun{ 
#' 
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
#'   trelli_foldchange_heatmap()
#'
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_heatmap <- function(trelliData,
                                      cognostics = c("n", "median", "mean", "sd"),
                                      p_value_test = FALSE,
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
                  acceptable_cognostics = c("n", "median", "mean", "sd"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example,
                  single_plot = single_plot,
                  p_value_test = p_value_test,
                  p_value_thresh = p_value_thresh)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is an emeta column
  if (attr(trelliData, "panel_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be paneled_by an e_meta column.")
  }
  
  # Make foldchange boxplot function--------------------------------------------
  
  fc_hm_plot_fun <- function(DF, title) {
    
    # Get edata cname
    edata_cname <- get_edata_cname(trelliData$statRes)
    
    if (p_value_test && p_value_thresh != 0) {
      
      # Get significant values
      DF <- determine_significance(DF, p_value_thresh)
      if (is.null(DF)) {return(NULL)}
      
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
    
    # Calculate stats and subset to selected choices 
    cog_to_trelli <- DF %>%
      dplyr::group_by(Comparison) %>%
      dplyr::summarise(
        "n" = sum(!is.nan(fold_change)), 
        "mean" = round(mean(fold_change, na.rm = TRUE), 4),
        "median" = round(median(fold_change, na.rm = TRUE), 4),
        "sd" = round(sd(fold_change, na.rm = TRUE), 4)
      ) %>%
      tidyr::pivot_longer(c(n, mean, median, sd)) %>%
      dplyr::filter(name %in% cognostics) %>%
      dplyr::mutate(
        name = paste(Comparison, lapply(name, function(x) {name_converter_foldchange[[x]]}) %>% unlist())
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Comparison) 
    
    # Add new cognostics 
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cog_to_trelli), function(row) {
      quick_cog(gsub("_", " ", cog_to_trelli$name[row]), cog_to_trelli$value[row])
    })) %>% tibble::tibble()
    
    return(cog_to_trelli)
    
  }
  
  # Build the trelliscope-------------------------------------------------------
  
  # Return a single plot if single_plot is TRUE
  if (single_plot) {
    
    singleData <- trelliData$trelliData.stat[test_example[1],]
    return(fc_hm_plot_fun(singleData$Nested_DF[[1]], unlist(singleData[1, 1])))
    
  } else {
    
    # Subset down to test example if applicable
    if (test_mode) {
      toBuild <- trelliData$trelliData.stat[test_example,]
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
                   ...)
    
  }

}
