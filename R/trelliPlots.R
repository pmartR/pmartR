# This function runs necessary checks for pmartR trelliscope plotting functions. 
# It cleans any parameters (rounding numerics to integers, etc.), and returns them.
trelli_precheck <- function(trelliData, trelliCheck,
                            cognostics, acceptable_cognostics,
                            ggplot_params,
                            interactive, 
                            test_mode,
                            test_example) {
  
  #######################
  ## trelliData checks ##
  #######################
  
  # trelliData object must be of the trelliData class
  if (any(class(trelliData) %in% c("trelliData")) == FALSE) {
    stop("trelliData must be of the class trelliData.")
  }
  
  # Check that trelliData has been passed to the "trelli_group_by" function.
  if (!attr(trelliData, "group_by")) {
    stop("trelliData must be grouped with trelli_group_by.")
  }
  
  # Check that omics data exists
  if ("omics" %in% trelliCheck) {
    
    # Assert that trelliData has omicsData
    if (is.null(trelliData$trelliData.omics)) {
      stop("trelliData must have omicsData for this plotting function.")
    }
    
  }
  
  # Check that statRes data exists 
  if ("stat" %in% trelliCheck) {
    
    # Assert that trelliData has statRes
    if (is.null(trelliData$trelliData.stat)) {
      stop("trelliData must have statRes for this plotting function.")
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
    stop("interactive must be a true or false.")
  }
  
  # test_mode must be a TRUE/FALSE
  if (!is.logical(test_mode) & !is.na(test_mode)) {
    stop("test_mode must be a true or false")
  }

  #########################
  ## TEST EXAMPLE CHECKS ##
  #########################
  
  # Ensure that test_example is an integer
  if (!is.numeric(test_example) | 0 %in% test_example) {
    "test_example should be a non-zero integer."
  }
  
  # Ensure that test_example is in the range of possibilities 
  if (max(test_example) > nrow(trelliData$trelliData.omics)) {
    stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.omics)))
  }
  
}

# A quick cognostic function 
quick_cog <- function(name, value) {
  dplyr::tibble(!!rlang::sym(name) := trelliscopejs::cog(value, desc = name))
}

# Create a list to convert from short name to long
name_converter <- list("n" = "Count", "mean" = "Mean Abundance", 
                       "median" = "Median Abundance", "sd" = "Standard Deviation Abundance", 
                       "skew" = "Skew Abundance", "p_value" = "P Value", "fold_change" = "Fold Change")

# This function builds all trelliscopes.
trelli_builder <- function(toBuild, cognostics, plotFUN, cogFUN, path, name, ...) {
  
  # Build trelliscope without cognostics if none are provided. Otherwise, build with cognostics.
  if (is.null(cognostics)) {
    
    toBuild %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), plotFUN)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = T, ...)
    
  } else {
    
    toBuild %>%
      dplyr::mutate(
        panel = trelliscopejs::map2_plot(Nested_DF, as.character(unlist(toBuild[,1])), plotFUN),
        cog = trelliscopejs::map2_cog(Nested_DF, as.character(unlist(toBuild[,1])), cogFUN)
      ) %>%
      trelliscopejs::trelliscope(path = path, name = name, nrow = 1, ncol = 1, thumb = T, ...)
    
  }
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
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "ylim(c(2,20))"). Default is NULL. 
#' @param interactive A logical argument indicating whether the plots should be interactive
#'    or not. Interactive plots are ggplots piped to ggplotly (for now). Default is FALSE.  
#' @param path The base directory of the trelliscope application. Default is Downloads. 
#' @param name The name of the display. Default is Trelliscope.
#' @param test_mode A logical to return a smaller trelliscope to confirm plot and design.
#'    Default is FALSE.
#' @param test_example The index number of the plot to return for test_mode. Default is 1. 
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
#' ## Other options include modifying the ggplot  
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName") %>% 
#'    trelli_abundance_boxplot(test_mode = T, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "ylim(c(2,20))"))
#' 
#' ## Or making the plot interactive 
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily") %>% trelli_abundance_boxplot(interactive = T)
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
                                     path = "~/Downloads/Trelliscope",
                                     name = "Trelliscope",
                                     test_mode = FALSE,
                                     test_example = 1,
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
                  test_example = test_example)
  
  
  # Remove stat specific options if no stats data was provided 
  if (is.null(trelliData$trelliData.stat)) {
    if (any(c("p_value", "fold_change") %in% cognostics) & is.null(trelliData$trelliData.stat)) {
      cognostics <- cognostics[-match(c("p_value", "fold_change"), cognostics, nomatch = 0)]
    }    
  }
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}

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
      "mean" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = T), 4), desc = "Mean Abundance")), 
      "median" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = T), 4), desc = "Median Abundance")), 
      "sd" = dplyr::tibble(`Standard Deviation Abundance` = trelliscopejs::cog(round(sd(DF$Abundance, na.rm = T), 4), desc = "Abundance Standard Deviation")), 
      "skew" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = T), 4), desc= "Abundance Skewness"))
    )
    
    # If cognostics are any of the cog, then add them 
    if (any(cognostics %in% c("n", "mean", "median", "sd", "skew"))) {
      cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% tibble::tibble()
    } else {cog_to_trelli <- NULL}
    
    # Get fdata cname and group_by selection
    fdata_cname <- pmartR::get_fdata_cname(trelliData$omicsData)
    group_by_choice <- attr(trelliData, "group_by_omics")
    
    # Additional group cognostics can be added only if group_designation was set and
    # trelli_group_by is not the fdata_cname
    if (!is.null(attributes(trelliData$omicsData)$group_DF) & fdata_cname != group_by_choice & !is.null(cog_to_trelli)) {
      
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
    
    # Add cognostics that only apply when stats data is grouped by edata_cname
    edata_cname <- pmartR::get_edata_cname(trelliData$omicsData)
    
    if (!is.null(trelliData$trelliData.stat) && edata_cname == attr(trelliData, "group_by_stat")) {
      
      # Downselect to only stats 
      stat_cogs <- cognostics[cognostics %in% c("fold_change", "p_value")]
      
      if (length(stat_cogs) != 0) {
        
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
            name = paste(Comparison, lapply(name, function(x) {name_converter[[x]]}) %>% unlist()),
            value = round(value, 4)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-Comparison)
        
        # Generate new cognostics 
        new_cogs <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
          quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
        })) %>% tibble::tibble()
        
        # Add new cognostics 
        if (!is.null(cog_to_trelli)) {
          cog_to_trelli <- cbind(cog_to_trelli, new_cogs) %>% tibble::tibble()
        } else {
          cog_to_trelli <- new_cogs
        }
        
      }
      
    }
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # If test_mode is on, then just build the required panels
  if (test_mode) {toBuild <- trelliData$trelliData.omics[test_example,]} else {toBuild <- trelliData$trelliData.omics}
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = box_plot_fun,
                 cogFUN = box_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}

#' @name trelli_abundance_histogram
#' 
#' @title Histogram trelliscope building function for abundance data   
#' 
#' @description Specify a plot design and cognostics for the abundance histogram trelliscope.
#'    Main_effects grouping are ignored. Data must be grouped by edata_cname. 
#' 
#' @param trelliData A trelliscope data object made by as.trelliData or as.trelliData.edata,
#'    and grouped by trelli_group_by. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    n, mean, median, sd, and skew.
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "ylim(c(1,2))"). Default is NULL. 
#' @param interactive A logical argument indicating whether the plots should be interactive
#'    or not. Interactive plots are ggplots piped to ggplotly (for now). Default is FALSE.  
#' @param path The base directory of the trelliscope application. Default is Downloads. 
#' @param name The name of the display. Default is Trelliscope.
#' @param test_mode A logical to return a smaller trelliscope to confirm plot and design.
#'    Default is FALSE.
#' @param test_example The index number of the plot to return for test_mode. Default is 1. 
#' 
#' @examples
#' \dontrun{
#' 
#' ## Build the abundance histogram with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName") %>% 
#'    trelli_abundance_histogram(test_mode = T, test_example = 1:10)
#' 
#' ## Build the abundance boxplot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'    trelli_abundance_histogram(test_mode = T, test_example = 1:10)
#'     
#' ## Build the abundance boxplot with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_group_by(trelliData = trelliData4, group = "LipidCommonName") %>%
#'    trelli_abundance_histogram(test_mode = T, test_example = 1:10)
#'    
#' ## Users can modify the plotting function with ggplot parameters and interactivity, 
#' ## and can also select certain cognostics.     
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName") %>% 
#'    trelli_abundance_histogram(test_mode = T, test_example = 1:10, 
#'      ggplot_params = c("ylab('')", "xlab('Abundance')"), interactive = TRUE,
#'      cognostics = c("mean", "median"))  
#'    
#' }
#' 
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_abundance_histogram <- function(trelliData,
                                     cognostics = c("n", "mean", "median", "sd", "skew", "p_value", "fold_change"),
                                     ggplot_params = NULL,
                                     interactive = FALSE,
                                     path = "~/Downloads/Trelliscope",
                                     name = "Trelliscope",
                                     test_mode = FALSE,
                                     test_example = 1,
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
                  test_example = test_example)

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
  if (edata_cname != attr(trelliData, "group_by_omics")) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Make histogram function-----------------------------------------------------
  
  # First, generate the boxplot function
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
      "mean" = dplyr::tibble(`Mean Abundance` = trelliscopejs::cog(round(mean(DF$Abundance, na.rm = T), 4), desc = "Mean Abundance")), 
      "median" = dplyr::tibble(`Median Abundance` = trelliscopejs::cog(round(median(DF$Abundance, na.rm = T), 4), desc = "Median Abundance")), 
      "sd" = dplyr::tibble(`Standard Deviation Abundance` = trelliscopejs::cog(round(sd(DF$Abundance, na.rm = T), 4), desc = "Abundance Standard Deviation")), 
      "skew" = dplyr::tibble(`Skew Abundance` = trelliscopejs::cog(round(e1071::skewness(DF$Abundance, na.rm = T), 4), desc= "Abundance Skewness"))
    )
    
    # If cognostics are any of the cog, then add them 
    if (any(cognostics %in% c("n", "mean", "median", "sd", "skew"))) {
      cog_to_trelli <- do.call(dplyr::bind_cols, lapply(cognostics, function(x) {cog[[x]]})) %>% tibble::tibble()
    } else {cog_to_trelli <- NULL}
    
    # Add statistics if applicable 
    if (!is.null(trelliData$trelliData.stat)) {
      
      # Downselect to only stats 
      stat_cogs <- cognostics[cognostics %in% c("fold_change", "p_value")]
      
      if (length(stat_cogs) != 0) {
        
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
            name = paste(Comparison, lapply(name, function(x) {name_converter[[x]]}) %>% unlist()),
            value = round(value, 4)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::select(-Comparison)
        
        # Generate new cognostics 
        new_cogs <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
          quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
        })) %>% tibble::tibble()
        
        # Add new cognostics 
        if (!is.null(cog_to_trelli)) {
          cog_to_trelli <- cbind(cog_to_trelli, new_cogs) %>% tibble::tibble()
        } else {
          cog_to_trelli <- new_cogs
        }
        
      }
      
    }
    
    return(cog_to_trelli)
    
  }
  
  # Build trelliscope display---------------------------------------------------
  
  # If test_mode is on, then just build the required panels
  if (test_mode) {toBuild <- trelliData$trelliData.omics[test_example,]} else {toBuild <- trelliData$trelliData.omics}
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = hist_plot_fun,
                 cogFUN = hist_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}


# trelli_abundance_heatmap (emeta only)

# trelli_missingness_bar 

# trelli_foldchange_bar (no emeta only)

# trelli_foldchange_boxplot (emeta only)

# trelli_foldchange_volcano (emeta only)

# trelli_foldchange_heatmap (emeta only)



