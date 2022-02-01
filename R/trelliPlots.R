# This function runs necessary checks for pmartR trelliscope plotting functions. 
# It cleans any parameters (rounding numerics to integers, etc.), and returns them.
trelli_precheck <- function(trelliData, trelliCheck,
                            cognostics, acceptable_cognostics,
                            ggplot_params,
                            interactive, 
                            test_mode,
                            test_example) {
  
  #########################
  ## TEST EXAMPLE CHECKS ##
  #########################
  
  # Ensure that test_example is an integer
  if (!is.numeric(test_example) | 0 %in% test_example) {
    "test_example should be a non-zero integer."
  }
  
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
    
    # Ensure that test_example is in the range of possibilities 
    if (max(test_example) > nrow(trelliData$trelliData.omics)) {
      stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.omics)))
    }
    
  }
  
  # Check that statRes data exists 
  if ("stat" %in% trelliCheck) {
    
    # Assert that trelliData has statRes
    if (is.null(trelliData$trelliData.stat)) {
      stop("trelliData must have statRes for this plotting function.")
    }
    
    # Ensure that test_example is in the range of possibilities 
    if (max(test_example) > nrow(trelliData$trelliData.stat)) {
      stop(paste("test_example must be in the range of possibilities, of 1 to", nrow(trelliData$trelliData.omics)))
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
  
}

# A quick cognostic function 
quick_cog <- function(name, value) {
  dplyr::tibble(!!rlang::sym(name) := trelliscopejs::cog(value, desc = name))
}

# Create a list to convert from short name to long
name_converter_abundance <- list("n" = "Count", "mean" = "Mean Abundance", 
                       "median" = "Median Abundance", "sd" = "Standard Deviation Abundance", 
                       "skew" = "Skew Abundance", "p_value" = "P Value", "fold_change" = "Fold Change")
name_converter_foldchange <- list("n" = "Count", "mean" = "Mean Fold Change", 
                                  "median" = "Median Fold Change", "sd" = "Standard Deviation Fold Change")

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

# Get downloads folder
getDownloadsFolder <- function() {
  if (Sys.info()['sysname'] == "Windows")
    folder <- dirname("~")
  else
    folder <- path.expand("~")
  folder <- file.path(folder, "Downloads")
  folder <- paste0(folder, .Platform$file.sep)
  return(folder)
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
                                     path = getDownloadsFolder(),
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
#'    and grouped by edata_cname in trelli_group_by. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    n, mean, median, sd, and skew. p_value and fold_change can be added if statRes
#'    is included.
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
#' ## Build the abundance histogram with an omicsData object. Generate trelliData in as.trelliData
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'    trelli_abundance_histogram(test_mode = T, test_example = 1:10)
#'     
#' ## Build the abundance histogram with an omicsData and statRes object. Generate trelliData in as.trelliData.
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
                                       path = getDownloadsFolder(),
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

#' @name trelli_abundance_heatmap
#' 
#' @title Heatmap trelliscope building function for abundance data   
#' 
#' @description Specify a plot design and cognostics for the abundance heatmap trelliscope.
#'    Data must be grouped by an emeta column. Main_effects order the y-variables. All 
#'    statRes data is ignored. 
#' 
#' @param trelliData A trelliscope data object made by as.trelliData,
#'    and grouped by an emeta variable. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    n, mean, median, sd, and skew per "main_effects" group designation. Otherwise, 
#'    no cognostics are returned. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "xlab('')"). Default is NULL. 
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
#' ## Build the abundance heatmap with an omicsData object with emeta variables. Generate trelliData in as.trelliData.
#' trelli_group_by(trelliData = trelliData2, group = "LipidFamily") %>% 
#'    trelli_abundance_heatmap(test_mode = T, test_example = 1:3)
#'    
#' ## Users can modify the plotting function with ggplot parameters and interactivity, 
#' ## and can also select certain cognostics.     
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily") %>% 
#'    trelli_abundance_heatmap(test_mode = T, test_example = 1:5, 
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
                                     path = getDownloadsFolder(),
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
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is grouped by an e_meta variable
  if (attr(trelliData, "group_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be grouped_by an e_meta column.")
  }
  
  # If no group designation, set cognostics to NULL. 
  if (is.null(attributes(trelliData2$omicsData)$group_DF)) {
    congnostics <- NULL
  }
  
  # Make heatmap function-------------------------------------------------------
  
  # First, generate the heatmap function
  hm_plot_fun <- function(DF, title) {
    
    # If group designation was set, then convert Group to a factor variable 
    if (is.null(attributes(trelliData$omicsData)$group_DF)) {
      DF$Group <- factor(DF$Sample_Name, levels = attributes(trelliData2$omicsData)$group_DF$Sample_Name)
    } 
    
    # Get edata cname
    edata_cname <- get_edata_cname(trelliData$omicsData)
    
    # Build plot: this should be edata_cname
    hm <- ggplot2::ggplot(DF, ggplot2::aes(x = .data[[edata_cname]], y = Sample_Name, fill = Abundance)) +
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
    cog_to_trelli <- do.call(cbind, lapply(1:nrow(cogs_to_add), function(row) {
      quick_cog(cogs_to_add$name[row], cogs_to_add$value[row])
    })) %>% tibble::tibble() 
    
    return(cog_to_trelli)
    
  } 
  
  # Build trelliscope display---------------------------------------------------
  
  # If test_mode is on, then just build the required panels
  if (test_mode) {toBuild <- trelliData$trelliData.omics[test_example,]} else {toBuild <- trelliData$trelliData.omics}
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = hm_plot_fun,
                 cogFUN = hm_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}

#' @name trelli_missingness_bar
#' 
#' @title Bar chart trelliscope building function for missing data   
#' 
#' @description Specify a plot design and cognostics for the missing barchart trelliscope.
#'    Missingness is displayed per group_by variable. Main_effects data is used to
#'    split samples when applicable. 
#' 
#' @param trelliData A trelliscope data object made by as.trelliData.edata or as.trelliData. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    is n and proportion.
#' @param proportion A logical to determine whether plots should display counts or proportions.
#'    Default is TRUE.     
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "xlab('')"). Default is NULL. 
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
#' ## Build the missingness bar plot with an edata file. Generate trelliData in as.trelliData.edata
#' trelli_group_by(trelliData = trelliData, group = "LipidCommonName") %>% 
#'   trelli_missingness_bar(test_mode = T, test_example = 1:10)
#' trelli_group_by(trelliData = trelliData, group = "Sample") %>% trelli_missingness_bar()
#' 
#' ## Build the missingness bar plot with an omicsData object. Generate trelliData in as.trelliData
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'   trelli_missingness_bar(test_mode = T, test_example = 1:10)
#' 
#' ## Build the missingness bar plot with a statRes object. Generate trelliData in as.trelliData
#' trelli_group_by(trelliData = trelliData3, group = "LipidCommonName") %>%
#'   trelli_missingness_bar(test_mode = T, test_example = 1:10)
#' 
#' ## Build the missingness bar plot with an omicsData and statRes object. Generate trelliData in as.trelliData.
#' trelli_group_by(trelliData = trelliData4, group = "LipidCommonName") %>%
#'   trelli_missingness_bar(test_mode = T, test_example = 1:10) 
#' 
#' ## Or making the plot interactive 
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'    trelli_missingness_bar(test_mode = T, test_example = 1:5)
#'    
#' ## Or visualize only count data 
#' trelli_group_by(trelliData = trelliData2, group = "LipidCommonName") %>% 
#'    trelli_missingness_bar(test_mode = T, test_example = 1:5, cognostics = "n", proportion = FALSE)
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
                                   path = getDownloadsFolder(),
                                   name = "Trelliscope",
                                   test_mode = FALSE,
                                   test_example = 1,
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
                  test_example = test_example)
  
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
    if (test_mode) {toBuild <- trelliData$trelliData.omics[test_example,]} else {toBuild <- trelliData$trelliData.omics}
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
    totalCounts <- attr(statRes, "group_DF")$Group %>% 
      table(dnn = "Group") %>% 
      data.frame() %>% 
      dplyr::rename(Total = Freq)
    
    if (test_mode) {toBuild <- toBuild[test_example,]} 
  }
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = missing_bar_plot_fun,
                 cogFUN = missing_bar_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}

#' @name trelli_foldchange_bar
#' 
#' @title Bar chart trelliscope building function for fold_change   
#' 
#' @description Specify a plot design and cognostics for the fold_change barchart trelliscope.
#'    Fold change must be grouped by edata_cname. 
#' 
#' @param trelliData A trelliscope data object with statRes results. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    are fold_change and p_value.
#' @param p_value_thresh A value between 0 and 1 to indicate a threshold to highlight
#'    significant biomolecules. Default is 0.05. Selecting 0 will remove this feature. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "xlab('')"). Default is NULL. 
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
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_group_by(trelliData = trelliData3, group = "LipidCommonName") %>% 
#'   trelli_foldchange_bar(test_mode = T, test_example = 1:10)
#'   
#' ## Or make the plot interactive  
#' trelli_group_by(trelliData = trelliData4, group = "LipidCommonName") %>% 
#'   trelli_foldchange_bar(test_mode = T, test_example = 1:10, interactive = T) 
#'    
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_bar <- function(trelliData,
                                  cognostics = c("fold_change", "p_value"),
                                  p_value_thresh = 0.05,
                                  ggplot_params = NULL,
                                  interactive = FALSE,
                                  path = getDownloadsFolder(),
                                  name = "Trelliscope",
                                  test_mode = FALSE,
                                  test_example = 1,
                                  ...) {
  
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = "stat",
                  cognostics = cognostics,
                  acceptable_cognostics = c("p_value", "fold_change"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is edata_cname
  edata_cname <- pmartR::get_edata_cname(trelliData$statRes)
  if (edata_cname != attr(trelliData, "group_by_stat")) {
    stop("trelliData must be grouped by edata_cname.")
  }
  
  # Check p_value threshold
  if (!is.numeric(p_value_thresh)) {
    stop("p_value_thresh must be a numeric.")
  }
  if (p_value_thresh < 0 | p_value_thresh > 1) {
    stop("p_value_thresh must be between 0 and 1.")
  }
  
  # Make foldchange bar function------------------------------------------------
  
  fc_bar_plot_fun <- function(DF, title) {
    
    # Change NaN to 0 just for the plotting functions
    DF$p_value[is.nan(DF$p_value)] <- 0
    DF$fold_change[is.nan(DF$fold_change)] <- 0
    
    # Indicate which comparisons should be highlighted
    Significant <- ifelse(DF$p_value <= p_value_thresh & DF$p_value != 0, "black", NA)
    
    # Make bar plot 
    bar <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
      ggplot2::geom_bar(stat = "identity", color = Significant) + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5), 
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none"
      ) + ggplot2::ylab("Fold Change")
    
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
  
  # Subset down to test example if applicable
  if (test_mode) {toBuild <- trelliData$trelliData.stat[test_example,]} else {toBuild <- trelliData$trelliData.stat}
  
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
#' @description Specify a plot design and cognostics for the fold_change boxplot trelliscope.
#'    Fold change must be grouped by an emeta column, which means both an omicsData
#'    object and statRes are required to make this plot. 
#' 
#' @param trelliData A trelliscope data object with omicsData and statRes results. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entries are
#'    are n ,mean, median, and sd. 
#' @param p_value_thresh A value between 0 and 1 to indicate a threshold to highlight
#'    significant biomolecules. Default is 0.05. Selecting 0 will remove this feature. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "xlab('')"). Default is NULL. 
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
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily") %>% 
#'   trelli_foldchange_boxplot()
#'
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_boxplot <- function(trelliData,
                                      cognostics = c("n", "median", "mean", "sd"),
                                      p_value_thresh = 0.05,
                                      ggplot_params = NULL,
                                      interactive = FALSE,
                                      path = getDownloadsFolder(),
                                      name = "Trelliscope",
                                      test_mode = FALSE,
                                      test_example = 1,
                                      ...) {
  
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("omics", "stat"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("n", "median", "mean", "sd"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is an emeta column
  if (attr(trelliData, "group_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be grouped_by an e_meta column.")
  }
  
  # Check p_value threshold
  if (!is.numeric(p_value_thresh)) {
    stop("p_value_thresh must be a numeric.")
  }
  if (p_value_thresh < 0 | p_value_thresh > 1) {
    stop("p_value_thresh must be between 0 and 1.")
  }
  
  # Make foldchange boxplot function--------------------------------------------
  
  fc_box_plot_fun <- function(DF, title) {
    
    # Change NaN to 0 just for the plotting functions
    DF$p_value[is.nan(DF$p_value)] <- 0
    DF$fold_change[is.nan(DF$fold_change)] <- 0
    
    # Indicate which comparisons should be highlighted
    Significant <- ifelse(DF$p_value <= p_value_thresh & DF$p_value != 0, "red", "black")
    
    # Make boxplot
    boxplot <- ggplot2::ggplot(DF, ggplot2::aes(x = Comparison, y = fold_change, fill = Comparison)) +
      ggplot2::geom_boxplot() + ggplot2::geom_point(color = Significant) + ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5), 
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none"
      ) + ggplot2::ylab("Fold Change") + ggplot2::ggtitle(title)
    
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
        "mean" = round(mean(fold_change, na.rm = T), 4),
        "median" = round(median(fold_change, na.rm = T), 4),
        "sd" = round(sd(fold_change, na.rm = T), 4)
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
  
  # Subset down to test example if applicable
  if (test_mode) {toBuild <- trelliData$trelliData.stat[test_example,]} else {toBuild <- trelliData$trelliData.stat}
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = fc_box_plot_fun,
                 cogFUN = fc_box_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}

#' @name trelli_foldchange_volcano
#' 
#' @title Volcano trelliscope building function for fold_change   
#' 
#' @description Specify a plot design and cognostics for the fold_change volcano trelliscope.
#'    Fold change must be grouped by an emeta column, which means both an omicsData
#'    object and statRes are required to make this plot. 
#' 
#' @param trelliData A trelliscope data object with omicsData and statRes results. Required. 
#' @param cognostics A vector of cognostic options for each plot. Valid entry is n. 
#' @param p_value_thresh A value between 0 and 1 to indicate a threshold to highlight
#'    significant biomolecules. Default is 0.05. Selecting 0 will remove this feature. 
#' @param ggplot_params An optional vector of strings of ggplot parameters to the backend ggplot
#'    function. For example, c("ylab('')", "xlab('')"). Default is NULL. 
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
#' ## Build fold_change bar plot with statRes data grouped by edata_colname.
#' trelli_group_by(trelliData = trelliData4, group = "LipidFamily") %>% 
#'   trelli_foldchange_volcano()
#'
#' }
#' 
#' @author David Degnan, Lisa Bramer
#' 
#' @export
trelli_foldchange_volcano <- function(trelliData,
                                      cognostics = c("n"),
                                      p_value_thresh = 0.05,
                                      ggplot_params = NULL,
                                      interactive = FALSE,
                                      path = getDownloadsFolder(),
                                      name = "Trelliscope",
                                      test_mode = FALSE,
                                      test_example = 1,
                                      ...) {
  
  # Run initial checks----------------------------------------------------------
  
  # Run generic checks 
  trelli_precheck(trelliData = trelliData, 
                  trelliCheck = c("omics", "stat"),
                  cognostics = cognostics,
                  acceptable_cognostics = c("n"),
                  ggplot_params = ggplot_params,
                  interactive = interactive,
                  test_mode = test_mode, 
                  test_example = test_example)
  
  # Round test example to integer 
  if (test_mode) {test_example <- unique(abs(round(test_example)))}
  
  # Check that group data is an emeta column
  if (attr(trelliData, "group_by_omics") %in% attr(trelliData, "emeta_col") == FALSE) {
    stop("trelliData must be grouped_by an e_meta column.")
  }
  
  # Check p_value threshold
  if (!is.numeric(p_value_thresh)) {
    stop("p_value_thresh must be a numeric.")
  }
  if (p_value_thresh < 0 | p_value_thresh > 1) {
    stop("p_value_thresh must be between 0 and 1.")
  }
  
  # Make foldchange volcano function--------------------------------------------
  
  fc_volcano_plot_fun <- function(DF, title) {
    
    # Indicate which comparisons should be highlighted
    DF$Significant <- lapply(1:nrow(DF), function(row) {
      if (!is.nan(DF$p_value[row]) && !is.nan(DF$fold_change[row]) && DF$p_value[row] <= p_value_thresh) {
        ifelse(DF$fold_change[row] > 0, "High", "Low")
      } else {return("Not Significant")}
    }) %>% unlist()
    
    # Make volcano plot
    volcano <- ggplot2::ggplot(DF, ggplot2::aes(x = fold_change, y = -log10(p_value), color = Significant)) +
      ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
      ggplot2::scale_color_manual(values = c("Low" = "blue", "Not Significant" = "black", "High" = "red")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::xlab("Fold Change") + ggplot2::ylab("-Log10 P Value")
    
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
    
    # Get count 
    counts <- lapply(1:nrow(DF), function(row) {
      if (!is.nan(DF$p_value[row]) && !is.nan(DF$fold_change[row]) && DF$p_value[row] <= p_value_thresh) {
        ifelse(DF$fold_change[row] > 0, "High", "Low")
      } else {return("Not Significant")}
    }) %>% 
      unlist() %>%
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
  
  # Subset down to test example if applicable
  if (test_mode) {toBuild <- trelliData$trelliData.stat[test_example,]} else {toBuild <- trelliData$trelliData.stat}
  
  # Pass parameters to trelli_builder function
  trelli_builder(toBuild = toBuild,
                 cognostics = cognostics, 
                 plotFUN = fc_volcano_plot_fun,
                 cogFUN = fc_volcano_cog_fun,
                 path = path,
                 name = name,
                 ...)
  
}


# trelli_foldchange_heatmap (emeta only)



