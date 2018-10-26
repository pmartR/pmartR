#' Produce a plot of a pmartR S3 Object
#'
#' This function will provide a plot for a \code{omicsData} object, any of the filter objects in pmartR, or a \code{corRes} object.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'lipidData', or 'metabData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param corRes_object an object of class corRes. A correlation matrix of all samples.
#' @param filter_object a filter object for the respective \code{omicsData} class.
#' @param ... further arguments
#'
#' @return a ggplot summarizing the pmartR object
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' cor_matrix <- cor_result(pep_object)
#' plot(cor_matrix)
#' plot(cor_matrix, interactive = TRUE)
#'
#' to_filter <- molecule_filter(pep_object)
#' plot(to_filter, min_num = 4)
#'
#' data("pro_object")
#' plot(pro_object)
#' plot(pro_object, order_by = "Condition")
#' plot(pro_object, color_by = "Condition")
#' plot(pro_object, facet_by = "Condition")
#' plot(pro_object, facet_by = "Condition", facet_cols = 1)
#'
#' pep_object2 <- group_designation(pep_object, main_effects = "Condition")
#' pep_object2 <- edata_transform(pep_object2, "log2")
#' pca_res <- dim_reduction(pep_object2, k = 2)
#' plot(pca_res)
#'}
#' @details
#' Various further arguments can be specified depending on the class of the object being plotted.
#'
#'
#'
#' @name plot-pmartR
#'
NULL


#' plot.corRes
#' 
#' For plotting an S3 object of type 'corRes':
#' 
#'@rdname plot-pmartR-corRes
#'@export
#'@param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#'@param interactive default value is FALSE. If TRUE, an interactive d3heatmap will be rendered, allowing you to zoom in on sections of the heatmap.
#'@param x_lab logical indicator of whether to label the x-axis with the sample names. Defaults to TRUE. If FALSE, no x-axis labels are included in the plot.
#'@param y_lab logical indicator of whether to label the y-axis with the sample names. Defaults to TRUE. If FALSE, no y-axis labels are included in the plot.
#'@param colorbar_lim numeric pair of numeric values specifying the minimum and maximum values to use in the heatmap color bar. Defaults to 'c(NA, NA)', in which case ggplot2 automatically sets the minimum and maximum values based on the correlation values in the data.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case "Sample Outlier Results /n p-value threshold = 0.xyz" is used, where 'xyz' is the pvalue_threshold supplied to the function. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{use_VizSampNames} \tab logical specifies whether to use custom sample names \cr
#' } 
plot.corRes <- function(corRes_object, omicsData = NULL, interactive = FALSE, x_lab=TRUE, y_lab=TRUE, colorbar_lim=c(NA, NA), ...){
  .plot.corRes(corRes_object, omicsData, interactive, x_lab, y_lab, colorbar_lim, ...)
}
.plot.corRes <- function(corRes_object, omicsData = NULL, interactive = FALSE, x_lab=TRUE, y_lab=TRUE, colorbar_lim=c(NA, NA), title_plot = NULL, title_size = 14, use_VizSampNames = FALSE){
  # check for a corRes object #
  if(!inherits(corRes_object, "corRes")) stop("object must be of class 'corRes'")
  
  # check title and colorbar options #
  if(!is.null(title_plot)) {
    if(!is.character(title_plot)) stop("title_plot must be a character vector")
  }
  if(!all(is.na(colorbar_lim))){
    if(!is.numeric(colorbar_lim) | length(colorbar_lim) != 2) stop("colorbar_lim must be a numeric vector of length 2")
  }
  if(!is.numeric(title_size)) message("Title size should be a numeric value, a default value will be used")
  
  # check that omicsData is of appropriate class #
  if(!is.null(omicsData)){
    if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  }
  # Workaround for certain "check" warnings
  Var1 <- Var2 <- value <- NULL
  
  if(is.null(title_plot)){
    if(!is.null(attributes(corRes_object)$data_norm)){
      if(attributes(corRes_object)$data_norm == TRUE) {
        plot_title <- "Correlations Among Samples (Normalized Data)"
      }else{
        if(attributes(corRes_object)$data_norm == FALSE){
          plot_title <- "Correlations Among Samples (Un-Normalized Data)"
        }
      }
      
    }else{
        plot_title <- "Correlations Among Samples (Un-Normalized Data)"
    }
    
  }else{
    plot_title <- title_plot
  }
  
  if(!interactive) {
    #pal <- colorRampPalette(c("blue","white","red"))
    pal <- colorRampPalette(c("red", "yellow")) # modified by KS 2/12/2018

    colnames(corRes_object) <- rownames(corRes_object)
    corRes_melt <- reshape2::melt(corRes_object)
    corRes_melt$Var1 <- abbreviate(corRes_melt$Var1, minlength=20)
    corRes_melt$Var2 <- abbreviate(corRes_melt$Var2, minlength=20)
    
    if(all(is.na(colorbar_lim))){
      heatmap <- ggplot2::ggplot(corRes_melt, ggplot2::aes(x = ordered(Var2, levels = rev(sort(unique(Var2)))),
                                                           y = ordered(Var1, levels = rev(sort(unique(Var1)))))) +
        ggplot2::geom_tile(ggplot2::aes(fill = value)) +
        ggplot2::scale_fill_gradientn("Correlation", colours = pal(50)) +
        ggplot2::xlab("") +  ggplot2::ylab("") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                       axis.ticks = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size))
    }else{
      heatmap <- ggplot2::ggplot(corRes_melt, ggplot2::aes(x = ordered(Var2, levels = rev(sort(unique(Var2)))),
                                                           y = ordered(Var1, levels = rev(sort(unique(Var1)))))) +
        ggplot2::geom_tile(ggplot2::aes(fill = value)) +
        ggplot2::scale_fill_gradientn("Correlation", limits = colorbar_lim, colours = pal(50)) +
        ggplot2::xlab("") +  ggplot2::ylab("") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90),
                       axis.ticks = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size))
    }
    
    
    if(nrow(corRes_object) > 40) heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    if(x_lab==FALSE) heatmap <- heatmap + ggplot2::theme(axis.text.x = ggplot2::element_blank()) # added by KS
    if(y_lab==FALSE) heatmap <- heatmap + ggplot2::theme(axis.text.y = ggplot2::element_blank()) # added by KS
    
    if(!is.null(attributes(corRes_object)$data_norm)){
      if(attributes(corRes_object)$data_norm == TRUE) {
        heatmap <- heatmap + ggplot2::ggtitle(plot_title)
      }else{heatmap <- heatmap + ggplot2::ggtitle(plot_title)}
    }else{
        heatmap <- heatmap + ggplot2::ggtitle(plot_title)
    }
    
    if(use_VizSampNames){
      if(!is.logical(as.logical(use_VizSampNames))) stop("Must specify a logical TRUE/FALSE value for use_VizSampNames")
      if(is.null(omicsData)) stop("If using custom sample names, specify the corresponding omicsData object")
      heatmap = heatmap + 
        ggplot2::scale_x_discrete(labels = omicsData$f_data$VizSampNames, breaks = levels(omicsData$f_data[,get_fdata_cname(omicsData)])) + 
        ggplot2::scale_y_discrete(labels = omicsData$f_data$VizSampNames, breaks = levels(omicsData$f_data[,get_fdata_cname(omicsData)]))
    }
    
  } else {
    if(!is.null(title_plot)) message("The ability to display a plot title is not available when interactive is set to TRUE")
    if(use_VizSampNames) message("Custom Sample names not available in interactive plots")
    heatmap <- d3heatmap::d3heatmap(corRes_object, dendrogram = 'none', reorderfun = function(x) ordered(x, levels = rev(sort(unique(x)))), title = plot_title)
  }
  
  return(heatmap)
}

#' plot.moleculeFilt
#' 
#' For plotting an S3 object of type 'moleculeFilt':
#' 
#'@export
#'@rdname plot-pmartR-moleculeFilt
#'@param min_num integer value specifying the minimum number of times each feature must be observed across all samples.  If a value is specified, a horizontal line will be drawn when \code{cumulative=TRUE}, and bars will be colored appropriatel if \code{cumulative=FALSE}.  Defaults to NULL.
#'@param cumulative logical indicating whether the number of biomolecules observed in \emph{at least} (TRUE) x number of samples or \emph{exactly} (FALSE) x number of samples should be plotted.  Defaults to TRUE.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
#' 
#' @examples
#' data(pep_object)
#' molfilt <- molecule_filter(pep_object)
#' plot(molfilt, min_num = 5)
#' plot(molfilt, min_num = 3, cumulative = FALSE)
#' 
plot.moleculeFilt <- function(filter_object, min_num = NULL, cumulative = TRUE,...) {
  .plot.moleculeFilt(filter_object, min_num, cumulative, ...)
}

.plot.moleculeFilt <- function(filter_object, min_num = NULL, cumulative = TRUE, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {
  
  ## initial checks ##
  if(!is.null(min_num)) {
    # check that min_num is numeric >= 0 #
    if(!inherits(min_num, c("numeric","integer")) | min_num < 0) stop("min_num must be an integer >= 0")
    # check that min_num is an integer #
    if(min_num %% 1 != 0) stop("min_num must be an integer >= 0")
    # check that min_num is less than the max number of observations #
    if(min_num > max(filter_object$Num_Observations)) stop("min_num cannot be greater than the number of samples")
  }
  if(!is.null(title_plot)) {
    if(!is.character(title_plot)) stop("title_plot must be a character vector")
  }
  if(!is.null(x_lab)) {
    if(!is.character(x_lab)) stop("x_lab must be a character vector")
  }
  if(!is.null(y_lab)) {
    if(!is.character(y_lab)) stop("y_lab must be a character vector")
  }
  
  if(!inherits(cumulative, "logical")) stop("argument cumulative must be either TRUE or FALSE")
  
  ## end of initial checks ##
  
  # make counts, colors, and plot shape based on value of cumulative
  
  # values that will be updated depending on input
  fill <- 0
  legend_scale <- NULL
  hline <- NULL
  
  if(cumulative) {
    # cumulative counts (>=)
    counts <- sapply(1:max(filter_object$Num_Observations), function(i){
        filter_object[filter_object$Num_Observations >= i,]$Num_Observations %>% length()
      })
    # append 1 to extend last step (looks awkward without this)
    counts <- c(counts, counts[length(counts)])
    # 1 appended step
    num_obs <- 1:(max(filter_object$Num_Observations)+1)
    
    # shape is a step function with fixed color
    shape <- ggplot2::geom_step(ggplot2::aes(x = num_observations, y = counts), color = "red")
    
    # draw a horizontal line if min_num is specified
    hline <- if(!is.null(min_num)) ggplot2::geom_hline(ggplot2::aes(color = "black"), yintercept = counts[min_num], linetype = "dashed") else NULL
    
    xlabel <- ifelse(is.null(x_lab), "Minimum Number of Times a Biomolecule Appears Across All Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Count of biomolecules observed in at least X number of samples", title_plot)
  }
  else if(!cumulative){
   # counts for a specific number of nonmissing biomolecules (==)
   counts <- sapply(1:max(filter_object$Num_Observations), function(i){
     filter_object[filter_object$Num_Observations == i,]$Num_Observations %>% length()
   })
   
   # color by which values are kept if min_num is specified
   if(!is.null(min_num)){
     fill <- ifelse(1:max(filter_object$Num_Observations) >= min_num, "retained", "dropped")
     shape <- ggplot2::geom_bar(ggplot2::aes(x = num_observations, y = counts, fill = fill), color = "black", stat = "identity")
     legend_scale <- ggplot2::scale_fill_manual(name = "", values = c("dropped" =  "red", "retained" =  "green"))
     }
   else{
     shape <- ggplot2::geom_bar(ggplot2::aes(x = num_observations, y = counts), fill = "royalblue1", color = "black", stat = "identity")
   }
   
   num_obs <- 1:max(filter_object$Num_Observations)

   xlabel <- ifelse(is.null(x_lab), "Number of Times a Biomolecule Appears Across All Samples", x_lab)
   ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)
   plot_title <- ifelse(is.null(title_plot), "Count of biomolecules observed in exactly X number of samples", title_plot)
  }
  
  # create plotting dataframe
  pep_observation_counts <- data.frame(num_observations = num_obs, 
                                       frequency_counts = counts,
                                       fill = fill)
  
  # value annotation - data argument is to shave off last row that is only there to extend the last step
  text <- ggplot2::geom_text(data = pep_observation_counts[1:max(filter_object$Num_Observations),], 
                             ggplot2::aes(x = num_observations, y = frequency_counts, label = frequency_counts), vjust = -1, nudge_x = ifelse(cumulative, 0.5, 0))
  
  # plot #
  p <- ggplot2::ggplot(pep_observation_counts) +
    shape + hline + text + legend_scale +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous(breaks = max(filter_object$Num_Observations):1) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size),
                   axis.title.x = ggplot2::element_text(size = x_lab_size),
                   axis.title.y = ggplot2::element_text(size = y_lab_size))
  
  if(bw_theme == TRUE){
    p <- p + ggplot2::theme_bw()  
  }
  
  return(p)
}

#' plot.proteomicsFilt
#' 
#' For plotting an S3 object of type 'proteomicsFilt':
#' 
#'@export
#'@rdname plot-pmartR-proteomicsFilt
#'@param mapping whether to display a histogram of the number of peptides mapping to a protein for all proteins (mapping = 'pep_to_pro'), a histogram of the number of proteins mapped to by each peptide (mapping = "pro_to_pep"), or 'both'.  Defaults to "both"
#'@param cumulative logical specifying whether the peptide to protein mapping should be cumulative or at each value i.e. number of proteins with \emph{at least} x (TRUE) peptides mapping to it, or \emph{exactly} x (FALSE) peptides mapping to it.  Defaults to TRUE
#'@param log_scale TRUE or FALSE specifying whether or not to plot the x-axis on the log scale, defaults to TRUE
#'@param min_num_peps an optional integer value between 1 and the maximum number of peptides that map to a protein in the data. The value specifies the minimum number of peptides that must map to a protein. Any protein with less than \code{min_num_peps} mapping to it will be returned as a protein that should be filtered. Default value is NULL.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{xlim} \tab length 2 numeric vector specifying x-axis limits for the peptide to protein plot.  Defaults to NULL, which gives default ggplot2 axis limits \cr
#' \code{ylim} \tab length 2 numeric vector specifying y-axis limits for the peptide to protein plot.  Defaults to NULL, which gives default ggplot2 axis limits \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
#' 
#' @examples
#' data(pep_object)
#' profilt <- proteomics_filter(pep_object)
#' plot(profilt, min_num_peps = 5)
#' plot(profilt, min_num_peps = 10, cumulative = FALSE)
#' 
plot.proteomicsFilt <- function(filter_object, mapping = "both", cumulative = TRUE, log_scale = TRUE, min_num_peps = NULL, ...) {
  .plot.proteomicsFilt(filter_object, mapping, cumulative, log_scale, min_num_peps, ...)
}

.plot.proteomicsFilt <- function(filter_object, mapping = "both", cumulative = TRUE, 
                                 log_scale = TRUE, min_num_peps = NULL,
                                 x_lab_pep = NULL, y_lab_pep = NULL, title.pep = NULL,
                                 x_lab_pro = NULL, y_lab_pro = NULL, title.pro = NULL,
                                 xlim = NULL, ylim = NULL,
                                 title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {
  # Error Checks
  if(!is.null(min_num_peps)) {
    # check that min_num_peps is numeric and >=1 #
    if(!inherits(min_num_peps, "numeric") | min_num_peps < 1) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if(min_num_peps %% 1 != 0) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if(length(min_num_peps) != 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if(min_num_peps > sum(filter_object$counts_by_pep$n)) stop("min_num_peps cannot be greater than the total number of peptides")
  }
  # check that cumulative and log_scale are logical #
  if(!inherits(cumulative, "logical")) stop("argument cumulative must be either TRUE or FALSE")
  if(!inherits(log_scale, "logical")) stop("argument log_scale must be either TRUE or FALSE")
  
  # check that a valid mapping argument is supplied
  if(!(mapping %in% c("both", "pep_to_pro", "pro_to_pep"))) stop("'mapping' argument must be one of 'both', 'pep_to_pro' or 'pro_to_pep'")
  
  # suppress ggplot2 warnings
  suppressMessages({
    
    #intantiate bin values for both mappings
    pep_bins = sort(unique(filter_object$counts_by_pep$n)) # bins for peptide COUNTS.  These bins represent a number of PROTEINS
    pro_bins = sort(unique(filter_object$counts_by_pro$n)) # bins for protein COUNTS.  These bins represent a number of PEPTIDES
    
    # fill value, fill labels, and horizontal line values that will be given real values if certain arguments are supplied
    fill <- 0
    fill_format <- NULL
    hline <- NULL
    
    # construct resurces for cumulative (step function) plot
    if(cumulative){
      # get counts of peptides that are mapped to by AT LEAST the number of proteins given in pep_bins
      pep_counts <- sapply(pep_bins, function(x){
        filter_object$counts_by_pep[filter_object$counts_by_pep$n >= x,] %>% nrow()
      })
      
      # get counts of proteins that are mapped to by AT LEAST the number of peptides given in pro_bins
      pro_counts <- sapply(pro_bins, function(x){
        filter_object$counts_by_pro[filter_object$counts_by_pro$n >= x,] %>% nrow()
      })
      
      # only use step function if there is more than 1 step
      if(length(pep_counts) > 1){
        # an extra count to extend the last step
        pep_counts <- c(pep_counts, pep_counts[length(pep_counts)])
        pep_bins <- c(pep_bins, max(pep_bins)+1)
        
        # store shape variable to be added to plot object later
        shape_pep <- ggplot2::geom_step(ggplot2::aes(x = bins, y = counts), color = "red")
      }
      else shape_pep <- ggplot2::geom_bar(ggplot2::aes(x = bins, y = counts), fill = "royalblue1", stat = "identity")
      
      if(length(pro_counts) > 1){
        pro_counts <- c(pro_counts, pro_counts[length(pro_counts)])
        pro_bins <- c(pro_bins, max(pro_bins)+1)
        shape_pro <- ggplot2::geom_step(ggplot2::aes(x = bins, y = counts), color = "red")
      }
      else shape_pro <- ggplot2::geom_bar(ggplot2::aes(x = bins, y = counts), fill = "royalblue1", stat = "identity")
      
      # if min_num_peps is specified, draw a horizontal line to indicate the number of peptides retained
      if(!is.null(min_num_peps)){
        hline <- ggplot2::geom_hline(ggplot2::aes(color = "black"), yintercept = pro_counts[which(pro_bins == min_num_peps)], linetype = "dashed")
      }
    
    }
    # get resources for bar chart (not cumulative)
    else if(!cumulative){
      # get counts of peptides that are mapped to by EXACTLY the number of proteins given in pep_bins
      pep_counts <- sapply(pep_bins, function(x){
        filter_object$counts_by_pep[filter_object$counts_by_pep$n == x,] %>% nrow()
      })
      # get counts of proteins that are mapped to by EXACTLY the number of peptides given in pro_bins
      pro_counts <- sapply(pro_bins, function(x){
        filter_object$counts_by_pro[filter_object$counts_by_pro$n == x,] %>% nrow()
      })
      
      # if min_num_peps is specified, add a coloring variable that is red for dropped values and green for retained values
      if(!is.null(min_num_peps)){
        fill <- ifelse(pro_bins >= min_num_peps, "retained", "dropped")
        shape_pro <- ggplot2::geom_bar(ggplot2::aes(x = bins, y = counts, fill = fill), stat = "identity")
        fill_format <- ggplot2::scale_fill_manual(name = "", values = c("dropped" =  "red", "retained" =  "green"))
      }
      else{
        shape_pro <- ggplot2::geom_bar(ggplot2::aes(x = bins, y = counts), fill = "royalblue1", stat = "identity")
      }
      
      shape_pep <- ggplot2::geom_bar(ggplot2::aes(x = bins, y = counts), fill = "royalblue1", stat = "identity") 
    }
    
    # log scaling
    if(log_scale){
      scale_pep <- ggplot2::scale_x_continuous(breaks = unique(pep_bins), trans = "log10")
      scale_pro <- ggplot2::scale_x_continuous(trans = "log10")
    }  
    else{
      scale_pro <- NULL
      scale_pep <- ggplot2::scale_x_continuous(breaks = seq(1:max(pep_counts)))
    }
  
    # make plotting df for peptides and proteins
    pep_counts_df <- data.frame(counts = pep_counts, bins = pep_bins)
    pro_counts_df <- data.frame(counts = pro_counts, bins = pro_bins, fill = fill)
    
    # make peptide counts plot
    if(mapping %in% c("both", "pep_to_pro")){
      ## make labels ##
      xlabel_pep <- ifelse(is.null(x_lab_pep), "Number of Peptides Mapped to a Protein", x_lab_pep)
      ylabel_pep <- ifelse(is.null(y_lab_pep), "Count of Proteins", y_lab_pep)
      plot_title_pep <- ifelse(is.null(title.pep), paste0("Number of proteins mapped to by ", ifelse(cumulative, "at least", "exactly")," X peptides"), title.pep)
      
      # main plot object that adds resources built in ifelse blocks above
      p1 <- ggplot2::ggplot(pro_counts_df) +
        shape_pro + scale_pro + hline + fill_format +
        ggplot2::ggtitle(plot_title_pep) +
        ggplot2::xlab(xlabel_pep) +
        ggplot2::ylab(ylabel_pep) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size))
      # other ggplot2 arguments
      if(bw_theme){
        p1 <- p1 + ggplot2::theme_bw()
      }
      if(!is.null(xlim)){
        p1 <- p1 + ggplot2::xlim(xlim)
      }
      if(!is.null(ylim)){
        p1 <- p1 + ggplot2::ylim(ylim)
      }
    }
    
    # make protein counts plot
    if(mapping %in% c("both","pro_to_pep")){
      xlabel_pro <- ifelse(is.null(x_lab_pro), "Number of Proteins Mapped to a Peptide", x_lab_pro)
      ylabel_pro <- ifelse(is.null(y_lab_pro), "Count of Peptides", y_lab_pro)
      plot_title_pro <- ifelse(is.null(title.pep), paste0("Number of peptides mapped to by ", ifelse(cumulative, "at least", "exactly")," X proteins"), title.pep)
      
      # text annotation (peptide count plot has too many bins to make this an option) 
      # data argument shaves off last row that is only used to extend last step in cumulative plot
      text <- ggplot2::geom_text(data = pep_counts_df[pep_counts_df$bins <= max(filter_object$counts_by_pep$n),], ggplot2::aes(x = bins, y = counts, label = counts), vjust = -1, hjust = ifelse(cumulative & (length(pep_counts) > 1), -0.5, 0.5))
      
      p2 <- ggplot2::ggplot(pep_counts_df) +
        shape_pep + text + scale_pep +
        ggplot2::ggtitle(plot_title_pro) +
        ggplot2::xlab(xlabel_pro) +
        ggplot2::ylab(ylabel_pro) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size))
      if(bw_theme){
        p2 <- p2 + ggplot2::theme_bw()
      }
    }
    
    # conditionally return plots depending on mapping input
    if(mapping == "pep_to_pro") return(p1)
    else if(mapping == "pro_to_pep") return(p2)
    else if(mapping == "both") return(Rmisc::multiplot(p1, p2, cols = 2))
    
  })
}

#' plot.imdanovaFilt
#'
#' For plotting an S3 object of type 'imdanovaFilt':
#' 
#'@export
#'@rdname plot-pmartR-imdanovaFilt
#'@param min_nonmiss_gtest the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Suggested value is 3.
#'@param min_nonmiss_anova the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Suggested value is 2.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
plot.imdanovaFilt <- function(filter_object, min_nonmiss_anova = NULL, min_nonmiss_gtest = NULL, ...) {
  .plot.imdanovaFilt(filter_object, min_nonmiss_anova, min_nonmiss_gtest, ...)
}

.plot.imdanovaFilt <- function(filter_object, min_nonmiss_anova = NULL, min_nonmiss_gtest = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11) {
  
  ## initial checks ##
  
  # check that at least one of min_nonmiss_anova and min_nonmiss_gtest are present #
  if(is.null(min_nonmiss_anova) & is.null(min_nonmiss_gtest)) stop("At least one of min_nonmiss_anova and min_nonmiss_gtest must be present")
  # check that if they aren't NULL, min_nonmiss_anova and min_nonmiss_gtest are numeric, >=2 and >=3, respectively, and neither are bigger than the minimum group size (group_sizes in an attribute of the filter_object, see below) #
  if(!is.null(min_nonmiss_anova)) {
    # check that min_nonmiss_anova is numeric >= 2 #
    if(!inherits(min_nonmiss_anova, c("numeric","integer")) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is an integer #
    if(min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is less than the minimum group size #
    if(min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_anova cannot be greater than the minimum group size")
  }
  if(!is.null(min_nonmiss_gtest)) {
    # check that min_nonmiss_gtest is numeric >= 3 #
    if(!inherits(min_nonmiss_gtest, c("numeric","integer")) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is an integer #
    if(min_nonmiss_gtest %% 1 != 0) stop("min_nonmiss_gtest must be an integer >= 3")
    # check that min_nonmiss_gtest is less than the minimum group size #
    if(min_nonmiss_gtest > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_gtest cannot be greater than the minimum group size")
  }
  
  ## end of initial checks ##
  
  # workaround for certain "check" warnings
  x_anova <-x_gtest <- Var2 <- Var1 <- value <- NULL
  
  # The smallest group size
  min_n_group <- min(attr(filter_object,"group_sizes")$n_group)
  
  # ANOVA only
  if(is.null(min_nonmiss_gtest)) {
    # labels for plot #
    xlabel <- ifelse(is.null(x_lab), "Minimum Number of Samples Per Group", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Number of Biomolecules to Test", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Number of Biomolecules After ANOVA Filter", title_plot)
    
    # plot histogram
    ## get values
    log_df <- lapply(2:min_n_group, function(x) filter_object[,-1] >= x)
    filtered_vec1 <- colSums(sapply(2:min_n_group, function(x) rowSums(log_df[[x-1]]) >= 2))
    plot_df_anova <- data.frame(x_anova = 2:min_n_group, filtered_vec1)
    
    ## generate plot
    p <- ggplot2::ggplot(plot_df_anova) +
      ggplot2::geom_rect(ggplot2::aes(xmin = x_anova - 0.5, xmax = x_anova + 0.5,
                                      ymin = 0, ymax = filtered_vec1), fill="royalblue1", col="black") +
      ggplot2::geom_text(ggplot2::aes(x = x_anova, y = filtered_vec1, label = filtered_vec1), vjust = -0.25) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::ylim(0, 1.1*max(filtered_vec1, na.rm=TRUE)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
    
    if(!is.null(min_nonmiss_anova)) {
      num_tested <- plot_df_anova$filtered_vec1[min_nonmiss_anova-1]
      p <- p + ggplot2::annotate(geom = "rect", xmin = min_nonmiss_anova - 0.5, xmax = min_nonmiss_anova + 0.5,
                                 ymin = 0, ymax = num_tested, fill="royalblue1", col="black", size=1.5)
    }
    
    
    # gtest only
  }else if(is.null(min_nonmiss_anova)) {
    # plot labels #
    xlabel <- ifelse(is.null(x_lab), "Minimum Number of Samples Per Group", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Number of Biomolecules to Test", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Number of Biomolecules After IMD Filter", title_plot)
    
    # plot histogram
    
    ## get values
    max_vec <- sapply(1:nrow(filter_object), function(i) max(filter_object[i,-1]))
    filtered_vec2 <- sapply(3:min_n_group, function(i) sum(max_vec >= i))
    plot_df_gtest <- data.frame(x_gtest = 3:min_n_group, filtered_vec2)
    
    ## generate plot
    p <- ggplot2::ggplot(plot_df_gtest) +
      ggplot2::geom_rect(ggplot2::aes(xmin = x_gtest - 0.5, xmax = x_gtest + 0.5,
                                      ymin = 0, ymax = filtered_vec2), fill="royalblue1", col="black") +
      ggplot2::geom_text(ggplot2::aes(x = x_gtest, y = filtered_vec2, label = filtered_vec2), vjust = -0.25) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::ylim(0, 1.1*max(filtered_vec2, na.rm=TRUE)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
    
    if(!is.null(min_nonmiss_gtest)) {
      num_tested <- plot_df_gtest$filtered_vec2[min_nonmiss_gtest-2]
      p <- p + ggplot2::annotate(geom = "rect", xmin = min_nonmiss_gtest - 0.5, xmax = min_nonmiss_gtest + 0.5,
                                 ymin = 0, ymax = num_tested, fill="royalblue1", col="black", size=1.5)
    }
    
    
    # combined
  }else if(!is.null(min_nonmiss_anova) & !is.null(min_nonmiss_gtest)) {
    # plot labels #
    xlabel <- ifelse(is.null(x_lab), "Minimum Non-Missing Biomolecules (IMD)", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Minimum Non-Missing Biomolecules (ANOVA)", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Number of Biomolecules After Combined Filter", title_plot)
    
    # generate heatmap
    # ## get values for anova
    
    heat_df <- matrix(rep(0,length(2:min_n_group)*length(3:min_n_group)), nrow=length(2:min_n_group), ncol = length(3:min_n_group))
    for(i in 2:min_n_group) {#min_nonmiss_anova
      for(j in 3:min_n_group) {#min_nonmiss_gtest
        anova_log <- filter_object[,-1] >= i
        gtest_log <- filter_object[,-1] >= j
        heat_df[i-1,j-2] <- sum(rowSums(anova_log)>=2 & rowSums(gtest_log)>=1)
      }
    }
    ## plot heatmap
    colnames(heat_df) <- 3:min_n_group
    rownames(heat_df) <- 2:min_n_group
    heat_df_melt <- reshape2::melt(heat_df)
    bold_label <- heat_df_melt[heat_df_melt$Var1 == min_nonmiss_anova & heat_df_melt$Var2 == min_nonmiss_gtest, "value"]
    
    p <- ggplot2::ggplot(heat_df_melt, ggplot2::aes(x = ordered(Var2, levels = sort(unique(Var2))), y = ordered(Var1, levels = sort(unique(Var1))))) +
      ggplot2::geom_tile(ggplot2::aes(fill = value)) +
      ggplot2::geom_text(ggplot2::aes(label = value), size=50/max(heat_df_melt$Var2)) +
      ggplot2::geom_segment(x=min_nonmiss_gtest-2.5, xend=min_nonmiss_gtest-2.5, y=0, yend=max(heat_df_melt$Var1)+1, linetype = 1, size = 1.5) +
      ggplot2::geom_segment(x=min_nonmiss_gtest-1.5, xend=min_nonmiss_gtest-1.5, y=0, yend=max(heat_df_melt$Var1)+1, linetype = 1, size = 1.5) +
      ggplot2::geom_segment(x=0, xend=max(heat_df_melt$Var2)+1, y=min_nonmiss_anova-1.5, yend=min_nonmiss_anova-1.5, linetype = 1, size = 1.5) +
      ggplot2::geom_segment(x=0, xend=max(heat_df_melt$Var2)+1, y=min_nonmiss_anova-0.5, yend=min_nonmiss_anova-0.5, linetype = 1, size = 1.5) +
      ggplot2::scale_fill_gradientn(colours = heat.colors(100)) +
      ggplot2::annotate("text", label = bold_label, x = min_nonmiss_gtest-2, y = min_nonmiss_anova-1, size = 50/max(heat_df_melt$Var2), fontface = "bold", col="blue") +
      ggplot2::xlab(xlabel) +  ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }
  
  
  return(p)
}

#' plot.rmdFilt
#'
#' For plotting an S3 object of type 'rmdFilt':
#'
#'
#'
#'@export
#'@rdname plot-pmartR-rmdFilt
#'@param pvalue_threshold A threshold for the Robust Mahalanobis Distance (RMD) p-value. If \code{sampleID} is NULL (see \code{sampleID} below), a horizontal line is plotted at the RMD value that corresponds with the threshold, and all samples above the line have a p-value below the threshold. If \code{sampleID} is not NULL, \code{pvalue_threshold} will do nothing. Default value is NULL.
#'@param sampleID If specified, the plot produces a boxplot instead of a scatterplot. The \code{sampleID} input will place an "x" at the value for each of the metrics on the boxplots. Default value is NULL.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value indicating the font size for the plot title. Defaults to 14.
#' \code{x_lab_size} \tab integer value indicating the font size for the plot title. Defaults to 11.
#' \code{y_lab_size} \tab integer value indicating the font size for the plot title. Defaults to 11.
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used.
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Deafults to "right".
#'} 
plot.rmdFilt <- function(filter_object, pvalue_threshold = NULL, sampleID = NULL, ...) {
  .plot.rmdFilt(filter_object, pvalue_threshold, sampleID, ...)
}

.plot.rmdFilt <- function(filter_object, pvalue_threshold = NULL, sampleID = NULL, x_lab = NULL, y_lab = NULL, legend_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", point_size = 4) {
  
  ## initial checks ##
  
  # check pvalue_threshold if not NULL
  if(!is.null(pvalue_threshold)) {
    # check that length is 1
    if(length(pvalue_threshold) > 1) stop("pvalue_threshold must be of length 1")
    # check that pvalue_threshold is numeric [0,1] #
    if(!is.numeric(pvalue_threshold)) stop("pvalue_threshold must be numeric between 0 and 1")
    if(pvalue_threshold < 0 | pvalue_threshold > 1) stop("pvalue_threshold must be numeric between 0 and 1")
  }
  # check sampleID if not NULL
  if(!is.null(sampleID)) {
    # check that sampleID is a character
    if(!is.character(sampleID)) stop("sampleID must be a character string")
    # check that length is 1
    if(length(sampleID) > 1) stop("sampleID must be of length 1")
  }
  #check point_size argument is numeric and >= to zero
  if(!is.numeric(point_size)) stop("point_size must be numeric")
  if(point_size < 0) stop("point_size must be greater than or equal to zero")
  
  ## end of initial checks ##
  
  samp_id <- names(attr(filter_object, "group_DF"))[1]
  metrics <- attributes(filter_object)$metrics
  
  # set text size parameters #
  # if(is.null(title_size)){
  #   title_size <- 14
  # }
  # if(is.null(x_lab_size)){
  #   x_lab_size <- 11
  # }
  # if(is.null(y_lab_size)){
  #   y_lab_size1 <- 11
  # }
  
  # determine how to melt based on the number of main effects
  group_df <- attributes(filter_object)$group_DF
  
  # determine the main effects, then melt the df #
  if(ncol(group_df) == 2) {
    main_eff_names <- "Group"
    dfmelt <- reshape2::melt(filter_object, id = c(samp_id, main_eff_names))
  } else if(ncol(group_df) > 2) {
    ## put main effect with more levels first, for plotting aesthetics ##
    temp <- droplevels(group_df[,-(1:2)])
    numlevels <- sapply(1:2, function(j) length(levels(as.factor(temp[,j]))))
    main_eff_names <- names(temp)[order(numlevels, decreasing = TRUE)]
    dfmelt <- reshape2::melt(filter_object[,-2], id = c(samp_id, main_eff_names))
  }
  
  dfsub <- dfmelt[dfmelt$variable %in% metrics,]
  
  
  # legend labels #
  ## function for shortening long main effect names ##
  abbrev_fun <- function(string) {
    if(nchar(string)>25) string=paste0(substr(string,1,23),"...")
    return(string)
  }
  display_names <- sapply(main_eff_names, abbrev_fun)
  
  legend_title_color <- ifelse(is.null(legend_lab), display_names[1], legend_lab[1])
  legend_title_shape <- NULL
  if(length(display_names) > 1) {
    if(length(legend_lab) > 1) {
      legend_title_shape <- legend_lab[2]
    } else {
      legend_title_shape <- display_names[2]
    }
  }
  
  
  # Make plot #
  if(!is.null(sampleID)) {
    levels(dfsub$variable)[levels(dfsub$variable)=="Fraction_Missing"] <- "Prop_missing"
    plot_title <- ifelse(is.null(title_plot), paste0("Summary of Sample ", sampleID, " and Metrics Used"), title_plot)
    xlabel <- ifelse(is.null(x_lab), " ", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
    
    if(bw_theme == FALSE){
      p <- ggplot2::ggplot(dfsub) +
        ggplot2::geom_boxplot(ggplot2::aes(x=rep(1,length(value)), y=value), fill=heat.colors(length(metrics))) +
        ggplot2::facet_wrap(~ variable, scales = "free", ncol=length(metrics)) +
        ggplot2::geom_point(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), size=point_size, pch=4) +
        ggplot2::geom_text(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), label = sampleID, vjust=1.5, size=3.5, fontface="bold") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title)
    }else{
      p <- ggplot2::ggplot(dfsub) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x=rep(1,length(value)), y=value), fill=heat.colors(length(metrics))) +
        ggplot2::facet_wrap(~ variable, scales = "free", ncol=length(metrics)) +
        ggplot2::geom_point(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), size=point_size, pch=4) +
        ggplot2::geom_text(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), label = sampleID, vjust=1.5, size=3.5, fontface="bold") +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title)
    }
    
  }else if(is.null(pvalue_threshold)) {
    if(length(main_eff_names)==1) {
      p <- ggplot2::ggplot(filter_object) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names[1]), size=point_size)
    } else {
      p <- ggplot2::ggplot(filter_object) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names[1], pch=main_eff_names[2]), size=point_size)
    }
    
    plot_title <- ifelse(is.null(title_plot), "Sample Outlier Results", title_plot)
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)
    
    if(bw_theme == FALSE){
      p <- p + ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                              plot.title = ggplot2::element_text(size=title_size),
                              axis.title.x = ggplot2::element_text(size=x_lab_size),
                              axis.title.y = ggplot2::element_text(size=y_lab_size),
                              legend.position = legend_position) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title) +
        ggplot2::scale_color_discrete(legend_title_color) +
        ggplot2::scale_shape_discrete(legend_title_shape)
    }else{
      p <- p + ggplot2::theme_bw() +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title) +
        ggplot2::scale_color_discrete(legend_title_color) +
        ggplot2::scale_shape_discrete(legend_title_shape)
    }
    
    
  }else {
    # get y-intercept for line
    df <- attributes(filter_object)$df
    yint <- log(qchisq(1 - pvalue_threshold, df = df), base = 2)
    
    # divide data by the threshold
    sub1 <- subset(filter_object, pvalue < pvalue_threshold)
    sub2 <- subset(filter_object, pvalue >= pvalue_threshold)
    
    # make title
    if(is.null(title_plot)) {
      plot_title <- bquote(atop("Sample Outlier Results", italic(paste("p-value threshold = ", .(pvalue_threshold)))))
    } else {
      plot_title <- title_plot
    }
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)
    
    if(length(main_eff_names)==1) {
      p <- ggplot2::ggplot(sub1) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names), size=point_size, bg="gray") +
        ggplot2::geom_point(data = sub2, ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names), alpha=0.5, size=point_size, bg="gray")
    } else {
      p <- ggplot2::ggplot(sub1) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", pch=main_eff_names[2], col=main_eff_names[1]), size=point_size, bg="gray") +
        ggplot2::geom_point(data = sub2, ggplot2::aes_string(x = samp_id, y="Log2.md", pch=main_eff_names[2], col=main_eff_names[1]), alpha=0.5, size=point_size, bg="gray")
    }
    
    if(bw_theme == FALSE){
      p <- p +
        ggplot2::geom_hline(yintercept = yint) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position) +
        ggplot2::guides(col = ggplot2::guide_legend(ncol = 1), pch = ggplot2::guide_legend(ncol = 1)) +
        ggplot2::scale_shape_manual(legend_title_shape, values=rep(c(16,9,22,13,11,3,4,15,0,5),6)) +
        ggplot2::scale_color_manual(legend_title_color, values=rep(c("blue","red","green","darkturquoise","goldenrod3","darkorchid2"),10)) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title)
    }else{
      p <- p +
        ggplot2::geom_hline(yintercept = yint) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position) +
        ggplot2::guides(col = ggplot2::guide_legend(ncol = 1), pch = ggplot2::guide_legend(ncol = 1)) +
        ggplot2::scale_shape_manual(legend_title_shape, values=rep(c(16,9,22,13,11,3,4,15,0,5),6)) +
        ggplot2::scale_color_manual(legend_title_color, values=rep(c("blue","red","green","darkturquoise","goldenrod3","darkorchid2"),10)) +
        ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + ggplot2::ggtitle(plot_title)
    }
    
  }
  
  return(p)
}

#' plot.cvFilt
#' 
#' For plotting an S3 object of type 'cvFilt':
#' 
#'@export
#'@rdname plot-pmartR-cvFilt
#'@param cv_threshold draws a vertical line at the pooled cv cutoff, values to the left of this cutoff are dropped if the filter is applied.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{log_scale} \tab logical indicating whether to use a log2 transformed x-axis. Defaults to TRUE.\cr
#' \code{n_bins} \tab integer value specifying the number of bins to draw in the histogram.  Defaults to 30. \cr
#' \bode{n_breaks} \tab integer value specifying the number of breaks to use.  You may get less breaks if rounding causes certain values to become non-unique.  Defaults to 15. \cr
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
#' 
#' @examples
#' data(pep_object)
#' pep_object <- group_designation(omicsData = pep_object, main_effects = "Condition")
#' cvfilt <- cv_filter(pep_object)
#' 
#' plot(cvfilt, cv_threshold = 20)
#' plot(cvfilt, cv_threshold = 10, log_scale = FALSE)
#' 
plot.cvFilt <- function(filter_object, cv_threshold = NULL, ...) {
  .plot.cvFilt(filter_object, cv_threshold, ...)
}

.plot.cvFilt <- function(filter_object, cv_threshold = NULL, log_scale = TRUE, n_breaks = 15, n_bins = 30, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {
  
  # checks for cv_threshold if not null
  if(!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if(!is.numeric(cv_threshold)) stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if(length(cv_threshold)>1) stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if(cv_threshold <= 1 | cv_threshold >= max(filter_object$CV_pooled, na.rm = TRUE)) stop("cv_threshold must be greater than 1 and less than the maximum CV_pooled value")
  }
  
  if(!is.logical(log_scale)) stop("log_scale must be logical: TRUE or FALSE")
  if(!(n_breaks%%1 == 0)) stop("n_breaks must be integer valued")
  if(!(n_bins%%1 == 0)) stop("n_bins must be integer valued")
  
  # plotting object
  new_object <- filter_object[!is.na(filter_object$CV_pooled),]
  max_x_val <- attributes(filter_object)$max_x_val
  
  # labels
  plot_title <- ifelse(is.null(title_plot), "Coefficient of Variation (CV)", title_plot)
  xlabel <- ifelse(is.null(x_lab), "Pooled CV", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Count", y_lab)
  
  ### Store ggplot2 layers
  
  # bw theme
  bw = if(bw_theme) ggplot2::theme_bw() else NULL
  
  # scale transform and breaks depending on x-axis scale
  if(log_scale){
    trans <- "log2"
    i <- log2(min(new_object$CV_pooled, na.rm = TRUE))
    breaks <- 0
    
    # define a step value that is evenly spaced in the log2 scale
    step = (max(log2(new_object$CV_pooled), na.rm = TRUE) - min(log2(new_object$CV_pooled), na.rm = TRUE))/n_breaks
    
    # create the normal scale labels that will be log2 transformed when passed to ggplot
    while(2^i < max(new_object$CV_pooled, na.rm = TRUE)){
      breaks <- c(breaks, 2^i)
      i <- i+step
    }
    # rounding for plot purposes
    breaks <- round(breaks, 2)
  }
  else{
    breaks <- scales::pretty_breaks(n = n_breaks)
    trans <- "identity" 
  }
 
  # change default title and draw a vertical line if cv_thresh specified
  if(!is.null(cv_threshold)) {
    if(is.null(title_plot)) {
      plot_title <- bquote(paste("Coefficient of Variation (CV):  ",italic(paste("CV Threshold = ",.(cv_threshold))),""))
    }
    cutoff <- ggplot2::geom_vline(xintercept = cv_threshold)
  }else{
    cutoff <- NULL
  }
  
  # main plot object
  p <- ggplot2::ggplot(new_object) +
    ggplot2::geom_histogram(ggplot2::aes(x=CV_pooled), bins = n_bins, fill = "steelblue1") + 
    cutoff + bw +
    ggplot2::scale_x_continuous(breaks = breaks, trans = trans) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) + 
    ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                   axis.title.x = ggplot2::element_text(size=x_lab_size),
                   axis.title.y = ggplot2::element_text(size=y_lab_size))
  
  return(p)
}

#' plot.customFilt
#'
#' For plotting an s3 object of type 'customFilt'
#' 
#'@export
#'@rdname plot-pmartR-customFilt
plot.customFilt <- function(filter_object, ...) {
  .plot.customFilt(filter_object, ...)
}

.plot.customFilt <- function(filter_object, x_lab = NULL, y_lab = NULL, title_plot = NULL) {
  
  warning("There is no plot method for objects of class 'customFilt'. See summary.customFilt instead.")
}

#' plot.pepData
#' 
#' For plotting an S3 object of type 'pepData'
#' 
#'@export
#'@rdname plot-pmartR-pepData
#'@param order_by a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If \code{order_by} is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data.
#'@param color_by a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If \code{color_by} is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color.
#'@param facet_by a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL.
#'@param facet_cols an optional integer specifying the number of columns to show in the facet plot.
#'
#'@param ... Additional arguments
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' \code{ylimit} \tab numeric vector of length 2 specifying y axis lower and upper limits. \cr
#' \code{use_VizSampNames} \tab logical specifies whether to use custom sample names \cr
#' }
plot.pepData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.pepData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.pepData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", ylimit = NULL, use_VizSampNames = FALSE) {
  
  ## initial checks ##
  if(!is.null(order_by)) {
    if(!is.character(order_by) | length(order_by) > 1) stop("order_by must be a character vector of length 1")
  }
  if(!is.null(color_by)) {
    if(!is.character(color_by) | length(color_by) > 1) stop("color_by must be a character vector of length 1")
  }
  if(!is.null(facet_by)) {
    if(!is.character(facet_by) | length(facet_by) > 1) stop("facet_by must be a character vector of length 1")
  }
  if(!is.null(facet_cols)) {
    if(is.null(facet_by)) stop("facet_by cannot be NULL when facet_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  
  ##checking that 'size' arguments are numeric
  if(!is.numeric(title_size) | !is.numeric(x_lab_size) | !is.numeric(y_lab_size)) stop("title_size, x_lab_size and y_lab_size must be integer values")
  
  ##checking that ylimit is numeric of length 2
  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) | length(ylimit)!= 2) stop("ylimit must be a numeric vector of length 2")
  }
  
  
  # add check for samples with all NAs and return message to user that these will not be plotted #
  sample_nas <- colSums(is.na(omicsData$e_data))
  if(any(sample_nas == nrow(omicsData$e_data))){
    empties <- names(omicsData$e_data)[which(sample_nas == nrow(omicsData$e_data))]
    message(paste("The following sample(s) are comprised entirely of missing data and will not be included in the plot: ", empties, sep = " "))
  }
  
  
  ## end of initial checks ##
  
  
  # organize data #
  e_data <- omicsData$e_data
  f_data = omicsData$f_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)
  
  if(inherits(omicsData, "isobaricpepData")){
    maintitle <- ifelse(attributes(omicsData)$data_info$isobaric_norm,
                        "Boxplots of Normalized Isobaric Peptide Data",
                        "Boxplots of Un-Normalized Isobaric Peptide Data")
  }
  else if(inherits(omicsData, "pepData")){
    maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                        "Boxplots of Normalized Peptide Data",
                        "Boxplots of Un-Normalized Peptide Data")
  }


  
  # get data and aesthetics for plots #
  
  ## if facet_by is not null and isn't the same as either order_by or color_by ##
  if(!is.null(facet_by)) {
    if(!(facet_by %in% c(order_by, color_by))) {
      facet_temp <- group_designation(omicsData, main_effects = facet_by)
      facetDF <- attributes(facet_temp)$group_DF
      colnames(facetDF) <- c("variable", facet_by)
      
      plot_data <- merge(plot_data, facetDF, by = "variable")
    }
  }
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data <- merge(plot_data, orderDF, by = "variable")
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    #plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- ggplot2::theme_bw()
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data <- merge(plot_data, colorDF, by = "variable")
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,color_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data <- merge(plot_data, tempdata, by = "variable")
    } else {
      plot_data <- merge(plot_data, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
  }
  
  
  # facet plot #
  if(!is.null(facet_by)) {
    if(is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
    }
  }
  
  # custom labels #
  if(!is.null(title_plot)) title <- title_plot
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
 
  if(is.null(y_lab)){
    if(attr(omicsData, "data_info")$data_scale == 'abundance'){
      ylabel<- "Abundance"}
    else if(attr(omicsData, "data_info")$data_scale == 'log'){
      ylabel<- "ln Abundance"
    }else ylabel <- paste(attr(omicsData, "data_info")$data_scale, "Abundance", sep = " ")
  }else ylabel = y_lab
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)
  
  if(!is.null(ylimit))
  {
    p <- p + ggplot2::scale_y_continuous(limits = ylimit)
  }
  
  
  return(p)
}

#' plot.proData
#' 
#' For plotting an S3 object of type 'proData'
#' 
#' \tabular{ll}{
#' \code{order_by} \tab a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If order_by is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data. \cr
#' \code{color_by} \tab a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If color_by is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color. \cr
#' \code{facet_by} \tab a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL. \cr
#' \code{facet_cols} \tab an optional integer specifying the number of columns to show in the facet plot. \cr
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' \code{ylimit} \tab numeric vector of length 2 specifying y axis lower and upper limits. \cr
#' \code{use_VizSampNames} \tab logical specifies whether to use custom sample names \cr
#' }
#'
#'@export
#'@rdname plot-pmartR-proData
plot.proData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.proData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.proData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", ylimit = NULL, use_VizSampNames = FALSE) {
  
  ## initial checks ##
  if(!is.null(order_by)) {
    if(!is.character(order_by) | length(order_by) > 1) stop("order_by must be a character vector of length 1")
  }
  if(!is.null(color_by)) {
    if(!is.character(color_by) | length(color_by) > 1) stop("color_by must be a character vector of length 1")
  }
  if(!is.null(facet_by)) {
    if(!is.character(facet_by) | length(facet_by) > 1) stop("facet_by must be a character vector of length 1")
  }
  if(!is.null(facet_cols)) {
    if(is.null(facet_by)) stop("facet_by cannot be NULL when facet_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  ##checking that 'size' arguments are numeric
  if(!is.numeric(title_size) | !is.numeric(x_lab_size) | !is.numeric(y_lab_size)) stop("title_size, x_lab_size and y_lab_size must be integer values")
  
  ##checking that ylimit is numeric of length 2
  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) | length(ylimit)!= 2) stop("ylimit must be a numeric vector of length 2")
  }  
  # add check for samples with all NAs and return message to user that these will not be plotted #
  sample_nas <- colSums(is.na(omicsData$e_data))
  if(any(sample_nas == nrow(omicsData$e_data))){
    empties <- names(omicsData$e_data)[which(sample_nas == nrow(omicsData$e_data))]
    message(paste("The following sample(s) are comprised entirely of missing data and will not be included in the plot: ", empties, sep = " "))
  }

  ## end of initial checks ##
  
  
  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  f_data = omicsData$f_data
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)
  
  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Protein Data",
                      "Boxplots of Un-Normalized Protein Data")
  
  # get data and aesthetics for plots #
  
  ## if facet_by is not null and isn't the same as either order_by or color_by ##
  if(!is.null(facet_by)) {
    # if(!use_VizSampNames) stop("if argument 'facet_by' is provided, argument 'use_VizSampNames' must be set to FALSE")
    if(!(facet_by %in% c(order_by, color_by))) {
      facet_temp <- group_designation(omicsData, main_effects = facet_by)
      facetDF <- attributes(facet_temp)$group_DF
      colnames(facetDF) <- c("variable", facet_by)
      
      plot_data <- merge(plot_data, facetDF, by = "variable")
    }
  }
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme == TRUE){
      p  <- p + ggplot2::theme_bw()
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data <- merge(plot_data, orderDF, by = "variable")
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    #plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data <- merge(plot_data, colorDF, by = "variable")
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,color_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data <- merge(plot_data, tempdata, by = "variable")
    } else {
      plot_data <- merge(plot_data, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
  }
  
  
  # facet plot #
  if(!is.null(facet_by)) {
    if(is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
    }
  }
  
  # custom labels #
  if(!is.null(title_plot)) title <- title_plot
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
  
  if(is.null(y_lab)){
    if(attr(omicsData, "data_info")$data_scale == 'abundance'){
      ylabel<- "Abundance"}
    else if(attr(omicsData, "data_info")$data_scale == 'log'){
      ylabel<- "ln Abundance"
    }else ylabel <- paste(attr(omicsData, "data_info")$data_scale, "Abundance", sep = " ")
  }else{ylabel = y_lab}
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)
  
  if(!is.null(ylimit))
  {
    p <- p + ggplot2::scale_y_continuous(limits = ylimit)
  }
  
  
  return(p)
}

#' plot.lipidData
#' 
#' For plotting an S3 object of type 'lipidData'
#' 
#' \tabular{ll}{
#' \code{order_by} \tab a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If order_by is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data. \cr
#' \code{color_by} \tab a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If color_by is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color. \cr
#' \code{facet_by} \tab a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL. \cr
#' \code{facet_cols} \tab an optional integer specifying the number of columns to show in the facet plot. \cr
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' \code{ylimit} \tab numeric vector of length 2 specifying y axis lower and upper limits. \cr
#' \code{use_VizSampNames} \tab logical specifies whether to use custom sample names \cr
#' }
#'@export
#'@rdname plot-pmartR-lipidData
plot.lipidData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.lipidData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.lipidData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", ylimit = NULL, use_VizSampNames = FALSE) {
  
  ## initial checks ##
  if(!is.null(order_by)) {
    if(!is.character(order_by) | length(order_by) > 1) stop("order_by must be a character vector of length 1")
  }
  if(!is.null(color_by)) {
    if(!is.character(color_by) | length(color_by) > 1) stop("color_by must be a character vector of length 1")
  }
  if(!is.null(facet_by)) {
    if(!is.character(facet_by) | length(facet_by) > 1) stop("facet_by must be a character vector of length 1")
  }
  if(!is.null(facet_cols)) {
    if(is.null(facet_by)) stop("facet_by cannot be NULL when facet_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  
  ##checking that 'size' arguments are numeric
  if(!is.numeric(title_size) | !is.numeric(x_lab_size) | !is.numeric(y_lab_size)) stop("title_size, x_lab_size and y_lab_size must be integer values")
  
  ##checking that ylimit is numeric of length 2
  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) | length(ylimit)!= 2) stop("ylimit must be a numeric vector of length 2")
  }  
  # add check for samples with all NAs and return message to user that these will not be plotted #
  sample_nas <- colSums(is.na(omicsData$e_data))
  if(any(sample_nas == nrow(omicsData$e_data))){
    empties <- names(omicsData$e_data)[which(sample_nas == nrow(omicsData$e_data))]
    message(paste("The following sample(s) are comprised entirely of missing data and will not be included in the plot: ", empties, sep = " "))
  }
  ## end of initial checks ##
  
  
  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  f_data = omicsData$f_data
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)
  
  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Lipid Data",
                      "Boxplots of Un-Normalized Lipid Data")
  
  # get data and aesthetics for plots #
  
  ## if facet_by is not null and isn't the same as either order_by or color_by ##
  if(!is.null(facet_by)) {
    # if(!use_VizSampNames) stop("if argument 'facet_by' is provided, argument 'use_VizSampNames' must be set to FALSE")
    if(!(facet_by %in% c(order_by, color_by))) {
      facet_temp <- group_designation(omicsData, main_effects = facet_by)
      facetDF <- attributes(facet_temp)$group_DF
      colnames(facetDF) <- c("variable", facet_by)
      
      plot_data <- merge(plot_data, facetDF, by = "variable")
    }
  }
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data <- merge(plot_data, orderDF, by = "variable")
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    #plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data <- merge(plot_data, colorDF, by = "variable")
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot::theme_bw()
      
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,color_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data <- merge(plot_data, tempdata, by = "variable")
    } else {
      plot_data <- merge(plot_data, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
  }
  
  
  # facet plot #
  if(!is.null(facet_by)) {
    if(is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
    }
  }
  
  # custom labels #
  if(!is.null(title_plot)) title <- title_plot
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
  
  if(is.null(y_lab)){
    if(attr(omicsData, "data_info")$data_scale == 'abundance'){
      ylabel<- "Abundance"}
    else if(attr(omicsData, "data_info")$data_scale == 'log'){
      ylabel<- "ln Abundance"
    }else ylabel <- paste(attr(omicsData, "data_info")$data_scale, "Abundance", sep = " ")
  }else{ylabel = y_lab}
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)
  
  if(!is.null(ylimit))
  {
    p <- p + ggplot2::scale_y_continuous(limits = ylimit)
  }
  
  
  return(p)
}

#' plot.metabData
#' 
#' For plotting an S3 object of type 'metabData'
#' 
#' \tabular{ll}{
#' \code{order_by} \tab a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If order_by is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data. \cr
#' \code{color_by} \tab a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If color_by is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color. \cr
#' \code{facet_by} \tab a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL. \cr
#' \code{facet_cols} \tab an optional integer specifying the number of columns to show in the facet plot. \cr
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' \code{ylimit} \tab numeric vector of length 2 specifying y axis lower and upper limits. \cr
#' \code{use_VizSampNames} \tab logical specifies whether to use custom sample names \cr
#' }
#'@export
#'@rdname plot-pmartR-metabData
plot.metabData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.metabData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.metabData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", ylimit = NULL, use_VizSampNames = FALSE) {
  
  ## initial checks ##
  if(!is.null(order_by)) {
    if(!is.character(order_by) | length(order_by) > 1) stop("order_by must be a character vector of length 1")
  }
  if(!is.null(color_by)) {
    if(!is.character(color_by) | length(color_by) > 1) stop("color_by must be a character vector of length 1")
  }
  if(!is.null(facet_by)) {
    if(!is.character(facet_by) | length(facet_by) > 1) stop("facet_by must be a character vector of length 1")
  }
  if(!is.null(facet_cols)) {
    if(is.null(facet_by)) stop("facet_by cannot be NULL when facet_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  
  ##checking that 'size' arguments are numeric
  if(!is.numeric(title_size) | !is.numeric(x_lab_size) | !is.numeric(y_lab_size)) stop("title_size, x_lab_size and y_lab_size must be integer values")
  
  ##checking that ylimit is numeric of length 2
  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) | length(ylimit)!= 2) stop("ylimit must be a numeric vector of length 2")
  }
  
  # add check for samples with all NAs and return message to user that these will not be plotted #
  sample_nas <- colSums(is.na(omicsData$e_data))
  if(any(sample_nas == nrow(omicsData$e_data))){
    empties <- names(omicsData$e_data)[which(sample_nas == nrow(omicsData$e_data))]
    message(paste("The following sample(s) are comprised entirely of missing data and will not be included in the plot: ", empties, sep = " "))
  }
  ## end of initial checks ##
  
  
  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  f_data = omicsData$f_data
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)
  
  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Metabolite Data",
                      "Boxplots of Un-Normalized Metabolite Data")
  
  # get data and aesthetics for plots #
  
  ## if facet_by is not null and isn't the same as either order_by or color_by ##
  if(!is.null(facet_by)) {
    if(!(facet_by %in% c(order_by, color_by))) {
      # if(!use_VizSampNames) stop("if argument 'facet_by' is provided, argument 'use_VizSampNames' must be set to FALSE")
      facet_temp <- group_designation(omicsData, main_effects = facet_by)
      facetDF <- attributes(facet_temp)$group_DF
      colnames(facetDF) <- c("variable", facet_by)
      
      plot_data <- merge(plot_data, facetDF, by = "variable")
    }
  }
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data <- merge(plot_data, orderDF, by = "variable")
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    #plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data <- merge(plot_data, colorDF, by = "variable")
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- maintitle
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,color_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data <- merge(plot_data, tempdata, by = "variable")
    } else {
      plot_data <- merge(plot_data, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data <- plot_data[order(plot_data[,order_by]),]
    plot_data$variable <- factor(plot_data$variable, levels=unique(plot_data$variable), ordered=TRUE)
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    
    if(bw_theme==TRUE){
      p <- p + ggplot2::theme_bw()
      
    }
    
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    if(use_VizSampNames == T){
      f_data = f_data[order(f_data[,order_by]),]
      p = p + ggplot2::scale_x_discrete(labels = f_data$VizSampNames, breaks = unique(f_data[,get_fdata_cname(omicsData)]))
    }
  }
  
  
  # facet plot #
  if(!is.null(facet_by)) {
    if(is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
    }
  }
  
  # custom labels #
  if(!is.null(title_plot)) title <- title_plot
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
  
  if(is.null(y_lab)){
    if(attr(omicsData, "data_info")$data_scale == 'abundance'){
      ylabel<- "Abundance"}
    else if(attr(omicsData, "data_info")$data_scale == 'log'){
      ylabel<- "ln Abundance"
    }else ylabel <- paste(attr(omicsData, "data_info")$data_scale, "Abundance", sep = " ")
  }else{ylabel = y_lab}
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                   axis.title.x = ggplot2::element_text(size=x_lab_size),
                   axis.title.y = ggplot2::element_text(size=y_lab_size),
                   legend.position = legend_position)
  
  if(!is.null(ylimit))
  {
    p <- p + ggplot2::scale_y_continuous(limits = ylimit)
  }
  
  
  return(p)
}

#' plot.dimRes
#' 
#' For plotting an S3 object of type 'dimRes':
#' 
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case a default label is used. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case a default title is used. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' }
#'@export
#'@rdname plot-pmartR-dimRes
plot.dimRes <- function(dimRes_object, ...) {
  .plot.dimRes(dimRes_object, ...)
}

.plot.dimRes <- function(dimRes_object, x_lab = NULL, y_lab = NULL, legend_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE, legend_position = "right", point_size = 4) {
  
  plotdata <- data.frame(SampleID = dimRes_object$SampleID, PC1 = dimRes_object$PC1, PC2 = dimRes_object$PC2)
  plotdata_name<-names(plotdata)[1]
  
  #check point_size argument is numeric and >= zero
  if(!is.numeric(point_size)) stop("point_size must be numeric")
  if(point_size < 0) stop("point_size must be greater than or equal to zero")
  
  # if there is a group designation #
  if(!is.null(attr(dimRes_object,"group_DF"))) {
    group_DF <- attr(dimRes_object,"group_DF")
    fdata_cname<- names(group_DF)[1]
    
    # if there's two main effects #
    if(ncol(group_DF) == 4) {
      main_eff_names <- names(group_DF[,3:4])
      # function for shortening long main effect names #
      abbrev_fun <- function(string) {
        if(nchar(string)>25) string=paste0(substr(string,1,23),"...")
        return(string)
      }

      # manage the length of legend titles #
      display_names <- sapply(main_eff_names, abbrev_fun)
      
      if(!identical(as.character(plotdata[[plotdata_name]]), as.character(group_DF[[fdata_cname]]))){
        group_DF<- group_DF[match(plotdata[[plotdata_name]], group_DF[[fdata_cname]]), ]
        plotdata<- merge.data.frame(plotdata, group_DF, by.x = plotdata_name, by.y = fdata_cname, sort = FALSE)
      }
      
      else plotdata<- merge.data.frame(plotdata, group_DF, by.x = plotdata_name, by.y = fdata_cname, sort = FALSE)
      
      color_var <- main_eff_names[1]
      pch_var <- main_eff_names[2]
      
      #plotdata[, colnames(plotdata) == pch_var] <- as.factor(plotdata[, colnames(plotdata) == pch_var])
      #plotdata[, colnames(plotdata) == color_var] <- as.factor(plotdata[, colnames(plotdata) == color_var])
      
      if(length(legend_lab) > length(main_eff_names)) warning("legend_lab length is greater than the number of main effects. Only the first two entries will be used.")
      
    } else {
      
      if(!identical(as.character(plotdata[[plotdata_name]]), as.character(group_DF[[fdata_cname]]))){
        group_DF<- group_DF[match(plotdata[[plotdata_name]], group_DF[[fdata_cname]]), ]
        plotdata<- merge.data.frame(plotdata, group_DF, by.x = plotdata_name, by.y = fdata_cname, sort = FALSE)
      }
      
      else plotdata<- merge.data.frame(plotdata, group_DF, by.x = plotdata_name, by.y = fdata_cname,sort = FALSE)
      
      color_var <- "Group"
      pch_var <- NULL
      display_names <- c("Group", NULL)
      
      if(length(legend_lab) > 1) warning("legend_lab length is greater than the number of main effects. Only the first entry will be used.")
    }
  } else {
    color_var <- NULL
    pch_var <- NULL
    display_names <- c(NULL, NULL)
    if(!is.null(legend_lab)) warning("There is no group designation, so legend_lab will go unused.")
  }
  
  # axis labels #
  xr2 <- paste(" = ", round(attr(dimRes_object, "R2")[1],3), ")", sep = "")
  yr2 <- paste(" = ", round(attr(dimRes_object, "R2")[2],3), ")", sep = "")
  pc1 <- "PC1 ("
  pc2 <- "PC2 ("
  
  # custom legend names #
  if(!is.null(legend_lab)) {
    # make the vector at least length 2 to avoid errors in the plot
    display_names[1:length(legend_lab)] <- legend_lab[1:min(2, length(legend_lab))]
  }
  
  # title #
  plot_title <- ifelse(is.null(title_plot), "Principal Components", title_plot)
  
  # plot #
  if(bw_theme==FALSE){
    p <- ggplot2::ggplot(plotdata, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(ggplot2::aes_string(col = color_var, pch = pch_var), size = point_size) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_color_discrete(display_names[1]) +
      ggplot2::scale_shape_discrete(display_names[2]) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
  }else{
    p <- ggplot2::ggplot(plotdata, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::theme_bw() +
      ggplot2::geom_point(ggplot2::aes_string(col = color_var, pch = pch_var), size = point_size) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_color_discrete(display_names[1]) +
      ggplot2::scale_shape_discrete(display_names[2]) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
  }
  
  
  if(is.null(x_lab) & is.null(y_lab)) {
    p <- p + ggplot2::xlab(substitute(paste(pc1, R^2, xr2))) +
      ggplot2::ylab(substitute(paste(pc2, R^2, yr2)))
  } else {
    if(!is.null(x_lab)) {
      p <- p + ggplot2::xlab(x_lab)
    }
    if(!is.null(y_lab)) {
      p <- p + ggplot2::ylab(y_lab)
    }
  }
  
  
  return(p)
}

#' plot.normRes
#'
#' For plotting an S3 object of type 'normRes'
#'
#'@export
#'@rdname plot-pmartR-normRes
#'
#'
plot.normRes <- function(normData, order_by = NULL, color_by = NULL, ...) {
  .plot.normRes(normData, order_by, color_by, ...)
}

.plot.normRes <- function(normData, order_by = NULL, color_by = NULL, x_lab = NULL, y_lab = NULL, title_plot_raw = NULL, title_plot_norm = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right", combined=FALSE) {
  
  ## initial checks ##
  if(!is.null(order_by)) {
    if(!is.character(order_by) | length(order_by) > 1) stop("order_by must be a character vector of length 1")
  }
  if(!is.null(color_by)) {
    if(!is.character(color_by) | length(color_by) > 1) stop("color_by must be a character vector of length 1")
  }
  # if(!is.null(facet_by)) {
  #   if(!is.character(facet_by) | length(facet_by) > 1) stop("facet_by must be a character vector of length 1")
  # }
  # if(!is.null(facet_cols)) {
  #   if(is.null(facet_by)) stop("facet_by cannot be NULL when facet_cols is specified")
  #   if(length(facet_cols)>1) stop("facet_cols must be of length 1")
  #   if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
  #   if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  # }
  ## end of initial checks ##
  
  
  # organize data #
  omicsData <- attributes(normData)$omicsData
  e_data_raw <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data_raw <- reshape2::melt(e_data_raw, id = e_data_cname, na.rm = TRUE)
  
  omicsDataNorm <- omicsData
  if(is.null(normData$parameters$normalization$scale)){
    # location param only
    if(is.null(normData$parameters$backtransform$location)){
      # no backtransform
      omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] <- omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] - matrix(rep(normData$parameters$normalization$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw))
    }else{
      # yes backtransform
      omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] <- omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] - matrix(rep(normData$parameters$normalization$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw)) + matrix(rep(normData$parameters$backtransform$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw))
    }
    
  }else{
    # location plus scale param
    if(is.null(normData$parameters$backtransform$location)){
      # no backtransform
      omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] <- (omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] - matrix(rep(normData$parameters$normalization$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw)))/matrix(rep(normData$parameters$normalization$scale, each=nrow(e_data_raw)), nrow=nrow(e_data_raw))
    }else{
      # yes backtransform
      omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] <- ( (omicsDataNorm$e_data[,-which(names(omicsDataNorm$e_data) == e_data_cname)] - matrix(rep(normData$parameters$normalization$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw)))/matrix(rep(normData$parameters$normalization$scale, each=nrow(e_data_raw)), nrow=nrow(e_data_raw)) ) * normData$parameters$backtransform$scale + matrix(rep(normData$parameters$backtransform$location, each=nrow(e_data_raw)), nrow=nrow(e_data_raw))
      
    }
  }
  e_data_norm <- omicsDataNorm$e_data
  plot_data_norm <- reshape2::melt(e_data_norm, id = e_data_cname, na.rm = TRUE)
  
  
  maintitle_raw <- "Boxplots of Un-Normalized Peptide Data"
  maintitle_norm <- "Boxplots of Normalized Peptide Data"
  
  # get data and aesthetics for plots #
  
  # ## if facet_by is not null and isn't the same as either order_by or color_by ##
  # if(!is.null(facet_by)) {
  #   if(!(facet_by %in% c(order_by, color_by))) {
  #     facet_temp <- group_designation(normData, main_effects = facet_by)
  #     facetDF <- attributes(facet_temp)$group_DF
  #     colnames(facetDF) <- c("variable", facet_by)
  #
  #     plot_data <- merge(plot_data, facetDF, by = "variable")
  #   }
  # }
  
  #### RAW DATA PLOT ####
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data_raw) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data_raw) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- maintitle_raw
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data_raw <- merge(plot_data_raw, orderDF, by = "variable")
    
    # reorder levels #
    plot_data_raw <- plot_data_raw[order(plot_data_raw[,order_by]),]
    plot_data_raw$variable <- factor(plot_data_raw$variable, levels=unique(plot_data_raw$variable), ordered=TRUE)
    plot_data_raw[[color_by]] <- factor(plot_data_raw[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data_raw) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data_raw) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- bquote(atop(.(maintitle_raw),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data_raw <- merge(plot_data_raw, colorDF, by = "variable")
    plot_data_raw[[color_by]] <- factor(plot_data_raw[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data_raw) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data_raw) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- maintitle_raw
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data_raw <- merge(plot_data_raw, tempdata, by = "variable")
    } else {
      plot_data_raw <- merge(plot_data_raw, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data_raw <- plot_data_raw[order(plot_data_raw[,order_by]),]
    plot_data_raw$variable <- factor(plot_data_raw$variable, levels=unique(plot_data_raw$variable), ordered=TRUE)
    plot_data_raw[[color_by]] <- factor(plot_data_raw[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data_raw) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data_raw) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- bquote(atop(.(maintitle_raw),atop(italic(paste("Ordered by ",.(order_by))),"")))
  }
  
  
  # facet plot #
  # if(!is.null(facet_by)) {
  #   if(is.null(facet_cols)) {
  #     p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
  #   } else {
  #     p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
  #   }
  # }
  
  # custom labels #
  if(!is.null(title_plot_raw)) title <- title_plot_raw
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)
  
  #### END OF RAW DATA PLOT ####
  
  #### ---------------------------------------------------------------------------- ####
  
  #### NORM DATA PLOT ####
  
  ## if both order_by and color_by are null ##
  if(is.null(order_by) & is.null(color_by)) {
    
    if(bw_theme==FALSE){
      p_norm <- ggplot2::ggplot(plot_data_norm) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p_norm <- ggplot2::ggplot(plot_data_norm) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- maintitle_norm
    
    ## if order_by is not null and color_by is ##
  } else if(!is.null(order_by) & is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    plot_data_norm <- merge(plot_data_norm, orderDF, by = "variable")
    
    # reorder levels #
    plot_data_norm <- plot_data_norm[order(plot_data_norm[,order_by]),]
    plot_data_norm$variable <- factor(plot_data_norm$variable, levels=unique(plot_data_norm$variable), ordered=TRUE)
    plot_data_norm[[color_by]] <- factor(plot_data_norm[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p_norm <- ggplot2::ggplot(plot_data_norm) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p_norm <- ggplot2::ggplot(plot_data_norm) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- bquote(atop(.(maintitle_norm),atop(italic(paste("Ordered by ",.(order_by))),"")))
    
    ## if color_by is not null and order_by is ##
  } else if(!is.null(color_by) & is.null(order_by)) {
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    plot_data_norm <- merge(plot_data_norm, colorDF, by = "variable")
    plot_data_norm[[color_by]] <- factor(plot_data_norm[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p_norm <- ggplot2::ggplot(plot_data_norm) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p_norm <- ggplot2::ggplot(plot_data_norm) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- maintitle_norm
    
    ## if neither order_by or color_by are null ##
  } else if(!is.null(order_by) & !is.null(color_by)) {
    if(order_by != "group_DF") {
      order_temp <- group_designation(omicsData, main_effects = order_by)
    } else {
      order_temp <- omicsData
    }
    orderDF <- attributes(order_temp)$group_DF
    colnames(orderDF)[1:2] <- c("variable", order_by)
    
    if(color_by != "group_DF") {
      color_temp <- group_designation(omicsData, main_effects = color_by)
    } else {
      color_temp <- omicsData
    }
    colorDF <- attributes(color_temp)$group_DF
    colnames(colorDF)[1:2] <- c("variable", color_by)
    
    # deal with case where both are equal #
    if(order_by != color_by) {
      tempdata <- merge(orderDF, colorDF, by = "variable")
      plot_data_norm <- merge(plot_data_norm, tempdata, by = "variable")
    } else {
      plot_data_norm <- merge(plot_data_norm, colorDF, by = "variable")
    }
    
    # reorder levels #
    plot_data_norm <- plot_data_norm[order(plot_data_norm[,order_by]),]
    plot_data_norm$variable <- factor(plot_data_norm$variable, levels=unique(plot_data_norm$variable), ordered=TRUE)
    plot_data_norm[[color_by]] <- factor(plot_data_norm[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))
    
    if(bw_theme==FALSE){
      p_norm <- ggplot2::ggplot(plot_data_norm) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p_norm <- ggplot2::ggplot(plot_data_norm) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }
    
    title <- bquote(atop(.(maintitle_norm),atop(italic(paste("Ordered by ",.(order_by))),"")))
  }
  
  
  # facet plot #
  # if(!is.null(facet_by)) {
  #   if(is.null(facet_cols)) {
  #     p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x")
  #   } else {
  #     p <- p + ggplot2::facet_wrap(formula(paste("~",facet_by)), scales = "free_x", ncol = facet_cols)
  #   }
  # }
  
  # custom labels #
  if(!is.null(title_plot_norm)) title <- title_plot_norm
  xlabel <- ifelse(is.null(x_lab), "Sample", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab
  
  # add additional features to plot #
  p_norm <- p_norm + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)
  
  #### END OF NORM DATA PLOT ####
  
  if(combined == FALSE){
    plot(p)
    par(ask=TRUE)
    plot(p_norm)
    par(ask=FALSE)
  }else{
    # the user wants the graphs in a single image
    #Rmisc:::multiplot(p, p_norm, cols=1)
    gridExtra:::grid.arrange(p, p_norm, ncol=1)
  }
  
  invisible(list(p, p_norm))
}

