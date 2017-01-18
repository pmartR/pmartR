#' Produce a plot of a MSomicsR S3 Object
#'
#' This function will provide a plot for a \code{omicsData} object, any of the filter objects in MSomicsR, or a \code{corRes} object.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'lipidData', or 'metabData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param corRes_object an object of class corRes. A correlation matrix of all samples.
#' @param filter_object a filter object for the respective \code{omicsData} class.
#' @param ... further arguments
#'
#' @return a ggplot summarizing the MSomicsR object
#'
#' @examples
#' donrun{
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
#' pca_res <- dim_reduction(pep_object2, k = 2)
#' plot(pca_res)
#'}
#' @details
#' Various further arguments can be specified depending on the class of the object being plotted.
#'
#' For plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData':
#' \tabular{ll}{
#' \code{order_by} \tab a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If order_by is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data. \cr
#' \code{color_by} \tab a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If color_by is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color. \cr
#' \code{facet_by} \tab a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL. \cr
#' \code{facet_cols} \tab an optional integer specifying the number of columns to show in the facet plot. \cr
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL, in which case "Samples" is used. \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL, in which case "log2(Robust Mahalanobis Distance)" is used. \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL, in which case "Sample Outlier Results /n p-value threshold = 0.xyz" is used, where 'xyz' is the pvalue_threshold supplied to the function. \cr
#' \code{legend_lab} \tab character string specifying the title label to use for the legend \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{legend_position} \tab character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Defaults to "right". \cr
#' }
#'
#' For plotting an S3 object of type 'corRes':
#' \tabular{ll}{
#' \code{interactive} \tab default value is FALSE. If TRUE, an interactive d3heatmap will be rendered, allowing you to zoom in on sections of the heatmap. \cr
#' \code{title_plot} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{x_lab} \tab logical indicator of whether to label the x-axis with the sample names. Defaults to TRUE. If FALSE, no x-axis labels are included in the plot. \cr
#' \code{y_lab} \tab logical indicator of whether to label the y-axis with the sample names. Defaults to TRUE. If FALSE, no y-axis labels are included in the plot. \cr
#' \code{colorbar_lim} \tab numeric pair of numeric values specifying the minimum and maximum values to use in the heatmap color bar. Defaults to 'c(NA, NA)', in which case ggplot2 automatically sets the minimum and maximum values based on the correlation values in the data. \cr
#' \code{title_size} See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \tab \cr
#' }
#'
#' For plotting an S3 object of type 'dimRes':
#' \tabular{ll}{
#' \code{x_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{y_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{legend_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{title_plot} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{title_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{x_lab_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{y_lab_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{bw_theme} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{legend_position} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' }
#'
#' The following arguments are available for plotting the \code{filter_object}s of type 'moleculeFilt', 'cvFilt', 'rmdFilt', 'proteomicsFilt', and 'imdanovaFilt':
#' \tabular{ll}{
#' \code{x_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{y_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{legend_lab} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{title_plot} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{title_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData'  \cr
#' \code{x_lab_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{y_lab_size} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{bw_theme} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' \code{legend_position} \tab See specification of extra arguments for plotting an S3 object of type 'pepData', 'proData', 'metabData', or 'lipidData' \cr
#' }
#'
#' The following give arguments specific to plotting individual \code{filter_objects}.
#'
#' For plotting an S3 object of type 'moleculeFilt':
#' \tabular{ll}{
#' \code{min_num} \tab an integer value specifying the minimum number of times each biomolecule must be observed across all samples in order to retain the biomolecule. Default value is 2. \cr
#' }
#'
#' For plotting an S3 object of type 'cvFilt':
#' \tabular{ll}{
#' \code{cv_threshold} \tab an integer value specifying the maximum coefficient of variation (CV) threshold for the biomolecules. Biomolecules with CV greater than this threshold will be filtered. Default is 150. \cr
#' }
#'
#' For plotting an S3 object of type 'rmdFilt':
#' \tabular{ll}{
#' \code{pvalue_threshold} \tab numeric value between 0 and 1, specifying the p-value, below which samples will be removed from the dataset. Default is 0.001. \cr
#' }
#'
#' For plotting an S3 object of type 'proteomicsFilt':
#' \tabular{ll}{
#' \code{min_num_peps} \tab an optional integer value between 1 and the maximum number of peptides that map to a protein in omicsData. The value specifies the minimum number of peptides that must map to a protein. Any protein with less than \code{min_num_peps} mapping to it will be returned as a protein that should be filtered. Default value is NULL. \cr
#' \code{degen_peps} \tab logical indicator of whether to filter out degenerate peptides (TRUE) or not (FALSE). Default value is FALSE.\cr
#' }
#'
#' For plotting an S3 object of type 'imdanovaFilt':
#' \tabular{ll}{
#' \code{min_nonmiss_anova} \tab integer value specifying the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Default value is 2. \cr
#' \code{min_nonmiss_gtest} \tab integer value specifying the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Default value is 3.\cr
#' }
#'
#' @export


#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param interactive default value is FALSE. If TRUE, an interactive d3heatmap will be rendered, allowing you to zoom in on sections of the heatmap.
#'@param x_lab logical indicator of whether to label the x-axis with the sample names. Defaults to TRUE. If FALSE, no x-axis labels are included in the plot.
#'@param y_lab logical indicator of whether to label the y-axis with the sample names. Defaults to TRUE. If FALSE, no y-axis labels are included in the plot.
#'@param colorbar_lim numeric pair of numeric values specifying the minimum and maximum values to use in the heatmap color bar. Defaults to 'c(NA, NA)', in which case ggplot2 automatically sets the minimum and maximum values based on the correlation values in the data.
plot.corRes <- function(corRes_object, ...){
  .plot.corRes(corRes_object, ...)
}
.plot.corRes <- function(corRes_object, interactive = FALSE, title_plot = NULL, x_lab=TRUE, y_lab=TRUE, colorbar_lim=c(NA, NA), title_size = 14){
  # check for a corRes object #
  if(class(corRes_object)[1] != "corRes") stop("object must be of class 'corRes'")
  # check plot title is a string #
  if(!is.null(title_plot)) {
    if(!is.character(title_plot)) stop("title_plot must be a character vector")
  }
  
  # Workaround for certain "check" warnings
  Var1 <- Var2 <- value <- NULL

  if(is.null(title_plot)){
    if(attributes(corRes_object)$data_norm == TRUE) {
      plot_title <- "Correlations Among Samples (Normalized Data)"
    }else {
      plot_title <- "Correlations Among Samples (Un-Normalized Data)"
    }
  }else{
    plot_title <- title_plot
  }

  if(!interactive) {
    pal <- colorRampPalette(c("blue","white","red"))

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

    if(attributes(corRes_object)$data_norm == TRUE) {
      heatmap <- heatmap + ggplot2::ggtitle(bquote(atop(plot_title)))
    }else {
      heatmap <- heatmap + ggplot2::ggtitle(plot_title)
    }

  } else {
    heatmap <- d3heatmap::d3heatmap(corRes_object, dendrogram = 'none', reorderfun = function(x) ordered(x, levels = rev(sort(unique(x)))))
  }

  return(heatmap)
}


#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param min_num an integer value specifying the minimum number of times each feature must be observed across all samples.
plot.moleculeFilt <- function(filter_object, min_num = NULL, ...) {
  .plot.moleculeFilt(filter_object, min_num, ...)
}

.plot.moleculeFilt <- function(filter_object, min_num = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {

  ## initial checks ##
  if(!is.null(min_num)) {
    # check that min_num is numeric >= 0 #
    if(!(class(min_num) %in% c("numeric","integer")) | min_num < 0) stop("min_num must be an integer >= 0")
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
  ## end of initial checks ##

  # get observation counts #
  cut_data <- table(cut(filter_object$Num_Observations, breaks = -1:max(filter_object$Num_Observations)))
  cumcounts <- cumsum(cut_data)
  pep_observation_counts <- data.frame(num_observations=0:(length(cumcounts)-1), frequency_counts=cumcounts)

  # make labels #
  xlabel <- ifelse(is.null(x_lab), "Minimum Number of Times a Biomolecule Appears Across All Samples", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Frequency of Biomolecules", y_lab)
  plot_title <- ifelse(is.null(title_plot), "Cumulative Frequency of Biomolecules in Samples", title_plot)

  # plot #
  if(bw_theme == FALSE){
    p <- ggplot2::ggplot(pep_observation_counts[-1, ]) +
      ggplot2::geom_rect(ggplot2::aes(xmin = num_observations - 0.5, xmax = num_observations + 0.5,
                                      ymin = 0, ymax = frequency_counts), fill="royalblue1", col="black") +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size),
                     axis.title.x = ggplot2::element_text(size = x_lab_size),
                     axis.title.y = ggplot2::element_text(size = y_lab_size))
  }else{
    p <- ggplot2::ggplot(pep_observation_counts[-1, ]) +
      ggplot2::theme_bw() +
      ggplot2::geom_rect(ggplot2::aes(xmin = num_observations - 0.5, xmax = num_observations + 0.5,
                                      ymin = 0, ymax = frequency_counts), fill="royalblue1", col="black") +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size),
                     axis.title.x = ggplot2::element_text(size = x_lab_size),
                     axis.title.y = ggplot2::element_text(size = y_lab_size))
  }


  if(!is.null(min_num)) {
    num_tested <- pep_observation_counts$frequency_counts[min_num+1]
    p <- p + ggplot2::annotate(geom = "rect", xmin = min_num - 0.5, xmax = min_num + 0.5,
                        ymin = 0, ymax = num_tested, fill="royalblue1", col="black", size=1.5)
  }

  return(p)
}



#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param min_num_peps an optional integer value between 1 and the maximum number of peptides that map to a protein in the data. The value specifies the minimum number of peptides that must map to a protein. Any protein with less than \code{min_num_peps} mapping to it will be returned as a protein that should be filtered. Default value is NULL.
#'@param degen_peps logical indicator of whether to filter out degenerate peptides (TRUE) or not (FALSE). Default value is FALSE.
plot.proteomicsFilt <- function(filter_object, min_num_peps = NULL, degen_peps = FALSE, ...) {
  .plot.proteomicsFilt(filter_object, min_num_peps, degen_peps, ...)
  }

.plot.proteomicsFilt <- function(filter_object, min_num_peps = NULL, degen_peps = FALSE,
                                 x_lab_pep = NULL, y_lab_pep = NULL, title.pep = NULL,
                                 x_lab_pro = NULL, y_lab_pro = NULL, title.pro = NULL,
                                 title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {

  # error checks for min_num_peps, if not NULL #
  if(!is.null(min_num_peps)) {
    # check that min_num_peps is numeric and >=1 #
    if(class(min_num_peps) != "numeric"| min_num_peps < 1) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if(min_num_peps %% 1 != 0) stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if(length(min_num_peps) != 1) stop("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if(min_num_peps > sum(filter_object$counts_by_pep$n)) stop("min_num_peps cannot be greater than the total number of peptides")
  }
  # check that degen_peps is logical #
  if(class(degen_peps) != "logical") stop("degen_peps must be either TRUE or FALSE")


  suppressMessages({
  # get counts
  pep_counts <- filter_object$counts_by_pep$n
  pro_counts <- filter_object$counts_by_pro$n

  # peptides plot
  ## make labels ##
  xlabel_pep <- ifelse(is.null(x_lab_pep), "Number of Peptides", x_lab_pep)
  ylabel_pep <- ifelse(is.null(y_lab_pep), "Number of Proteins", y_lab_pep)
  plot_title_pep <- ifelse(is.null(title.pep), "Number of Peptides \nPer Protein", title.pep)

  ## plot histogram
  if(bw_theme==FALSE){
    p1 <- ggplot2::ggplot(as.data.frame(pep_counts), ggplot2::aes(x=pep_counts)) +
      ggplot2::geom_histogram(fill="royalblue4", binwidth = 1) +
      ggplot2::ggtitle(plot_title_pep) +
      ggplot2::xlab(xlabel_pep) +
      ggplot2::ylab(ylabel_pep)  +
      ggplot2::scale_x_continuous(breaks = 1, labels = "1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }else{
    p1 <- ggplot2::ggplot(as.data.frame(pep_counts), ggplot2::aes(x=pep_counts)) +
      ggplot2::theme_bw() +
      ggplot2::geom_histogram(fill="royalblue4", binwidth = 1) +
      ggplot2::ggtitle(plot_title_pep) +
      ggplot2::xlab(xlabel_pep) +
      ggplot2::ylab(ylabel_pep)  +
      ggplot2::scale_x_continuous(breaks = 1, labels = "1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }


  ## get max count for line segment, draw segment on plot
  if(degen_peps) {
    max_counts <- ggplot2::ggplot_build(p1)$data[[1]]$count[1]
    p1 <- p1 + ggplot2::annotate(geom="rect", xmin=0.5, xmax=1.5, ymin=0, ymax=max_counts, col="black", alpha = 0, size = 1)
  }

  # proteins plot
  ## make labels ##
  xlabel_pro <- ifelse(is.null(x_lab_pro), "Number of Proteins", x_lab_pro)
  ylabel_pro <- ifelse(is.null(y_lab_pro), "Number of Peptides", y_lab_pro)
  plot_title_pro <- ifelse(is.null(title.pro), "Number of Proteins \nPer Peptide", title.pro)

  ## plot histogram
  if(bw_theme == FALSE){
    p2 <- ggplot2::ggplot(as.data.frame(pro_counts), ggplot2::aes(x=pro_counts)) +
      ggplot2::geom_histogram(fill="springgreen4", binwidth = 1) +
      ggplot2::ggtitle(plot_title_pro) +
      ggplot2::xlab(xlabel_pro) +
      ggplot2::ylab(ylabel_pro) +
      ggplot2::scale_x_continuous(breaks = 1, labels = "1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }else{
    p2 <- ggplot2::ggplot(as.data.frame(pro_counts), ggplot2::aes(x=pro_counts)) +
      ggplot2::theme_bw() +
      ggplot2::geom_histogram(fill="springgreen4", binwidth = 1) +
      ggplot2::ggtitle(plot_title_pro) +
      ggplot2::xlab(xlabel_pro) +
      ggplot2::ylab(ylabel_pro) +
      ggplot2::scale_x_continuous(breaks = 1, labels = "1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }



  ## get max count for line segment, draw segment on plot
  if(!is.null(min_num_peps)) {
    max_counts <- max(ggplot2::ggplot_build(p2)$data[[1]]$count)
    p2 <- p2 + ggplot2::annotate(geom="rect", xmin=0.5, xmax=min_num_peps + 0.5, ymin=0, ymax=max_counts, col="red", fill="red", alpha = 0.2, size = 0.5)
  }

  return(Rmisc::multiplot(p1, p2, cols = 2))
  })
}



#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param min_nonmiss_gtest the minimum number of non-missing feature values allowed per group for \code{gtest_filter}. Suggested value is 3.
#'@param min_nonmiss_anova the minimum number of non-missing feature values allowed per group for \code{anova_filter}. Suggested value is 2.
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
    if(!(class(min_nonmiss_anova) %in% c("numeric","integer")) | min_nonmiss_anova < 2) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is an integer #
    if(min_nonmiss_anova %% 1 != 0) stop("min_nonmiss_anova must be an integer >= 2")
    # check that min_nonmiss_anova is less than the minimum group size #
    if(min_nonmiss_anova > min(attributes(filter_object)$group_sizes$n_group)) stop("min_nonmiss_anova cannot be greater than the minimum group size")
  }
  if(!is.null(min_nonmiss_gtest)) {
    # check that min_nonmiss_gtest is numeric >= 3 #
    if(!(class(min_nonmiss_gtest) %in% c("numeric","integer")) | min_nonmiss_gtest < 3) stop("min_nonmiss_gtest must be an integer >= 3")
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
      ggplot2::geom_text(ggplot2::aes(fill = value, label = value), size=50/max(heat_df_melt$Var2)) +
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



#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param pvalue_threshold A threshold for the Robust Mahalanobis Distance (RMD) p-value. If \code{sampleID} is NULL (see \code{sampleID} below), a horizontal line is plotted at the RMD value that corresponds with the threshold, and all samples above the line have a p-value below the threshold. If \code{sampleID} is not NULL, \code{pvalue_threshold} will do nothing. Default value is NULL.
#'@param sampleID If specified, the plot produces a boxplot instead of a scatterplot. The \code{sampleID} input will place an "x" at the value for each of the metrics on the boxplots. Default value is NULL.
#'@param x_lab character string to be used for x-axis label. Defaults to NULL, in which case "Samples" is used.
#'@param y_lab character string to be used for y-axis label. Defaults to NULL, in which case "log2(Robust Mahalanobis Distance)" is used.
#'@param legend_lab character string to be used for the legend title. Defaults to NULL, in which case "Group" is used.
#'@param title_plot character string to be used for the plot title. Defaults to NULL, in which case "Sample Outlier Results /n p-value threshold = 0.xyz" is used, where 'xyz' is the pvalue_threshold supplied to the function.
#'@param title_size integer value indicating the font size for the plot title. Defaults to 14.
#'@param x_lab_size integer value indicating the font size for the plot title. Defaults to 11.
#'@param y_lab_size integer value indicating the font size for the plot title. Defaults to 11.
#'@param bw_theme logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used.
#'@param legend_position character string specifying one of "right", "left", "top", or "bottom" for the location of the legend. Deafults to "right".
plot.rmdFilt <- function(filter_object, pvalue_threshold = NULL, sampleID = NULL, ...) {
  .plot.rmdFilt(filter_object, pvalue_threshold, sampleID, ...)
}

.plot.rmdFilt <- function(filter_object, pvalue_threshold = NULL, sampleID = NULL, x_lab = NULL, y_lab = NULL, legend_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right") {

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
      ggplot2::geom_point(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), size=4, pch=4) +
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
        ggplot2::geom_point(data = dfsub[dfsub[,samp_id]==sampleID,], ggplot2::aes(x=rep(1,length(value)), y=value), size=4, pch=4) +
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
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names[1]), size=4)
    } else {
      p <- ggplot2::ggplot(filter_object) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names[1], pch=main_eff_names[2]), size=4)
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
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names), size=4, bg="gray") +
        ggplot2::geom_point(data = sub2, ggplot2::aes_string(x = samp_id, y="Log2.md", col=main_eff_names), alpha=0.5, size=4, bg="gray")
    } else {
      p <- ggplot2::ggplot(sub1) +
        ggplot2::geom_point(ggplot2::aes_string(x = samp_id, y="Log2.md", pch=main_eff_names[2], col=main_eff_names[1]), size=4, bg="gray") +
        ggplot2::geom_point(data = sub2, ggplot2::aes_string(x = samp_id, y="Log2.md", pch=main_eff_names[2], col=main_eff_names[1]), alpha=0.5, size=4, bg="gray")
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




#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param cv_threshold shades the area on the histogram below the given threshold. Default value is NULL.
plot.cvFilt <- function(filter_object, cv_threshold = NULL, ...) {
  .plot.cvFilt(filter_object, cv_threshold, ...)
}

.plot.cvFilt <- function(filter_object, cv_threshold = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE) {

  # checks for cv_threshold if not null
  if(!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if(!is.numeric(cv_threshold)) stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if(length(cv_threshold)>1) stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if(cv_threshold <= 1 | cv_threshold >= max(filter_object$CV_pooled, na.rm = TRUE)) stop("cv_threshold must be greater than 1 and less than the maximum CV_pooled value")
  }

  new_object <- filter_object[!is.na(filter_object$CV_pooled),]
  max_x_val <- attributes(filter_object)$max_x_val

  # get number of biomolecules with CV > max_x_val & display a warning #
  n_not_displayed <- sum(new_object$CV_pooled > max_x_val)
  warning(paste("For display purposes, biomolecules with pooled CV greater than ", round(max_x_val, 2), " are not displayed in the graph. This corresponds to ", n_not_displayed, " biomolecules.", sep=""))

  # labels
  plot_title <- ifelse(is.null(title_plot), "Coefficient of Variation (CV)", title_plot)
  xlabel <- ifelse(is.null(x_lab), "Pooled CV", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Count", y_lab)

  if(bw_theme==FALSE){
    p <- ggplot2::ggplot(new_object) +
      ggplot2::geom_histogram(ggplot2::aes(x=CV_pooled), fill = "steelblue1", breaks = seq(0,max_x_val,length.out = 20)) +
      ggplot2::scale_x_continuous(breaks = pretty(new_object$CV_pooled[new_object$CV_pooled < max_x_val], n=10), limits = c(0, max_x_val)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }else{
    p <- ggplot2::ggplot(new_object) +
      ggplot2::theme_bw() +
      ggplot2::geom_histogram(ggplot2::aes(x=CV_pooled), fill = "steelblue1", breaks = seq(0,max_x_val,length.out = 20)) +
      ggplot2::scale_x_continuous(breaks = pretty(new_object$CV_pooled[new_object$CV_pooled < max_x_val], n=10), limits = c(0, max_x_val)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=5)) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size))
  }


  if(!is.null(cv_threshold)) {
    if(is.null(title_plot)) {
      plot_title <- bquote(atop("Coefficient of Variation (CV)",atop(italic(paste("CV Threshold = ",.(cv_threshold))),"")))
    }
    p <- p + ggplot2::annotate("rect", xmin = 0, xmax = cv_threshold, ymin = 0, ymax = Inf, alpha = 0.3, fill = "steelblue1") +
      ggplot2::ggtitle(plot_title)
  }

  return(p)
}


#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
plot.customFilt <- function(filter_object, ...) {
  .plot.customFilt(filter_object, ...)
}

.plot.customFilt <- function(filter_object, x_lab = NULL, y_lab = NULL, title_plot = NULL) {

  warning("There is no plot method for objects of class 'customFilt'. See summary.customFilt instead.")
}


#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
#'@param order_by a character string specifying a main effect by which to order the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If \code{order_by} is "group_DF", the boxplots will be ordered by the group variable from the group_designation function. If NULL (default), the boxplots will be displayed in the order they appear in the data.
#'@param color_by a character string specifying a main effect by which to color the boxplots. This main effect must be found in the column names of f_data in the omicsData object. If \code{color_by} is "group_DF", the boxplots will be colored by the group variable from the group_designation function. If NULL (default), the boxplots will have one default color.
#'@param facet_by a character string specifying a main effect with which to create a facet plot. This main effect must be found in the column names of f_data in the omicsData object. Default value is NULL.
#'@param facet_cols an optional integer specifying the number of columns to show in the facet plot.
plot.pepData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.pepData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.pepData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right") {

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
    if(is.null(facet_by)) stop("facet_by cannot be NULL when fact_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  ## end of initial checks ##


  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)

  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Peptide Data",
                      "Boxplots of Un-Normalized Peptide Data")

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
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
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab

  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)

  return(p)
}




#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
plot.proData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.proData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.proData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right") {

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
    if(is.null(facet_by)) stop("facet_by cannot be NULL when fact_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  ## end of initial checks ##


  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)

  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Protein Data",
                      "Boxplots of Un-Normalized Protein Data")

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

    if(bw_theme == FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        + ggplot2::theme_bw()
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
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
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab

  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)

  return(p)
}




#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
plot.lipidData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.lipidData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.lipidData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right") {

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
    if(is.null(facet_by)) stop("facet_by cannot be NULL when fact_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  ## end of initial checks ##


  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)

  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Lipid Data",
                      "Boxplots of Un-Normalized Lipid Data")

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))

    p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
      ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                     axis.title.x = ggplot2::element_text(size=x_lab_size),
                     axis.title.y = ggplot2::element_text(size=y_lab_size),
                     legend.position = legend_position)
    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }


    title <- maintitle

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
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
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
  legend_title <- color_by
  if(!is.null(legend_lab)) legend_title <- legend_lab

  # add additional features to plot #
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::xlab("Sample") +
    ggplot2::ggtitle(title) + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_discrete(legend_title)

  return(p)
}





#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
plot.metabData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, ...) {
  .plot.metabData(omicsData, order_by, color_by, facet_by, facet_cols, ...)
}

.plot.metabData <- function(omicsData, order_by = NULL, color_by = NULL, facet_by = NULL, facet_cols = NULL, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_lab = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme=FALSE, legend_position = "right") {

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
    if(is.null(facet_by)) stop("facet_by cannot be NULL when fact_cols is specified")
    if(length(facet_cols)>1) stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols)) stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 | facet_cols <= 0) stop("facet_cols must be an integer greater than zero")
  }
  ## end of initial checks ##


  # organize data #
  e_data <- omicsData$e_data
  e_data_cname <- attributes(omicsData)$cnames$edata_cname
  plot_data <- reshape2::melt(e_data, id = e_data_cname, na.rm = TRUE)

  maintitle <- ifelse(attributes(omicsData)$data_info$data_norm,
                      "Boxplots of Normalized Metabolite Data",
                      "Boxplots of Un-Normalized Metabolite Data")

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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
    plot_data[[color_by]] <- factor(plot_data[[color_by]], levels = unique(factor(omicsData$f_data[[color_by]])))

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = value), fill = "deepskyblue1") +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- maintitle

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

    if(bw_theme==FALSE){
      p <- ggplot2::ggplot(plot_data) + ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }else{
      p <- ggplot2::ggplot(plot_data) +
        ggplot2::theme_bw() +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "variable", y = "value", fill = color_by)) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=title_size),
                       axis.title.x = ggplot2::element_text(size=x_lab_size),
                       axis.title.y = ggplot2::element_text(size=y_lab_size),
                       legend.position = legend_position)
    }

    title <- bquote(atop(.(maintitle),atop(italic(paste("Ordered by ",.(order_by))),"")))
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
  ylabel <- ifelse(is.null(y_lab), "Value", y_lab)
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

  return(p)
}



#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
plot.dimRes <- function(dimRes_object, ...) {
  .plot.dimRes(dimRes_object, ...)
}

.plot.dimRes <- function(dimRes_object, x_lab = NULL, y_lab = NULL, legend_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE, legend_position = "right") {

  plotdata <- data.frame(SampleID = dimRes_object$SampleID, PC1 = dimRes_object$PC1, PC2 = dimRes_object$PC2)

  # if there is a group designation #
  if(!is.null(attr(dimRes_object,"group_DF"))) {
    group_DF <- attr(dimRes_object,"group_DF")

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

      plotdata[,4:5] <- group_DF[,3:4]
      names(plotdata[,4:5]) <- main_eff_names

      color_var <- main_eff_names[1]
      pch_var <- main_eff_names[2]

      #plotdata[, colnames(plotdata) == pch_var] <- as.factor(plotdata[, colnames(plotdata) == pch_var])
      #plotdata[, colnames(plotdata) == color_var] <- as.factor(plotdata[, colnames(plotdata) == color_var])

      if(length(legend_lab) > length(main_eff_names)) warning("legend_lab length is greater than the number of main effects. Only the first two entries will be used.")

    } else {
      plotdata$Group <- group_DF$Group

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
      ggplot2::geom_point(ggplot2::aes_string(col = color_var, pch = pch_var)) +
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
      ggplot2::geom_point(ggplot2::aes_string(col = color_var, pch = pch_var)) +
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



#'@export
#'@rdname plot-MSomicsR
#'@name plot-MSomicsR
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
  #   if(is.null(facet_by)) stop("facet_by cannot be NULL when fact_cols is specified")
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


  #return(list(p, p_norm))
}

