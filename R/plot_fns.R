#' Plots an object of class dataRes
#'
#' For plotting an S3 object of type dataRes
#'
#' @param dataRes_obj An object of class dataRes. A dataRes object is created by
#'   the \code{edata_summary} function.
#' @param metric A character string indicating which metric to use in plot, one
#'   of 'mean', 'median', 'sd, 'pct_obs', 'min', or 'max'
#' @param density logical default to FALSE, if set to true a density plot of the
#'   specified metric is returned
#' @param ncols An integer specifying the number columns for the histogram
#'   facet_wrap. This argument is used when \code{metric} is not null. The
#'   default is NULL.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label when the metric
#'   argument is NULL. The default is NULL in which case the x-axis label will
#'   be "count".
#' @param x_lab_sd A character string used for the x-axis label for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param x_lab_median A character string used for the x-axis label for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param y_lab_sd A character string used for the y-axis label for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param y_lab_median A character string used for the y-axis label for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#' @param title_lab A character string specifying the plot title when the
#'   \code{metric} argument is NULL.
#' @param title_lab_sd A character string used for the plot title for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param title_lab_median A character string used for the plot title for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", or "bottom". The default is
#'   "right".
#' @param point_size An integer specifying the size of the points. The default
#'   is 2.
#' @param bin_width An integer indicating the bin width in a histogram. The
#'   default is 0.5.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#'
#' @details This function can only create plots of dataRes objects whose 'by' ==
#'   'molecule' and 'groupvar' attribute is non NULL
#'
#' @return plots ggplot2 object
#'
#' @examples
#' \dontrun{
#' library (pmartRdata)
#'
#' data(lipid_object)
#'
#' lipid_object = edata_transform(lipid_object, "log2")
#'
#' result = edata_summary(omicsData = lipid_object,
#'                        by = "molecule",
#'                        groupvar = "Condition")
#'
#' plot(result)
#' }
#'
#' @rdname plot-dataRes
#'
#' @export
#'
plot.dataRes <- function (dataRes_obj, metric = NULL, density = FALSE,
                          ncols = NULL, interactive = FALSE, x_lab = NULL,
                          x_lab_sd = NULL, x_lab_median = NULL, y_lab = NULL,
                          y_lab_sd = NULL, y_lab_median = NULL, x_lab_size = 11,
                          y_lab_size = 11, x_lab_angle = NULL, title_lab = NULL,
                          title_lab_sd = NULL, title_lab_median = NULL,
                          title_lab_size = 14, legend_lab = NULL,
                          legend_position = "right", point_size = 2,
                          bin_width = 1, bw_theme = TRUE, palette = NULL) {

  # Preliminaries --------------------------------------------------------------

  #check that attr(dataRes_obj, "by") == "molecule"
  if(attr(dataRes_obj, "by") != "molecule") {

    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop (paste("can only plot a dataRes object if its 'by' attribute is equal",
                "to 'molecule'",
                sep = " "))

  }

  #check that attr(dataRes_obj, "groupvar") is not NULL
  if(is.null(attr(dataRes_obj, "groupvar"))) {

    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop (paste("can only plot a dataRes object if its 'groupvar' attribute is",
                "not NULL",
                sep = " "))

  }

  # Check if the data is on a log scale. If it is not the histograms will not
  # work properly when metric = mean, sd, ... because of the wide range of
  # abundance values.
  if (attr(dataRes_obj, "data_scale") == "abundance") {

    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop ("Data must be on the log scale to plot.")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # extract molecule name
  edata_cname = attr(dataRes_obj, "cnames")$edata_cname

  # Make labels for various aspects of the plot. (Super tedious!!!)
  xLabelSd <- if (is.null(x_lab_sd)) "mean" else x_lab_sd
  yLabelSd <- if (is.null(y_lab_sd)) "sd" else y_lab_sd
  xLabelMedian <- if (is.null(x_lab_median)) "mean" else x_lab_median
  yLabelMedian <- if (is.null(y_lab_median)) "median" else y_lab_median
  plotTitleSd <- if (is.null(title_lab_sd))
    "Mean x Standard Deviation" else
      title_lab_sd
  plotTitleMedian <- if (is.null(title_lab_median))
    "Mean x Median" else
      title_lab_median
  legendLabel <- if (is.null(legend_lab)) "Group" else legend_lab

  # Build beautiful plots ------------------------------------------------------

  if (is.null(metric)) {

    #scatterplots
    #subseting dataRes_obj object
    mean <- dataRes_obj$mean
    median <- dataRes_obj$median
    sd <- dataRes_obj$sd

    #melting data frames from dataRes object
    mean_melt <- reshape2::melt(mean, id.vars = edata_cname)
    names(mean_melt)[3] <- "mean"
    sd_melt <- reshape2::melt(sd, id.vars = edata_cname)
    names(sd_melt)[3] <- "sd"
    median_melt <- reshape2::melt(median, id.vars = edata_cname)
    names(median_melt)[3] <- "median"

    data_mean_sd <- merge(mean_melt,
                          sd_melt,
                          by = c(edata_cname, "variable"))
    data_mean_median <- merge(mean_melt,
                              median_melt,
                              by = c(edata_cname, "variable"))

    q <- ggplot2::ggplot(data_mean_sd,
                         ggplot2::aes(x = mean,
                                      y = sd,
                                      color = variable)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::xlab(xLabelSd) +
      ggplot2::ylab(yLabelSd) +
      ggplot2::ggtitle(plotTitleSd)

    p <- ggplot2::ggplot(data_mean_median,
                         ggplot2::aes(x = mean,
                                      y = median,
                                      color = variable)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::xlab(xLabelMedian) +
      ggplot2::ylab(yLabelMedian) +
      ggplot2::ggtitle(plotTitleMedian)

    # Want the black and white theme? As you wish.
    if (bw_theme) {

      q <- q + ggplot2::theme_bw()
      p <- p + ggplot2::theme_bw()

    }

    # Create a generic theme that will be used for both plot.
    the_theme <- ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

    # Add information to the plot from theme. Must be done after theme_bw()
    # otherwise the black and white theme overwrites everything that was
    # specified.
    q <- q +
      ggplot2::labs(color = legendLabel) +
      the_theme
    p <- p +
      ggplot2::labs(color = legendLabel) +
      the_theme

    # Want a plot with beautiful colors? As you wish.
    if (!is.null(palette)) {

      # Use the ColorBrewer color and create the legend title
      p <- p +
        ggplot2::scale_color_brewer(palette = palette,
                                    name = legendLabel)
      q <- q +
        ggplot2::scale_color_brewer(palette = palette,
                                    name = legendLabel)

    }

    # Want an interactive plot? As you wish.
    if (interactive) {

      q <- plotly::ggplotly(q)
      p <- plotly::ggplotly(p)

      plotly::subplot(q, p, nrows = 1)

    } else {

      # Combine the plots into one plot when interactive is FALSE.

      q + p


    }


  } else if (!is.null(metric)) {

    if (!(metric %in% c('mean', 'median','sd', 'pct_obs', 'min', 'max'))) {

      # There is a shortage of perfect plots in the world. It is a shame to ruin
      # this one.
      stop ("metric must be one of mean, median, sd, pct_obs, min or max")

    }

    if (!is.logical(density)) stop ("density must be either TRUE or FALSE")

    # if density == F, will plot faceted histograms.
    if (density == FALSE) {

      # More tedious label creating .... (deep sigh).
      xlabel <- if (is.null(x_lab)) metric else x_lab
      ylabel <- if (is.null(y_lab)) "count" else y_lab
      plotTitle <- if (is.null(title_lab))
        paste("Histograms for ", metric, sep = "") else
          title_lab

      #subsetting dataRes object
      data = dataRes_obj[[metric]]
      data_melt = reshape2::melt(data, id.vars = edata_cname)

      r <- ggplot2::ggplot(data_melt,
                           ggplot2::aes(x = value,
                                        fill = variable)) +
        ggplot2::geom_histogram(binwidth = bin_width, colour = "white") +
        ggplot2::facet_wrap(~variable, ncol = ncols)

      # Create a density plot. Following code runs when density = TRUE.
    } else {

      # More tedious label creating .... (defeated sigh).
      xlabel <- if (is.null(x_lab)) metric else x_lab
      ylabel <- if (is.null(y_lab)) "density" else y_lab
      plotTitle <- if (is.null(title_lab))
        paste("Density plots for ", metric, sep = "") else
          title_lab

      #if density == T, will plot geom_density
      data = dataRes_obj[[metric]]
      data_melt = reshape2::melt(data, id.vars = edata_cname)

      r <- ggplot2::ggplot(data_melt,
                           ggplot2::aes(x = value,
                                        colour = variable)) +
        ggplot2::geom_density()

    }

    # Want the black and white theme? As you wish.
    if (bw_theme) r <- r + ggplot2::theme_bw() +
        ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    # Add the remaining labels and text sizes/orientation to the plot.
    r <- r +
      ggplot2::ggtitle(plotTitle) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::labs(color = legendLabel) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = title_lab_size),
        axis.title.x = ggplot2::element_text(size = x_lab_size),
        axis.title.y = ggplot2::element_text(size = y_lab_size),
        axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
        legend.position = legend_position
      )

    # Want a plot with gorgeous colors? As you wish.
    if (!is.null(palette)) {

      # Use the ColorBrewer color and create the legend title
      r <- r +
        ggplot2::scale_fill_brewer(palette = palette,
                                   name = legendLabel)

    }

    # Want an interactive plot? As you wish.
    if (interactive) r <- plotly::ggplotly(r)

    return (r)

  }

}

#' Plots an object of class isobaricnormRes
#'
#' Creates box plots for an S3 object of type 'isobaricnormRes'
#'
#' @param isobaricnormRes_obj an object of type isobaricnormRes, created by
#'   \code{\link{normalize_isobaric}}
#' @param order Logical. If TRUE the samples will be ordered by the column of
#'   f_data containing the experiment/plate information.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#'
#' @return plots ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' data(isobaric_object)
#'
#' isobaric_object = edata_transform(isobaric_object, "log2")
#'
#' result = normalize_isobaric(isobaric_object, exp_cname = "Set",
#'                             apply_norm = FALSE,
#'                             channel_cname = "iTRAQ.Channel",
#'                             refpool_channel = "116")
#'
#' plot(result)
#' }
#'
#' @importFrom rlang .data
#'
#' @rdname plot-isobaricnormRes
#'
#' @export
#'
plot.isobaricnormRes <- function (isobaricnormRes_obj, order = FALSE,
                                  interactive = FALSE, x_lab = NULL,
                                  y_lab = NULL, x_lab_size = 11,
                                  y_lab_size = 11, x_lab_angle = NULL,
                                  title_lab = NULL, title_lab_size = 14,
                                  legend_lab = NULL, legend_position = "none",
                                  bw_theme = TRUE, palette = NULL) {

  # Preliminaries --------------------------------------------------------------

  # Check if input is an isobaricnormRes class object.
  if (!inherits(isobaricnormRes_obj, "isobaricnormRes")) {

    stop ("object must be of class 'isobaricnormRes'")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # extracting attributes from isobaricnormRes_obj
  exp_cname = attr(isobaricnormRes_obj, "isobaric_info")$exp_cname
  fdata_cname = attr(isobaricnormRes_obj, "cnames")$fdata_cname
  edata_cname <- attr(isobaricnormRes_obj, "cnames")$edata_cname

  # Transform the isobaricnormRes data frames into a format usable by ggplot2.
  tall_data <- prime_iso(isonormRes = isobaricnormRes_obj,
                         exp_cname = exp_cname,
                         fdata_cname = fdata_cname,
                         edata_cname = edata_cname)

  # Do all the tedious plot label crap.
  xlabel <- if (is.null(x_lab)) "Reference Sample" else x_lab
  ylabel <- if (is.null(y_lab)) "Log Abundance" else y_lab
  plot_title <- if (is.null(title_lab))
    "Reference Sample Profile" else
      title_lab
  legendLabel <- if (is.null(legend_lab)) "Sample" else legend_lab

  # Create pretty plots --------------------------------------------------------

  # If order is TRUE order the box plots by experiment name/value.
  if (order == TRUE) {

    xlabel <- if (is.null(x_lab)) exp_cname else x_lab

    p <- ggplot2::ggplot(data = tall_data,
                         ggplot2::aes(x = .data[[exp_cname]],
                                      y = values,
                                      fill = Sample))

    # Otherwise separate the box plots by sample name.
  } else {

    p <- ggplot2::ggplot(data = tall_data,
                         ggplot2::aes(x = .data[[fdata_cname]],
                                      y = values,
                                      fill = Sample))

  }

  # Want your plot to have a black and white background? As you wish.
  # The black and white theme needs to come before the code that creates the box
  # plot. If the theme_bw code follows the code that creates the box plot it
  # will add the legend back to the graph (if legend_position = "none").
  if (bw_theme == TRUE) p <- p + ggplot2::theme_bw()

  # Create the box plot with all of the users input.
  p <- p +
    ggplot2::geom_boxplot() +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::labs(color = legendLabel) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Want to use beautiful non-default colors? As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(palette = palette,
                                 name = legendLabel)

  }

  # Want an interactive plot? As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

# Takes the isobaricnormRes object and returns a data frame in the correct
# format for plotting in ggplot2.
prime_iso <- function (isonormRes, exp_cname,
                       fdata_cname, edata_cname) {

  # Find the column in e_data pertaining to the peptide IDs.
  id_col <- which(names(isonormRes$e_data) == edata_cname)

  # Stack the columns of e_data. This will include the sample names as an
  # additional column.
  tall_data <- stack(isonormRes$e_data[, -id_col])

  # Change the name of the column containing sample names to match the
  # corresponding column name in f_data.
  names(tall_data)[2] <- fdata_cname

  # Include the experiment column in tall_data.
  tall_data <- dplyr::inner_join(tall_data,
                                 isonormRes$f_data[, c(fdata_cname, exp_cname)],
                                 by = fdata_cname)

  # Extract the indices for each column in tall_data.
  exp_cname_ind <- which(names(tall_data) == exp_cname)
  fdata_cname_ind <- which(names(tall_data) == fdata_cname)
  values_col_ind <- which(names(tall_data) == "values")

  # Convert the experiment column into a factor.
  tall_data[[exp_cname]] <- factor(tall_data[[exp_cname]])

  # Return the tall data!
  return (tall_data)

}

#' Plots an object of class nmrnormRes
#'
#' Creates a scatter plot for an S3 object of type 'nmrnormRes'
#'
#' @param nmrnormRes_obj an object of type nmrnormRes, created by
#'   \code{\link{normalize_nmr}}
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 2.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#'
#' @return plots ggplot2 object
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(nmr_object_identified)
#'
#' nmr_object = edata_transform(nmr_object_identified, "log2")
#' nmr_norm = normalize_nmr(nmr_object,
#'                          apply_norm = FALSE,
#'                          metabolite_name = "unkm1.53")
#' plot(nmr_norm)
#'
#' # alternate specification: #
#' data(nmr_object_identified)
#'
#' nmr_object = edata_transform(nmr_object, "log2")
#' nmr_norm = normalize_nmr(nmr_object,
#'                          apply_norm = FALSE,
#'                          sample_property_cname = "Concentration")
#' plot(nmr_norm)
#' }
#'
#' @rdname plot-nmrnormRes
#'
#' @export
#'
plot.nmrnormRes <- function (nmrnormRes_obj, interactive = FALSE,
                             x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                             y_lab_size = 11, x_lab_angle = 90,
                             title_lab = NULL, title_lab_size = 14,
                             legend_lab = NULL, legend_position = "none",
                             point_size = 2, bw_theme = TRUE) {

  # Preliminaries --------------------------------------------------------------

  if (!inherits(nmrnormRes_obj, "nmrnormRes")) {

    # My name is Evan Martin. You defiled pmart. Prepare to die.
    stop("object must be of class 'nmrnormRes'")

  }

  #extracting attributes from nmrnormRes_obj
  sample_property_cname <- attr(nmrnormRes_obj,
                                "nmr_info")$sample_property_cname
  metabolite_name <- attr(nmrnormRes_obj, "nmr_info")$metabolite_name
  fdata_cname <- attr(nmrnormRes_obj, "cnames")$fdata_cname

  # organize nmrnormRes_obj
  data <- data.frame(Sample = nmrnormRes_obj$Sample,
                     value = nmrnormRes_obj$value)

  # Do all the tedious plot label crap.
  xlabel <- if (is.null(x_lab)) "Sample ID" else x_lab
  legendLabel <- if (is.null(legend_lab)) "Sample" else legend_lab
  # Create y-axis and plot title labels when a reference metabolite is used.
  if (!is.null(metabolite_name)) {

    ylabel <- if (is.null(y_lab)) metabolite_name else y_lab
    plot_title <- if (is.null(title_lab))
      "Reference Metabolite Profile" else
        title_lab

    # Create y-axis and plot title labels when the sample property is used.
  } else if (!is.null(sample_property_cname)) {

    ylabel <- if (is.null(y_lab)) sample_property_cname else y_lab
    plot_title <- if (is.null(title_lab))
      "Sample Property Profile" else
        title_lab

  }

  # Generate dazzling plots ----------------------------------------------------

  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes(x = Sample,
                                    y = value))

  # Want your plot to have a black and white background? As you wish. The black
  # and white theme needs to come before the code that creates the scatter plot.
  # If the theme_bw code follows the code that creates the scatter plot it will
  # add the legend back to the graph (if legend_position = "none").
  if (bw_theme == TRUE) p <- p + ggplot2::theme_bw()

  # Add the points and labels to the dazzling plot.
  p <- p +
    ggplot2::geom_point(size = point_size) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::labs(color = legendLabel) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Want an interactive plot? As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

#' Plots an object of class SPANSRes
#'
#' For plotting an S3 object of type 'SPANSRes'
#'
#' @param SPANSRes_obj an object of the class 'SPANSRes', usually created by
#'   \code{\link{spans_procedure()}}.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is NULL.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param color_low A character string specifying the color of the gradient for
#'   low values.
#' @param color_high A character string specifying the color of the gradient for
#'   high values.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#'
#' @return plots a plotly object
#'
#' @rdname plot-SPANSRes
#'
#' @export
#'
plot.SPANSRes <- function (SPANSRes_obj, interactive = FALSE,
                           x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                           y_lab_size = 11, x_lab_angle = NULL,
                           title_lab = NULL, title_lab_size = 14,
                           legend_lab = NULL, legend_position = "right",
                           color_low = NULL, color_high = NULL,
                           bw_theme = TRUE) {

  # Preliminaries --------------------------------------------------------------

  if (!inherits(SPANSRes_obj, "SPANSRes")) {

    # Suffer the wrath of Dread Pirate Roberts!!!!!
    stop("object must be of class 'SPANSRes'")

  }

  # plotting object with numeric SPANS_score, the normalization method, and a
  # modified string specifying the subset method + parameters for that method
  SPANSRes_obj <- SPANSRes_obj %>%
    dplyr::mutate(ss_par = paste(subset_method, parameters, sep = " | "),
                  SPANS_score = as.numeric(SPANS_score)) %>%
    dplyr::left_join(attr(SPANSRes_obj, "method_selection_pvals"),
                     by = c("subset_method",
                            "normalization_method",
                            "parameters"))

  # get subset/normalization names for the best scored methods
  best_ss <- SPANSRes_obj %>%
    dplyr::top_n(1, wt = SPANS_score) %>%
    {.$ss_par}

  best_norm <- SPANSRes_obj %>%
    dplyr::top_n(1, wt = SPANS_score) %>%
    {.$normalization_method}

  # Do all the tedious plot label crap.
  xlabel <- if (is.null(x_lab)) "Normalization Method" else x_lab
  ylabel <- if (is.null(y_lab)) "Subset Parameters" else y_lab
  plot_title <- if (is.null(title_lab)) NULL else title_lab
  legendLabel <- if (is.null(legend_lab)) "Score" else legend_lab

  # Produce magnificent plots --------------------------------------------------

  p <- ggplot2::ggplot(data = SPANSRes_obj) +
    ggplot2::geom_tile(ggplot2::aes(x = normalization_method,
                                    y = ss_par,
                                    alpha = 1),
                       color = 'black') +
    ggplot2::geom_tile(ggplot2::aes(x = normalization_method,
                                    y = ss_par,
                                    fill = SPANS_score),
                       color = 'black') +
    ggplot2::geom_point(data = SPANSRes_obj %>%
                          dplyr::filter(ss_par %in% best_ss,
                                        normalization_method %in% best_norm),
                        ggplot2::aes(x = normalization_method,
                                     y = ss_par, shape = '1')) +
    ggplot2::scale_alpha_continuous(name = 'Not Scored',
                                    labels = '') +
    ggplot2::scale_shape_discrete(name = 'Best Scores',
                                  labels = '') +
    ggplot2::scale_fill_gradient(
      low = if (is.null(color_low)) "#132B43" else color_low,
      high = if (is.null(color_high)) "#56B1F7" else color_high,
    ) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel)

  # Want the black and white theme? As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  p <- p +
    ggplot2::labs(fill = legendLabel) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position,
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  # Want an interactive plot? As you wish.
  if (interactive) {

    p <- plotly::plot_ly(
      SPANSRes_obj,
      x = ~normalization_method,
      y = ~ss_par,
      z = ~SPANS_score,
      hoverinfo = 'text',
      text = ~paste('</br> F(-Log10(HSmPV)):',
                    F_log_HSmPV,
                    '</br> F(Log10(NSmPV)): ',
                    F_log_NSmPV,
                    '</br> Scale p-value: ',
                    scale_p_value,
                    '</br> Location p-value',
                    location_p_value),
      colors = grDevices::colorRamp(
        c(if (is.null(color_low)) "#132B43" else color_low,
          if (is.null(color_high)) "#56B1F7" else color_high)
      ),
      type = "heatmap") %>%
      plotly::add_trace(x = best_norm,
                        y = best_ss,
                        type = 'scatter',
                        mode = "markers",
                        marker = list(color = "black"),
                        name = "Top SPANS scores",
                        inherit = FALSE) %>%
      plotly::colorbar(title = "SPANS score") %>%
      plotly::layout(plot_bgcolor = 'black',
                     xaxis = list(title = "Normalization Method"),
                     yaxis = list(title = "Subset Method"),
                     showlegend = TRUE
    )

  }

  return (p)

}

#' Plots an object of class naRes
#'
#' For plotting an S3 object of type 'naRes'
#'
#' @param naRes_obj A list of two data frames, one contains NA values by sample,
#'   the second contains NA values by molecule
#' @param omicsData an object of class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param plot_type A character string specifying which type of plot to produce.
#'   The two options are 'bar' or 'scatter'.
#' @param order_by A character string specifying a main effect by which to order
#'   the bar plot. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the bar plot
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the bar plot will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the bar plot. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the bar plot
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the bar plot will have one default color.
#'
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab_bar A character string used for the x-axis label for the bar
#'   plot.
#' @param x_lab_scatter A character string used for the x-axis label for the
#'   scatter plot.
#' @param y_lab_bar A character string used for the y-axis label for the bar
#'   plot.
#' @param y_lab_scatter A character string used for the y-axis label for the
#'   scatter plot.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#' @param title_lab_bar A character string used for the plot title when
#'   \code{plot_type} is 'bar'.
#' @param title_lab_scatter A character string used for the plot title when
#'   \code{plot_type} is 'scatter'.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab_bar A character string specifying the legend title when
#'   creating a bar plot.
#' @param legend_lab_scatter A character string specifying the legend title when
#'   creating a scatter plot.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", or "bottom". The default is
#'   "right".
#' @param point_size An integer specifying the size of the points. The default
#'   is 2.
#' @param text_size An integer specifying the size of the text (number of
#'   missing values by sample) within the bar plot. The default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count Logical. Indicates whether the missing value counts by
#'   sample will be displayed on the bar plot. The default is TRUE.
#' @param coordinate_flip Logical. Indicates whether the x and y axes will be
#'   flipped. The default is FALSE.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @return plots ggplot2 object
#'
#' @details This function takes in an object of class naRes and creates either a
#'   bar or scatter plot of missing values. When plot_type = 'bar', a sample
#'   name by missing values count bar chart is returned. When plot_type =
#'   'scatter' a mean intensity vs number of missing values (per molecule)
#'   scatter plot is returned. Note: If the omicsData object has had
#'   \code{\link{group_designation}} applied to it, the points in the plot will
#'   be colored by group.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' data("lipid_object")
#'
#' result<- missingval_result(lipid_object)
#'
#' plot(result, plog_type = "bar", x_lab_angle = 50)
#' }
#'
#' @rdname plot-naRes
#'
#' @export
#'
plot.naRes <- function (naRes_obj, omicsData, plot_type = "bar",
                        order_by = NULL, color_by = NULL,
                        interactive = FALSE, x_lab_bar = NULL,
                        x_lab_scatter = NULL, y_lab_bar = NULL,
                        y_lab_scatter = NULL, x_lab_size = 11, y_lab_size = 11,
                        x_lab_angle = 60, title_lab_bar = NULL,
                        title_lab_scatter = NULL, title_lab_size = 14,
                        legend_lab_bar = NULL, legend_lab_scatter = NULL,
                        legend_position = "right", point_size = 2,
                        text_size = 3, bar_width = 0.8, bw_theme = TRUE,
                        palette = NULL, display_count = TRUE,
                        coordinate_flip = FALSE, use_VizSampNames = FALSE) {

  # Preliminaries --------------------------------------------------------------

  # check for a naRes object #
  if(!inherits(naRes_obj, "naRes")) stop("object must be of class 'naRes'")

  # Check that omicsData is the correct class.
  if (!inherits(omicsData, c("proData", "pepData", "lipidData",
                                 "metabData", "nmrData"))) {

    # Fezzik, tear his arms off.
    stop ("omicsData is not an appropriate class")

  }

  # Check that type is either bar or scatter.
  if(!(plot_type %in% c("bar", "scatter"))) {

    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop ("plot_type must be either 'bar' or 'scatter'")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # Check if order_by is present. If it is additional checks need to be
  # performed to make sure it is a valid input.
  if (!is.null(order_by)) {

    # Make sure order_by is either "Group" or is in f_data.
    if (order_by != "Group" && !(order_by %in% names(omicsData$f_data))) {

      stop (paste("order_by must either be 'Group' or the name of a column",
                  "from f_data.",
                  sep = " "))

    }

    # If order_by is "Group" the group_DF attribute must not be NULL.
    if (order_by == "Group" && is.null(attr(omicsData, "group_DF"))) {

      # Welcome to the pit of despair!! You will never escape!!!
      stop (paste("group_DF must not be NULL. Run the group_designation",
                  "function prior to running the missingval_result function.",
                  sep = " "))

    }

  }

  # Check if color_by is present. If it is additional checks need to be
  # performed to make sure it is a valid input.
  if (!is.null(color_by)) {

    # Make sure color_by is either "Group" or is in f_data.
    if (color_by != "Group" && !(color_by %in% names(omicsData$f_data))) {

      stop (paste("color_by must either be 'Group' or the name of a column",
                  "from f_data.",
                  sep = " "))

    }

    # If color_by is "Group" the group_DF attribute must not be NULL.
    if (color_by == "Group" && is.null(attr(omicsData, "group_DF"))) {

      # Welcome to the pit of despair!! You will never escape!!!
      stop (paste("group_DF must not be NULL. Run the group_designation",
                  "function prior to running the missingval_result function.",
                  sep = " "))

    }

  }

  # Extract info from naRes_obj
  na.by.sample <- naRes_obj$na.by.sample
  na.by.molecule <- naRes_obj$na.by.molecule
  num_missing_vals <- na.by.molecule$num_NA
  edata_cname <- attr(naRes_obj, "cnames")$edata_cname
  fdata_cname <- attr(naRes_obj, "cnames")$fdata_cname

  # Extract info from omicsData
  edata <- omicsData$e_data
  edata_cname_id <- which(names(edata) == edata_cname)
  group_df <- attr(omicsData, "group_DF")

  # Bar plot order_by and group_by crap ---------------

  # Check if order_by is NULL and update the plot_data object accordingly.
  if (!is.null(order_by)) {

    # Reorder the rows of na.by.sample so the bar plot will be displayed in the
    # correct order.
    na.by.sample <- na.by.sample[order(na.by.sample[, order_by]), ]
    na.by.sample[[fdata_cname]] <- factor(
      na.by.sample[[fdata_cname]],
      levels = na.by.sample[[fdata_cname]],
      ordered = TRUE
    )

  }

  # Check if color_by is NULL and update na.by.sample accordingly.
  if (!is.null(color_by)) {

    # Create factors to color by according to the input of color_by.
    color_levels <- if (color_by != "Group")
      unique(factor(omicsData$f_data[[color_by]])) else
        unique(factor(attr(omicsData, "group_DF")[["Group"]]))
    na.by.sample[[color_by]] <- factor(na.by.sample[[color_by]],
                                       levels = color_levels)

  }

  # Fashion astonishing plots --------------------------------------------------

  # Make me a bar plot. As you wish.
  if (plot_type == "bar") {

    p <- na_bar(na.by.sample = na.by.sample, x_lab_bar = x_lab_bar,
                y_lab_bar = y_lab_bar, x_lab_size = x_lab_size,
                y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                title_lab_bar = title_lab_bar, title_lab_size = title_lab_size,
                legend_lab_bar = legend_lab_bar,
                legend_position = legend_position, text_size = text_size,
                bar_width = bar_width, bw_theme = bw_theme, palette = palette,
                display_count = display_count,
                coordinate_flip = coordinate_flip,
                use_VizSampNames = use_VizSampNames, interactive = interactive,
                fdata_cname = fdata_cname, color_by = color_by)

  }

  # Make me a scatter plot. As you wish.
  if (plot_type == "scatter") {

    p <- na_scatter(edata = edata, group_df = group_df,
                    na.by.sample = na.by.sample,
                    num_missing_vals = num_missing_vals,
                    edata_cname = edata_cname, edata_cname_id = edata_cname_id,
                    fdata_cname = fdata_cname, x_lab_scatter = x_lab_scatter,
                    y_lab_scatter = y_lab_scatter,
                    title_lab_scatter = title_lab_scatter,
                    legend_lab_scatter = legend_lab_scatter,
                    legend_position = legend_position,
                    title_lab_size = title_lab_size, x_lab_size = x_lab_size,
                    y_lab_size = y_lab_size, text_size = text_size,
                    bw_theme = bw_theme, palette = palette,
                    x_lab_angle = x_lab_angle,
                    coordinate_flip = coordinate_flip,
                    interactive = interactive, point_size = point_size)

  }

  return (p)

}

na_bar <- function (na.by.sample, x_lab_bar, y_lab_bar, x_lab_size, y_lab_size,
                    x_lab_angle, title_lab_bar, title_lab_size, legend_lab_bar,
                    legend_position, text_size, bar_width, bw_theme, palette,
                    display_count, coordinate_flip, use_VizSampNames,
                    interactive, fdata_cname, color_by) {



  # Farm boy, color the plots based on the input. As you wish.
  if (!is.null(color_by)) {

    # Forge the basic sample bar plot with group info. More details will be
    # added according to the users input later.
    samp <- ggplot2::ggplot(data = na.by.sample,
                            ggplot2::aes(x = .data[[fdata_cname]],
                                         y = num_NA,
                                         fill = !!rlang::sym(color_by))) +
      ggplot2::geom_bar(stat = "identity", width = bar_width)

  } else {

    # Check if palette is NULL or not. Hopefully it isn't so the plot will be
    # created with colors other than the super hideous default ggplot2 colors.
    if (!is.null(palette)) {

      # Create a color from the color brewer package if a palette is provided.
      colas <- RColorBrewer::brewer.pal(5, palette)

    }

    # Create an object for the first default ggplot2 color.
    hideous <- grDevices::hcl(h = 15,
                              c = 100,
                              l = 65)

    # Forge the basic sample bar plot without group info. More details will be
    # added according to the users input later.
    samp <- ggplot2::ggplot(data = na.by.sample,
                            ggplot2::aes(x = .data[[fdata_cname]],
                                         y = num_NA)) +
      ggplot2::geom_bar(stat = "identity",
                        width = bar_width,
                        fill = if (is.null(palette))
                          hideous else
                            colas[[3]])

  }

  # Add the counts to the bar plot and histogram if the user so desires.
  if (display_count) {

    samp <- samp +
      ggplot2::geom_text(ggplot2::aes(label = num_NA),
                         vjust = 2,
                         color = "black",
                         size = text_size)

  }

  # Want the black and white theme? As you wish.
  if (bw_theme) samp <- samp + ggplot2::theme_bw()

  # Make labels for both the plot by sample and plot by molecule.
  xLabelBar <- if (is.null(x_lab_bar))
    "Sample Name" else
      x_lab_bar
  yLabelBar <- if (is.null(y_lab_bar))
    "Number of missing values" else
      y_lab_bar
  plotTitleBar <- if (is.null(title_lab_bar))
    "Missing values by sample" else
      title_lab_bar
  legendLabelBar <- if (is.null(legend_lab_bar))
    "Group" else
      legend_lab_bar

  # Add labels to the plots. Even more tedious work to do. Will it ever end!?!
  samp <- samp +
    ggplot2::xlab(xLabelBar) +
    ggplot2::ylab(yLabelBar) +
    ggplot2::ggtitle(plotTitleBar) +
    ggplot2::labs(fill = legendLabelBar) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Want to use beautiful non-default colors? As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    samp <- samp +
      ggplot2::scale_fill_brewer(palette = palette,
                                 name = legendLabelBar)

  }

  # Want to use short sample names? As you wish.
  if (use_VizSampNames) {

    samp <- samp +
      ggplot2::scale_x_discrete(labels = na.by.sample$VizSampNames)

  }

  # Want to be bass ackwards? As you wish. (This is for the crazies.)
  if (coordinate_flip) samp <- samp + ggplot2::coord_flip()

  # Want an interactive plot? As you wish.
  if (interactive) samp <- plotly::ggplotly(samp)

  return (samp)

}

na_scatter <- function (edata, group_df, na.by.sample, num_missing_vals,
                        edata_cname, edata_cname_id, fdata_cname, x_lab_scatter,
                        y_lab_scatter, title_lab_scatter, legend_lab_scatter,
                        legend_position, title_lab_size, x_lab_size, y_lab_size,
                        text_size, bw_theme, palette, x_lab_angle,
                        coordinate_flip, interactive, point_size) {

  # More tedious label making.
  xLabelScatter <- if (is.null(x_lab_scatter))
    "Mean Intensity" else
      x_lab_scatter
  yLabelScatter <- if (is.null(y_lab_scatter))
    "Number of Missing Values (per molecule)" else
      y_lab_scatter
  plotTitleScatter <- if (is.null(title_lab_scatter))
    "Mean Intensity vs NA per Molecule" else
      title_lab_scatter
  legendLabelScatter <- if (is.null(legend_lab_scatter))
    "Group" else
      legend_lab_scatter

  if (is.null(group_df)) {

    # Calculate the mean intensity for each molecule when the group_DF attribute
    # is NULL.
    mean_intensity <- rowMeans(edata[, -edata_cname_id], na.rm = TRUE)

    plot_data <- as.data.frame(cbind(mean_intensity, num_missing_vals))

    # Start the scatter plot when the group_DF attribute is NULL.
    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(mean_intensity,
                                      num_missing_vals)) +
      ggplot2::geom_point(size = point_size)

  } else {

    # Extract group information to calculate the group-wise mean for each
    # molecule.
    levels <- unique(group_df$Group)
    indices_list <- vector(mode = "list",
                           length = length(levels))

    for (i in 1:length(levels)) {

      # Extract the column indices for the ith group. This will be used to
      # subset e_data.
      inds <- which(group_df$Group == levels[i])
      indices_list[[i]] <- inds

    }

    # Calculate the mean intensity for each molecule by group. NaN can appear if
    # an entire row has all NA values.
    mean_by_group <- lapply(indices_list,
                            function (x, temp_edata) rowMeans(temp_edata[, x],
                                                              na.rm = TRUE),
                            temp_edata = edata[, -edata_cname_id])

    names(mean_by_group) <- levels

    mean_intensity <- do.call(cbind, mean_by_group)

    plot_data <- cbind(num_missing_vals, mean_intensity)
    plot_data <- as.data.frame(plot_data)
    plot_data <- reshape2::melt(plot_data, id.vars = "num_missing_vals")

    # Start the scatter plot when the group_DF attribute is present.
    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(value,
                                      num_missing_vals)) +
      ggplot2::geom_point(ggplot2::aes(color = variable),
                          size = point_size)

  }

  # Add the x, y, and title labels.
  p <- p +
    ggplot2::xlab(xLabelScatter) +
    ggplot2::ylab(yLabelScatter) +
    ggplot2::ggtitle(plotTitleScatter)

  # Want the black and white theme? As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Want to be bass ackwards? As you wish.
  if (coordinate_flip) p <- p + ggplot2::coord_flip()

  # Add all of the theme crap now that the theme_bw crap is out of the way.
  p <- p +
    ggplot2::labs(color = legendLabelScatter) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Want to use beautiful non-default colors? As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_color_brewer(palette = palette,
                                  name = legendLabelScatter)

  }

  # Want an interactive plot? As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

#' plot.corRes
#'
#' For plotting an S3 object of type 'corRes'
#'
#' @param corRes_obj An object of class "corRes" created by calling
#'   \code{cor_result} on an omicsData object.
#' @param omicsData an object of the class 'pepData', 'isobaricpepData',
#'   'proData', 'lipidData', 'metabData', or 'nmrData' usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}},
#'   \code{\link{as.proData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, or\code{\link{as.nmrData}}, respectively.
#' @param x_text Logical. Indicates whether the x-axis will be labeled with the
#'   sample names. The default is TRUE.
#' @param y_text Logical. Indicates whether the y-axis will be labeled with the
#'   sample names. The default is TRUE.
#' @param colorbar_lim A pair of numeric values specifying the minimum and
#'   maximum values to use in the heat map color bar. Defaults to 'c(NA, NA)',
#'   in which case ggplot2 automatically sets the minimum and maximum values
#'   based on the correlation values in the data.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 90.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param color_low A character string specifying the color of the gradient for
#'   low values.
#' @param color_high A character string specifying the color of the gradient for
#'   high values.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-corRes
#'
#' @export
#'
plot.corRes <- function (corRes_obj, omicsData = NULL, colorbar_lim = c(NA, NA),
                         x_text = TRUE, y_text = TRUE, interactive = FALSE,
                         x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                         y_lab_size = 11, x_lab_angle = 90, title_lab = NULL,
                         title_lab_size = 14, legend_lab = NULL,
                         legend_position = "right", color_low = NULL,
                         color_high = NULL, bw_theme = TRUE,
                         use_VizSampNames = FALSE) {

  # Preliminaries --------------------------------------------------------------

  # check for a corRes object #
  if (!inherits(corRes_obj, "corRes")) {

    # You will be given iocane powder due to your carelessness.
    stop ("corRes_obj must be of class 'corRes'")

  }

  if (!all(is.na(colorbar_lim))) {

    if (!is.numeric(colorbar_lim) || length(colorbar_lim) != 2) {

      # I can clearly not choose the arguments in front of you!
      stop ("colorbar_lim must be a numeric vector of length 2")

    }

  }

  # check that omicsData is of appropriate class #
  if(!is.null(omicsData)){

    if (!inherits(omicsData, c("pepData", "proData", "metabData",
                               "lipidData", "nmrData"))) {

      stop ("omicsData is not an appropriate class")

    }

  }

  # Workaround for certain "check" warnings
  Var1 <- Var2 <- value <- NULL

  # Create the data frame that will be used to produce the correlation heatmap.
  corRes_melt <- reshape2::melt(corRes_obj)
  corRes_melt$Var1 <- abbreviate(corRes_melt$Var1, minlength = 20)
  corRes_melt$Var2 <- abbreviate(corRes_melt$Var2, minlength = 20)
  sampleIDx = ordered(corRes_melt$Var2,
                      levels = rev(sort(unique(corRes_melt$Var2))))
  sampleIDy = ordered(corRes_melt$Var1,
                      levels = rev(sort(unique(corRes_melt$Var1))))

  # Create all the plot labels. Life is pain!!!
  xlabel <- if (is.null(x_lab)) "" else x_lab
  ylabel <- if (is.null(y_lab)) "" else y_lab
  legendLabel <- if (is.null(legend_lab)) "Correlation" else legend_lab
  if (is.null(title_lab)) {

    # Determine the plot's title based on whether the data have been normalized.
    if (attributes(corRes_obj)$is_normalized) {

      plotTitle <- "Correlations Among Samples (Normalized Data)"

    } else {

      plotTitle <- "Correlations Among Samples (Un-Normalized Data)"

    }

    # Runs when title_lab is not NULL (the user specified title).
  } else {

    plotTitle <- title_lab

  }

  # Design spectacular plots ---------------------------------------------------

  # Start the skeleton of the heat map. Other aspects are forthcoming.
  hm <- ggplot2::ggplot(corRes_melt,
                        ggplot2::aes(x = sampleIDx,
                                     y = sampleIDy)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value))

  # Farm boy, make me a plot with the black and white theme. As you wish.
  if (bw_theme) hm <- hm + ggplot2::theme_bw()

  # Add plot and axis labels and theme elements.
  hm <- hm +
    ggplot2::ggtitle(plotTitle) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position,
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )


  # Farm boy, make me a plot with a custom correlation range. As you wish.
  if (all(is.na(colorbar_lim))) {

    # Use default range for the color gradient.
    hm <- hm +
      ggplot2::scale_fill_gradient(
        legendLabel,
        low = if (is.null(color_low)) "#132B43" else color_low,
        high = if (is.null(color_high)) "#56B1F7" else color_high,
      )

  } else {

    # Use the user specified range for the color gradient.
    hm <- hm +
      ggplot2::scale_fill_gradient(
        legendLabel,
        limits = colorbar_lim,
        low = if (is.null(color_low)) "#132B43" else color_low,
        high = if (is.null(color_high)) "#56B1F7" else color_high,
      )

  }

  # If there are too many samples don't display the sample names on the plot.
  if (nrow(corRes_obj) > 40) {

    hm <- hm +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank())

  }

  # Farm boy, make me a plot without the x-axis sample labels. As you wish.
  if (!x_text) hm <- hm + ggplot2::theme(axis.text.x = ggplot2::element_blank())

  # Farm boy, make me a plot without the y-axis sample labels. As you wish.
  if (!y_text) hm <- hm + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  # Farm boy, make me a plot with short sample names. As you wish.
  if (use_VizSampNames) {

    if (is.null(omicsData))
      stop ("If using custom sample names omicsData must not be NULL.")

    hm <- hm +
      ggplot2::scale_x_discrete(
        labels = omicsData$f_data$VizSampNames,
        breaks = levels(omicsData$f_data[, get_fdata_cname(omicsData)])
      ) +
      ggplot2::scale_y_discrete(
        labels = omicsData$f_data$VizSampNames,
        breaks = levels(omicsData$f_data[, get_fdata_cname(omicsData)])
      )

  }

  # Farm boy, make me an interactive plot. As you wish.
  if (interactive) hm <- plotly::ggplotly(hm)

  return (hm)

}

#' plot.dimRes
#'
#' For plotting an S3 object of type 'dimRes'
#'
#' @param dimRes_obj An object of class dimRes. A dimRes object is created by
#'   the \code{dim_reduction} function.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 4.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#'
#' @rdname plot-dimRes
#'
#' @export
#'
plot.dimRes <- function (dimRes_obj, interactive = FALSE, x_lab = NULL,
                         y_lab = NULL, x_lab_size = 11, y_lab_size = 11,
                         x_lab_angle = 0, title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         point_size = 4, bw_theme = TRUE, palette = NULL) {

  # Preliminaries --------------------------------------------------------------

  # Evan, make sure the input is the correct class. As you wish.
  if (!inherits(dimRes_obj, "dimRes")) {

    # Welcome to the pit of despair!
    stop ("dimRes_obj must be of class 'dimRes'")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  plotdata <- data.frame(SampleID = dimRes_obj$SampleID,
                         PC1 = dimRes_obj$PC1,
                         PC2 = dimRes_obj$PC2)
  plotdata_name <- names(plotdata)[1]

  # if there is a group designation #
  if(!is.null(attr(dimRes_obj,"group_DF"))) {
    group_DF <- attr(dimRes_obj,"group_DF")
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

      if (!identical(as.character(plotdata[[plotdata_name]]),
                    as.character(group_DF[[fdata_cname]]))) {
        group_DF <- group_DF[match(plotdata[[plotdata_name]],
                                   group_DF[[fdata_cname]]), ]
        plotdata<- merge.data.frame(plotdata,
                                    group_DF,
                                    by.x = plotdata_name,
                                    by.y = fdata_cname,
                                    sort = FALSE)
      } else {

        plotdata<- merge.data.frame(plotdata,
                                    group_DF,
                                    by.x = plotdata_name,
                                    by.y = fdata_cname,
                                    sort = FALSE)

      }

      color_var <- main_eff_names[1]
      pch_var <- main_eff_names[2]

      if(length(legend_lab) > length(main_eff_names)) {

        warning (paste("legend_lab length is greater than the number of main",
                       "effects. Only the first two entries will be used.",
                       sep = " "))

      }

    } else {

      if (!identical(as.character(plotdata[[plotdata_name]]),
                     as.character(group_DF[[fdata_cname]]))) {
        group_DF<- group_DF[match(plotdata[[plotdata_name]],
                                  group_DF[[fdata_cname]]), ]
        plotdata<- merge.data.frame(plotdata,
                                    group_DF,
                                    by.x = plotdata_name,
                                    by.y = fdata_cname,
                                    sort = FALSE)
      } else {

        plotdata<- merge.data.frame(plotdata,
                                    group_DF,
                                    by.x = plotdata_name,
                                    by.y = fdata_cname,
                                    sort = FALSE)

      }

      color_var <- "Group"
      pch_var <- NULL
      display_names <- c("Group", NULL)

      if (length(legend_lab) > 1) {

        warning (paste("legend_lab length is greater than the number of main",
                       "effects. Only the first entry will be used.",
                       sep = " "))

      }

    }

    # Runs when there is no group information.
  } else {

    color_var <- NULL
    pch_var <- NULL
    display_names <- c(NULL, NULL)

    if (!is.null(legend_lab)) {

      warning("There is no group designation, so legend_lab will go unused.")

    }

  }

  # axis labels #
  xr2 <- paste(" = ", round(attr(dimRes_obj, "R2")[1],3), ")", sep = "")
  yr2 <- paste(" = ", round(attr(dimRes_obj, "R2")[2],3), ")", sep = "")
  pc1 <- "PC1 ("
  pc2 <- "PC2 ("

  # custom legend names #
  if(!is.null(legend_lab)) {
    # make the vector at least length 2 to avoid errors in the plot
    display_names[1:length(legend_lab)] <- legend_lab[1:min(2,
                                                            length(legend_lab))]
  }

  # title #
  plot_title <- ifelse(is.null(title_lab), "Principal Components", title_lab)

  # Construct impressive plots -------------------------------------------------

  # Create the bare bones plot.
  p <- ggplot2::ggplot(plotdata,
                       ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(ggplot2::aes_string(col = color_var,
                                            pch = pch_var),
                        size = point_size)

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Create objects for the percent variation for PC1 and PC2
  pc1 <- round(attr(dimRes_obj, "R2")[1], 3)
  pc2 <- round(attr(dimRes_obj, "R2")[2], 3)

  # Tedious label making crap.
  # These two lines don't currently work as intended. They display R^2 without
  # the 2 as a superscript. However, the label cannot be created with substitute
  # because it throws an error when converting to an interactive plot.
  x_lab <- if (is.null(x_lab)) paste0("PC1 (",
                                      expression(R^2),
                                      " = ",
                                      pc1,
                                      ")") else x_lab
  y_lab <- if (is.null(y_lab)) paste0("PC2 (",
                                      expression(R^2),
                                      " = ",
                                      pc2,
                                      ")") else y_lab

  # Add labels and thematic elements.
  p <- p +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::scale_shape_discrete(display_names[2])
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me a plot with dazzling colors. As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_color_brewer(palette = palette,
                                  name = display_names[1])

    # Evan, make me a plot with hideous default ggplot2 colors. As you wish.
  } else {

    p <- p +
      ggplot2::scale_color_discrete(display_names[1])

  }

  # Evan, make me an interactive plot. As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

#' plot.moleculeFilt
#'
#' For plotting an S3 object of type 'moleculeFilt':
#'
#' @param filter_obj An object of class moleculeFilt that contains the molecule
#'   identifier and the number of samples for which the molecule was measured
#'   (not NA).
#' @param min_num An integer specifying the minimum number of samples in which a
#'   biomolecule must appear. If a value is specified, a horizontal line will be
#'   drawn when \code{cumulative=TRUE}, and bars will be colored appropriately
#'   if \code{cumulative=FALSE}.  Defaults to NULL.
#' @param cumulative logical indicating whether the number of biomolecules
#'   observed in \emph{at least} (TRUE) x number of samples or \emph{exactly}
#'   (FALSE) x number of samples should be plotted.  Defaults to TRUE.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size An integer specifying the size of the text (number of
#'   biomolecules by sample) within the bar plot. The default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count Logical. Indicates whether the missing value counts by
#'   sample will be displayed on the bar plot. The default is TRUE.
#'
#' @examples
#' \dontrun{
#' data(pep_object)
#' molfilt <- molecule_filter(pep_object)
#' plot(molfilt, min_num = 5)
#' plot(molfilt, min_num = 3, cumulative = FALSE)
#' }
#'
#' @rdname plot-moleculeFilt
#'
#' @export
#'
plot.moleculeFilt <- function (filter_obj, min_num = NULL, cumulative = TRUE,
                               interactive = FALSE, x_lab = NULL, y_lab = NULL,
                               x_lab_size = 11, y_lab_size = 11,
                               x_lab_angle = 0, title_lab = NULL,
                               title_lab_size = 14, legend_lab = NULL,
                               legend_position = "right", text_size = 3,
                               bar_width = 0.8, bw_theme = TRUE,
                               palette = NULL, display_count = TRUE) {

  # Preliminaries --------------------------------------------------------------

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_obj, "moleculeFilt")) {

    # Fezzik, tear his arms off.
    stop ("filter_obj must be of class moleculeFilt")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # Forge sensational plots ----------------------------------------------------

  # make counts, colors, and plot shape based on value of cumulative

  # values that will be updated depending on input
  fill <- 0
  hline <- NULL

  if(cumulative) {
    # cumulative counts (>=)
    counts <- sapply(1:max(filter_obj$Num_Observations), function(i){
      filter_obj[filter_obj$Num_Observations >= i,]$Num_Observations %>% length()
    })

    # append 1 to extend last step (looks awkward without this)
    counts <- c(counts, counts[length(counts)])
    # 1 appended step
    num_obs <- 1:(max(filter_obj$Num_Observations)+1)

    # shape is a step function with fixed color
    shape <- ggplot2::geom_step(ggplot2::aes(x = num_observations,
                                             y = counts),
                                color = "red")

    # draw a horizontal line if min_num is specified
    hline <- if (!is.null(min_num)) ggplot2::geom_hline(
      ggplot2::aes(color = "black"),
      yintercept = counts[min_num],
      linetype = "dashed"
    ) else NULL

    xlabel <- if (is.null(x_lab)) "Number of samples" else x_lab

    ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)

    plot_title <- if (is.null(title_lab))
      "Count of biomolecules observed in at least X number of samples" else
        title_lab

  } else if (!cumulative) {

    # counts for a specific number of nonmissing biomolecules (==)
    counts <- sapply(1:max(filter_obj$Num_Observations), function (i) {
      filter_obj[filter_obj$Num_Observations == i, ]$Num_Observations %>%
        length()
    })

    # color by which values are kept if min_num is specified
    if (!is.null(min_num)){

      fill <- ifelse(1:max(filter_obj$Num_Observations) >= min_num,
                     "retained",
                     "dropped")
      shape <- ggplot2::geom_bar(ggplot2::aes(x = num_observations,
                                              y = counts,
                                              fill = fill),
                                 stat = "identity",
                                 width = bar_width)

    } else {

      # Check if palette is NULL or not. Hopefully it isn't so the plot will be
      # created with colors other than the super hideous default ggplot2 colors.
      if (!is.null(palette)) {

        # Create a color from the color brewer package if a palette is provided.
        colas <- RColorBrewer::brewer.pal(5, palette)

      }

      # Create an object for the first default ggplot2 color.
      hideous <- grDevices::hcl(h = 15,
                                c = 100,
                                l = 65)

      shape <- ggplot2::geom_bar(ggplot2::aes(x = num_observations,
                                              y = counts),
                                 fill = if (is.null(palette))
                                   hideous else
                                     colas[[3]],
                                 stat = "identity",
                                 width = bar_width)

    }

    num_obs <- 1:max(filter_obj$Num_Observations)

    xlabel <- if (is.null(x_lab)) "Number of samples" else x_lab
    ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)
    plot_title <- ifelse(
      is.null(title_lab),
      "Count of biomolecules observed in exactly X number of samples",
      title_lab
    )

  }

  # create plotting dataframe
  pep_observation_counts <- data.frame(num_observations = num_obs,
                                       frequency_counts = counts,
                                       fill = fill)

  # plot #
  p <- ggplot2::ggplot(pep_observation_counts) +
    shape +
    hline +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous(breaks = max(filter_obj$Num_Observations):1)

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Evan, display the biomolecule counts on the graph. As you wish.
  if (display_count) p <- p + ggplot2::geom_text(
    data = pep_observation_counts[1:max(filter_obj$Num_Observations),],
    ggplot2::aes(x = num_observations,
                 y = frequency_counts,
                 label = frequency_counts),
    vjust = -0.5,
    nudge_x = ifelse(cumulative, 0.5, 0),
    size = text_size
  )

  # Add the theme elements to the plot.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(palette = palette,
                                 name = legend_lab)

  } else {

    p <- p + ggplot2::scale_fill_manual(
      name = legend_lab,
      values = c("dropped" =  "red",
                 "retained" =  "green")
    )

  }

  # Evan, make me an interactive plot. As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

#' plot.imdanovaFilt
#'
#' For plotting an S3 object of type 'imdanovaFilt'
#'
#' @param filter_obj Object of class imdanovaFilt (also a data frame) containing
#'   the molecule identifier and number of samples in each group with
#'   non-missing values for that molecule.
#' @param min_nonmiss_gtest An integer indicating the minimum number of
#'   non-missing feature values allowed per group for \code{gtest_filter}.
#'   Suggested value is 3.
#' @param min_nonmiss_anova An integer indicating the minimum number of
#'   non-missing feature values allowed per group for \code{anova_filter}.
#'   Suggested value is 2.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 3.
#' @param line_size An integer specifying the thickness of the line. The default
#'   is 0.75.
#' @param text_size An integer specifying the size of the text (number of
#'   biomolecules per group). The default is 3.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count Logical. Indicates whether the missing value counts by
#'   sample will be displayed on the bar plot. The default is TRUE.
#'
#' @rdname plot-imdanovaFilt
#'
#' @export
#'
plot.imdanovaFilt <- function (filter_obj, min_nonmiss_anova = NULL,
                               min_nonmiss_gtest = NULL, interactive = FALSE,
                               x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                               y_lab_size = 11, x_lab_angle = 0,
                               title_lab = NULL, title_lab_size = 14,
                               legend_lab = NULL, legend_position = "right",
                               point_size = 3, line_size = 0.75, text_size = 3,
                               bw_theme = TRUE, palette = NULL,
                               display_count = TRUE) {

  # Preliminaries --------------------------------------------------------------

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_obj, "imdanovaFilt")) {

    # You used the wrong input for filter_obj. Savvy?
    stop ("filter_obj must be of class imdanovaFilt")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # Nab group size minimum that is not a singleton
  group_sizes <- attr(filter_obj, "group_sizes")$n_group
  group_sizes_valid <- group_sizes[group_sizes > 1]
  max_x <- min(group_sizes_valid)
  obs <- as.data.frame(filter_obj)[-1]

  n_biomolecules_anova <- purrr::map_int(0:max_x, function (num) {

    # Number of groups greater than num
    n_greater_obs <- apply(obs, 1, function (row) sum(!(row < num)))
    # Must be in at least two groups for ANOVA
    sum(n_greater_obs > 1)

  })

  n_biomolecules_gtest <- purrr::map_int(0:max_x, function (num) {

    # biomolecules with 1+ group > num
    n_greater_obs <- sum(apply(obs, 1, function (row) any(!(row < num))))

  })

  # ANOVA table
  plotter1 <- data.frame(
    Min_obs = c(0:max_x),
    Count_biomolecules = c(n_biomolecules_anova),
    Statistic = c(rep("Within 2+ groups (ANOVA)", length(0:max_x)))
  )

  # G-test table
  plotter2 <- data.frame(
    Min_obs = c(0:max_x),
    Count_biomolecules = c(n_biomolecules_gtest),
    Statistic = c(rep("Within 1+ groups (G-Test)", length(0:max_x)))
  )

  # Evan, get me the colors of my choosing. As you wish.
  if (!is.null(palette)) colas <- RColorBrewer::brewer.pal(5, palette)

  # More monotonous label making code.
  xlabel <- if (is.null(x_lab)) "Number of biomolecules" else x_lab
  ylabel <- if (is.null(y_lab))
    "Minimum number of observations per group" else
        y_lab
  titleLabel <- if (is.null(title_lab)) "IMD-ANOVA filter" else title_lab

  # Fabricate glorious plots ---------------------------------------------------

  # Construct the ggplot skeleton.
  p <- ggplot2::ggplot() +
    ggplot2::xlim(0, ceiling(max(n_biomolecules_gtest) * 1.25)) +
    ggplot2::ylim(0, max_x)

  # Evan, make the background black and white. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Evan, add title and axis labels. As you wish.
  p <- p +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(titleLabel)

  # Evan, add the gtest points to the plot. As you wish.
  p <- p +
    ggplot2::geom_point(data = plotter2,
                        ggplot2::aes(x = Count_biomolecules,
                                     y = Min_obs,
                                     color = Statistic),
                        size = point_size)

  # Evan, add the anova points to the plot. As you wish.
  p <- p +
    ggplot2::geom_point(data = plotter1,
                        ggplot2::aes(x = Count_biomolecules,
                                     y = Min_obs,
                                     color = Statistic),
                        size = point_size)

  # Evan, display the counts on the plot. As you wish.
  if (display_count) p <- p +
    ggplot2::geom_text(
      data = plotter1,
      ggplot2::aes(x = Count_biomolecules,
                   y = Min_obs,
                   label = Count_biomolecules),
      size = text_size,
      hjust = -0.5
    ) +
    ggplot2::geom_text(
      data = plotter2,
      ggplot2::aes(x = Count_biomolecules,
                   y = Min_obs,
                   label = Count_biomolecules),
      size = text_size,
      hjust = -0.5
    )

  # Evan, add gtest info to the plot. As you wish.
  if (!is.null(min_nonmiss_gtest)) {

    p <- p +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = n_biomolecules_gtest[min_nonmiss_gtest + 1],
                     color = "G-test applied filter"),
        linetype = "dashed",
        size = if (is.null(line_size)) 1 else line_size
      )

  }

  # Evan, add anova info to the plot. As you wish.
  if (!is.null(min_nonmiss_anova)) {

    p <- p +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = n_biomolecules_anova[min_nonmiss_anova + 1],
                     color = "ANOVA applied filter"),
        linetype = "dashed",
        size = if (is.null(line_size)) 1 else line_size
      )

  }

  # Evan, add a custom legend according to anova and gtest inputs. As you wish.
  # A minimum for gtest IS supplied and a minimum for anova IS NOT supplied.
  if (!is.null(min_nonmiss_gtest) && is.null(min_nonmiss_anova)) {

    # Add a customized legend when min_nonmiss_gtest is not NULL and
    # min_nonmiss_anova is NULL.
    p <- p +
      ggplot2::scale_color_manual(
        name = if (is.null(legend_lab)) "" else legend_lab,
        values = c(
          `Within 1+ groups (G-Test)` = if (is.null(palette))
            "#FFC107" else
              colas[[2]],
          `G-test applied filter` = if (is.null(palette))
            "#FFC107" else
              colas[[2]],
          `Within 2+ groups (ANOVA)` = if (is.null(palette))
            "#004D40" else
              colas[[3]]
        ),
        breaks = c(
          plotter2$Statistic[[1]], "G-test applied filter",
          plotter1$Statistic[[1]], "ANOVA applied filter"
        ),
        guide = ggplot2::guide_legend(
          override.aes = list(
            linetype = c(0, 2, 0),
            shape = c(16, NA, 16),
            color = c(
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#004D40" else colas[[3]]
            )
          )
        )
      )

    # A minimum for gtest IS NOT supplied and a minimum for anova IS supplied.
  } else if (is.null(min_nonmiss_gtest) && !is.null(min_nonmiss_anova)) {

    # Add a customized legend when min_nonmiss_gtest is NULL and
    # min_nonmiss_anova is not NULL.
    p <- p +
      ggplot2::scale_color_manual(
        name = if (is.null(legend_lab)) "" else legend_lab,
        values = c(
          `Within 1+ groups (G-Test)` = if (is.null(palette))
            "#FFC107" else
              colas[[2]],
          `Within 2+ groups (ANOVA)` = if (is.null(palette))
            "#004D40" else
              colas[[3]],
          `ANOVA applied filter` = if (is.null(palette))
            "#004D40" else
              colas[[3]]
        ),
        breaks = c(
          plotter2$Statistic[[1]], "G-test applied filter",
          plotter1$Statistic[[1]], "ANOVA applied filter"
        ),
        guide = ggplot2::guide_legend(
          override.aes = list(
            linetype = c(0, 0, 2),
            shape = c(16, 16, NA),
            color = c(
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#004D40" else colas[[3]],
              if (is.null(palette)) "#004D40" else colas[[3]]
            )
          )
        )
      )

    # A minimum for gtest AND anova IS supplied.
  } else if (!is.null(min_nonmiss_gtest) && !is.null(min_nonmiss_anova)) {

    # Add a customized legend when both min_nonmiss_gtest and min_nonmiss_anova
    # are not NULL.
    p <- p +
      ggplot2::scale_color_manual(
        name = if (is.null(legend_lab)) "" else legend_lab,
        values = c(
          if (is.null(palette)) "#FFC107" else colas[[2]],
          if (is.null(palette)) "#FFC107" else colas[[2]],
          if (is.null(palette)) "#004D40" else colas[[3]],
          if (is.null(palette)) "#004D40" else colas[[3]]
        ),
        breaks = c(
          plotter2$Statistic[[1]], "G-test applied filter",
          plotter1$Statistic[[1]], "ANOVA applied filter"
        ),
        guide = ggplot2::guide_legend(
          override.aes = list(
            linetype = c(0, 2, 0, 2),
            shape = c(16, NA, 16, NA),
            color = c(
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#004D40" else colas[[3]],
              if (is.null(palette)) "#004D40" else colas[[3]]
            )
          )
        )
      )

    # Neither a minimum for gtest nor anova is supplied.
  } else {

    # Add a customized legend when both min_nonmiss_gtest and min_nonmiss_anova
    # are NULL.
    p <- p +
      ggplot2::scale_color_manual(
        name = if (is.null(legend_lab)) "" else legend_lab,
        values = c(
          if (is.null(palette)) "#FFC107" else colas[[2]],
          if (is.null(palette)) "#004D40" else colas[[3]]
        ),
        guide = ggplot2::guide_legend(
          override.aes = list(
            shape = c(16, 16),
            color = c(
              if (is.null(palette)) "#FFC107" else colas[[2]],
              if (is.null(palette)) "#004D40" else colas[[3]]
            )
          )
        )
      )

  }

  # Evan, add thematic elements to the plot. As you wish.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me an interactive plot. As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  return (p)

}

#' plot.proteomicsFilt
#'
#' For plotting an S3 object of type 'proteomicsFilt':
#'
#' @param filter_obj Object of class proteomicsFilt which is a list containing
#'   two elements. The first element is a data frame with the counts of proteins
#'   mapping to each peptide. The second element is also a data frame with the
#'   counts of peptides mapping to each protein.
#' @param min_num_peps an optional integer value between 1 and the maximum
#'   number of peptides that map to a protein in the data. The value specifies
#'   the minimum number of peptides that must map to a protein. Any protein with
#'   less than \code{min_num_peps} mapping to it will be returned as a protein
#'   that should be filtered. Default value is NULL.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab_pep A character string used for the x-axis label for the
#'   peptide-to-protein plot. The default is NULL in which case the default
#'   x-axis label will be used.
#' @param x_lab_pro A character string used for the x-axis label for the
#'   protein-to-peptide plot. The default is NULL in which case the default
#'   x-axis label will be used.
#' @param y_lab_pep A character string used for the y-axis label for the
#'   peptide-to-protein plot. The default is NULL in which case the default
#'   y-axis label will be used.
#' @param y_lab_pro A character string used for the y-axis label for the
#'   protein-to-peptide plot. The default is NULL in which case the default
#'   y-axis label will be used.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab_pep A character string specifying the peptide-to-protein
#'   plot title. The default is NULL in which case the default title will be
#'   used.
#' @param title_lab_pro A character string specifying the protein-to-peptide
#'   plot title. The default is NULL in which case the default title will be
#'   used.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size An integer specifying the size of the text (number of
#'   peptides or proteins depending on the plot) within the bar plot. The
#'   default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count Logical. Indicates whether the peptide or protein counts
#'   will be displayed on the bar plot. The default is TRUE.
#'
#' @examples
#' \dontrun{
#' data(pep_object)
#' profilt <- proteomics_filter(pep_object)
#' plot(profilt, min_num_peps = 5)
#' }
#'
#' @rdname plot-proteomicsFilt
#'
#' @export
#'
plot.proteomicsFilt <- function (filter_obj, min_num_peps = NULL,
                                 interactive = FALSE, x_lab_pep = NULL,
                                 x_lab_pro = NULL, y_lab_pep = NULL,
                                 y_lab_pro = NULL, x_lab_size = 11,
                                 y_lab_size = 11, x_lab_angle = 0,
                                 title_lab_pep = NULL, title_lab_pro = NULL,
                                 title_lab_size = 14, legend_lab = NULL,
                                 legend_position = "right", text_size = 3,
                                 bar_width = 0.8, bw_theme = TRUE,
                                 palette = NULL, display_count = TRUE) {

  # Preliminaries --------------------------------------------------------------

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_obj, "proteomicsFilt")) {

    # The filter object is the wrong class. Savvy?
    stop ("filter_obj must be of class proteomicsFilt")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # Error Checks
  if (!is.null(min_num_peps)) {
    # check that min_num_peps is numeric and >=1 #
    if (!inherits(min_num_peps, "numeric") || min_num_peps < 1)
      stop ("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if (min_num_peps %% 1 != 0)
      stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if (length(min_num_peps) != 1)
      stop ("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if (min_num_peps > max(filter_obj$counts_by_pro$n))
      stop (paste("min_num_peps cannot be greater than the maximum number of",
                  "peptides that map to a protein.",
                  sep = " "))
  }

  # Seize unique values for peptide to protein and protein to peptide counts:
  # These bins represent the number of PROTEINS each peptide maps to. Unless
  # there are degenerate peptides this will be a vector with one value: 1.
  pep_bins <- sort(unique(filter_obj$counts_by_pep$n))
  # These bins represent the number of PEPTIDES each protein maps to.
  pro_bins <- sort(unique(filter_obj$counts_by_pro$n))

  # get counts of peptides that are mapped to by EXACTLY the number of
  # proteins given in pep_bins
  pep_counts <- sapply(pep_bins, function(x){
    filter_obj$counts_by_pep[filter_obj$counts_by_pep$n == x,] %>% nrow()
  })

  # get counts of proteins that are mapped to by EXACTLY the number of
  # peptides given in pro_bins
  pro_counts <- sapply(pro_bins, function(x){
    filter_obj$counts_by_pro[filter_obj$counts_by_pro$n == x, ] %>% nrow()
  })

  # fill value and fill labels that will be given real values if certain
  # arguments are supplied.
  fill <- 0
  fill_format <- NULL

  # Mind numbing label making bit.
  xlabel_pep <- if (is.null(x_lab_pep))
    "Number of Peptides" else
      x_lab_pep
  ylabel_pep <- if (is.null(y_lab_pep)) "Count of Proteins" else y_lab_pep
  titleLabelPep <- if (is.null(title_lab_pep))
    "Y peptides mapped to by exactly X proteins" else
      title_lab_pep
  xlabel_pro <- if (is.null(x_lab_pro))
    "Number of Peptides" else
      x_lab_pro
  ylabel_pro <- if (is.null(y_lab_pro)) "Count of Proteins" else y_lab_pro
  titleLabelPro <- if (is.null(title_lab_pro))
    "Y proteins mapped to by exactly X peptides" else
      title_lab_pro

  # create plotting dataframe
  pro_counts_df <- data.frame(counts = pro_counts,
                              bins = pro_bins)
  pep_counts_df <- data.frame(counts = pep_counts,
                              bins = pep_bins)

  # Manufacture phenomenal plots -----------------------------------------------

  #!#!#!#! p represents the protein plot #!#!#!#!
  #!#!#!#! q represents the pepe plot #!#!#!#!

  # Create the bare bones protein and peptide plots.
  p <- ggplot2::ggplot(pro_counts_df)
  q <- ggplot2::ggplot(pep_counts_df)

  # if min_num_peps is specified, add a coloring variable that is red for
  # dropped values and green for retained values
  if (!is.null(min_num_peps)) {

    fill <- ifelse(pro_bins >= min_num_peps, "retained", "dropped")

    # Add the fill column in pro_counts_df to reflect the retained and
    # dropped proteins.
    pro_counts_df$fill <- fill

    # We change bins to a factor so the x-axis tick labels are the value in bins
    # but the width of the bars and x-axis does not change according to the
    # numeric value of bins.
    p <- p +
      ggplot2::geom_bar(
        ggplot2::aes(x = as.factor(bins),
                     y = counts,
                     fill = fill),
        stat = "identity",
        width = bar_width
      )

    # Creates bar charts when there is only one group (all bars will have the
    # same color). These bar charts are created when min_num_peps is NOT
    # specified. Color palettes must be defined here because there are no groups
    # to fill by later in the script (when color palettes are defined for bar
    # charts with groups).
  } else {

    # Check if palette is NULL or not. Hopefully it isn't so the plot will be
    # created with colors other than the super hideous default ggplot2 colors.
    if (!is.null(palette)) {

      # Create a color from the color brewer package if a palette is provided.
      colas <- RColorBrewer::brewer.pal(5, palette)

    }

    # Create an object for the first default ggplot2 color.
    hideous <- grDevices::hcl(h = 15,
                              c = 100,
                              l = 65)

    # We change bins to a factor so the x-axis tick labels are the value in bins
    # but the width of the bars and x-axis does not change according to the
    # numeric value of bins.
    p <- p +
      ggplot2::geom_bar(
        ggplot2::aes(x = as.factor(bins),
                     y = counts),
        fill = if (is.null(palette))
          hideous else
            colas[[3]],
        stat = "identity",
        width = bar_width
      )
    q <- q +
      ggplot2::geom_bar(
        ggplot2::aes(x = as.factor(bins),
                     y = counts),
        fill = if (is.null(palette))
          hideous else
            colas[[3]],
        stat = "identity",
        width = bar_width
      )

  }

  # Evan, add plot labels for me. As you wish.
  p <- p +
    ggplot2::xlab(xlabel_pro) +
    ggplot2::ylab(ylabel_pro) +
    ggplot2::ggtitle(titleLabelPro)
  q <- q +
    ggplot2::xlab(xlabel_pep) +
    ggplot2::ylab(ylabel_pep) +
    ggplot2::ggtitle(titleLabelPep)

  # Evan, make the plot theme black and white. As you wish.
  if (bw_theme) {

    p <- p + ggplot2::theme_bw()
    q <- q + ggplot2::theme_bw()

  }

  # Evan, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) {

    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(palette = palette,
                                 name = legend_lab)

  } else {

    p <- p + ggplot2::scale_fill_manual(name = legend_lab,
                                        values = c("dropped" = "red",
                                                   "retained" = "green"))

  }

  # Evan, add the protein counts to the plot. As you wish.
  if (display_count) p <- p + ggplot2::geom_text(
    data = pro_counts_df,
    ggplot2::aes(x = as.factor(bins),
                 y = counts,
                 label = counts),
    vjust = -0.5,
    size = text_size
  )

  # Create generic theme settings.
  axs <- ggplot2::theme(
    plot.title = ggplot2::element_text(size = title_lab_size),
    axis.title.x = ggplot2::element_text(size = x_lab_size),
    axis.title.y = ggplot2::element_text(size = y_lab_size),
    axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
    legend.position = legend_position
  )

  # Evan, add thematic elements to the plot. As you wish.
  p <- p + axs
  q <- q + axs

  # Return both peptide and protein plots if there are degenerate peptides.
  if (length(pep_bins) > 1) {

    # Evan, make me an interactive plot. As you wish.
    if (interactive) {

      p <- plotly::ggplotly(p)
      q <- plotly::ggplotly(q)

      plotly::subplot(p, q, nrows = 1)

    } else {

      p + q


    }

  } else {

    # Evan, make me an interactive plot. As you wish.
    if (interactive) p <- plotly::ggplotly(p)

    return (p)

  }

}

#' plot.rmdFilt
#'
#' For plotting an S3 object of type 'rmdFilt'
#'
#' @param filter_obj An object of class rmdFilt which is a data frame containing
#'   columns for the sample ID, group, log2 robust Mahalanobis distance,
#'   p-values, and metrics used to calcualte the robust Mahalanobis distance.
#' @param pvalue_threshold A threshold for the Robust Mahalanobis Distance (RMD)
#'   p-value. If \code{sampleID} is NULL (see \code{sampleID} below), a
#'   horizontal line is plotted at the RMD value that corresponds with the
#'   threshold, and all samples above the line have a p-value below the
#'   threshold. If \code{sampleID} is not NULL, \code{pvalue_threshold} will do
#'   nothing. Default value is NULL.
#' @param sampleID A character vector specifying the sample names to be plotted.
#'   If specified, the plot function produces a boxplot instead of a
#'   scatterplot. A point, colored by sample, will be placed on each boxplot for
#'   that sample's value for the given metric. The default value is NULL.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 90.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 3.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-rmdFilt
#'
#' @export
#'
plot.rmdFilt <- function (filter_obj, pvalue_threshold = NULL, sampleID = NULL,
                          interactive = FALSE, x_lab = NULL, y_lab = NULL,
                          x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",

                          point_size = 3, bw_theme = TRUE, palette = NULL,
                          use_VizSampNames = FALSE) {


  # Preliminaries --------------------------------------------------------------

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_obj, "rmdFilt")) {

    # The filter object is the wrong class. Savvy?
    stop ("filter_obj must be of class rmdFilt")

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # Have a looksie at the sampleID argument.
  if (!is.null(sampleID)) {

    # Make sure the sample IDs provided actually exist.
    if (!all(sampleID %in% attr(filter_obj, "sample_names")))
      stop ("The sample IDs provided do not match the sample IDs in the data.")

  }

  samp_id <- names(attr(filter_obj, "group_DF"))[1]
  metrics <- attributes(filter_obj)$metrics

  # determine how to melt based on the number of main effects
  group_df <- attributes(filter_obj)$group_DF

  # determine the main effects, then melt the df #
  if (ncol(group_df) == 2) {

    main_eff_names <- "Group"
    dfmelt <- reshape2::melt(filter_obj, id = c(samp_id, main_eff_names))

  } else if (ncol(group_df) > 2) {

    ## put main effect with more levels first, for plotting aesthetics ##
    temp <- droplevels(group_df[,-(1:2)])
    numlevels <- sapply(1:2, function (j) length(levels(as.factor(temp[, j]))))
    main_eff_names <- names(temp)[order(numlevels, decreasing = TRUE)]
    dfmelt <- reshape2::melt(filter_obj[,-2], id = c(samp_id, main_eff_names))

  }

  # Data frame that has information to create rmd box plots.
  dfsub <- dfmelt[dfmelt$variable %in% metrics, ]

  # legend labels #
  ## function for shortening long main effect names ##
  abbrev_fun <- function(string) {
    if(nchar(string)>25) string=paste0(substr(string,1,23),"...")
    return(string)
  }

  # Shorten the names of the main effects (if the are over 25 characters). This
  # vector will have more than one element if there is more than one main
  # effect.
  display_names <- sapply(main_eff_names, abbrev_fun)

  # Create the title for the legend that describes differences by color.
  legend_title_color <- ifelse(is.null(legend_lab),
                               display_names[1],
                               legend_lab[1])

  # Initialize the title for the legend that describes differences by point
  # shape to NULL. It will only be used if there is more than one main effect.
  legend_title_shape <- NULL

  # Use the legend label (if it is provided in the input) for the legend
  # describing the differences by point shape. The user must specify a legend
  # label with more than one element. It is not created within the function.
  if(length(display_names) > 1) {
    if(length(legend_lab) > 1) {
      legend_title_shape <- legend_lab[2]
    } else {
      legend_title_shape <- display_names[2]
    }
  }

  # Assemble captivating plots -------------------------------------------------

  # Create plot skeleton that will be filled in depending on the input to the
  # pvalue_threshold argument.

  # Beautiful box plots ---------------

  # Check is sampleID is NULL. Generate box plots if it is not NULL.
  if (!is.null(sampleID)) {

    levels(dfsub$variable)[levels(dfsub$variable) == "Fraction_Missing"] <-
      "Prop_missing"

    plot_title <- ifelse(is.null(title_lab),
                         "Summary of RMD metrics used",
                         title_lab)
    xlabel <- ifelse(is.null(x_lab), " ", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Value", y_lab)

    p <- ggplot2::ggplot(dfsub) +
      ggplot2::geom_boxplot(ggplot2::aes(x = rep(1,length(value)),
                                         y = value)) +
      ggplot2::facet_wrap(~ variable,
                          scales = "free",
                          ncol = length(metrics)) +
      ggplot2::geom_point(
        data = dfsub[dfsub[, samp_id] %in% sampleID, ],
        ggplot2::aes(x = rep(1, length(value)),
                     y = value,
                     color = !!rlang::sym(samp_id)),
        size = point_size
      ) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title)

    # Stunning scatter plots: no threshold ---------------

  } else if (is.null(pvalue_threshold)) {

    # Start scatter plot skeleton when there is no p-value threshold.
    p <- ggplot2::ggplot(filter_obj)

    # Start plot when there is only one main effect.
    if (length(main_eff_names) == 1) {

      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = forcats::fct_inorder(!!rlang::sym(samp_id)),
                       y = Log2.md,
                       color = !!rlang::sym(main_eff_names[1])),
          size = point_size
        )


      # Start plot when there are two main effects.
    } else {

      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = forcats::fct_inorder(!!rlang::sym(samp_id)),
                       y = Log2.md,
                       color = !!rlang::sym(main_eff_names[1]),
                       shape = !!rlang::sym(main_eff_names[2])),
          size = point_size
        )


    }

    plot_title <- ifelse(is.null(title_lab),
                         "Sample Outlier Results",
                         title_lab)
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)

    p <- p  +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::labs(color = legend_title_color,
                    shape = legend_title_shape)

    # Stunning scatter plots: with threshold ---------------

  } else {

    # get y-intercept for line
    df <- attributes(filter_obj)$df
    yint <- log(qchisq(1 - pvalue_threshold, df = df), base = 2)

    # make title
    if(is.null(title_lab)) {
      plot_title <- "Sample Outlier Results"
      subtitle <- paste("p-value threshold = ", pvalue_threshold)
    } else {
      plot_title <- title_lab
      subtitle <- ggplot2::waiver()
    }
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)

    # Start scatter plot skeleton when a p-value threshold is specified.
    p <- ggplot2::ggplot(filter_obj)

    # Start scatter plot skeleton when a p-value threshold has been provided.
    if (length(main_eff_names) == 1) {

      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = forcats::fct_inorder(!!rlang::sym(samp_id)),
                       y = Log2.md,
                       col = !!rlang::sym(main_eff_names)),
          alpha = ifelse(filter_obj$pvalue < pvalue_threshold, 1, 0.5),
          size = point_size
        )

    } else {

      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = forcats::fct_inorder(!!rlang::sym(samp_id)),
                       y = Log2.md,
                       col = !!rlang::sym(main_eff_names[1]),
                       shape = !!rlang::sym(main_eff_names[2])),
          alpha = ifelse(filter_obj$pvalue < pvalue_threshold, 1, 0.5),
          size = point_size
        )

    }

    # Add title, axis labels, and other crap.
    p <- p +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title,
                       subtitle = subtitle) +
      ggplot2::geom_hline(yintercept = yint) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 1),
                      shape = ggplot2::guide_legend(ncol = 1)) +
      ggplot2::labs(color = legend_title_color,
                    shape = legend_title_shape)

  }

  # Farm boy, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) p <- p +
    ggplot2::scale_color_brewer(name = legend_title_color,
                                palette = palette)

  # Farm boy, create thematic elements for the plot. As you wish.
  mytheme <- ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
    plot.title = ggplot2::element_text(size = title_lab_size),
    plot.subtitle = ggplot2::element_text(face = 'italic'),
    axis.title.x = ggplot2::element_text(size = x_lab_size),
    axis.title.y = ggplot2::element_text(size = y_lab_size),
    legend.position = legend_position
  )

  # Farm boy, remove useless box plot x-axis labels. As you wish.
  if (!is.null(sampleID)) {

    mytheme <- mytheme +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())

  }

  # Farm boy, add the black and white theme to the plot. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

  # Farm boy, add the thematic elements you created to the plot. As you wish.
  p <- p + mytheme

  # Farm boy, use my custom sample names in the plot. As you wish.
  if (use_VizSampNames) {

    # Change the sample names of the scatter plot.
    p <- p +
      ggplot2::scale_x_discrete(labels = attr(filter_obj, "VizSampNames"))

    # We only want to change the legend if the box plots are created.
    if (!is.null(sampleID)) {

      # Nab the indices of the samples that will be highlighted in the box
      # plots.
      idx <- which(attr(filter_obj, "sample_names") %in% sampleID)

      # Change the names in the legend of the box plots.
      p <- p +
        ggplot2::scale_color_hue(
          labels = attr(filter_obj, "VizSampNames")[idx]
        )

    }

  }


  # Farm boy, make the plot interactive. As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  # Farm boy, return the plot so the entire world can enjoy it. As you wish.
  return (p)

}

#' plot.cvFilt
#'
#' For plotting an S3 object of type 'cvFilt'
#'
#' @param filter_obj An object of class cvFilt which is a data frame containing
#'   a column of biomolecule IDs and a column with the coefficient of variation
#'   (CV) for each biomolecule. If the group_designation function was run
#'   previously the CV column will contain pooled CV values.
#' @param cv_threshold draws a vertical line at the pooled cv cutoff, values to
#'   the left of this cutoff are dropped if the filter is applied.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param log_scale Logical. Indicates whether to use a log2 transformed x-axis.
#'   The default is TRUE.
#' @param n_bins An integer value specifying the number of bins to draw in the
#'   histogram.  The default is 30.
#' @param n_breaks An integer value specifying the number of breaks to use. You
#'   may get less breaks if rounding causes certain values to become non-unique.
#'   The default is 15.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#'
#'
#' @examples
#' \dontrun{
#' data(pep_object)
#'
#' pep_object <- group_designation(omicsData = pep_object,
#'                                 main_effects = "Condition")
#'
#' cvfilt <- cv_filter(pep_object)
#'
#' plot(cvfilt, cv_threshold = 20)
#' plot(cvfilt, cv_threshold = 10, log_scale = FALSE)
#' }
#'
#' @rdname plot-cvFilt
#'
#' @export
#'
plot.cvFilt <- function (filter_obj, cv_threshold = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 0,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         log_scale = TRUE, n_breaks = 15, n_bins = 30,
                         bw_theme = TRUE, palette = NULL) {

  # Preliminaries --------------------------------------------------------------

  # checks for cv_threshold if not null
  if(!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if(!is.numeric(cv_threshold))
      stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if(length(cv_threshold)>1)
      stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if(cv_threshold <= 1 || cv_threshold >= max(filter_obj$CV, na.rm = TRUE))
      stop("cv_threshold must be greater than 1 and less than the max CV value")
  }

  if(!is.logical(log_scale)) stop("log_scale must be logical: TRUE or FALSE")
  if(!(n_breaks%%1 == 0)) stop("n_breaks must be integer valued")
  if(!(n_bins%%1 == 0)) stop("n_bins must be integer valued")

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {

    if (!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu",
                         "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges",
                         "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues",
                         "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
                         "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
                         "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) {

      # INCONCEIVABLE!!!
      stop ("palette must be an RColorBrewer palette")

    }

  }

  # plotting object
  new_object <- filter_obj[!is.na(filter_obj$CV), ]
  max_x_val <- attributes(filter_obj)$max_x_val

  # labels
  plot_title <- ifelse(is.null(title_lab),
                       "Coefficient of Variation (CV)",
                       title_lab)
  ylabel <- ifelse(is.null(y_lab), "Count", y_lab)

  # Label x-axis according to whether group_designation() was run.
  if (is.null(x_lab)) {
    xlabel <- if (attributes(filter_obj)$pooled) "Pooled CV" else "CV"
  }

  ### Store ggplot2 layers

  # scale transform and breaks depending on x-axis scale
  if (log_scale) {
    trans <- "log2"
    i <- log2(min(new_object$CV, na.rm = TRUE))
    breaks <- 0

    # define a step value that is evenly spaced in the log2 scale
    step = (max(log2(new_object$CV), na.rm = TRUE) -
              min(log2(new_object$CV), na.rm = TRUE)) / n_breaks

    # create the normal scale labels that will be log2 transformed when passed
    # to ggplot
    while(2^i < max(new_object$CV, na.rm = TRUE)){
      breaks <- c(breaks, 2^i)
      i <- i+step
    }
    # rounding for plot purposes
    breaks <- round(breaks, 2)
  } else {
    breaks <- scales::pretty_breaks(n = n_breaks)
    trans <- "identity"
  }

  # change default title and draw a vertical line if cv_thresh specified
  if(!is.null(cv_threshold)) {

    if(is.null(title_lab)) {
      plot_title <- bquote(paste("Coefficient of Variation (CV):  ",
                                 italic(paste("CV Threshold = ",
                                              .(cv_threshold))),
                                 ""))
    }

    cutoff <- ggplot2::geom_vline(xintercept = cv_threshold)

  } else {
    cutoff <- NULL
  }

  # Make magnificent plots -----------------------------------------------------

  # Farm boy, start the plot for me. As you wish.
  p <- ggplot2::ggplot(new_object)

  # Check if palette is NULL or not. Hopefully it isn't so the plot will be
  # created with colors other than the super hideous default ggplot2 colors.
  if (!is.null(palette)) {

    # Farm boy, make colors to replace the default color. As you wish.
    colas <- RColorBrewer::brewer.pal(5, palette)

    # Farm boy, add the beautiful color to the histogram. As you wish.
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(x = CV),
                              bins = n_bins,
                              fill = colas[[3]])

    # Runs when palette is NULL.
  } else {

    # Farm boy, make the histogram with default colors. As you wish.
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(x = CV),
                              bins = n_bins)

  }

  # Farm boy, add the cutoff and make the histogram pretty. As you wish.
  p <- p +
    cutoff +
    ggplot2::scale_x_continuous(breaks = breaks,
                                trans = trans) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

  # Farm boy, add labels to the plot. As you wish.
  p <- p +
    ggplot2::ggtitle(plot_title) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel)

  # Farm boy, make the plot black and white. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Farm boy, add thematic elements to the plot. As you wish.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1)
    )

  # Farm boy, send my plot into the world. As you wish.
  return (p)

}

#' plot.customFilt
#'
#' Currently plotting a customFilt object is not supported
#'
#' @param filter_obj An object of class customFilt.
#'
#' @rdname plot-customFilt
#'
#' @export
#'
plot.customFilt <- function (filter_obj) {

  message (paste("There is no plot method for objects of class 'customFilt'.",
                 "See summary.customFilt instead.",
                 sep = " "))

}

#' plot.normRes
#'
#' For plotting an S3 object of type 'normRes'
#'
#' @param normRes A normRes object created by the normalize_global function.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-normRes
#'
#' @export
#'
plot.normRes <- function (normRes_obj, order_by = NULL, color_by = NULL,
                          facet_by = NULL, facet_cols = NULL,
                          interactive = FALSE, x_lab = NULL, y_lab = NULL,
                          x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",
                          ylimit = NULL, bw_theme = TRUE, palette = NULL,
                          use_VizSampNames = FALSE) {

  # Preliminaries --------------------------------------------------------------

  # Create the raw and normalized omicsData objects.
  omicsRaw <- attr(normRes_obj, "omicsData")

  # e_data will be updated after applying the normalization.
  omicsNorm <- attr(normRes_obj, "omicsData")

  # Extricate the ID column index.
  id_col <- which(names(omicsRaw$e_data) == get_edata_cname(omicsRaw))

  # Apply normalization depending on the location, scale, and backtransform
  # values.
  if (is.null(normRes_obj$parameters$normalization$scale)) {

    # Normalize with just the location parameter (no backtransform). There isn't
    # a section where backtransform is TRUE because the normalize_global
    # function will not calculate backtransform parameters if the data are not
    # normalized.
    omicsNorm$e_data[, -id_col] <- (
      omicsNorm$e_data[, -id_col] -
        rep(normRes_obj$parameters$normalization$location,
            each = nrow(omicsNorm$e_data))
    )

    # Backtransmute if there is a location backtransform parameter.
    if (!is.null(normRes_obj$parameters$backtransform$location)) {

      omicsNorm$e_data[, -id_col] <- (
        omicsNorm$e_data[, -id_col] +
          normRes_obj$parameters$backtransform$location
      )

    }

  } else {

    # Normalize with the location and scale parameters (no backtransform). There
    # isn't a section where backtransform is TRUE because the normalize_global
    # function will not calculate backtransform parameters if the data are not
    # normalized.
    omicsNorm$e_data[, -id_col] <- (
      omicsNorm$e_data[, -id_col] -
        rep(normRes_obj$parameters$normalization$location,
            each = nrow(omicsNorm$e_data)) /
        rep(normRes_obj$parameters$normalization$scale,
            each = nrow(omicsNorm$e_data))
    )

    # Backtransmute if there are location and scale backtransform parameters.
    if (!is.null(normRes_obj$parameters$backtransform$location) &&
        !is.null(normRes_obj$parameters$backtransform$scale)) {

      omicsNorm$e_data[, -id_col] <- (
        omicsNorm$e_data[, -id_col] *
          normRes_obj$parameters$backtransform$scale +
          normRes_obj$parameters$backtransform$location
      )

    }

  }

  # Update the norm_info element of the data_info attribute for the omicsNorm
  # object.
  attributes(omicsNorm)$data_info$norm_info$is_normalized <- TRUE

  # Generate majestic plots ----------------------------------------------------

  # Farm boy, make me a plot with a raw object. As you wish.
  p_raw <- plot_omicsData(omicsData = omicsRaw, order_by = order_by,
                          color_by = color_by, facet_by = facet_by,
                          facet_cols = facet_cols, interactive = interactive,
                          x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                          y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                          title_lab = title_lab,
                          title_lab_size = title_lab_size,
                          legend_lab = legend_lab,
                          legend_position = legend_position, ylimit = ylimit,
                          bw_theme = bw_theme, palette = palette,
                          use_VizSampNames = use_VizSampNames)

  # Farm boy, make me a plot with a normalized object. As you wish.
  p_norm <- plot_omicsData(omicsData = omicsNorm, order_by = order_by,
                           color_by = color_by, facet_by = facet_by,
                           facet_cols = facet_cols, interactive = interactive,
                           x_lab = x_lab, y_lab = y_lab,
                           x_lab_size = x_lab_size, y_lab_size = y_lab_size,
                           x_lab_angle = x_lab_angle, title_lab = title_lab,
                           title_lab_size = title_lab_size,
                           legend_lab = legend_lab,
                           legend_position = legend_position, ylimit = ylimit,
                           bw_theme = bw_theme, palette = palette,
                           use_VizSampNames = use_VizSampNames)

  # Farm boy, combine the plots and send them into the world. As you wish.
  if (!interactive) {

    # Return the regular plots side-by-side.

    p_raw + p_norm

  } else {

    # Return the interactive plots side-by-side.
    plotly::subplot(p_raw, p_norm, nrows = 1)

  }

}

#' plot.isobaricpepData
#'
#' For plotting isobaricpepData S3 objects
#'
#' @param omicsData An isobaricpepData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-isobaricpepData
#'
#' @export
#'
plot.isobaricpepData <- function (omicsData, order_by = NULL, color_by = NULL,
                                  facet_by = NULL, facet_cols = NULL,
                                  interactive = FALSE, x_lab = NULL,
                                  y_lab = NULL, x_lab_size = 11,
                                  y_lab_size = 11, x_lab_angle = 90,
                                  title_lab = NULL, title_lab_size = 14,
                                  legend_lab = NULL, legend_position = "right",
                                  ylimit = NULL, bw_theme = TRUE,
                                  palette = NULL, use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with an isobaricpepData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

#' plot.lipidData
#'
#' For plotting lipidData S3 objects
#'
#' @param omicsData A lipidData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-lipidData
#'
#' @export
#'
plot.lipidData <- function (omicsData, order_by = NULL, color_by = NULL,
                            facet_by = NULL, facet_cols = NULL,
                            interactive = FALSE, x_lab = NULL, y_lab = NULL,
                            x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                            title_lab = NULL, title_lab_size = 14,
                            legend_lab = NULL, legend_position = "right",
                            ylimit = NULL, bw_theme = TRUE, palette = NULL,
                            use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with a lipidData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

#' plot.metabData
#'
#' For plotting metabData S3 objects
#'
#' @param omicsData A metabData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-metabData
#'
#' @export
#'
plot.metabData <- function (omicsData, order_by = NULL, color_by = NULL,
                            facet_by = NULL, facet_cols = NULL,
                            interactive = FALSE, x_lab = NULL, y_lab = NULL,
                            x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                            title_lab = NULL, title_lab_size = 14,
                            legend_lab = NULL, legend_position = "right",
                            ylimit = NULL, bw_theme = TRUE, palette = NULL,
                            use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with a metabData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

#' plot.nmrData
#'
#' For plotting nmrData S3 objects
#'
#' @param omicsData An nmrData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-nmrData
#'
#' @export
#'
plot.nmrData <- function (omicsData, order_by = NULL, color_by = NULL,
                          facet_by = NULL, facet_cols = NULL,
                          interactive = FALSE, x_lab = NULL, y_lab = NULL,
                          x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",
                          ylimit = NULL, bw_theme = TRUE, palette = NULL,
                          use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with a nmrData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

#' plot.pepData
#'
#' For plotting pepData S3 objects
#'
#' @param omicsData A pepData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-pepData
#'
#' @export
#'
plot.pepData <- function (omicsData, order_by = NULL, color_by = NULL,
                          facet_by = NULL, facet_cols = NULL,
                          interactive = FALSE, x_lab = NULL, y_lab = NULL,
                          x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",
                          ylimit = NULL, bw_theme = TRUE, palette = NULL,
                          use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with a pepData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

#' plot.proData
#'
#' For plotting proData S3 objects
#'
#' @param omicsData A proData object.
#' @param order_by A character string specifying a main effect by which to order
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by A character string specifying a main effect by which to color
#'   the boxplots. This main effect must be found in the column names of f_data
#'   in the omicsData object. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by A character string specifying a main effect with which to
#'   create a facet plot. This main effect must be found in the column names of
#'   f_data in the omicsData object. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#'
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label.
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab A character string specifying the plot title.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames Logical. Indicates whether to use custom sample
#'   names. The default is FALSE.
#'
#' @rdname plot-proData
#'
#' @export
#'
plot.proData <- function (omicsData, order_by = NULL, color_by = NULL,
                          facet_by = NULL, facet_cols = NULL,
                          interactive = FALSE, x_lab = NULL, y_lab = NULL,
                          x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",
                          ylimit = NULL, bw_theme = TRUE, palette = NULL,
                          use_VizSampNames = FALSE) {

  # Farm boy, make me a plot with a proData object. As you wish.
  plot_omicsData(omicsData = omicsData, order_by = order_by,
                 color_by = color_by, facet_by = facet_by,
                 facet_cols = facet_cols, interactive = interactive,
                 x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
                 y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
                 title_lab = title_lab, title_lab_size = title_lab_size,
                 legend_lab = legend_lab, legend_position = legend_position,
                 ylimit = ylimit, bw_theme = bw_theme, palette = palette,
                 use_VizSampNames = use_VizSampNames)

}

# The function that does all the heavy lifting for the isobaricpepData,
# lipidData, metabData, nmrData, pepData, proData, and normRes plot methods.
plot_omicsData <- function (omicsData, order_by, color_by, facet_by, facet_cols,
                            interactive, x_lab, y_lab, x_lab_size, y_lab_size,
                            x_lab_angle, title_lab, title_lab_size, legend_lab,
                            legend_position, ylimit, bw_theme, palette,
                            use_VizSampNames) {

  # Preliminaries --------------------------------------------------------------

  # Keeping the user honest ---------------

  # Farm boy, make sure the data is the correct class. As you wish.
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData",
                             "lipidData", "nmrData"))) {

    # INCONCEIVABLE!!!
    stop (paste("omicsData must be of class 'isobaricpepData', 'lipidData'",
                "'metabData', 'nmrData', 'pepData', or 'proData'.",
                sep = " "))

  }

  if(!is.null(order_by)) {
    if(!is.character(order_by) || length(order_by) > 1)
      stop("order_by must be a character vector of length 1")
  }

  if(!is.null(color_by)) {
    if(!is.character(color_by) || length(color_by) > 1)
      stop("color_by must be a character vector of length 1")
  }

  if(!is.null(facet_by)) {
    if(!is.character(facet_by) || length(facet_by) > 1)
      stop("facet_by must be a character vector of length 1")
  }

  if (!is.null(order_by) || !is.null(color_by) || !is.null(facet_by)) {

    # Make sure the group designation function has been run.
    if (is.null(attr(omicsData, "group_DF"))) {

      # Welcome to the pit of despair!! You will never escape!!!
      stop (paste("group_DF must not be NULL. Run the group_designation",
                  "function prior to plotting with any of the order_by,",
                  "color_by, or facet_by arguments.",
                  sep = " "))

    }

  }

  if(!is.null(facet_cols)) {
    if(is.null(facet_by))
      stop("facet_by cannot be NULL when facet_cols is specified")
    if(length(facet_cols) > 1)
      stop("facet_cols must be of length 1")
    if(!is.numeric(facet_cols))
      stop("facet_cols must be an integer greater than zero")
    if(facet_cols %% 1 != 0 || facet_cols <= 0)
      stop("facet_cols must be an integer greater than zero")
  }

  if(!is.null(ylimit)){
    if(!is.numeric(ylimit) || length(ylimit) != 2)
      stop("ylimit must be a numeric vector of length 2")
  }

  # Checking for NAs ---------------

  # Check for samples with all NAs and return message to user that these will
  # not be plotted.
  sample_nas <- colSums(is.na(omicsData$e_data))

  if (any(sample_nas == nrow(omicsData$e_data))) {

    empties <- names(omicsData$e_data)[which(
      sample_nas == nrow(omicsData$e_data)
    )]
    message(paste("The following sample(s) are comprised entirely of missing",
                  "data and will not be included in the plot: ",
                  empties,
                  sep = " "))

  }

  # Label crap ---------------

  # Farm boy, make me a title depending on data type and norm_info. As you wish.
  if (inherits(omicsData, "isobaricpepData")) {
    maintitle <- if (attr(omicsData, "isobaric_info")$norm_info$is_normalized)
      "Boxplots of Normalized Isobaric Peptide Data" else
        "Boxplots of Un-Normalized Isobaric Peptide Data"
  } else if(inherits(omicsData, "pepData")){
    maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
      "Boxplots of Normalized Peptide Data" else
        "Boxplots of Un-Normalized Peptide Data"
  } else if(inherits(omicsData, "proData")){
    maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
      "Boxplots of Normalized Protein Data" else
        "Boxplots of Un-Normalized Protein Data"
  } else if(inherits(omicsData, "lipidData")){
    maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
      "Boxplots of Normalized Lipid Data" else
        "Boxplots of Un-Normalized Lipid Data"
  } else if(inherits(omicsData, "metabData")){
    maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
      "Boxplots of Normalized Metabolite Data" else
        "Boxplots of Un-Normalized Metabolite Data"
  } else if(inherits(omicsData, "nmrData")){
    maintitle <- if (attr(omicsData, "nmr_info")$norm_info$is_normalized)
      "Boxplots of Normalized NMR Data" else
        "Boxplots of Un-Normalized NMR Data"
  }

  # Farm boy, create an  plot subtitle object. As you wish.
  subtitle <- NULL

  # Farm boy, fashion plot, axis, and legend title objects. As you wish.
  title <- if (is.null(title_lab)) maintitle else title_lab
  xlabel <- if (is.null(x_lab)) "Sample" else x_lab
  ylabel <- if (is.null(y_lab)) {
    if (get_data_scale(omicsData) == "abundance")
      "Abundance" else
        paste(get_data_scale(omicsData), "Abundance", sep = " ")
  } else y_lab
  legend_title <- if (is.null(legend_lab)) color_by else legend_lab

  # Data crap ---------------

  # Farm boy, melt the data for me. As you wish.
  e_data_cname <- get_edata_cname(omicsData)
  plot_data <- reshape2::melt(omicsData$e_data,
                              id = e_data_cname,
                              na.rm = TRUE)

  # Farm boy, extract the group_DF attribute. As you wish.
  groupDF <- attr(omicsData, "group_DF")

  # If facet_by is not null and isn't the same as either order_by or color_by.
  if (!is.null(facet_by)) {

    if (!(facet_by %in% c(order_by, color_by))) {

      # Extract the group_DF attribute. This will be used to facet the plots
      # later in the function.
      facetDF <- attr(
        group_designation(omicsData = omicsData, main_effects = facet_by),
        "group_DF"
      )

      # Rename the columns so they can be merged with the plot_data object and
      # referred to by name when plotting.
      colnames(facetDF) <- c("variable", facet_by)

      plot_data <- merge(plot_data, facetDF, by = "variable")

    }

  }

  # Check if order_by is NULL and update the plot_data object accordingly.
  if (!is.null(order_by)) {

    # If order_by is not "group_DF" then run group designation with the
    # specified variable as the main effect.
    if (order_by != "Group") {

      # Select the column in f_data containing the sample name as well as the
      # column corresponding to the order_by input.
      orderDF <- dplyr::select(omicsData$f_data,
                               !!rlang::sym(get_fdata_cname(omicsData)),
                               !!rlang::sym(order_by))

    } else {

      # Use the original group_DF attribute for ordering the samples. This
      # occurs when order_by = "group_DF".
      orderDF <- groupDF

    }

    # Rename the columns so they can be merged with the plot_data object and
    # referred to by name when plotting.
    colnames(orderDF)[1:2] <- c("variable", order_by)

    # Add the variables to plot_data that will be used to order the box plots.
    plot_data <- merge(plot_data, orderDF, by = "variable")

    # Reorder the rows of plot_data so the box plots will be displayed in the
    # correct order.
    plot_data <- plot_data[order(plot_data[, order_by]), ]
    plot_data$variable <- factor(plot_data$variable,
                                 levels = unique(plot_data$variable),
                                 ordered = TRUE)

    # Farm boy, update the subtitle object. As you wish.
    subtitle <- paste('Ordered by ', order_by)

  }

  # Check if color_by is NULL and update plot_data accordingly.
  if (!is.null(color_by)) {

    # If color_by is not "group_DF" then run group designation with the
    # specified variable as the main effect.
    if (color_by != "Group") {

      # Extricate the group_DF data frame from omicsData after creating the
      # group_DF attribute with the color_by input. This will be combined with
      # the plot_data object so the samples can be colored by the main effect.
      colorDF <- attr(
        group_designation(omicsData = omicsData, main_effects = color_by),
        "group_DF"
      )

    } else {

      # Use the original group_DF attribute for coloring the samples. This
      # occurs when color_by = "group_DF".
      colorDF <- groupDF

    }

    # Rename the columns so they can be merged with the plot_data object and
    # referred to by name when plotting.
    colnames(colorDF)[1:2] <- c("variable", color_by)

    # Check the status of order_by. If order_by is not NULL AND it is different
    # from color_by OR if order_by is NULL AND color_by is not NULL the colorDF
    # object needs to be merged with plot_data.
    if ((!is.null(order_by) && order_by != color_by) ||
        (is.null(order_by) && !is.null(color_by))) {

      # Merge the colorDF object with plot_data so the samples can be colored
      # correctly.
      plot_data <- merge(plot_data, colorDF, by = "variable")

    }

    # Create factors to color by according to the input of color_by.
    color_levels <- if (color_by != "Group")
      unique(factor(omicsData$f_data[[color_by]])) else
        unique(factor(attr(omicsData, "group_DF")[["Group"]]))
    plot_data[[color_by]] <- factor(plot_data[[color_by]],
                                    levels = color_levels)

  }

  # Form awe-inspiring plots ---------------------------------------------------

  # Farm boy, make the graph skeleton. As you wish.
  p <- ggplot2::ggplot(plot_data)

  # Farm boy, color the plots based on the input. As you wish.
  if (is.null(color_by)) {

    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes(x = variable,
                                         y = value))

  } else {

    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes(x = variable,
                                         y = value,
                                         fill = !!rlang::sym(color_by)))

  }

  # Farm boy, use my custom sample names in the plot. As you wish.
  if (use_VizSampNames) p <- p + ggplot2::scale_x_discrete(
    labels = omicsData$f_data$VizSampNames,
    breaks = unique(omicsData$f_data[, get_fdata_cname(omicsData)])
  )

  # Farm boy, add the black and white theme to my plot. As you wish.
  if (bw_theme) {

    p <- p + ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

  }

  # Farm boy, group my plots for me. As you wish.
  if(!is.null(facet_by)) {
    if(is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~", facet_by)),
                                   scales = "free_x")
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~", facet_by)),
                                   scales = "free_x",
                                   ncol = facet_cols)
    }
  }

  # Farm boy, add thematic elements to my plot. As you wish.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      plot.subtitle = ggplot2::element_text(face = 'italic')
    )

  # Farm boy, add labels to my plot. As you wish.
  p <- p +
    ggplot2::ggtitle(title, subtitle) +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::scale_color_discrete(legend_title)


  # Farm boy, add limits to my plot. As you wish.
  if (!is.null(ylimit)) {
    p <- p + ggplot2::scale_y_continuous(limits = ylimit)
  }

  # Farm boy, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) p <- p +
    ggplot2::scale_fill_brewer(name = legend_title,
                               palette = palette)

  # Farm boy, make me an interactive plot. As you wish.
  if (interactive) p <- plotly::ggplotly(p)

  # Farm boy, send my plot into the world to be adored by all. As you wish.
  return (p)

}
