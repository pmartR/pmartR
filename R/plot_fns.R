#' Plot dataRes object
#'
#' For plotting an S3 object of type dataRes
#'
#' @param x object of class dataRes, created by the
#'   \code{\link{edata_summary}} function
#' @param metric character string indicating which metric to use in plot:
#'   'mean', 'median', 'sd, 'pct_obs', 'min', or 'max'
#' @param density logical value, defaults to FALSE. If TRUE, a density plot of
#'   the specified metric is returned.
#' @param ncols integer value specifying the number columns for the histogram
#'   facet_wrap. This argument is used when \code{metric} is not null. The
#'   default is NULL.
#' @param interactive logical value. If TRUE, produces an interactive plot.
#' @param x_lab character string specifying the x-axis label when the metric
#'   argument is NULL. The default is NULL in which case the x-axis label will
#'   be "count".
#' @param x_lab_sd character string used for the x-axis label for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param x_lab_median character string used for the x-axis label for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param y_lab_sd character string used for the y-axis label for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param y_lab_median character string used for the y-axis label for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param x_lab_size integer value indicating the font size for the x-axis. The
#'   default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis. The
#'   default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels
#' @param title_lab character string specifying the plot title when the
#'   \code{metric} argument is NULL.
#' @param title_lab_sd character string used for the plot title for the
#'   mean/standard deviation plot when the \code{metric} argument is not NULL.
#' @param title_lab_median character string used for the plot title for the
#'   mean/median plot when the \code{metric} argument is not NULL.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", or "bottom". The default is
#'   "right".
#' @param point_size integer specifying the size of the points. The default is
#'   2.
#' @param bin_width integer indicating the bin width in a histogram. The default
#'   is 0.5.
#' @param bw_theme logical value. If TRUE, uses the ggplot2 black and white
#'   theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details This function can only create plots for dataRes objects whose 'by' =
#'   'molecule' and 'groupvar' attribute is non NULL
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
#' result <- edata_summary(
#'   omicsData = mylipid,
#'   by = "molecule",
#'   groupvar = "Virus"
#' )
#' plot(result)
#'
#' @rdname plot-dataRes
#'
#' @export
#'
plot.dataRes <- function(x, metric = NULL, density = FALSE,
                         ncols = NULL, interactive = FALSE, x_lab = NULL,
                         x_lab_sd = NULL, x_lab_median = NULL, y_lab = NULL,
                         y_lab_sd = NULL, y_lab_median = NULL, x_lab_size = 11,
                         y_lab_size = 11, x_lab_angle = NULL, title_lab = NULL,
                         title_lab_sd = NULL, title_lab_median = NULL,
                         title_lab_size = 14, legend_lab = NULL,
                         legend_position = "right", point_size = 2,
                         bin_width = 1, bw_theme = TRUE, palette = NULL, ...) {
  dataRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # check that attr(dataRes_obj, "by") == "molecule"
  if (attr(dataRes_obj, "by") != "molecule") {
    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop(paste("can only plot a dataRes object if its 'by' attribute is equal",
      "to 'molecule'",
      sep = " "
    ))
  }

  # check that attr(dataRes_obj, "groupvar") is not NULL
  if (is.null(attr(dataRes_obj, "groupvar"))) {
    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop(paste("can only plot a dataRes object if its 'groupvar' attribute is",
      "not NULL",
      sep = " "
    ))
  }

  # Check if the data is on a log scale. If it is not the histograms will not
  # work properly when metric = mean, sd, ... because of the wide range of
  # abundance values.
  if (attr(dataRes_obj, "data_scale") == "abundance") {
    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop("Data must be on the log scale to plot.")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
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
    # scatterplots
    # subseting dataRes_obj object
    mean <- dataRes_obj$mean
    median <- dataRes_obj$median
    sd <- dataRes_obj$sd

    # melting data frames from dataRes object
    mean_melt <- tidyr::pivot_longer(
      mean, -!!edata_cname,
      cols_vary = "slowest",
      names_to = "variable",
      values_to = "mean"
    )
    sd_melt <- tidyr::pivot_longer(
      sd, -!!edata_cname,
      cols_vary = "slowest",
      names_to = "variable",
      values_to = "sd"
    )
    median_melt <- tidyr::pivot_longer(
      median, -!!edata_cname,
      cols_vary = "slowest",
      names_to = "variable",
      values_to = "median"
    )

    data_mean_sd <- merge(mean_melt,
      sd_melt,
      by = c(edata_cname, "variable")
    )
    data_mean_median <- merge(mean_melt,
      median_melt,
      by = c(edata_cname, "variable")
    )

    q <- ggplot2::ggplot(
      dplyr::filter(
        data_mean_sd,
        dplyr::if_all(dplyr::one_of(c("mean", "sd")), ~!is.na(.x))
      ),
      ggplot2::aes(
        x = mean,
        y = sd,
        color = variable
      )
    ) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::xlab(xLabelSd) +
      ggplot2::ylab(yLabelSd) +
      ggplot2::ggtitle(plotTitleSd)

    p <- ggplot2::ggplot(
      dplyr::filter(
        data_mean_median,
        dplyr::if_all(dplyr::one_of(c("mean", "median")), ~!is.na(.x))
      ),
      ggplot2::aes(
        x = mean,
        y = median,
        color = variable
      )
    ) +
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
        ggplot2::scale_color_brewer(
          palette = palette,
          name = legendLabel
        )
      q <- q +
        ggplot2::scale_color_brewer(
          palette = palette,
          name = legendLabel
        )
    }

    # Want an interactive plot? As you wish.
    if (interactive && requirePlotly()) {
      q <- plotly::ggplotly(q)
      p <- plotly::ggplotly(p)

      plotly::subplot(q, p, nrows = 1)
    } else {
      # Combine the plots into one plot when interactive is FALSE.

      patchwork::wrap_plots(q, p)
    }
  } else if (!is.null(metric)) {
    if (!(metric %in% c('mean', 'median', 'sd', 'pct_obs', 'min', 'max'))) {
      # There is a shortage of perfect plots in the world. It is a shame to ruin
      # this one.
      stop("metric must be one of mean, median, sd, pct_obs, min or max")
    }

    if (!is.logical(density)) stop("density must be either TRUE or FALSE")

    # if density == F, will plot faceted histograms.
    if (density == FALSE) {
      # More tedious label creating .... (deep sigh).
      xlabel <- if (is.null(x_lab)) metric else x_lab
      ylabel <- if (is.null(y_lab)) "count" else y_lab
      plotTitle <- if (is.null(title_lab))
        paste("Histograms for ", metric, sep = "") else
        title_lab

      # subsetting dataRes object
      data = dataRes_obj[[metric]]
      data_melt = tidyr::pivot_longer(
        data, -!!edata_cname,
        cols_vary = "slowest",
        names_to = "variable",
        values_to = "value"
      )

      r <- ggplot2::ggplot(
        data_melt,
        ggplot2::aes(
          x = value,
          fill = variable
        )
      ) +
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

      # if density == TRUE, will plot geom_density
      data = dataRes_obj[[metric]]
      data_melt = tidyr::pivot_longer(
        data, -!!edata_cname,
        cols_vary = "slowest",
        names_to = "variable",
        values_to = "value"
      )

      r <- ggplot2::ggplot(
        data_melt,
        ggplot2::aes(
          x = value,
          colour = variable
        )
      ) +
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
        ggplot2::scale_fill_brewer(
          palette = palette,
          name = legendLabel
        )
    }

    # Want an interactive plot? As you wish.
    if (interactive && requirePlotly()) r <- plotly::ggplotly(r)

    return(r)
  }
}

#' Plot isobaricnormRes object
#'
#' Creates box plots for an S3 object of type 'isobaricnormRes'
#'
#' @param x an object of type isobaricnormRes, created by
#'   \code{\link{normalize_isobaric}} with apply_norm = FALSE
#' @param order logical value. If TRUE the samples will be ordered by the column
#'   of f_data containing the experiment/plate information.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label
#' @param x_lab_size integer value indicating the font size for the x-axis. The
#'   default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis. The
#'   default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white
#'   theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
#' result <- normalize_isobaric(myiso,
#'   exp_cname = "Plex",
#'   apply_norm = FALSE,
#'   refpool_cname = "Virus",
#'   refpool_notation = "Pool"
#' )
#' plot(result)
#'
#' @importFrom dplyr .data
#'
#' @rdname plot-isobaricnormRes
#'
#' @export
#'
plot.isobaricnormRes <- function(x, order = FALSE,
                                 interactive = FALSE, x_lab = NULL,
                                 y_lab = NULL, x_lab_size = 11,
                                 y_lab_size = 11, x_lab_angle = NULL,
                                 title_lab = NULL, title_lab_size = 14,
                                 legend_lab = NULL, legend_position = "none",
                                 bw_theme = TRUE, palette = NULL, ...) {
  isobaricnormRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Check if input is an isobaricnormRes class object.
  if (!inherits(isobaricnormRes_obj, "isobaricnormRes")) {
    stop("object must be of class 'isobaricnormRes'")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # extracting attributes from isobaricnormRes_obj
  exp_cname = attr(isobaricnormRes_obj, "isobaric_info")$exp_cname
  fdata_cname = attr(isobaricnormRes_obj, "cnames")$fdata_cname
  edata_cname <- attr(isobaricnormRes_obj, "cnames")$edata_cname

  # Transform the isobaricnormRes data frames into a format usable by ggplot2.
  tall_data <- prime_iso(
    isonormRes = isobaricnormRes_obj,
    exp_cname = exp_cname,
    fdata_cname = fdata_cname,
    edata_cname = edata_cname
  )

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

    p <- ggplot2::ggplot(
      data = dplyr::filter(tall_data, !is.na(values)),
      ggplot2::aes(
        x = .data[[exp_cname]],
        y = values,
        fill = .data[[exp_cname]]
      )
    )

    # Otherwise separate the box plots by sample name.
  } else {
    p <- ggplot2::ggplot(
      data = dplyr::filter(tall_data, !is.na(values)),
      ggplot2::aes(
        x = .data[[fdata_cname]],
        y = values,
        fill = .data[[exp_cname]]
      )
    )
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
      axis.text.x = ggplot2::element_text(angle = x_lab_angle),
      legend.position = legend_position
    )

  # Want to use beautiful non-default colors? As you wish.
  if (!is.null(palette)) {
    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(
        palette = palette,
        name = legendLabel
      )
  }

  # Want an interactive plot? As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

# Takes the isobaricnormRes object and returns a data frame in the correct
# format for plotting in ggplot2.
prime_iso <- function(isonormRes, exp_cname,
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
    by = fdata_cname
  )

  # Extract the indices for each column in tall_data.
  exp_cname_ind <- which(names(tall_data) == exp_cname)
  fdata_cname_ind <- which(names(tall_data) == fdata_cname)
  values_col_ind <- which(names(tall_data) == "values")

  # Convert the experiment column into a factor.
  tall_data[[exp_cname]] <- factor(tall_data[[exp_cname]])

  # Return the tall data!
  return(tall_data)
}

#' Plot nmrnormRes Object
#'
#' Creates a scatter plot for an S3 object of type 'nmrnormRes'
#'
#' @param x an object of type nmrnormRes, created by
#'   \code{\link{normalize_nmr}}
#' @param nmrData An nmrData object.
#' @param order_by A character string specifying a column in f_data by which to
#'   order the samples.
#' @param color_by A character string specifying a column in f_data by which to
#'   color the points.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label
#' @param x_lab_size integer value indicating the font size for the x-axis. The
#'   default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis. The
#'   default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size integer specifying the size of the points. The default is
#'   2.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white
#'   theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
#' mynmrnorm <- normalize_nmr(
#'   omicsData = mynmr,
#'   apply_norm = FALSE,
#'   metabolite_name = "unkm1.53"
#' )
#' plot(mynmrnorm)
#'
#' mynmrnorm2 <- normalize_nmr(
#'   omicsData = mynmr,
#'   apply_norm = FALSE,
#'   sample_property_cname = "Concentration"
#' )
#' plot(mynmrnorm2)
#'
#' @rdname plot-nmrnormRes
#' @export
#'
plot.nmrnormRes <- function(x, nmrData = NULL, order_by = NULL,
                            color_by = NULL, interactive = FALSE,
                            x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                            y_lab_size = 11, x_lab_angle = 90,
                            title_lab = NULL, title_lab_size = 14,
                            legend_lab = NULL, legend_position = "none",
                            point_size = 2, bw_theme = TRUE, palette = NULL,
                            ...) {
  nmrnormRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  if (!inherits(nmrnormRes_obj, "nmrnormRes")) {
    # My name is Evan Martin. You defiled pmart. Prepare to die.
    stop("object must be of class 'nmrnormRes'")
  }

  if (!is.null(order_by)) {
    # Farm boy, make sure an nmrData object is provided. As you wish.
    if (is.null(nmrData)) {
      # To the pain!
      stop("An nmrData object must be provided if order_by is not NULL.")
    }

    # Farm boy, make sure order_by exists in f_data. As you wish.
    if (!order_by %in% names(nmrData$f_data)) {
      # Do you hear that user? Those are the shrieking eels!
      stop("order_by must be a column in f_data.")
    }
  }

  if (!is.null(color_by)) {
    # Farm boy, make sure an nmrData object is provided. As you wish.
    if (is.null(nmrData)) {
      # To the pain!
      stop("An nmrData object must be provided if color_by is not NULL.")
    }

    # Farm boy, make sure color_by exists in f_data. As you wish.
    if (!color_by %in% names(nmrData$f_data)) {
      # Fezzik rip his arms off. Oh, you mean these column names.
      stop("color_by must be a column in f_data.")
    }
  }

  # extracting attributes from nmrnormRes_obj
  sample_property_cname <- attr(
    nmrnormRes_obj,
    "nmr_info"
  )$sample_property_cname
  metabolite_name <- attr(nmrnormRes_obj, "nmr_info")$metabolite_name
  fdata_cname <- attr(nmrnormRes_obj, "cnames")$fdata_cname

  # organize nmrnormRes_obj
  data <- data.frame(
    Sample = nmrnormRes_obj$Sample,
    value = nmrnormRes_obj$value
  )

  # Check if order_by is NULL and update the plot data object accordingly.
  if (!is.null(order_by)) {
    # Reorder the rows of data so the plot will be displayed in the correct
    # order. Dread Pirate Roberts only likes plots that are ordered. What Dread
    # Pirate Roberts doesn't like doesn't get plotted. Savvy?
    #
    # ggplot thinks it knows what is best for everyone in every situation. I
    # disagree. However, we must follow ggplot convention. ggplot orders the
    # data how they want to unless you explicitly give them the order you want.
    # The code below is telling ggplot to shove it and to plot it in the order
    # we want. BAM!!
    data$Sample <- factor(
      data$Sample,
      levels = data$Sample[order(nmrData$f_data[, order_by])],
      ordered = TRUE
    )
  }

  # Check if color_by is NULL and update the plot data accordingly.
  if (!is.null(color_by)) {
    # Create factors to color by according to the input of color_by.
    color_levels <- unique(factor(nmrData$f_data[[color_by]]))

    # Farm boy, add a variable specifying the color for each point.
    data$Color <- factor(nmrData$f_data[[color_by]],
                         levels = sort(color_levels))

  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

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

  # Create the plot skeleton according to the color_by argument.
  if (!is.null(color_by)) {
    p <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(
        x = Sample,
        y = value,
        color = Color
      )
    )

    # Check if palette is NULL or not. Hopefully it isn't so the plot will be
    # created with colors other than the super hideous default ggplot2 colors.
    if (!is.null(palette)) {
      # Create a color from the color brewer package if a palette is provided.
      colas <- RColorBrewer::brewer.pal(5, palette)
    }
  } else {
    p <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(
        x = Sample,
        y = value
      )
    )
  }

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

  # Farm boy, make me a plot with beautiful colors? As you wish.
  if (!is.null(palette)) {
    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_color_brewer(
        palette = palette,
        name = legendLabel
      )
  }

  # Farm boy, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot SPANSRes Object
#'
#' For plotting an S3 object of type 'SPANSRes'
#'
#' @param x an object of the class 'SPANSRes', created by
#'   \code{\link{spans_procedure}}
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is NULL.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param color_low character string specifying the color of the gradient for
#'   low values.
#' @param color_high character string specifying the color of the gradient for
#'   high values
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' library(pmartRdata)
#' data(pep_object)
#' mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' mypep <- group_designation(omicsData = mypep, main_effects = "Phenotype")
#' myspans <- spans_procedure(omicsData = mypep)
#' plot(myspans)
#' }
#' @rdname plot-SPANSRes
#'
#' @export
#'
plot.SPANSRes <- function(x, interactive = FALSE,
                          x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                          y_lab_size = 11, x_lab_angle = NULL,
                          title_lab = NULL, title_lab_size = 14,
                          legend_lab = NULL, legend_position = "right",
                          color_low = NULL, color_high = NULL,
                          bw_theme = TRUE, ...) {
  SPANSRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  if (!inherits(SPANSRes_obj, "SPANSRes")) {
    # Suffer the wrath of Dread Pirate Roberts!!!!!
    stop("object must be of class 'SPANSRes'")
  }

  # plotting object with numeric SPANS_score, the normalization method, and a
  # modified string specifying the subset method + parameters for that method
  SPANSRes_obj <- SPANSRes_obj %>%
    dplyr::mutate(
      ss_par = paste(subset_method, parameters, sep = " | "),
      SPANS_score = as.numeric(SPANS_score)
    ) %>%
    dplyr::left_join(attr(SPANSRes_obj, "method_selection_pvals"),
      by = c(
        "subset_method",
        "normalization_method",
        "parameters"
      )
    )

  # Farm boy, fix all the problems. As you wish. Filter rows with the highest
  # SPANS score. This subsetted/filtered data frame will be used to add points
  # to the plot for the best scoring methods. The best methods are those with
  # the highest SPANS_score. The slice_max function will select ALL rows where
  # the highest score occurs.
  da_best <- SPANSRes_obj %>%
    dplyr::slice_max(SPANS_score)

  # Do all the tedious plot label crap.
  xlabel <- if (is.null(x_lab)) "Normalization Method" else x_lab
  ylabel <- if (is.null(y_lab)) "Subset Parameters" else y_lab
  plot_title <- if (is.null(title_lab)) NULL else title_lab
  legendLabel <- if (is.null(legend_lab)) "Score" else legend_lab

  # Produce magnificent plots --------------------------------------------------

  p <- ggplot2::ggplot(data = SPANSRes_obj) +
    ggplot2::geom_tile(
      ggplot2::aes(
        x = normalization_method,
        y = ss_par,
        alpha = 1
      ),
      color = 'black'
    ) +
    ggplot2::geom_tile(
      ggplot2::aes(
        x = normalization_method,
        y = ss_par,
        fill = SPANS_score
      ),
      color = 'black'
    ) +
    ggplot2::geom_point(
      data = da_best,
      ggplot2::aes(
        x = normalization_method,
        y = ss_par,
        shape = '1'
      )
    ) +
    ggplot2::scale_alpha_continuous(
      name = 'Not Scored',
      labels = ''
    ) +
    ggplot2::scale_shape_discrete(
      name = 'Best Scores',
      labels = ''
    ) +
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
      axis.text.x = ggplot2::element_text(angle = x_lab_angle),
      legend.position = legend_position,
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  # Want an interactive plot? As you wish.
  if (interactive && requirePlotly()) {
    p <- plotly::plot_ly(
      SPANSRes_obj,
      x = ~normalization_method,
      y = ~ss_par,
      z = ~SPANS_score,
      hoverinfo = 'text',
      text = ~ paste(
        '</br> F(-Log10(HSmPV)):',
        F_log_HSmPV,
        '</br> F(Log10(NSmPV)): ',
        F_log_NSmPV,
        '</br> Scale p-value: ',
        scale_p_value,
        '</br> Location p-value',
        location_p_value
      ),
      colors = grDevices::colorRamp(
        c(
          if (is.null(color_low)) "#132B43" else color_low,
          if (is.null(color_high)) "#56B1F7" else color_high
        )
      ),
      type = "heatmap"
    ) %>%
      plotly::add_trace(
        x = da_best$normalization_method,
        y = da_best$ss_par,
        type = 'scatter',
        mode = "markers",
        marker = list(color = "black"),
        name = "Top SPANS scores",
        inherit = FALSE
      ) %>%
      plotly::colorbar(title = "SPANS score") %>%
      plotly::layout(
        plot_bgcolor = 'black',
        xaxis = list(title = "Normalization Method"),
        yaxis = list(title = "Subset Method"),
        showlegend = TRUE
      )
  }

  return(p)
}

#' Plot naRes Object
#'
#' For plotting an S3 object of type 'naRes'
#'
#' @param x list of two data frames, one containing the number of
#'   missing values by sample, and the other containing missing values by
#'   molecule
#' @param omicsData object of class 'pepData', 'proData', 'metabData',
#'   'lipidData', nmrData', or 'seqData', created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, or \code{\link{as.seqData}}, respectively.
#' @param plot_type character string specifying which type of plot to produce.
#'   The two options are 'bar' or 'scatter'.
#' @param nonmissing logical value. When FALSE, plots missing values. When TRUE,
#'   plots non-missing values.
#' @param proportion logical value. When TRUE, plots the proportion of missing 
#'   (or non-missing if \code{nonmissing} is TRUE) to the total number of 
#'   values. Only works with a plot type of 'bar'.
#' @param order_by A character string specifying a column in f_data by which to
#'   order the samples. Specifying "Group" will use the "Group" column of the
#'   object's \code{group_DF} attribute to order the samples. Only works with a 
#'   plot type of 'bar'.
#' @param color_by A character string specifying a column in f_data by which to
#'   color the bars or the points depending on the \code{plot_type}. Specifying
#'   "Group" will use the "Group" column of the object's \code{group_DF}
#'   attribute to color the samples. Only works with a plot type of 'bar'.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab_bar character string used for the x-axis label for the bar
#'   plot
#' @param x_lab_scatter character string used for the x-axis label for the
#'   scatter plot
#' @param y_lab_bar character string used for the y-axis label for the bar
#'   plot
#' @param y_lab_scatter character string used for the y-axis label for the
#'   scatter plot
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#' @param title_lab_bar character string used for the plot title when
#'   \code{plot_type} is 'bar'.
#' @param title_lab_scatter character string used for the plot title when
#'   \code{plot_type} is 'scatter'.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab_bar character string specifying the legend title when
#'   creating a bar plot.
#' @param legend_lab_scatter character string specifying the legend title when
#'   creating a scatter plot.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", or "bottom". The default is
#'   "right".
#' @param point_size An integer specifying the size of the points. The default
#'   is 2.
#' @param text_size An integer specifying the size of the text (number of
#'   missing values by sample) within the bar plot. The default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count logical value. Indicates whether the missing value
#'   counts by sample will be displayed on the bar plot. The default is TRUE.
#' @param coordinate_flip logical value. Indicates whether the x and y axes will
#'   be flipped. The default is FALSE.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @details This function takes in an object of class naRes and creates either a
#'   bar or scatter plot of missing values. When plot_type = 'bar', a sample
#'   name by missing values count bar chart is returned. When plot_type =
#'   'scatter' a mean intensity vs number of missing values (per molecule)
#'   scatter plot is returned. Note: If the omicsData object has had
#'   \code{\link{group_designation}} applied to it, the points in the plot will
#'   be colored by group.
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mylipid <- group_designation(omicsData = lipid_neg_object, main_effects = "Virus")
#' result <- missingval_result(omicsData = mylipid)
#' plot(result, omicsData = mylipid, plot_type = "bar",
#'      x_lab_angle = 50, order_by = "Virus", color_by = "Virus")
#' plot(result, omicsData = mylipid, plot_type = "scatter",
#'      x_lab_angle = 50, color_by = "Virus")
#'
#' result <- missingval_result(omicsData = rnaseq_object)
#' plot(result, omicsData = rnaseq_object, plot_type = "bar")
#'
#' @rdname plot-naRes
#'
#' @export
#'
plot.naRes <- function (x, omicsData, plot_type = "bar", 
                        nonmissing = FALSE, proportion = FALSE,
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
                        coordinate_flip = FALSE, use_VizSampNames = FALSE,
                        ...) {

  naRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # check for a naRes object #
  if (!inherits(naRes_obj, "naRes")) stop("object must be of class 'naRes'")

  # Check that omicsData is the correct class.
  if (!inherits(omicsData, c(
    "proData", "pepData", "lipidData",
    "metabData", "nmrData", "seqData"
  ))) {
    # Fezzik, tear his arms off.
    stop("omicsData is not an appropriate class")
  }

  # Check that type is either bar or scatter.
  if (!(plot_type %in% c("bar", "scatter"))) {
    # My name is Evan Martin. You killed my plot. Prepare to die.
    stop("plot_type must be either 'bar' or 'scatter'")
  }
  
  # Check to make sure proportion is only used with bar plot
  if (proportion && plot_type != "bar") {
    
    stop("plot_type must be 'bar' if proportion is TRUE")
    
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # Farm boy, make sure order_by exists in f_data. As you wish.
  if (!is.null(order_by) && order_by != "Group") {

    if (!order_by %in% names(omicsData$f_data)) {
      # I'm a pmartR developer. You killed my plot. Prepare to receive an error.
      stop ("order_by: column '", order_by, "' not found in f_data.")

    }
  }

  # Farm boy, make sure color_by exists in f_data. As you wish.
  if (!is.null(color_by) && color_by != "Group") {
    
    if (!color_by %in% names(omicsData$f_data)) {
      # Clearly you cannot choose a column name in f_data!
      stop ("color_by: column '", color_by, "' not found in f_data.")

    }
  }

  # Extract info from naRes_obj
  if (inherits(omicsData, "seqData")) {
    na.by.sample <- naRes_obj$zeros.by.sample
    na.by.molecule <- naRes_obj$zeros.by.molecule
    num_missing_vals <- na.by.molecule$num_zeros
    num_nonmissing_vals <- na.by.molecule$num_nonzeros
    missing_proportion <- na.by.molecule$zeros_proportion
    names(na.by.sample)[which(names(na.by.sample) == "num_zeros")] <- "num_NA"
    names(na.by.molecule)[which(names(na.by.molecule) == "num_zeros")] <- "num_NA"
    names(na.by.sample)[which(names(na.by.sample) == "num_nonzeros")] <- "num_non_NA"
    names(na.by.molecule)[which(names(na.by.molecule) == "num_nonzeros")] <- "num_non_NA"
  } else {
    na.by.sample <- naRes_obj$na.by.sample
    na.by.molecule <- naRes_obj$na.by.molecule
    num_missing_vals <- na.by.molecule$num_NA
    num_nonmissing_vals <- na.by.molecule$num_non_NA
  }
  
  # Get the proportion of missing/nonmissing to total (missing + nonmissing)
  target_NA_column <- ifelse(nonmissing, "num_non_NA", "num_NA")
  na.by.molecule$NA_proportion <- na.by.molecule[[target_NA_column]] / 
    (na.by.molecule$num_NA + na.by.molecule$num_non_NA)
  na.by.sample$NA_proportion <- na.by.sample[[target_NA_column]] / 
    (na.by.sample$num_NA + na.by.sample$num_non_NA)

  edata_cname <- attr(naRes_obj, "cnames")$edata_cname
  fdata_cname <- attr(naRes_obj, "cnames")$fdata_cname

  # Extract info from omicsData
  edata <- omicsData$e_data
  edata_cname_id <- which(names(edata) == edata_cname)
  group_df <- get_group_DF(omicsData)

  # Bar plot order_by and group_by crap ---------------

  # Check if order_by is NULL and update the plot_data object accordingly.
  if (!is.null(order_by)) {

    if (order_by == "Group") {
      na.by.sample <- na.by.sample[order(na.by.sample$Group), ]
      na.by.sample[[fdata_cname]] <- factor(
        na.by.sample[[fdata_cname]],
        levels = na.by.sample[[fdata_cname]],
        ordered = TRUE
      )
      
    } else {
    
      # Reorder the rows of na.by.sample so the bar plot will be displayed in the
      # correct order.
      na.by.sample <- na.by.sample[order(na.by.sample[, order_by]), ]
      na.by.sample[[fdata_cname]] <- factor(
        na.by.sample[[fdata_cname]],
        levels = na.by.sample[[fdata_cname]],
        ordered = TRUE
      )
      
    }

  }

  # Check if color_by is NULL and update na.by.sample accordingly.
  if (!is.null(color_by)) {
    # Create factors to color by according to the input of color_by.
    color_levels <- if (color_by != "Group")
      unique(factor(omicsData$f_data[[color_by]])) else
      unique(factor(attr(omicsData, "group_DF")[["Group"]]))
    na.by.sample[[color_by]] <- factor(na.by.sample[[color_by]],
                                       levels = sort(color_levels))
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
                fdata_cname = fdata_cname, color_by = color_by,
                nonmissing = nonmissing, proportion = proportion)

  }

  # Make me a scatter plot. As you wish.
  if (plot_type == "scatter") {

    p <- na_scatter(edata = edata, group_df = group_df,
                    na.by.molecule = na.by.molecule,
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
                    interactive = interactive, point_size = point_size,
                    nonmissing = nonmissing)

  }

  return(p)
}

na_bar <- function (na.by.sample, x_lab_bar, y_lab_bar, x_lab_size, y_lab_size,
                    x_lab_angle, title_lab_bar, title_lab_size, legend_lab_bar,
                    legend_position, text_size, bar_width, bw_theme, palette,
                    display_count, coordinate_flip, use_VizSampNames,
                    interactive, fdata_cname, color_by, nonmissing,
                    proportion) {

  # Select which column of na.by.sample will be used
  if (proportion) {
    y_axis = "NA_proportion"
  } else if (nonmissing) {
    y_axis = "num_non_NA"
  } else {
    y_axis = "num_NA"
  }

  # Farm boy, color the plots based on the input. As you wish.
  if (!is.null(color_by)) {
    
    # Forge the basic sample bar plot with group info. More details will be
    # added according to the users input later.
    samp <- ggplot2::ggplot(
      data = na.by.sample,
      ggplot2::aes(
        x = .data[[fdata_cname]],
        y = .data[[y_axis]],
        fill = !!dplyr::sym(color_by)
      )
    ) +
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
                                         y = .data[[y_axis]])) +
      ggplot2::geom_bar(stat = "identity",
                        width = bar_width,
                        fill = if (is.null(palette))
                          hideous else
                            colas[[3]])

  }

  # Add the counts to the bar plot and histogram if the user so desires.
  if (display_count) {
    samp <- samp +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[[y_axis]], digits=4)),
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
    sprintf(
      "%s of %s values",
      ifelse(proportion, "Proportion", "Number"),
      ifelse(nonmissing, "non-missing", "missing")
    ) else
      y_lab_bar
  plotTitleBar <- if (is.null(title_lab_bar))
    sprintf(
      "%s values by sample",
      ifelse(nonmissing, "Non-missing", "Missing")
    ) else
      title_lab_bar
  legendLabelBar <- if (is.null(legend_lab_bar))
    color_by else
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
      ggplot2::scale_fill_brewer(
        palette = palette,
        name = legendLabelBar
      )
  }

  # Want to use short sample names? As you wish.
  if (use_VizSampNames) {
    samp <- samp +
      ggplot2::scale_x_discrete(labels = na.by.sample$VizSampNames)
  }

  # Want to be bass ackwards? As you wish. (This is for the crazies.)
  if (coordinate_flip) samp <- samp + ggplot2::coord_flip()

  # Want an interactive plot? As you wish.
  if (interactive && requirePlotly()) samp <- plotly::ggplotly(samp)

  return(samp)
}

na_scatter <- function (edata, group_df, na.by.molecule, edata_cname,
                        edata_cname_id, fdata_cname, x_lab_scatter,
                        y_lab_scatter, title_lab_scatter, legend_lab_scatter,
                        legend_position, title_lab_size, x_lab_size, y_lab_size,
                        text_size, bw_theme, palette, x_lab_angle,
                        coordinate_flip, interactive, point_size,
                        nonmissing) {

  # Select missing/nonmissing
  if (nonmissing)
    num_missing_vals <- na.by.molecule$num_non_NA
  else
    num_missing_vals <- na.by.molecule$num_NA
  
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
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        mean_intensity,
        num_missing_vals
      )
    ) +
      ggplot2::geom_point(size = point_size)
  } else {
    # Extract group information to calculate the group-wise mean for each
    # molecule.
    levels <- unique(group_df$Group)
    indices_list <- vector(
      mode = "list",
      length = length(levels)
    )

    for (i in 1:length(levels)) {
      # Extract the column indices for the ith group. This will be used to
      # subset e_data.
      inds <- which(group_df$Group == levels[i])
      indices_list[[i]] <- inds
    }

    # Calculate the mean intensity for each molecule by group. NaN can appear if
    # an entire row has all NA values.
    mean_by_group <- lapply(indices_list,
      function(x, temp_edata) rowMeans(temp_edata[x],
        na.rm = TRUE
      ),
      temp_edata = edata[, -edata_cname_id]
    )

    names(mean_by_group) <- levels

    mean_intensity <- do.call(cbind, mean_by_group)

    plot_data <- cbind(num_missing_vals, mean_intensity)
    plot_data <- as.data.frame(plot_data)
    plot_data <- plot_data %>%
      tidyr::pivot_longer(
        -num_missing_vals,
        cols_vary = "slowest",
        names_to = "variable",
        values_to = "value"
      ) %>% dplyr::mutate(
        variable = factor(
          variable,
          levels = names(plot_data)[
            -which(names(plot_data) == "num_missing_vals")
          ]
        )
      ) %>% data.frame

    # Start the scatter plot when the group_DF attribute is present.
    p <- ggplot2::ggplot(
      dplyr::filter(plot_data, !is.na(value)),
      ggplot2::aes(
        value,
        num_missing_vals
      )
    ) +
      ggplot2::geom_point(ggplot2::aes(color = variable),
        size = point_size
      )
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
      ggplot2::scale_color_brewer(
        palette = palette,
        name = legendLabelScatter
      )
  }

  # Want an interactive plot? As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot corRes Object
#'
#' For plotting an S3 object of type 'corRes'
#'
#' @param x An object of class "corRes" created via \code{cor_result}
#' @param omicsData an object of the class 'pepData', 'isobaricpepData',
#'   'proData', 'lipidData', 'metabData', 'nmrData' or 'seqData' created via
#'   \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}},
#'   \code{\link{as.proData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, \code{\link{as.nmrData}}, or
#'   \code{\link{as.seqData}}, respectively.
#' @param order_by A character string specifying a column in f_data by which to
#'   order the samples.
#' @param x_text logical value. Indicates whether the x-axis will be labeled
#'   with the sample names. The default is TRUE.
#' @param y_text logical value. Indicates whether the y-axis will be labeled
#'   with the sample names. The default is TRUE.
#' @param colorbar_lim A pair of numeric values specifying the minimum and
#'   maximum values to use in the heat map color bar. Defaults to 'c(NA, NA)',
#'   in which case ggplot2 automatically sets the minimum and maximum values
#'   based on the correlation values in the data.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 90.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param color_low character string specifying the color of the gradient for
#'   low values
#' @param color_high character string specifying the color of the gradient for
#'   high values
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @rdname plot-corRes
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' mymetab <- group_designation(omicsData = mymetab, main_effects = "Phenotype")
#' my_correlation <- cor_result(omicsData = mymetab)
#' plot(my_correlation, omicsData = mymetab, order_by = "Phenotype")
#'
#' \donttest{
#' myseq_correlation <- cor_result(omicsData = rnaseq_object)
#' plot(myseq_correlation)
#' }
#'
#' @export
#'
plot.corRes <- function(x, omicsData = NULL, order_by = NULL,
                        colorbar_lim = c(NA, NA), x_text = TRUE,
                        y_text = TRUE, interactive = FALSE, x_lab = NULL,
                        y_lab = NULL, x_lab_size = 11, y_lab_size = 11,
                        x_lab_angle = 90, title_lab = NULL,
                        title_lab_size = 14, legend_lab = NULL,
                        legend_position = "right", color_low = NULL,
                        color_high = NULL, bw_theme = TRUE,
                        use_VizSampNames = FALSE, ...) {
  corRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # check for a corRes object #
  if (!inherits(corRes_obj, "corRes")) {
    # You will be given iocane powder due to your carelessness.
    stop("corRes_obj must be of class 'corRes'")
  }

  if (!all(is.na(colorbar_lim))) {
    if (!is.numeric(colorbar_lim) || length(colorbar_lim) != 2) {
      # I can clearly not choose the arguments in front of you!
      stop("colorbar_lim must be a numeric vector of length 2")
    }
  }

  # check that omicsData is of appropriate class #
  if (!is.null(omicsData)) {
    if (!inherits(omicsData, c(
      "pepData", "proData", "metabData",
      "lipidData", "nmrData", "seqData"
    ))) {
      stop("omicsData is not an appropriate class")
    }
  }

  if (!is.null(order_by)) {
    # Farm boy, make sure order_by exists in f_data. As you wish.
    if (!order_by %in% names(omicsData$f_data)) {
      # I do not think that input means what you think it means.
      stop("order_by must be a column in f_data.")
    }
  }

  # Create plot data matrix ----------------------------------------------------

  # Workaround for certain "check" warnings
  Var1 <- Var2 <- value <- NULL

  # Nab groupDF when omicsData is available and order_by is specified. BOOM!
  if (!is.null(order_by)) {
    # Farm boy, extract the group_DF attribute. As you wish.
    groupDF <- attr(omicsData, "group_DF")

    # If order_by is not "Group" then run group designation with the specified
    # variable as the main effect.
    if (order_by != "Group") {
      # Fish out the group_DF data frame from omicsData after creating the
      # group_DF attribute with the order_by input. This will be combined with
      # the plot_data object so the samples can be ordered by the main effect.
      orderDF <- attr(
        group_designation(omicsData = omicsData, main_effects = order_by),
        "group_DF"
      )
    } else {
      # Use the original group_DF attribute for ordering the samples. This
      # occurs when order_by = "group_DF".
      orderDF <- groupDF
    }

    # Reorder the orderDF according to the order_by input. We can subset using
    # "Group" because it will always be the column we are after. If the input to
    # order_by = "Group" we will use the Group column from the group_DF
    # attribute from the original input. If order_by is a variable other than
    # group, we will run the group_designation function and the name of this
    # variable will be changed to Group.
    orderDF <- orderDF[order(orderDF$Group), ]

    # Save the original is_normalized attribute. This attribute is erased when
    # the rows and columns are reordered based on the orderDF object. This
    # attribute is needed later in the function to correctly label the plot.
    normal_attr <- attr(corRes_obj, "is_normalized")

    # Reorder the rows and columns of the corRes object to match the order of
    # the groups in orderDF. We can hard-code the first column of orderDF
    # because this will always be the column containing the sample names.
    # 1. Reorder the columns of corRes_obj to match the order of orderDF.
    corRes_obj <- corRes_obj[, match(orderDF[, 1], colnames(corRes_obj))]
    # 2. Reorder the rows of corRes_obj to match the order of orderDF.
    corRes_obj <- corRes_obj[match(orderDF[, 1], rownames(corRes_obj)), ]

    # Add the lost attributes back to the corRes_obj because they are needed
    # later on.
    attr(corRes_obj, "is_normalized") <- normal_attr
  }

  # Create the data frame that will be used to produce the correlation heat map.
  corRes_obj_df <- data.frame(corRes_obj, check.names = FALSE)
  corRes_obj_df <- cbind(Var1 = rownames(corRes_obj_df), corRes_obj_df)
  rownames(corRes_obj_df) <- 1:nrow(corRes_obj_df)
  corRes_melt <- corRes_obj_df %>%
    tidyr::pivot_longer(
      -Var1,
      names_to = "Var2",
      cols_vary = "slowest"
    ) %>%
    dplyr::mutate(
      Var1 = factor(
        abbreviate(Var1, minlength = 20),
        levels = abbreviate(rownames(corRes_obj), minlength = 20)
      ),
      Var2 = factor(
        abbreviate(Var2, minlength = 20),
        levels = abbreviate(rownames(corRes_obj), minlength = 20)
      )
    ) %>%
    data.frame(check.names = FALSE)

  # Create all the plot labels. Life is pain!!!
  xlabel <- if (is.null(x_lab)) "" else x_lab
  ylabel <- if (is.null(y_lab)) "" else y_lab
  legendLabel <- if (is.null(legend_lab)) "Correlation" else legend_lab
  if (is.null(title_lab)) {
    # include correlation method in title
    plotTitle <- paste0(
      "Correlations Among Samples (",
      stringr::str_to_title(attr(corRes_obj, "cor_method")),
      ")"
    )

    # Runs when title_lab is not NULL (the user specified title).
  } else {
    plotTitle <- title_lab
  }

  # Design spectacular plots ---------------------------------------------------

  # Start the skeleton of the heat map. Other aspects are forthcoming.
  hm <- ggplot2::ggplot(
    corRes_melt,
    ggplot2::aes(
      x = Var2,
      y = Var1
    )
  ) +
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
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank()
      )
  }

  # Farm boy, make me a plot without the x-axis sample labels. As you wish.
  if (!x_text) hm <- hm + ggplot2::theme(axis.text.x = ggplot2::element_blank())

  # Farm boy, make me a plot without the y-axis sample labels. As you wish.
  if (!y_text) hm <- hm + ggplot2::theme(axis.text.y = ggplot2::element_blank())

  # Farm boy, make me a plot with short sample names. As you wish.
  if (use_VizSampNames) {
    if (is.null(omicsData))
      stop("If using custom sample names omicsData must not be NULL.")

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
  if (interactive && requirePlotly()) hm <- plotly::ggplotly(hm)

  return(hm)
}

#' Plot dimRes Object
#'
#' For plotting an S3 object of type 'dimRes'
#'
#' @param x object of class dimRes created by the \code{dim_reduction}
#'   function
#' @param omicsData optional omicsData for use in specifying a column name in
#'   fdata when using \code{color_by} or \code{shape_by}.
#' @param color_by character string specifying which column to use to control 
#'   the color for plotting. NULL indicates the default value of the main effect 
#'   levels (if present). "Group" uses the "Group" column of group_DF. NA 
#'   indicates no column will be used, and will use the default theme color. If 
#'   an omicsData object is passed, any other value will use the specified 
#'   column of f_data.
#' @param shape_by character string specifying which column to use to control 
#'   the shape for plotting. NULL indicates the default value of the second main
#'   effect levels (if present). "Group" uses the "Group" column of group_DF. NA 
#'   indicates no column will be used, and will use the default theme shape. If 
#'   an omicsData object is passed, any other value will use the specified 
#'   column of f_data.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 4.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#'
#' mylipid <- edata_transform(omicsData = lipid_neg_object, data_scale = "log2")
#' mylipid <- group_designation(omicsData = mylipid, main_effects = "Virus")
#' pca_lipids <- dim_reduction(omicsData = mylipid)
#' plot(pca_lipids)
#'
#' \donttest{
#' myseq <- group_designation(omicsData = rnaseq_object, main_effects = "Virus")
#' pca_seq <- dim_reduction(omicsData = myseq)
#' plot(pca_seq)
#' }
#'
#' @rdname plot-dimRes
#'
#' @export
#'
plot.dimRes <- function (x, omicsData = NULL,
                         color_by = NULL, shape_by = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL, 
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 0, 
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         point_size = 4, bw_theme = TRUE, palette = NULL, ...) {

  dimRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Evan, make sure the input is the correct class. As you wish.
  if (!inherits(dimRes_obj, "dimRes")) {
    # Welcome to the pit of despair!
    stop("dimRes_obj must be of class 'dimRes'")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }
  plotdata <- data.frame(SampleID = dimRes_obj$SampleID,
                         PC1 = dimRes_obj$PC1,
                         PC2 = dimRes_obj$PC2)
  plotdata_name <- names(plotdata)[1]

  # if there is a group designation #
  if (!is.null(attr(dimRes_obj, "group_DF"))) {
    group_DF <- attr(dimRes_obj, "group_DF")
    fdata_cname <- names(group_DF)[1]

    # if there's two main effects #
    if (ncol(group_DF) == 4) {
      main_eff_names <- names(group_DF[, 3:4])
      # function for shortening long main effect names #
      abbrev_fun <- function(string) {
        if (nchar(string) > 25) string = paste0(substr(string, 1, 23), "...")
        return(string)
      }

      # manage the length of legend titles #
      display_names <- sapply(main_eff_names, abbrev_fun)

      if (!identical(
        as.character(plotdata[[plotdata_name]]),
        as.character(group_DF[[fdata_cname]])
      )) {
        group_DF <- group_DF[match(
          plotdata[[plotdata_name]],
          group_DF[[fdata_cname]]
        ), ]
        plotdata <- merge.data.frame(plotdata,
          group_DF,
          by.x = plotdata_name,
          by.y = fdata_cname,
          sort = FALSE
        )
      } else {
        plotdata <- merge.data.frame(plotdata,
          group_DF,
          by.x = plotdata_name,
          by.y = fdata_cname,
          sort = FALSE
        )
        plotdata[[fdata_cname]] <- plotdata[[plotdata_name]]
      }

      color_var <- main_eff_names[1]
      pch_var <- main_eff_names[2]

      if (length(legend_lab) > length(main_eff_names)) {
        warning(paste("legend_lab length is greater than the number of main",
          "effects. Only the first two entries will be used.",
          sep = " "
        ))
      }
    } else {
      if (!identical(
        as.character(plotdata[[plotdata_name]]),
        as.character(group_DF[[fdata_cname]])
      )) {
        group_DF <- group_DF[match(
          plotdata[[plotdata_name]],
          group_DF[[fdata_cname]]
        ), ]
        plotdata <- merge.data.frame(plotdata,
          group_DF,
          by.x = plotdata_name,
          by.y = fdata_cname,
          sort = FALSE
        )
      } else {
        plotdata <- merge.data.frame(plotdata,
          group_DF,
          by.x = plotdata_name,
          by.y = fdata_cname,
          sort = FALSE
        )
        plotdata[[fdata_cname]] <- plotdata[[plotdata_name]]
      }

      color_var <- "Group"
      pch_var <- NULL
      display_names <- c("Group", NULL)

      if (length(legend_lab) > 1) {
        warning(paste("legend_lab length is greater than the number of main",
          "effects. Only the first entry will be used.",
          sep = " "
        ))
      }
    }
    
    # Add any columns from f_data to the plotdata if present (for color_by
    # and shape_by)
    if (!is.null(omicsData)) {
      if (is.null(omicsData$f_data)) {
        stop("omicsData does not have f_data")
      }
      
      fdata_concat <- omicsData$f_data[
        c(
          which(colnames(omicsData$f_data) == fdata_cname),
          which(!colnames(omicsData$f_data) %in% colnames(plotdata))
        )
      ]
      
      plotdata <- merge.data.frame(plotdata,
                                   fdata_concat,
                                   by.x = plotdata_name,
                                   by.y = fdata_cname,
                                   sort = FALSE)
      plotdata[[fdata_cname]] <- plotdata[[plotdata_name]]
      
    }

    # Runs when there is no group information.
  } else {
    color_var <- NULL
    pch_var <- NULL
    display_names <- c(NULL, NULL)

    if (!is.null(legend_lab)) {
      warning("There is no group designation, so legend_lab will go unused.")
    }
    
    # Add any columns from f_data to the plotdata if present (for color_by
    # and shape_by)
    if (!is.null(omicsData)) {
      if (is.null(omicsData$f_data)) {
        stop("omicsData does not have f_data")
      }
      
      fdata_cname <- get_fdata_cname(omicsData)
      
      fdata_concat <- omicsData$f_data[
        c(
          which(colnames(omicsData$f_data) == fdata_cname),
          which(!colnames(omicsData$f_data) %in% colnames(plotdata))
        )
      ]
      
      plotdata <- merge.data.frame(plotdata,
                                   fdata_concat,
                                   by.x = plotdata_name,
                                   by.y = fdata_cname,
                                   sort = FALSE)
      plotdata[[fdata_cname]] <- plotdata[[plotdata_name]]
      
    }

  }
  
  if (!is.null(color_by)) {
    if (is.na(color_by)) {
      color_var <- NULL
    } else if (color_by == "Group") {
      color_var <- "Group"
      display_names[1] <- "Group"
    } else {
      if(is.null(omicsData)) {
          stop("color_by value is invalid. Did you mean to add an omicsData?")
      }
      
      if (!color_by %in% colnames(omicsData$f_data)) {
        if (is.null(omicsData)) {
        } else {
          stop("color_by: column '", color_by, "' not found in f_data.")
        }
      }
      
      color_var <- color_by
      display_names[1] <- color_by
      plotdata[[color_by]] <- as.factor(plotdata[[color_by]])
    }
  }
  
  if (!is.null(shape_by)) {
    if (is.na(shape_by)) {
      pch_var <- NULL
    } else if (shape_by == "Group") {
      pch_var <- "Group"
      display_names[2] <- "Group"
    } else {
      if(is.null(omicsData)) {
        stop("shape_by value is invalid. Did you mean to add an omicsData?")
      }
      
      if (!shape_by %in% colnames(omicsData$f_data)) {
        if (is.null(omicsData)) {
        } else {
          stop("shape_by: column '", shape_by, "' not found in f_data.")
        }
      }
      
      pch_var <- shape_by
      display_names[2] <- shape_by
      plotdata[[shape_by]] <- as.factor(plotdata[[shape_by]])
    }
  }

  # custom legend names #
  if (!is.null(legend_lab)) {
    # make the vector at least length 2 to avoid errors in the plot
    display_names[1:length(legend_lab)] <- legend_lab[1:min(
      2,
      length(legend_lab)
    )]
  }

  # title #
  title_default <- ifelse(is.null(attr(dimRes_obj, "R2")),
    "Principal Components (GLM-PCA)",
    "Principal Components"
  )

  plot_title <- ifelse(is.null(title_lab), title_default, title_lab)

  # Construct impressive plots -------------------------------------------------

  # Create the bare bones plot.
  p <- ggplot2::ggplot(
    plotdata,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      text = paste("Sample name: ", SampleID)
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        col = if(!is.null(color_var)) .data[[color_var]] else NULL,
        pch = if(!is.null(pch_var)) .data[[pch_var]] else NULL
      ),
      size = point_size
    ) + ggplot2::labs(col = color_var, pch = pch_var)

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  if (is.null(attr(dimRes_obj, "R2"))) {
    x_lab <- ifelse(is.null(x_lab), "PC1", x_lab)
    y_lab <- ifelse(is.null(y_lab), "PC2", y_lab)
  } else {
    # Create objects for the percent variation for PC1 and PC2
    pc1 <- round(attr(dimRes_obj, "R2")[1], 3)
    pc2 <- round(attr(dimRes_obj, "R2")[2], 3)

    # Tedious label making crap.
    # These two lines don't currently work as intended. They display R^2 without
    # the 2 as a superscript. However, the label cannot be created with substitute
    # because it throws an error when converting to an interactive plot.
    x_lab <- if (is.null(x_lab)) paste0(
      "PC1 (",
      expression(R^2),
      " = ",
      pc1,
      ")"
    ) else x_lab
    y_lab <- if (is.null(y_lab)) paste0(
      "PC2 (",
      expression(R^2),
      " = ",
      pc2,
      ")"
    ) else y_lab
  }

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
      ggplot2::scale_color_brewer(
        palette = palette,
        name = display_names[1]
      )

    # Evan, make me a plot with hideous default ggplot2 colors. As you wish.
  } else {
    p <- p +
      ggplot2::scale_color_discrete(display_names[1])
  }

  # Evan, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot moleculeFilt Object
#'
#' For plotting an S3 object of type 'moleculeFilt':
#'
#' @param x object of class moleculeFilt that contains the molecule
#'   identifier and the number of samples for which the molecule was measured
#'   (not NA)
#' @param min_num An integer specifying the minimum number of samples in which a
#'   biomolecule must appear. If a value is specified, a horizontal line will be
#'   drawn when \code{cumulative=TRUE}, and bars will be colored appropriately
#'   if \code{cumulative=FALSE}. Defaults to NULL.
#' @param cumulative logical indicating whether the number of biomolecules
#'   observed in \emph{at least} (TRUE) x number of samples or \emph{exactly}
#'   (FALSE) x number of samples should be plotted.  Defaults to TRUE.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size integer specifying the size of the text (number of
#'   biomolecules by sample) within the bar plot. The default is 3.
#' @param bar_width integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count logical value. Indicates whether the missing value counts by
#'   sample will be displayed on the bar plot. The default is TRUE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' data(pep_object)
#' molfilt <- molecule_filter(omicsData = pep_object)
#' plot(molfilt, min_num = 5)
#' plot(molfilt, min_num = 3, cumulative = FALSE)
#'
#' @rdname plot-moleculeFilt
#'
#' @export
#'
plot.moleculeFilt <- function(x, min_num = NULL, cumulative = TRUE,
                              interactive = FALSE, x_lab = NULL, y_lab = NULL,
                              x_lab_size = 11, y_lab_size = 11,
                              x_lab_angle = 0, title_lab = NULL,
                              title_lab_size = 14, legend_lab = NULL,
                              legend_position = "right", text_size = 3,
                              bar_width = 0.8, bw_theme = TRUE,
                              palette = NULL, display_count = TRUE, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "moleculeFilt")) {
    # Fezzik, tear his arms off.
    stop("filter_object must be of class moleculeFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # Forge sensational plots ----------------------------------------------------

  # make counts, colors, and plot shape based on value of cumulative

  # values that will be updated depending on input
  fill <- 0
  hline <- NULL
  # create a TRUE/FALSE for if we are working with batch or not
  use_batch <- attributes(filter_object)$use_batch
  use_groups <- attributes(filter_object)$use_groups

  if (cumulative) {
    # cumulative counts (>=)
    counts <- sapply(1:max(filter_object$Num_Observations), function(i) {
      filter_object[filter_object$Num_Observations >= i, ]$Num_Observations %>% length()
    })

    # append 1 to extend last step (looks awkward without this)
    counts <- c(counts, counts[length(counts)])
    # 1 appended step
    num_obs <- 1:(max(filter_object$Num_Observations) + 1)

    # shape is a step function with fixed color
    shape <- ggplot2::geom_step(
      ggplot2::aes(
        x = num_observations,
        y = counts
      ),
      color = "red"
    )

    # draw a horizontal line if min_num is specified
    hline <- if (!is.null(min_num)) ggplot2::geom_hline(
      yintercept = counts[min_num],
      linetype = "dashed"
    ) else NULL

    xlabel <- if (is.null(x_lab)) {
      if (!use_batch & !use_groups) "Number of Samples"
      else if (use_batch & !use_groups) "Number of Samples per Batch"
      else if (!use_batch & use_groups) "Number of Samples per Group"
      else "Number of Samples per Batch per Group"
    }

    ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)

    plot_title <- if (is.null(title_lab)) {
      if (!use_batch & !use_groups) "Count of biomolecules observed in at least X number of samples"
      else if (use_batch & !use_groups) "Minimum count of biomolecules observed in at least X number of samples per batch"
      else if (!use_batch & use_groups) "Minimum count of biomolecules observed in at least X number of samples per group"
      else "Minimum count of biomolecules observed in at least X number of samples per batch per group"
    } else
      title_lab
  } else if (!cumulative) {
    # counts for a specific number of nonmissing biomolecules (==)
    counts <- sapply(1:max(filter_object$Num_Observations), function(i) {
      filter_object[filter_object$Num_Observations == i, ]$Num_Observations %>%
        length()
    })

    # color by which values are kept if min_num is specified
    if (!is.null(min_num)) {
      fill <- ifelse(1:max(filter_object$Num_Observations) >= min_num,
        "retained",
        "dropped"
      )
      shape <- ggplot2::geom_bar(
        ggplot2::aes(
          x = num_observations,
          y = counts,
          fill = fill
        ),
        stat = "identity",
        width = bar_width
      )
    } else {
      # Check if palette is NULL or not. Hopefully it isn't so the plot will be
      # created with colors other than the super hideous default ggplot2 colors.
      if (!is.null(palette)) {
        # Create a color from the color brewer package if a palette is provided.
        colas <- RColorBrewer::brewer.pal(5, palette)
      }

      # Create an object for the first default ggplot2 color.
      hideous <- grDevices::hcl(
        h = 15,
        c = 100,
        l = 65
      )

      shape <- ggplot2::geom_bar(
        ggplot2::aes(
          x = num_observations,
          y = counts
        ),
        fill = if (is.null(palette))
          hideous else
          colas[[3]],
        stat = "identity",
        width = bar_width
      )
    }

    num_obs <- 1:max(filter_object$Num_Observations)

    xlabel <- if (is.null(x_lab)) "Number of samples" else x_lab
    ylabel <- ifelse(is.null(y_lab), "Count of Biomolecules", y_lab)
    plot_title <- ifelse(
      is.null(title_lab),
      "Count of biomolecules observed in exactly X number of samples",
      title_lab
    )
  }

  # create plotting dataframe
  pep_observation_counts <- data.frame(
    num_observations = num_obs,
    frequency_counts = counts,
    fill = fill
  )

  # plot #
  p <- ggplot2::ggplot(pep_observation_counts) +
    shape +
    hline +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous(breaks = max(filter_object$Num_Observations):1)

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Evan, display the biomolecule counts on the graph. As you wish.
  if (display_count) p <- p + ggplot2::geom_text(
    data = pep_observation_counts[1:max(filter_object$Num_Observations), ],
    ggplot2::aes(
      x = num_observations,
      y = frequency_counts,
      label = frequency_counts
    ),
    vjust = -0.5,
    nudge_x = ifelse(cumulative, 0.5, 0),
    size = text_size
  )

  # Add the theme elements to the plot.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ifelse(use_batch | use_groups, 11, title_lab_size)),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) {
    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(
        palette = palette,
        name = legend_lab
      )
  } else {
    p <- p + ggplot2::scale_fill_manual(
      name = legend_lab,
      values = c(
        "dropped" = "red",
        "retained" = "green"
      )
    )
  }

  # Evan, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}


#' Plot totalCountFilt Object
#'
#' For plotting an S3 object of type 'totalCountFilt':
#'
#' @param x object of class totalCountFilt that contains the molecule
#'   identifier and the number of total counts for which the molecule was measured
#'   (not NA).
#' @param min_count integer specifying the minimum number of samples in which a
#'   biomolecule must appear. Defaults to NULL.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size integer specifying the size of the text (number of
#'   biomolecules by sample) within the bar plot. The default is 3.
#' @param bar_width integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' \donttest{
#' library(pmartRdata)
#' seqfilt <- total_count_filter(omicsData = rnaseq_object)
#' plot(seqfilt, min_count = 15)
#' }
#'
#' @rdname plot-totalCountFilt
#'
#' @export
#'
plot.totalCountFilt <- function(x, min_count = NULL,
                                interactive = FALSE, x_lab = NULL, y_lab = NULL,
                                x_lab_size = 11, y_lab_size = 11,
                                x_lab_angle = 0, title_lab = NULL,
                                title_lab_size = 14, legend_lab = "",
                                legend_position = "right", text_size = 3,
                                bar_width = 0.8, bw_theme = TRUE,
                                palette = NULL, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "totalCountFilt")) {
    # Fezzik, tear his arms off.
    stop("filter_object must be of class totalCountFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  ## Checks for min_count as single numeric
  if (!is.null(min_count) &&
    (length(min_count) > 1 ||
      !is.numeric(min_count))) stop("min_count must be numeric of length 1")



  # Forge sensational plots ----------------------------------------------------

  plot_data <- attr(filter_object, "e_data_lcpm")
  title_default <- "Observation Density by LCPM"
  if (!is.null(min_count)) {
    biomols <- filter_object[[1]][filter_object$Total_Counts >= min_count]
    if (length(biomols) == 0) stop(
      "min_count exceeds maximum total count (",
      max(filter_object$Total_Counts), ")"
    )
    plot_data <- plot_data[plot_data[[1]] %in% biomols, ]
    subtitle <- paste0(min_count, "+ total counts per transcript")
  } else subtitle <- NULL

  xlabel <- if (is.null(x_lab)) "Log Counts per Million" else x_lab
  ylabel <- ifelse(is.null(y_lab), "Observation Density", y_lab)
  plot_title <- ifelse(
    is.null(title_lab),
    title_default,
    title_lab
  )

  # Use the ColorBrewer color
  values <- if (!is.null(palette)) {
    c(
      "Samples" = RColorBrewer::brewer.pal(5, palette)[[3]],
      "Average Density" = RColorBrewer::brewer.pal(5, palette)[[5]]
    )
  } else {
    c("Samples" = "grey", "Average Density" = "black")
  }

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = lcpm)) +
    ggplot2::geom_density(ggplot2::aes(
      group = !!dplyr::sym(colnames(plot_data)[2]),
      color = !!dplyr::sym(colnames(plot_data)[2])
    )) +
    ggplot2::geom_density(ggplot2::aes(color = "Average Density")) +
    ggplot2::scale_color_manual(
      name = legend_lab,
      values = values,
      aesthetics = c("color")
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = subtitle,
      x = xlabel, y = ylabel
    )

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Add the theme elements to the plot.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot RNAFilt Object
#'
#' For plotting an S3 object of type 'RNAFilt'
#'
#' @param x object of class RNAFilt that contains the sample
#'   identifier, library size, number of non-zero biomolecules, and proportion
#'   of non-zero biomolecules
#' @param plot_type character string, specified as "library" or "biomolecule".
#' "library" displays library size for each sample, "biomolecule" displays the
#' number of unique biomolecules with non-zero counts per sample.
#' @param min_nonzero integer or float between 0 and 1. Cut-off for number of
#' unique biomolecules with non-zero counts or as a proportion of total
#' biomolecules. Defaults to NULL.
#' @param size_library integer cut-off for sample library size (i.e. number
#' of reads). Defaults to NULL.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size An integer specifying the size of the text (number of
#'   biomolecules by sample) within the bar plot. The default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' seqfilt <- RNA_filter(omicsData = rnaseq_object)
#' plot(seqfilt)
#'
#' @rdname plot-RNAFilt
#'
#' @export
#'
plot.RNAFilt <- function(x, plot_type = "library",
                         size_library = NULL, min_nonzero = NULL,
                         interactive = FALSE,
                         x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11,
                         x_lab_angle = 90, title_lab = NULL,
                         title_lab_size = 14, legend_lab = "",
                         legend_position = "right", text_size = 3,
                         bar_width = 0.8, bw_theme = TRUE,
                         palette = NULL, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "RNAFilt")) {
    # Fezzik, tear his arms off.
    stop("filter_object must be of class RNAFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  ## Checks for plot_type as single string
  if (is.null(plot_type) ||
    length(plot_type) > 1 ||
    !(plot_type %in% c('library', 'biomolecule'))) stop(
    "plot_type must be either 'library' or 'biomolecule'"
  )

  ## Checks for size_library as single integer
  if (!is.null(size_library) &&
    (length(size_library) > 1 ||
      size_library %% 1 != 0 ||
      size_library > max(filter_object$LibrarySize)
    )
  ) stop(
    paste0(
      "size_library must be integer of length 1 less than max library size (",
      max(filter_object$LibrarySize),
      ")"
    )
  )

  ## Checks for min_nonzero as single numeric
  if (!is.null(min_nonzero)) {
    ## Length
    if (length(min_nonzero) > 1) stop("min_nonzero must be length 1")

    ## proportion or int
    if (!(min_nonzero %% 1 == 0 || (min_nonzero > 0 && min_nonzero < 1))) stop(
      "min_nonzero must be integer or numeric between 0 and 1."
    )

    ## Within appropriate bounds
    if (min_nonzero %% 1 == 0 && min_nonzero > max(filter_object$NonZero)) stop(
      paste0(
        "min_nonzero exceeds maximum number of non-zero biomolecules (",
        max(filter_object$NonZero),
        ")"
      )
    )

    if (min_nonzero %% 1 != 0 &&
      min_nonzero > max(filter_object$ProportionNonZero)) stop(
      paste0(
        "min_nonzero exceeds maximum proportion of non-zero biomolecules (",
        signif(max(filter_object$ProportionNonZero)),
        ")"
      )
    )
  }

  # Forge sensational plots ----------------------------------------------------

  ## Set labels
  xlabel <- if (is.null(x_lab)) "SampleID" else x_lab
  plot_title <- ifelse(is.null(title_lab), "Library Size by Sample", title_lab)

  if (plot_type == "library") {
    ylabel <- ifelse(is.null(y_lab), "Library Size (Total Reads)", y_lab)
  } else {
    ylabel <- ifelse(is.null(y_lab), "N Non-zero Biomolecules", y_lab)
  }

  subtitle <- ""

  if (!is.null(min_nonzero)) {
    subtitle <- paste0(
      subtitle,
      ifelse(min_nonzero %% 1 == 0,
        paste0("At least ", min_nonzero, " non-zero biomolecules"),
        paste0("At least ", min_nonzero * 100, "% non-zero biomolecules")
      )
    )
  }

  if (!is.null(size_library)) {
    subtitle <- paste0(
      subtitle,
      ifelse(subtitle == "",
        "",
        "\n"
      ),
      paste0("At least ", size_library, " reads in sample library")
    )
  }

  # Use the ColorBrewer color
  values <- if (!is.null(palette)) {
    RColorBrewer::brewer.pal(5, palette)[2, 4]
  } else {
    c(Keep = "deepskyblue3", Remove = "tomato")
  }

  temp_obj <- filter_object
  temp_obj$keep_nz <- T
  temp_obj$keep_sz <- T

  if (!is.null(min_nonzero)) {
    column_use <- ifelse(min_nonzero %% 1 == 0,
      "NonZero", "ProportionNonZero"
    )
    temp_obj$keep_nz <- temp_obj[[column_use]] >= min_nonzero
  }

  if (!is.null(size_library)) {
    temp_obj$keep_sz <- temp_obj[["LibrarySize"]] >= size_library
  }

  temp_obj$color <- as.character(Reduce("&", list(temp_obj$keep_nz, temp_obj$keep_sz)))
  temp_obj$color <- gsub("FALSE", "Remove", gsub("TRUE", "Keep", temp_obj$color))

  # Plot
  if (plot_type == "library") {
    p <- ggplot2::ggplot(
      temp_obj, ggplot2::aes(x = !!dplyr::sym(colnames(temp_obj)[1]), y = LibrarySize, fill = color)
    ) +
      ggplot2::geom_col(show.legend = F) +
      ggplot2::scale_fill_manual(
        name = legend_lab,
        values = values,
        aesthetics = c("fill")
      ) +
      ggplot2::labs(
        title = plot_title,
        subtitle = subtitle,
        x = xlabel, y = ylabel
      )

    if (!is.null(size_library)) {
      p <- p + ggplot2::geom_segment(
        y = size_library, yend = size_library,
        xend = length(unique(temp_obj[[colnames(temp_obj)[1]]])) + .5,
        x = 0, linetype = "dashed"
      )
    }
  } else {
    mt <- round(filter_object$NonZero[[1]] / filter_object$ProportionNonZero[[1]])

    p <- ggplot2::ggplot(
      temp_obj, ggplot2::aes(x = !!dplyr::sym(colnames(temp_obj)[1]), y = NonZero, fill = color)
    ) +
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(
          trans = ~ . / mt,
          name = "Proportion of all biomolecules"
        )
      ) +
      ggplot2::geom_col(show.legend = F) +
      ggplot2::scale_fill_manual(
        name = legend_lab,
        values = values,
        aesthetics = c("fill")
      ) +
      ggplot2::labs(
        title = plot_title,
        subtitle = subtitle,
        x = xlabel, y = ylabel
      )

    if (!is.null(min_nonzero)) {
      p <- p + ggplot2::geom_segment(
        y = min_nonzero, yend = min_nonzero,
        xend = length(unique(temp_obj[[colnames(temp_obj)[1]]])) + .5,
        x = 0, linetype = "dashed"
      )
    }
  }

  # Evan, make me a plot with the black and white theme. As you wish.
  if (bw_theme) p <- p + ggplot2::theme_bw()

  # Add the theme elements to the plot.
  p <- p +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_lab_size),
      axis.title.x = ggplot2::element_text(size = x_lab_size),
      axis.title.y = ggplot2::element_text(size = y_lab_size),
      axis.text.x = ggplot2::element_text(angle = x_lab_angle, hjust = 1),
      legend.position = legend_position
    )

  # Evan, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}


#' Plot imdanovaFilt Object
#'
#' For plotting an S3 object of type 'imdanovaFilt'
#'
#' @param x Object of class imdanovaFilt (also a data frame) containing
#'   the molecule identifier and number of samples in each group with
#'   non-missing values for that molecule
#' @param min_nonmiss_gtest An integer indicating the minimum number of
#'   non-missing feature values allowed per group for \code{gtest_filter}.
#'   Suggested value is 3.
#' @param min_nonmiss_anova An integer indicating the minimum number of
#'   non-missing feature values allowed per group for \code{anova_filter}.
#'   Suggested value is 2.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size integer specifying the size of the points. The default
#'   is 3.
#' @param line_size integer specifying the thickness of the line. The default
#'   is 0.75.
#' @param text_size integer specifying the size of the text (number of
#'   biomolecules per group). The default is 3.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count logical value. Indicates whether the missing value counts by
#'   sample will be displayed on the bar plot. The default is TRUE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' data(pep_object)
#' mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
#' to_filter <- imdanova_filter(omicsData = mypep)
#' plot(to_filter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)
#'
#' @rdname plot-imdanovaFilt
#'
#' @export
#'
plot.imdanovaFilt <- function(x, min_nonmiss_anova = NULL,
                              min_nonmiss_gtest = NULL, interactive = FALSE,
                              x_lab = NULL, y_lab = NULL, x_lab_size = 11,
                              y_lab_size = 11, x_lab_angle = 0,
                              title_lab = NULL, title_lab_size = 14,
                              legend_lab = NULL, legend_position = "right",
                              point_size = 3, line_size = 0.75, text_size = 3,
                              bw_theme = TRUE, palette = NULL,
                              display_count = TRUE, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "imdanovaFilt")) {
    # You used the wrong input for filter_object. Savvy?
    stop("filter_object must be of class imdanovaFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # Nab group size minimum that is not a singleton
  group_sizes <- attr(filter_object, "group_sizes")$n_group
  group_sizes_valid <- group_sizes[group_sizes > 1]
  max_x <- min(group_sizes_valid)
  obs <- as.data.frame(filter_object)[-1]

  # Count number of groups in each row of the obs data frame that have > 0, > 1,
  # ..., > max_x non-missing values. For ANOVA we need at least two groups that
  # meet this criteria unless the data are paired and there are no main effects.
  n_biomolecules_anova <- purrr::map_int(0:max_x, function(num) {
    # Number of groups greater than num
    n_greater_obs <- apply(obs, 1, function(row) sum(!(row < num)))

    # If there are no main effects we need to count differently (in this case we
    # do not need at least two groups to run ANOVA).
    if ("paired_diff" %in% names(obs)) {
      # biomolecules with 1+ group > num
      n_greater_obs <- sum(apply(obs, 1, function(row) any(!(row < num))))
    } else {
      # Biomolecules with 2+ group > num.
      sum(n_greater_obs > 1)
    }
  })

  # Count number of groups in each row of the obs data frame that have > 0, > 1,
  # ..., > max_x non-missing values. For G-Test we need at least one group that
  # meets this criteria.
  n_biomolecules_gtest <- purrr::map_int(0:max_x, function(num) {
    # biomolecules with 1+ group > num
    n_greater_obs <- sum(apply(obs, 1, function(row) any(!(row < num))))
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
    ggplot2::geom_point(
      data = plotter2,
      ggplot2::aes(
        x = Count_biomolecules,
        y = Min_obs,
        color = Statistic
      ),
      size = point_size
    )

  # Evan, display the counts on the plot. As you wish.
  if (display_count) p <- p +
    ggplot2::geom_text(
      data = plotter2,
      ggplot2::aes(
        x = Count_biomolecules,
        y = Min_obs,
        label = Count_biomolecules
      ),
      size = text_size,
      hjust = -0.5
    )

  # Evan, add the anova points to the plot. As you wish.
  p <- p +
    ggplot2::geom_point(
      data = plotter1,
      ggplot2::aes(
        x = Count_biomolecules,
        y = Min_obs,
        color = Statistic
      ),
      size = point_size
    )

  # Evan, display the counts on the plot. As you wish.
  if (display_count) p <- p +
    ggplot2::geom_text(
      data = plotter1,
      ggplot2::aes(
        x = Count_biomolecules,
        y = Min_obs,
        label = Count_biomolecules
      ),
      size = text_size,
      hjust = -0.5
    )


  # Evan, add gtest threshold to the plot. As you wish.
  if (!is.null(min_nonmiss_gtest)) {
    # Evan, add a vertical line for the gtest threshold. As you wish.
    p <- p +
      ggplot2::geom_hline(
        ggplot2::aes(
          yintercept = min_nonmiss_gtest,
          color = "G-test applied filter"
        ),
        linetype = "dashed",
        linewidth = if (is.null(line_size)) 1 else line_size
      )
  }

  # Evan, add anova threshold to the plot. As you wish.
  if (!is.null(min_nonmiss_anova)) {
    # Evan, add a vertical line for the anova threshold. As you wish.
    p <- p +
      ggplot2::geom_hline(
        ggplot2::aes(
          yintercept = min_nonmiss_anova,
          color = "ANOVA applied filter"
        ),
        linetype = "dotted",
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
            linetype = c(0, 0, 3),
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
            linetype = c(0, 2, 0, 3),
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
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  return(p)
}

#' Plot proteomicsFilt Object
#'
#' For plotting an S3 object of type 'proteomicsFilt':
#'
#' @param x object of class proteomicsFilt, which is a list with two
#'   elements. The first element is a data frame of counts for each unique
#'   peptide. The second element is a data frame with the counts for the number
#'   of peptides that map to each unique protein.
#' @param plot_type character string specifying the type of plot to be
#'   displayed. The available options are "num_peps" or "redundancy". If
#'   "num_peps" the plot is displayed that shows the counts of proteins that
#'   have a specific number of peptides mapping to them. If "redundancy" the
#'   plot showing the counts of peptides that map to a specific number of
#'   proteins is displayed.
#' @param min_num_peps an optional integer value between 1 and the maximum
#'   number of peptides that map to a protein in the data. The value specifies
#'   the minimum number of peptides that must map to a protein. Any protein with
#'   less than \code{min_num_peps} mapping to it will be returned as a protein
#'   that should be filtered. Default value is NULL.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab_pep character string used for the x-axis label for the
#'   num_peps plot. The default is NULL in which case the default x-axis label
#'   will be used.
#' @param x_lab_pro character string used for the x-axis label for the
#'   redundancy plot. The default is NULL in which case the default x-axis label
#'   will be used.
#' @param y_lab_pep character string used for the y-axis label for the
#'   num_peps plot. The default is NULL in which case the default y-axis label
#'   will be used.
#' @param y_lab_pro character string used for the y-axis label for the
#'   redundancy plot. The default is NULL in which case the default y-axis label
#'   will be used.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab_pep character string specifying the num_peps plot title.
#'   The default is NULL in which case the default title will be used.
#' @param title_lab_pro character string specifying the redundancy plot title.
#'   The default is NULL in which case the default title will be used.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size An integer specifying the size of the text (number of
#'   peptides or proteins depending on the plot) within the bar plot. The
#'   default is 3.
#' @param bar_width An integer indicating the width of the bars in the bar plot.
#'   The default is 0.8.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param display_count logical value. Indicates whether the peptide or protein counts
#'   will be displayed on the bar plot. The default is TRUE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' data(pep_object)
#' my_filter <- proteomics_filter(omicsData = pep_object)
#' plot(my_filter, min_num_peps = 3)
#' plot(my_filter, plot_type = "redundancy")
#'
#' @rdname plot-proteomicsFilt
#'
#' @export
#'
plot.proteomicsFilt <- function(x,
                                plot_type = "num_peps",
                                min_num_peps = NULL,
                                interactive = FALSE, x_lab_pep = NULL,
                                x_lab_pro = NULL, y_lab_pep = NULL,
                                y_lab_pro = NULL, x_lab_size = 11,
                                y_lab_size = 11, x_lab_angle = 0,
                                title_lab_pep = NULL, title_lab_pro = NULL,
                                title_lab_size = 14, legend_lab = NULL,
                                legend_position = "right", text_size = 3,
                                bar_width = 0.8, bw_theme = TRUE,
                                palette = NULL, display_count = TRUE,
                                ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "proteomicsFilt")) {
    # The filter object is the wrong class. Savvy?
    stop("filter_object must be of class proteomicsFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # Error Checks
  if (!is.null(min_num_peps)) {
    # check that min_num_peps is numeric and >=1 #
    if (!inherits(min_num_peps, "numeric") || min_num_peps < 1)
      stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is an integer #
    if (min_num_peps %% 1 != 0)
      stop("min_num_peps must be an integer greater than or equal to 1")
    # check that min_num_peps is of length 1 #
    if (length(min_num_peps) != 1)
      stop("min_num_peps must be of length 1")
    # check that min_num_peps is less than the total number of peptides #
    if (min_num_peps > max(filter_object$counts_by_pro$n))
      stop(paste("min_num_peps cannot be greater than the maximum number of",
        "peptides that map to a protein.",
        sep = " "
      ))
  }

  if (!plot_type %in% c("num_peps", "redundancy")) {
    stop("plot_type must be either 'num_peps' or 'redundancy'.")
  }

  # Seize unique values for peptide to protein and protein to peptide counts:
  # These bins represent the number of PROTEINS each peptide maps to. Unless
  # there are degenerate peptides this will be a vector with one value: 1.
  pep_bins <- sort(unique(filter_object$counts_by_pep$n))
  # These bins represent the number of PEPTIDES each protein maps to.
  pro_bins <- sort(unique(filter_object$counts_by_pro$n))

  # get counts of peptides that are mapped to by EXACTLY the number of
  # proteins given in pep_bins
  pep_counts <- sapply(pep_bins, function(x) {
    filter_object$counts_by_pep[filter_object$counts_by_pep$n == x, ] %>% nrow()
  })

  # get counts of proteins that are mapped to by EXACTLY the number of
  # peptides given in pro_bins
  pro_counts <- sapply(pro_bins, function(x) {
    filter_object$counts_by_pro[filter_object$counts_by_pro$n == x, ] %>% nrow()
  })

  # fill value and fill labels that will be given real values if certain
  # arguments are supplied.
  fill <- 0
  fill_format <- NULL

  # Mind numbing label making bit.
  xlabel_pep <- if (is.null(x_lab_pep))
    "Number of Proteins" else
    x_lab_pep
  ylabel_pep <- if (is.null(y_lab_pep)) "Count of Peptides" else y_lab_pep
  titleLabelPep <- if (is.null(title_lab_pep))
    "Y peptides map to exactly X proteins" else
    title_lab_pep
  xlabel_pro <- if (is.null(x_lab_pro))
    "Number of Peptides" else
    x_lab_pro
  ylabel_pro <- if (is.null(y_lab_pro)) "Count of Proteins" else y_lab_pro
  titleLabelPro <- if (is.null(title_lab_pro))
    "Y proteins mapped to by exactly X peptides" else
    title_lab_pro

  # create plotting dataframe
  pro_counts_df <- data.frame(
    counts = pro_counts,
    bins = pro_bins
  )
  pep_counts_df <- data.frame(
    counts = pep_counts,
    bins = pep_bins
  )

  # Manufacture phenomenal plots -----------------------------------------------

  # !#!#!#! p represents the protein plot (num_peps) #!#!#!#!
  # !#!#!#! q represents the pepe plot (redundancy) #!#!#!#!

  # Create the bare bones protein and peptide plots.
  p <- ggplot2::ggplot(pro_counts_df)
  q <- ggplot2::ggplot(pep_counts_df)

  # Create an object for the first default ggplot2 color.
  hideous <- grDevices::hcl(
    h = 15,
    c = 100,
    l = 65
  )

  # Check if palette is NULL or not. Hopefully it isn't so the plot will be
  # created with colors other than the super hideous default ggplot2 colors.
  if (!is.null(palette)) {
    # Create a color from the color brewer package if a palette is provided.
    # This color will be used if the plot only has one color. For example, when
    # redundancy is selected or if num_peps is selected but
    # min_num_peps is not specified.
    colas <- RColorBrewer::brewer.pal(5, palette)
  }

  # If min_num_peps is specified, add a coloring variable that shows whether a
  # peptide will be retained or dropped based on the input provided.

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
        ggplot2::aes(
          x = as.factor(bins),
          y = counts,
          fill = fill
        ),
        stat = "identity",
        width = bar_width
      )

    # Creates bar charts when there is only one group (all bars will have the
    # same color). These bar charts are created when min_num_peps is NOT
    # specified. Color palettes must be defined here because there are no groups
    # to fill by later in the script (when color palettes are defined for bar
    # charts with groups).
  } else {
    # We change bins to a factor so the x-axis tick labels are the value in bins
    # but the width of the bars and x-axis does not change according to the
    # numeric value of bins.
    p <- p +
      ggplot2::geom_bar(
        ggplot2::aes(
          x = as.factor(bins),
          y = counts
        ),
        fill = if (is.null(palette))
          hideous else
          colas[[3]],
        stat = "identity",
        width = bar_width
      )
  }

  # Evan, deal with all the issues when we ask you to make seemingly small
  # changes to functions. AS YOU WISH.
  # The q plot (redundant plot) has to be created outside the if else statement
  # above because this plot does not change based on the input to min_num_peps.
  # I tried to be smooth and leave the function mostly the same and add changes
  # to the q plot together with the p plot when min_num_peps is specified.
  # However, that turned out to be a nightmare. This attempt better work without
  # any issues. Maybe I was too wishful: This attempt better not have the same
  # completely ridiculous and unsolvable errors as the last attempt.
  q <- q +
    ggplot2::geom_bar(
      ggplot2::aes(
        x = as.factor(bins),
        y = counts
      ),
      fill = if (is.null(palette)) hideous else colas[[3]],
      stat = "identity",
      width = bar_width
    )

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
  # NOTE: This code only changes colors if not all bars will have the same color
  # in the num_peps plot. If a plot has bars with all the same color, the
  # coloring (according to whether palette was specified) was taken care of
  # previously.
  if (!is.null(palette)) {
    # Use the ColorBrewer color and create the legend title
    p <- p +
      ggplot2::scale_fill_brewer(
        palette = palette,
        name = legend_lab
      )
  } else {
    p <- p + ggplot2::scale_fill_manual(
      name = legend_lab,
      values = c(
        "dropped" = "red",
        "retained" = "green"
      )
    )
  }

  # Evan, add the protein counts to the plot. As you wish.
  if (display_count) p <- p + ggplot2::geom_text(
    data = pro_counts_df,
    ggplot2::aes(
      x = as.factor(bins),
      y = counts,
      label = counts
    ),
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

  # Evan, just display the num_peps plot. As you wish.
  if (plot_type == "num_peps") {
    # Evan, make me an interactive num_peps plot. As you wish.
    if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

    return(p)

    # Evan, just display the redundancy plot. As you wish.
  } else if (plot_type == "redundancy") {
    # Evan, make me an interactive redundancy plot. As you wish.
    if (interactive && requirePlotly()) q <- plotly::ggplotly(q)

    return(q)
  }
}

#' Plot rmdFilt Object
#'
#' For plotting an S3 object of type 'rmdFilt'
#'
#' @param x object of class rmdFilt created via \code{\link{rmd_filter}}
#' @param pvalue_threshold numeric threshold for the Robust Mahalanobis Distance (RMD)
#'   p-value. If \code{sampleID} is NULL (see \code{sampleID} below), a
#'   horizontal line is plotted at the RMD value that corresponds with the
#'   threshold, and all samples above the line have a p-value below the
#'   threshold. If \code{sampleID} is not NULL, \code{pvalue_threshold} will do
#'   nothing. Default value is NULL.
#' @param sampleID character vector specifying the sample names to be plotted.
#'   If specified, the plot function produces a boxplot instead of a
#'   scatterplot. A point, colored by sample, will be placed on each boxplot for
#'   that sample's value for the given metric. The default value is NULL.
#' @param order_by character string specifying a variable by which to order
#'   the samples in the plot. This variable must be found in the column names of
#'   fdata from the rmdFilt object. If \code{order_by} is "Group", the plot will
#'   be ordered by the group variable from the group_designation function. If
#'   NULL (default), the samples will be displayed in the order in which they
#'   first appear.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 90.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param point_size An integer specifying the size of the points. The default
#'   is 3.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' mymetab <- group_designation(omicsData = mymetab, main_effects = "Phenotype")
#' rmd_results <- rmd_filter(omicsData = mymetab, metrics = c("MAD", "Skewness", "Correlation"))
#' plot(rmd_results, pvalue_threshold = 0.0001, order_by = "Phenotype")
#'
#' @rdname plot-rmdFilt
#'
#' @export
#'
plot.rmdFilt <- function(x, pvalue_threshold = NULL, sampleID = NULL,
                         order_by = NULL, interactive = FALSE, x_lab = NULL,
                         y_lab = NULL, x_lab_size = 11, y_lab_size = 11,
                         x_lab_angle = 90, title_lab = NULL,
                         title_lab_size = 14, legend_lab = NULL,
                         legend_position = "right", point_size = 3,
                         bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Have a looksie at the class of the filter object.
  if (!inherits(filter_object, "rmdFilt")) {
    # The filter object is the wrong class. Savvy?
    stop("filter_object must be of class rmdFilt")
  }

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # Check if order_by is present. If it is additional checks need to be
  # performed to make sure it is a valid input.
  if (!is.null(order_by)) {
    # Make sure order_by is either "Group" or is in f_data.
    if (order_by != "Group" &&
      !(order_by %in% names(attr(filter_object, "fdata")))) {
      # Your plot is only mostly dead! There is still a slim hope that you can
      # save it.
      stop(paste("order_by must either be 'Group' or the name of a column",
        "from f_data.",
        sep = " "
      ))
    }
  }

  # Have a looksie at the sampleID argument.
  if (!is.null(sampleID)) {
    # Make sure the sample IDs provided actually exist.
    if (!all(sampleID %in% attr(filter_object, "sample_names")))
      stop("The sample IDs provided do not match the sample IDs in the data.")
  }

  samp_id <- names(attr(filter_object, "group_DF"))[1]
  metrics <- attributes(filter_object)$metrics

  # determine how to melt based on the number of main effects
  group_df <- attributes(filter_object)$group_DF

  # determine the main effects, then melt the df #
  if (ncol(group_df) == 2) {
    main_eff_names <- "Group"
    dfmelt <- filter_object %>%
      tidyr::pivot_longer(
        -dplyr::all_of(c(!!samp_id, !!main_eff_names)),
        names_to = "variable",
        cols_vary = "slowest"
      ) %>% data.frame
  } else if (ncol(group_df) > 2) {
    ## put main effect with more levels first, for plotting aesthetics ##
    temp <- droplevels(group_df[, -(1:2)])
    numlevels <- sapply(1:2, function(j) length(levels(as.factor(temp[, j]))))
    main_eff_names <- names(temp)[order(numlevels, decreasing = TRUE)]
    dfmelt <- tidyr::pivot_longer(
      filter_object[, -2],
      -dplyr::all_of(c(!!samp_id, !!main_eff_names)),
      names_to = "variable",
      cols_vary = "slowest"
    ) %>% data.frame
  }

  # Data frame that has information to create rmd box plots.
  dfsub <- dfmelt[dfmelt$variable %in% metrics, ]

  # legend labels #
  ## function for shortening long main effect names ##
  abbrev_fun <- function(string) {
    if (nchar(string) > 25) string = paste0(substr(string, 1, 23), "...")
    return(string)
  }

  # Shorten the names of the main effects (if they are over 25 characters). This
  # vector will have more than one element if there is more than one main
  # effect.
  display_names <- sapply(main_eff_names, abbrev_fun)

  # Create the title for the legend that describes differences by color.
  legend_title_color <- ifelse(is.null(legend_lab),
    display_names[1],
    legend_lab[1]
  )

  # Initialize the title for the legend that describes differences by point
  # shape to NULL. It will only be used if there is more than one main effect.
  legend_title_shape <- NULL

  # Use the legend label (if it is provided in the input) for the legend
  # describing the differences by point shape. The user must specify a legend
  # label with more than one element. It is not created within the function.
  if (length(display_names) > 1) {
    if (length(legend_lab) > 1) {
      legend_title_shape <- legend_lab[2]
    } else {
      legend_title_shape <- display_names[2]
    }
  }

  # Assemble captivating plots -------------------------------------------------

  # Order the samples according to the order_argument. If order_by is not
  # provided order the samples according to the order they appear in the filter
  # object.
  if (!is.null(order_by)) {
    if (order_by == "Group") {
      # Order the samples by the Group column from group_DF.
      filter_object[[samp_id]] <- factor(
        filter_object[[samp_id]],
        levels = filter_object[[samp_id]][order(attr(
          filter_object,
          "group_DF"
        )$Group)],
        ordered = TRUE
      )
    } else {
      # Order the samples by the specified column from f_data.
      filter_object[[samp_id]] <- factor(
        filter_object[[samp_id]],
        levels = filter_object[[samp_id]][order(attr(
          filter_object,
          "fdata"
        )[[order_by]])],
        ordered = TRUE
      )
    }
  } else {
    # Order the data in the order the samples appear in the samp_id column.
    filter_object[[samp_id]] <- factor(
      filter_object[[samp_id]],
      levels = unique(filter_object[[samp_id]]),
    )
  }

  # Beautiful box plots ---------------

  # Check is sampleID is NULL. Generate box plots if it is not NULL.
  if (!is.null(sampleID)) {
    levels(dfsub$variable)[levels(dfsub$variable) == "Fraction_Missing"] <-
      "Prop_missing"

    plot_title <- ifelse(is.null(title_lab),
      "Summary of RMD metrics used",
      title_lab
    )
    xlabel <- ifelse(is.null(x_lab), " ", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Value", y_lab)

    p <- ggplot2::ggplot(dfsub) +
      ggplot2::geom_boxplot(ggplot2::aes(
        x = rep(1, length(value)),
        y = value
      )) +
      ggplot2::facet_wrap(~variable,
        scales = "free",
        ncol = length(metrics)
      ) +
      ggplot2::geom_point(
        data = dfsub[dfsub[, samp_id] %in% sampleID, ],
        ggplot2::aes(
          x = rep(1, length(value)),
          y = value,
          color = !!dplyr::sym(samp_id)
        ),
        size = point_size
      ) +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title)

    # Stunning scatter plots: no threshold ---------------
  } else if (is.null(pvalue_threshold)) {
    # Start scatter plot skeleton when there is no p-value threshold.
    p <- ggplot2::ggplot(filter_object)

    # Start plot when there is only one main effect.
    if (length(main_eff_names) == 1) {
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(
            x = !!dplyr::sym(samp_id),
            y = Log2.md,
            color = !!dplyr::sym(main_eff_names[1])
          ),
          size = point_size
        )


      # Start plot when there are two main effects.
    } else {
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(
            x = !!dplyr::sym(samp_id),
            y = Log2.md,
            color = !!dplyr::sym(main_eff_names[1]),
            shape = !!dplyr::sym(main_eff_names[2])
          ),
          size = point_size
        )
    }

    plot_title <- ifelse(is.null(title_lab),
      "Sample Outlier Results",
      title_lab
    )
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)

    p <- p +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::labs(
        color = legend_title_color,
        shape = legend_title_shape
      )

    # Stunning scatter plots: with threshold ---------------
  } else {
    # Determine which samples fall below the p-value threshold. This is a
    # logical vector that will be used to find pairs that are split (one sample
    # is filtered and the other is not) and to change the transparency of the
    # point in the plot.
    goodies_alpha <- filter_object$pvalue < pvalue_threshold

    # If data are paired make sure no sample is left behind!
    if (!is.null(attr(attr(filter_object, "group_DF"), "pair_id"))) {
      # Find the corresponding sample(s) if one sample from a pair is below the
      # p-value threshold and the other sample is not.
      if (sum(goodies_alpha) > 0) {
        # Snatch the sample name and pair name as they will be used in multiple
        # places. It will save some typing. ... However, all the typing I just
        # saved has probably been undone by writing this comment.
        sample_name <- attr(filter_object, "fdata_cname")
        pair_name <- attr(attr(filter_object, "group_DF"), "pair_id")

        # Grab the names of filtered samples.
        filtad <- as.character(filter_object[goodies_alpha, sample_name])

        # !#!#!#!#!#!#!#!#!#!
        # The following code checks if the samples in a pair will be split. For
        # example, one sample in a pair will be filtered and another sample in
        # the pair will not be filtered. If a pair is split remove ALL samples
        # associated with that pair.
        # !#!#!#!#!#!#!#!#!#!

        # Snag the associated pair IDs for the samples that will be filtered.
        filtad_pairs <- attr(filter_object, "fdata") %>%
          dplyr::filter(!!dplyr::sym(sample_name) %in% filtad) %>%
          dplyr::pull(!!dplyr::sym(pair_name))

        # Go back to f_data and nab all the sample names corresponding to the
        # pair IDs associated with the original samples that were selected for
        # removal. These sample names will be used to change the point shape.
        # The points that will be different are the ones corresponding to the
        # samples that belong to pair where only one sample falls below the
        # threshold.
        filtad_too <- attr(filter_object, "fdata") %>%
          dplyr::filter(!!dplyr::sym(pair_name) %in% filtad_pairs) %>%
          dplyr::pull(!!dplyr::sym(sample_name)) %>%
          as.character()

        # Update the goodies_alpha vector to reflect the additional samples that
        # will be filtered.
        goodies_alpha <- filter_object[, sample_name] %in% filtad_too

        # Create a vector of samples who are guilty by association. (They do not
        # fall below the threshold but their partners do.)
        condemned <- setdiff(filtad_too, filtad)

        # Create a logical vector that will determine the point shape.
        goodies_pch <- filter_object[, sample_name] %in% condemned

      } else goodies_pch <- rep(FALSE, nrow(filter_object))

    } else {
      # The data are not paired so all points should be a solid circle.
      goodies_pch <- rep(FALSE, nrow(filter_object))
    }

    # get y-intercept for line
    df <- attributes(filter_object)$df
    yint <- log(qchisq(1 - pvalue_threshold, df = df), base = 2)

    # make title
    if (is.null(title_lab)) {
      plot_title <- "Sample Outlier Results"
      subtitle <- paste("p-value threshold = ", pvalue_threshold)
    } else {
      plot_title <- title_lab
      subtitle <- ggplot2::waiver()
    }
    xlabel <- ifelse(is.null(x_lab), "Samples", x_lab)
    ylabel <- ifelse(is.null(y_lab), "log2(Robust Mahalanobis Distance)", y_lab)

    # Start scatter plot skeleton when a p-value threshold is specified.
    p <- ggplot2::ggplot(filter_object)

    # Start scatter plot skeleton when a p-value threshold has been provided.
    if (length(main_eff_names) == 1) {
      # we need separate plots for if we have pair id or not

      if (!is.null(attr(attr(filter_object, "group_DF"), "pair_id"))) {
        # if !null then we need to include pair fill information
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names),
              # Add a fill layer that will not affect how the plot looks. This
              # layer is used to create a legend when one sample from a pair is
              # an outlier but the other sample belonging to the pair is not.
              fill = "Removed because paired with outlier"
            ),
            alpha = ifelse(goodies_alpha, 1, 0.5),
            size = point_size,
            shape = ifelse(goodies_pch, 15, 16)
          )
      } else {
        # if null pair id then we do not need to account for fill
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names)
            ),
            alpha = ifelse(goodies_alpha, 1, 0.5),
            size = point_size,
            shape = ifelse(goodies_pch, 15, 16)
          )
      }
    } else {
      if (!is.null(attr(attr(filter_object, "group_DF"), "pair_id"))) {
        # if !null then we need to include pair fill information
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names[1]),
              shape = !!dplyr::sym(main_eff_names[2]),
              # Add a fill layer that will not affect how the plot looks. This
              # layer is used to create a legend when one sample from a pair is
              # an outlier but the other sample belonging to the pair is not.
              fill = "Removed because paired with outlier"
            ),
            alpha = ifelse(goodies_alpha, 1, 0.5),
            # Make guilty-by-association points really small because we will
            # plot a square point over them with the next lines of code. We
            # don't want the edges of a circle or triangle peeking out from
            # behind a square. That would make a hideous and confusing plot if
            # that happened.
            size = ifelse(goodies_pch, 0, point_size)
          ) +
          # Add another layer of points with guilty-by-association samples
          # plotted as a square.
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names[1])
            ),
            size = ifelse(goodies_pch, point_size, 0),
            shape = ifelse(goodies_pch, 15, 16)
          )
      } else {
        # if null pair id then we do not need to account for fill
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names[1]),
              shape = !!dplyr::sym(main_eff_names[2])
            ),
            alpha = ifelse(goodies_alpha, 1, 0.5),
            # Make guilty-by-association points really small because we will
            # plot a square point over them with the next lines of code. We
            # don't want the edges of a circle or triangle peeking out from
            # behind a square. That would make a hideous and confusing plot if
            # that happened.
            size = ifelse(goodies_pch, 0, point_size)
          ) +
          # Add another layer of points with guilty-by-association samples
          # plotted as a square.
          ggplot2::geom_point(
            ggplot2::aes(
              x = !!dplyr::sym(samp_id),
              y = Log2.md,
              col = !!dplyr::sym(main_eff_names[1])
            ),
            size = ifelse(goodies_pch, point_size, 0),
            shape = ifelse(goodies_pch, 15, 16)
          )
      }
    }

    # Add title, axis labels, and other crap.
    p <- p +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title,
        subtitle = subtitle
      ) +
      ggplot2::geom_hline(yintercept = yint) +
      ggplot2::guides(
        color = ggplot2::guide_legend(ncol = 1),
        shape = ggplot2::guide_legend(ncol = 1)
      ) +
      ggplot2::labs(
        color = legend_title_color,
        shape = legend_title_shape
      )

    # Add a manual legend that describes the shape of the points that are
    # removed because they belong to a pair with one sample being an outlier.
    if (sum(goodies_pch > 0)) {
      p <- p +
        ggplot2::scale_fill_manual(
          name = NULL,
          values = 1,
          breaks = "Removed because paired with outlier",
          guide = ggplot2::guide_legend(
            override.aes = list(
              linetype = 0,
              shape = 15,
              color = "black"
            )
          )
        )
    }
  }

  # Farm boy, make me a plot with beautiful colors. As you wish.
  if (!is.null(palette)) p <- p +
    ggplot2::scale_color_brewer(
      name = legend_title_color,
      palette = palette
    )

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
      ggplot2::scale_x_discrete(labels = attr(filter_object, "VizSampNames"))

    # We only want to change the legend if the box plots are created.
    if (!is.null(sampleID)) {
      # Nab the indices of the samples that will be highlighted in the box
      # plots.
      idx <- which(attr(filter_object, "sample_names") %in% sampleID)

      # Change the names in the legend of the box plots.
      p <- p +
        ggplot2::scale_color_hue(
          labels = attr(filter_object, "VizSampNames")[idx]
        )
    }
  }

  # Farm boy, remove the default legend(s) if the data are paired and there are
  # no main effects. As you wish.
  # If there are no main effects do not print a legend.
  if ("paired_diff" %in% attr(filter_object, "group_DF")$Group) {
    # Remove the color legend. This legend is automatically created because we
    # specified a variable to color by in aes().
    p <- p +
      ggplot2::guides(color = "none")
  }

  # Farm boy, make the plot interactive. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  # Farm boy, return the plot so the entire world can enjoy it. As you wish.
  return(p)
}

#' Plot cvFilt Object
#'
#' For plotting an S3 object of type 'cvFilt'
#'
#' @param x object of class cvFilt generated via
#'   \code{\link{cv_filter}}
#' @param cv_threshold numeric value for cv threshold above which to remove the
#'   corresponding biomolecules
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param log_scale logical value. Indicates whether to use a log2 transformed x-axis.
#'   The default is TRUE.
#' @param n_bins integer value specifying the number of bins to draw in the
#'   histogram.  The default is 30.
#' @param n_breaks integer value specifying the number of breaks to use. You
#'   may get less breaks if rounding causes certain values to become non-unique.
#'   The default is 15.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' data(pep_object)
#' mypep <- group_designation(
#'   omicsData = pep_object,
#'   main_effects = "Phenotype"
#' )
#'
#' cvfilt <- cv_filter(omicsData = mypep)
#'
#' plot(cvfilt, cv_threshold = 20)
#' plot(cvfilt, cv_threshold = 10, log_scale = FALSE)
#'
#' @rdname plot-cvFilt
#'
#' @export
#'
plot.cvFilt <- function(x, cv_threshold = NULL,
                        interactive = FALSE, x_lab = NULL, y_lab = NULL,
                        x_lab_size = 11, y_lab_size = 11, x_lab_angle = 0,
                        title_lab = NULL, title_lab_size = 14,
                        legend_lab = NULL, legend_position = "right",
                        log_scale = TRUE, n_breaks = 15, n_bins = 30,
                        bw_theme = TRUE, palette = NULL, ...) {
  filter_object <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # checks for cv_threshold if not null
  if (!is.null(cv_threshold)) {
    # check that cv_threshold is numeric
    if (!is.numeric(cv_threshold))
      stop("cv_threshold must be numeric of length 1")
    # chack that cv_threshold is of length 1
    if (length(cv_threshold) > 1)
      stop("cv_threshold must be numeric of length 1")
    # check that cv_threshold is more than 1 and less than max CV value
    if (cv_threshold <= 1 || cv_threshold >= max(filter_object$CV, na.rm = TRUE))
      stop("cv_threshold must be greater than 1 and less than the max CV value")
  }

  if (!is.logical(log_scale)) stop("log_scale must be logical: TRUE or FALSE")
  if (!(n_breaks %% 1 == 0)) stop("n_breaks must be integer valued")
  if (!(n_bins %% 1 == 0)) stop("n_bins must be integer valued")

  # Make sure palette is one of the RColorBrewer options if it is not NULL.
  if (!is.null(palette)) {
    if (!(palette %in% c(
      "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds", "RdPu",
      "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd", "Oranges",
      "Greys", "Greens", "GnBu", "BuPu", "BuGn", "Blues",
      "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired",
      "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu",
      "RdGy", "RdBu", "PuOr", "PRGn", "PiYG", "BrBG"
    ))) {
      # INCONCEIVABLE!!!
      stop("palette must be an RColorBrewer palette")
    }
  }

  # plotting object
  new_object <- filter_object[!is.na(filter_object$CV), ]
  max_x_val <- attributes(filter_object)$max_x_val

  # Check if any CV values are zero. If there are any zero values the log_scale
  # input must be changed to FALSE.
  if (any(new_object$CV == 0)) {
    log_scale <- FALSE

    message(paste("Because there are CV values = 0 the x-axis cannot be",
      "converted to the log2 scale.",
      "The original scale will be used.",
      sep = " "
    ))
  }

  # labels
  plot_title <- ifelse(is.null(title_lab),
    "Coefficient of Variation (CV)",
    title_lab
  )
  ylabel <- ifelse(is.null(y_lab), "Count", y_lab)

  # Label x-axis according to whether group_designation() was run.
  if (is.null(x_lab)) {
    xlabel <- if (attributes(filter_object)$pooled) "Pooled CV" else "CV"
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
    while (2^i < max(new_object$CV, na.rm = TRUE)) {
      breaks <- c(breaks, 2^i)
      i <- i + step
    }
    # rounding for plot purposes
    breaks <- round(breaks, 2)
  } else {
    breaks <- scales::pretty_breaks(n = n_breaks)
    trans <- "identity"
  }

  # change default title and draw a vertical line if cv_thresh specified
  if (!is.null(cv_threshold)) {
    if (is.null(title_lab)) {
      plot_title <- bquote(paste(
        "Coefficient of Variation (CV):  ",
        italic(paste(
          "CV Threshold = ",
          .(cv_threshold)
        )),
        ""
      ))
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
        fill = colas[[3]]
      )

    # Runs when palette is NULL.
  } else {
    # Farm boy, make the histogram with default colors. As you wish.
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(x = CV),
        bins = n_bins
      )
  }

  # Farm boy, add the cutoff and make the histogram pretty. As you wish.
  p <- p +
    cutoff +
    ggplot2::scale_x_continuous(
      breaks = breaks,
      trans = trans
    ) +
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
  return(p)
}

#' Plot customFilt Object
#'
#' Currently plotting a customFilt object is not supported
#'
#' @param x An object of class customFilt.
#' @param ... further arguments passed to or from other methods.
#'
#' @return No return value, implemented to provide information to user.
#' 
#' @rdname plot-customFilt
#'
#' @export
#'
plot.customFilt <- function(x, ...) {
  message(paste("There is no plot method for objects of class 'customFilt'.",
    "See summary.customFilt instead.",
    sep = " "
  ))
}

#' Plot normRes Object
#'
#' For plotting an S3 object of type 'normRes'
#'
#' @param x normRes object created by the normalize_global function
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mymetab <- edata_transform(
#'   omicsData = metab_object,
#'   data_scale = "log2"
#' )
#' mymetab <- group_designation(
#'   omicsData = mymetab,
#'   main_effects = "Phenotype"
#' )
#' norm_object <- normalize_global(
#'   omicsData = mymetab,
#'   subset_fn = "all",
#'   norm_fn = "median"
#' )
#' plot(norm_object, order_by = "Phenotype", color_by = "Phenotype")
#'
#' @rdname plot-normRes
#'
#' @export
#'
plot.normRes <- function(x, order_by = NULL, color_by = NULL,
                         facet_by = NULL, facet_cols = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         ylimit = NULL, bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, ...) {
  normRes_obj <- x

  # Preliminaries --------------------------------------------------------------

  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

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
          each = nrow(omicsNorm$e_data)
        )
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
          each = nrow(omicsNorm$e_data)
        ) /
          rep(normRes_obj$parameters$normalization$scale,
            each = nrow(omicsNorm$e_data)
          )
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
  p_raw <- plot_omicsData(
    omicsData = omicsRaw, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab,
    title_lab_size = title_lab_size,
    legend_lab = legend_lab,
    legend_position = legend_position, ylimit = ylimit,
    bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )

  # Farm boy, make me a plot with a normalized object. As you wish.
  p_norm <- plot_omicsData(
    omicsData = omicsNorm, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab,
    x_lab_size = x_lab_size, y_lab_size = y_lab_size,
    x_lab_angle = x_lab_angle, title_lab = title_lab,
    title_lab_size = title_lab_size,
    legend_lab = legend_lab,
    legend_position = legend_position, ylimit = ylimit,
    bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )

  # Farm boy, combine the plots and send them into the world. As you wish.
  if (interactive && requirePlotly()) {
    # Return the interactive plots side-by-side.
    plotly::subplot(p_raw, p_norm, nrows = 1)
  } else {
    # Return the regular plots side-by-side.
    p_raw + p_norm
  }
}

#' Plot isobaricpepData Object
#'
#' For plotting isobaricpepData S3 objects
#'
#' @param x An isobaricpepData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
#' plot(myiso)
#'
#' @rdname plot-isobaricpepData
#'
#' @export
#'
plot.isobaricpepData <- function(x, order_by = NULL, color_by = NULL,
                                 facet_by = NULL, facet_cols = NULL,
                                 interactive = FALSE, x_lab = NULL,
                                 y_lab = NULL, x_lab_size = 11,
                                 y_lab_size = 11, x_lab_angle = 90,
                                 title_lab = NULL, title_lab_size = 14,
                                 legend_lab = NULL, legend_position = "right",
                                 ylimit = NULL, bw_theme = TRUE,
                                 palette = NULL, use_VizSampNames = FALSE,
                                 ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with an isobaricpepData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

#' Plot lipidData Object
#'
#' For plotting lipidData S3 objects
#'
#' @param x lipidData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
#' plot(mylipid, order_by = "Virus", color_by = "Virus")
#'
#' @rdname plot-lipidData
#'
#' @export
#'
plot.lipidData <- function(x, order_by = NULL, color_by = NULL,
                           facet_by = NULL, facet_cols = NULL,
                           interactive = FALSE, x_lab = NULL, y_lab = NULL,
                           x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                           title_lab = NULL, title_lab_size = 14,
                           legend_lab = NULL, legend_position = "right",
                           ylimit = NULL, bw_theme = TRUE, palette = NULL,
                           use_VizSampNames = FALSE, ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a lipidData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

#' Plot metabData Object
#'
#' For plotting metabData S3 objects
#'
#' @param x metabData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' plot(mymetab, order_by = "Phenotype", color_by = "Phenotype")
#'
#' @rdname plot-metabData
#'
#' @export
#'
plot.metabData <- function(x, order_by = NULL, color_by = NULL,
                           facet_by = NULL, facet_cols = NULL,
                           interactive = FALSE, x_lab = NULL, y_lab = NULL,
                           x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                           title_lab = NULL, title_lab_size = 14,
                           legend_lab = NULL, legend_position = "right",
                           ylimit = NULL, bw_theme = TRUE, palette = NULL,
                           use_VizSampNames = FALSE, ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a metabData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

#' Plot nmrData Object
#'
#' For plotting nmrData S3 objects
#'
#' @param x nmrData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit A numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
#' plot(mynmr)
#'
#' @rdname plot-nmrData
#'
#' @export
#'
plot.nmrData <- function(x, order_by = NULL, color_by = NULL,
                         facet_by = NULL, facet_cols = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         ylimit = NULL, bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a nmrData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

#' Plot seqData Object
#'
#' For plotting seqData S3 objects
#'
#' @param x seqData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param transformation character string. String of length 1 defining a
#'   transformation for visualizing count data. Valid options are 'lcpm',
#'   'upper', and 'median'. 'lcpm' - For each column: scale column intensities by
#'   (total column sum/1 million), then log2 transform. 'median' - For each
#'   column: scale column intensities by median column intensities, then
#'   back-transform to original scale. 'upper' - For each column: scale column
#'   intensities by 75th quantile column intensities, then back-transform to
#'   original scale. For 'median' and 'upper' options, all zeros are converted
#'   to NAs.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' plot(rnaseq_object, transformation = "lcpm")
#'
#' @rdname plot-seqData
#'
#' @export
#'
plot.seqData <- function(x, order_by = NULL, color_by = NULL,
                         facet_by = NULL, facet_cols = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         ylimit = NULL, bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, transformation = NULL,
                         ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a seqData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames,
    transformation = transformation
  )
}

#' Plot pepData Object
#'
#' For plotting pepData S3 objects
#'
#' @param x pepData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols An optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' data(pep_object)
#' mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' plot(mypep, order_by = "Phenotype", color_by = "Phenotype")
#'
#' @rdname plot-pepData
#'
#' @export
#'
plot.pepData <- function(x, order_by = NULL, color_by = NULL,
                         facet_by = NULL, facet_cols = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         ylimit = NULL, bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a pepData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

#' Plot proData Object
#'
#' For plotting proData S3 objects
#'
#' @param x proData object
#' @param order_by character string specifying the column name of f_data by
#'   which to order the boxplots. If \code{order_by} is "Group", the boxplots
#'   will be ordered by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will be displayed in the order they appear
#'   in the data.
#' @param color_by character string specifying the column name of f_data by
#'   which to color the boxplots. If \code{color_by} is "Group", the boxplots
#'   will be colored by the group variable from the group_designation function.
#'   If NULL (default), the boxplots will have one default color.
#' @param facet_by character string specifying the column name of f_data with
#'   which to create a facet plot. Default value is NULL.
#' @param facet_cols optional integer specifying the number of columns to
#'   show in the facet plot.
#' @param interactive logical value. If TRUE produces an interactive plot.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#'   The default is 0.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param ylimit numeric vector of length 2 specifying y-axis lower and upper
#'   limits.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white theme.
#' @param palette character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param use_VizSampNames logical value. Indicates whether to use custom sample
#'   names. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @examplesIf requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' plot(pro_object, order_by = "Phenotype", color_by = "Phenotype")
#'
#' @rdname plot-proData
#'
#' @export
#'
plot.proData <- function(x, order_by = NULL, color_by = NULL,
                         facet_by = NULL, facet_cols = NULL,
                         interactive = FALSE, x_lab = NULL, y_lab = NULL,
                         x_lab_size = 11, y_lab_size = 11, x_lab_angle = 90,
                         title_lab = NULL, title_lab_size = 14,
                         legend_lab = NULL, legend_position = "right",
                         ylimit = NULL, bw_theme = TRUE, palette = NULL,
                         use_VizSampNames = FALSE, ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, make me a plot with a proData object. As you wish.
  plot_omicsData(
    omicsData = x, order_by = order_by,
    color_by = color_by, facet_by = facet_by,
    facet_cols = facet_cols, interactive = interactive,
    x_lab = x_lab, y_lab = y_lab, x_lab_size = x_lab_size,
    y_lab_size = y_lab_size, x_lab_angle = x_lab_angle,
    title_lab = title_lab, title_lab_size = title_lab_size,
    legend_lab = legend_lab, legend_position = legend_position,
    ylimit = ylimit, bw_theme = bw_theme, palette = palette,
    use_VizSampNames = use_VizSampNames
  )
}

# The function that does all the heavy lifting for the isobaricpepData,
# lipidData, metabData, nmrData, pepData, proData, and normRes plot methods.
plot_omicsData <- function(omicsData, order_by, color_by, facet_by, facet_cols,
                           interactive, x_lab, y_lab, x_lab_size, y_lab_size,
                           x_lab_angle, title_lab, title_lab_size, legend_lab,
                           legend_position, ylimit, bw_theme, palette,
                           use_VizSampNames, transformation = NULL) {
  # Preliminaries --------------------------------------------------------------

  # Keeping the user honest ---------------

  # Farm boy, make sure the data is the correct class. As you wish.
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c(
    "pepData", "proData", "metabData",
    "lipidData", "nmrData", "seqData"
  ))) {
    # INCONCEIVABLE!!!
    stop(paste("omicsData must be of class 'isobaricpepData', 'lipidData'",
      "'metabData', 'nmrData', 'pepData', 'proData', or 'seqData'.",
      sep = " "
    ))
  }

  if (inherits(omicsData, "seqData")) {
    if (!is.null(transformation)) {
      if (!(transformation %in% c('lcpm', 'upper', 'median'))) {
        # Tell the user that the input to data_scale is an abomination!
        stop(paste(transformation, "is not a valid option for 'transformation'.",
          "Refer to ?seqData.plot for specific seqData options.",
          sep = " "
        ))
      }
    } else {
      warning(paste0(
        "Using transformation argument for plotting on counts",
        " is recommended for seqData visualization. See ?seqData.plot for details."
      ))
    }
  }

  if (!inherits(omicsData, "seqData") && !is.null(transformation)) {
    warning("Transformation agrument is only applicable for seqData -- this argument is otherwise ignored.")
  }


  if (!is.null(order_by)) {
    if (!is.character(order_by) || length(order_by) > 1)
      stop("order_by must be a character vector of length 1")
  }

  if (!is.null(color_by)) {
    if (!is.character(color_by) || length(color_by) > 1)
      stop("color_by must be a character vector of length 1")
  }

  if (!is.null(facet_by)) {
    if (!is.character(facet_by) || length(facet_by) > 1)
      stop("facet_by must be a character vector of length 1")
  }

  if (!is.null(facet_cols)) {
    if (is.null(facet_by))
      stop("facet_by cannot be NULL when facet_cols is specified")
    if (length(facet_cols) > 1)
      stop("facet_cols must be of length 1")
    if (!is.numeric(facet_cols))
      stop("facet_cols must be an integer greater than zero")
    if (facet_cols %% 1 != 0 || facet_cols <= 0)
      stop("facet_cols must be an integer greater than zero")
  }

  if (!is.null(ylimit)) {
    if (!is.numeric(ylimit) || length(ylimit) != 2)
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
      sep = " "
    ))
  }

  # Label crap ---------------

  datatype_text <- switch(class(omicsData)[1],
    isobaricpepData = "Isobaric Peptide ",
    pepData = "Peptide ",
    proData = "Protein ",
    lipidData = "Lipid ",
    metabData = "Metabolite ",
    nmrData = "NMR ",
    seqData = "Transcript "
  )

  norm_info <- get_data_norm(omicsData)
  norm_text <- ifelse(norm_info, "Normalized ", "Un-Normalized ")

  ref_info <- if (inherits(omicsData, c("isobaricpepData", "nmrData"))) {
    infos <- c("isobaric_info", "nmr_info")
    res <- sapply(infos, function(info) {
      norm <- attr(omicsData, info)$norm_info$is_normalized
      !is.null(norm) && norm
    })
    any(res)
  } else FALSE
  ref_text <- if(ref_info) "Reference Standardized " else NULL

  maintitle <- paste0(
    "Boxplots of ", ref_text, norm_text, datatype_text, "Data"
  )

  ######## Is this right? This should be separate from statistical normalization? #########
  # Farm boy, make me a title depending on data type and norm_info. As you wish.
  # if (inherits(omicsData, "isobaricpepData")) {
  #   maintitle <- if (attr(omicsData, "isobaric_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized Isobaric Peptide Data" else
  #       "Boxplots of Un-Normalized Isobaric Peptide Data"
  # } else if(inherits(omicsData, "pepData")){
  #   maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized Peptide Data" else
  #       "Boxplots of Un-Normalized Peptide Data"
  # } else if(inherits(omicsData, "proData")){
  #   maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized Protein Data" else
  #       "Boxplots of Un-Normalized Protein Data"
  # } else if(inherits(omicsData, "lipidData")){
  #   maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized Lipid Data" else
  #       "Boxplots of Un-Normalized Lipid Data"
  # } else if(inherits(omicsData, "metabData")){
  #   maintitle <- if (attr(omicsData, "data_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized Metabolite Data" else
  #       "Boxplots of Un-Normalized Metabolite Data"
  # } else if(inherits(omicsData, "nmrData")){
  #   maintitle <- if (attr(omicsData, "nmr_info")$norm_info$is_normalized)
  #     "Boxplots of Normalized NMR Data" else
  #       "Boxplots of Un-Normalized NMR Data"
  # }

  # Farm boy, create an  plot subtitle object. As you wish.
  subtitle <- NULL

  # Farm boy, fashion plot, axis, and legend title objects. As you wish.
  title <- if (is.null(title_lab)) maintitle else title_lab
  xlabel <- if (is.null(x_lab)) "Sample" else x_lab
  ylabel <- if (is.null(y_lab)) {
    # Abundance based
    if (get_data_scale(omicsData) == "abundance") {
      out <- "Abundance"
    } else if (get_data_scale(omicsData) %in% c("log", "log2", "log10")) {
      out <- paste(get_data_scale(omicsData), "Abundance", sep = " ")
    } else if (get_data_scale(omicsData) == "counts") {
      if (is.null(transformation)) {
        out <- "Counts"
      } else {
        out <- switch(transformation,
          lcpm = "Log Counts per Million",
          upper = "Upper-quantile transformed Counts",
          median = "Median Counts"
        )
      }
    }
    out
  } else y_lab
  legend_title <- if (is.null(legend_lab)) color_by else legend_lab

  # Data crap ---------------

  # Farm boy, melt the data for me. As you wish.
  e_data_cname <- get_edata_cname(omicsData)

  temp_data <- omicsData$e_data
  iCol <- which(names(omicsData$e_data) == get_edata_cname(omicsData))

  if (!is.null(transformation)) {
    transform_data <- temp_data[, -iCol]

    if (transformation == 'lcpm') {
      ## log cpm, limma voom method and a similar method used for visualizations in edgeR

      ## EdgeR
      # First scales the prior.count/pseudo-count and adds 2x the scaled prior count to the libsize
      # prior.count.scaled <- lib.size/mean(lib.size)*prior.count
      # lib.size <- lib.size+2*prior.count.scaled
      # lib.size <- 1e-6*lib.size
      # Calculates log2 log2(t( (t(x)+prior.count.scaled) / lib.size ))

      # sum library size
      samp_sum <- apply(transform_data,
        2,
        sum,
        na.rm = TRUE
      ) + 1

      # divide adjusted (ensure non-zero) counts by library size
      div_sum <- sweep((transform_data + .5), 2, samp_sum, `/`)

      # Apply per million multiplier and log2
      temp_data[, -iCol] <- log2(div_sum * 10^6)
    } else if (transformation == 'upper') {
      warning("Zeros will be regarded as NA for 'upper' transformation")

      transform_data[transform_data == 0] <- NA

      # Grab non-zero upper quantile of data
      samp_upper <- apply(transform_data,
        2,
        quantile,
        na.rm = TRUE,
        probs = .75
      )

      g.q <- quantile(unlist(transform_data), probs = .75, na.rm = TRUE)

      # Divide each count by the upper quantile in respective columns
      div_75 <- sweep(transform_data, 2, samp_upper, `/`)

      # Set new data
      temp_data[, -iCol] <- div_75 * g.q # back transform
    } else if (transformation == 'median') {
      warning("Zeros will be regarded as NA for 'median' transformation")

      transform_data[transform_data == 0] <- NA

      # Grab non-zero median of data
      samp_med <- apply(transform_data,
        2,
        median,
        na.rm = TRUE
      )

      # Divide each count by the upper quantile in respective columns
      div_med <- sweep(transform_data, 2, samp_med, `/`)

      g.q <- median(unlist(transform_data), na.rm = TRUE)

      # Set new data
      temp_data[, -iCol] <- div_med * g.q
    }
  }

  plot_data <- temp_data %>%
    tidyr::pivot_longer(
      -!!e_data_cname,
      names_to = "variable",
      values_drop_na = TRUE,
      cols_vary = "slowest"
    ) %>% dplyr::mutate(
      variable = factor(
        variable,
        levels = names(temp_data)[
          -which(names(temp_data) == e_data_cname)
        ]
      )
    ) %>% data.frame

  # Farm boy, extract the group_DF attribute. As you wish.
  # groupDF <- attr(omicsData, "group_DF")
  groupDF <- get_group_DF(omicsData)

  # If facet_by is not null and isn't the same as either order_by or color_by.
  if (!is.null(facet_by)) {
    if (!(facet_by %in% c(order_by, color_by))) {
      # Extract the group_DF attribute. This will be used to facet the plots
      # later in the function.
      # facetDF <- attr(
      #   group_designation(omicsData = omicsData, main_effects = facet_by),
      #   "group_DF"
      # )
      facetDF <- get_group_DF(
        group_designation(omicsData = omicsData, main_effects = facet_by)
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
      orderDF <- dplyr::select(
        omicsData$f_data,
        !!dplyr::sym(get_fdata_cname(omicsData)),
        !!dplyr::sym(order_by)
      )
    } else {
      # Use the original group_DF attribute for ordering the samples. This
      # occurs when order_by = "Group".
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
      ordered = TRUE
    )

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
      # colorDF <- attr(
      #   group_designation(omicsData = omicsData, main_effects = color_by),
      #   "group_DF"
      # )
      colorDF <- get_group_DF(
        group_designation(omicsData = omicsData, main_effects = color_by)
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
                                    levels = sort(color_levels))

  }

  # Form awe-inspiring plots ---------------------------------------------------

  # Farm boy, make the graph skeleton. As you wish.
  p <- ggplot2::ggplot(plot_data)

  # Farm boy, color the plots based on the input. As you wish.
  if (is.null(color_by)) {
    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes(
        x = variable,
        y = value
      ))
  } else {
    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes(
        x = variable,
        y = value,
        fill = !!dplyr::sym(color_by)
      ))
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
  if (!is.null(facet_by)) {
    if (is.null(facet_cols)) {
      p <- p + ggplot2::facet_wrap(formula(paste("~", facet_by)),
        scales = "free_x"
      )
    } else {
      p <- p + ggplot2::facet_wrap(formula(paste("~", facet_by)),
        scales = "free_x",
        ncol = facet_cols
      )
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
    ggplot2::scale_fill_brewer(
      name = legend_title,
      palette = palette
    )

  # Farm boy, make me an interactive plot. As you wish.
  if (interactive && requirePlotly()) p <- plotly::ggplotly(p)

  # Farm boy, send my plot into the world to be adored by all. As you wish.
  return(p)
}

#' Plot statRes Object
#'
#' Produces plots that summarize the results contained in a `statRes` object.
#'
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#' @param plot_type defines which plots to be produced, options are "bar",
#'   "volcano", "gheatmap", "fcheatmap"; defaults to "bar".  See details for
#'   plot descriptions.
#' @param fc_threshold optional threshold value for fold change estimates.
#'   Modifies the volcano plot as follows:  Vertical lines are added at
#'   (+/-)\code{fc_threshold} and all observations that have absolute fold
#'   change less than \code{abs(fc_threshold)} are colored as 'non-significant'
#'   (as specified by \code{fc_colors}).
#' @param fc_colors vector of length three with character color values
#'   interpretable by ggplot. i.e. c("orange", "black", "blue") with the values
#'   being used to color negative, non-significant, and positive fold changes
#'   respectively
#' @param stacked TRUE/FALSE for whether to stack positive and negative fold
#'   change sections in the barplot, defaults to TRUE
#' @param show_sig This input is used when \code{plot_type = "gheatmap"}. A
#'   logical value. If TRUE a visual indicator that a certain bin combination is
#'   significant by the g-test is shown.
#' @param color_low This input is used when \code{plot_type = "gheatmap"}. A
#'   character string specifying the color of the gradient for low count values.
#' @param color_high This input is used when \code{plot_type = "gheatmap"}. A
#'   character string specifying the color of the gradient for high count
#'   values.
#' @param plotly_layout This input is used when \code{plot_type = "gheatmap"}. A
#'   list of arguments, not including the plot, to be passed to
#'   \code{plotly::layout} if \code{interactive = TRUE}.
#' @param interactive TRUE/FALSE for whether to create an interactive plot using
#'   plotly. Not valid for all plots.
#' @param x_lab character string specifying the x-axis label.
#' @param x_lab_size integer value indicating the font size for the x-axis. The
#'   default is 11.
#' @param x_lab_angle integer value indicating the angle of x-axis labels.
#' @param y_lab character string specifying the y-axis label.
#' @param y_lab_size integer value indicating the font size for the y-axis. The
#'   default is 11.
#' @param title_lab character string specifying the plot title.
#' @param title_lab_size integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab character string specifying the legend title.
#' @param legend_position character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", "bottom", or "none". The
#'   default is "none".
#' @param text_size integer specifying the size of the text (number of
#'   non-missing values) within the plot. The default is 3.
#' @param bw_theme logical value. If TRUE uses the ggplot2 black and white
#'   theme.
#' @param display_count logical value. Indicates whether the non-missing counts
#'   will be displayed on the bar plot. The default is TRUE.
#' @param custom_theme a ggplot `theme` object to be applied to non-interactive
#'   plots, or those converted by plotly::ggplotly().
#' @param cluster logical for heatmaps; TRUE will cluster biomolecules on X
#'   axis. defaults to TRUE for seqData statistics and FALSE for all others.
#' @param free_y_axis Logical. If TRUE the y axis for each bar plot can have its
#'   own range. The default is FALSE.
#' @param ... further arguments passed to or from other methods.
#'
#' @details Plot types:
#' \itemize{
#'  \item{"bar"} \code{?pmartR::statres_barplot} Bar-chart with bar heights
#'  indicating the number of significant biomolecules, grouped by test type and
#'  fold change direction.
#'  \item{"volcano"} \code{?pmartR::statres_volcano_plot} Scatter plot showing
#'  negative-log-pvalues against fold change. Colored by statistical
#'  significance and fold change.
#'  \item{"gheatmap"} \code{?pmartR::gtest_heatmap} Heatmap with x and y axes
#'  indicating the number of nonmissing values for two groups. Colored by
#'  number of biomolecules that fall into that combination of nonmissing values.
#'  \item{"fcheatmap"} Heatmap showing all biomolecules across comparisons,
#'  colored by fold change.
#' }
#'
#' @return ggplot2 plot object if interactive is FALSE, or plotly plot object if
#'   interactive is TRUE
#'
#' @export
#' @method plot statRes
#' @examplesIf identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)
#' library(pmartRdata)
#' # Group the data by condition
#' mypro <- group_designation(
#'   omicsData = pro_object,
#'   main_effects = c("Phenotype")
#' )
#'
#' # Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = mypro)
#' mypro <- applyFilt(
#'   filter_object = imdanova_Filt,
#'   omicsData = mypro,
#'   min_nonmiss_anova = 2
#' )
#'
#' # Implement the IMD ANOVA method and compuate all pairwise comparisons
#' # (i.e. leave the `comparisons` argument NULL)
#' anova_res <- imd_anova(omicsData = mypro, test_method = 'anova')
#' plot(anova_res)
#' plot(anova_res, plot_type = "volcano")
#'
#' imd_res <- imd_anova(omicsData = mypro, test_method = 'gtest')
#' plot(imd_res)
#'
#' imd_anova_res <- imd_anova(
#'   omicsData = mypro,
#'   test_method = 'comb',
#'   pval_adjust_a_multcomp = 'bon',
#'   pval_adjust_g_multcomp = 'bon'
#' )
#' plot(imd_anova_res, bw_theme = TRUE)
#' plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE)
#'
plot.statRes <- function(x,
                         plot_type = "bar",
                         fc_threshold = NULL,
                         fc_colors = c("red", "black", "green"),
                         stacked = TRUE,
                         show_sig = TRUE,
                         color_low = NULL,
                         color_high = NULL,
                         plotly_layout = NULL,
                         interactive = FALSE,
                         x_lab = NULL,
                         x_lab_size = 11,
                         x_lab_angle = NULL,
                         y_lab = NULL,
                         y_lab_size = 11,
                         title_lab = NULL,
                         title_lab_size = 14,
                         legend_lab = NULL,
                         legend_position = "right",
                         text_size = 3,
                         bw_theme = TRUE,
                         display_count = TRUE,
                         custom_theme = NULL,
                         cluster = FALSE,
                         free_y_axis = FALSE,
                         ...) {
  # Make sure we only have valid arguments
  if (length(list(...)) > 0) {
    warning("unused argument(s): ",
             toString(as.list(tail(match.call(), length(list(...))))))
  }

  # Farm boy, fix all the problems. As you wish.

  # Most plots are based on "number_significant" data frame so pull it out
  comp_df <- attr(x, "number_significant")

  ## --------##
  # Go through the given plot_types and remove any that aren't currently avilable
  plt_tyj <- try(
    match.arg(
      tolower(plot_type),
      c(
        "bar",
        "volcano",
        "gheatmap",
        "fcheatmap",
        "ma"
      )
    ),
    silent = TRUE
  )
  if (inherits(plt_tyj, 'try-error')) {
    warning(paste0(
      "Plot type '",
      plot_type,
      "' is not currently available, defaulting to bar plot."
    ))
    plot_type <- "bar"
  } else {
    plot_type <- plt_tyj
  }

  # Don't make biomolecule heatmaps if there's only one comparison
  if (plot_type %in% c("fcheatmap") & nrow(comp_df) == 1) {
    stop(paste("Fold change heatmaps not supported when only one comparison",
      "is being made.",
      sep = " "
    ))
  }

  if (plot_type %in% c("ma") && attr(x, "data_class") != "seqData") {
    stop("MA plots are only supported for transcriptomic data.",
      sep = " "
    )
  }

  # specified theme parameters
  if (!is.null(custom_theme)) {
    if (bw_theme)
      warning(paste("Setting both bw_theme to TRUE and specifying a custom",
        "theme may cause undesirable results",
        sep = " "
      ))
    if (!inherits(custom_theme, c("theme", "gg")))
      stop("custom_theme must be a valid 'theme' object as used in ggplot")
    mytheme = custom_theme
  } else mytheme = ggplot2::theme(
    plot.title = ggplot2::element_text(size = title_lab_size),
    axis.title.x = ggplot2::element_text(size = x_lab_size),
    axis.title.y = ggplot2::element_text(size = y_lab_size),
    axis.text.x = ggplot2::element_text(angle = x_lab_angle),
    legend.position = legend_position
  )

  # Both the volcano plot and heatmaps need a dataframe of fold changes by
  # comparison/biomolecule
  if (plot_type %in% c("volcano", "gheatmap", "fcheatmap")) {
    volcano <- make_volcano_plot_df(x)
  }

  # Bar plot -------------------------------------------------------------------

  if ("bar" %in% plot_type) {
    p <- statres_barplot(
      x = x,
      stacked = stacked,
      fc_colors = fc_colors,
      text_size = text_size,
      display_count = display_count,
      x_lab = x_lab,
      y_lab = y_lab,
      title_lab = title_lab,
      legend_lab = legend_lab,
      free_y_axis = free_y_axis
    )

    if (bw_theme) p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    return(p + mytheme)
  }

  # Volcano plot
  else if ("volcano" %in% plot_type) {
    if (!attr(x, "statistical_test") %in% c(
      "anova", "combined", "EdgeR_F", "Voom_T", "DESeq_Wald", "DESeq_LRT"
    )
    ) {
      stop(paste("imd_anova must have been run with test_method = 'anova' or",
        "'combined' to make the volcano plot. For seqData,",
        "DE_wrapper must have been run.",
        sep = " "
      ))
    }

    # still returns a ggplot, even if interactive = T

    p <-
      statres_volcano_plot(
        volcano = volcano,
        data_scale = attr(x, "data_info")$data_scale,
        pval_thresh = attr(x, "pval_thresh"),
        fc_colors = fc_colors,
        fc_threshold = fc_threshold,
        interactive = interactive,
        x_lab = x_lab,
        y_lab = y_lab,
        title_lab = title_lab,
        legend_lab = legend_lab
      )

    if (bw_theme) p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    p <- p + mytheme

    if (interactive && requirePlotly())
      return(plotly::ggplotly(p, tooltip = c("text"))) else
      return(p)
  }

  # g-test heat map ------------------------------------------------------------

  else if ("gheatmap" %in% plot_type) {
    if (!attr(x, "statistical_test") %in% c("gtest", "combined")) {
      stop(paste("imd_anova must have been run with test_method = 'gtest' or",
        "'combined' to make the g-test heatmap. Not valid for seqData.",
        sep = " "
      ))
    }

    p <-
      gtest_heatmap(
        volcano = volcano,
        pval_thresh = attr(x, "pval_thresh"),
        show_sig = show_sig,
        interactive = interactive,
        color_low = color_low,
        color_high = color_high,
        plotly_layout = plotly_layout,
        text_size = text_size,
        display_count = display_count,
        x_lab = x_lab,
        y_lab = y_lab,
        title_lab = title_lab,
        legend_lab = legend_lab
      )

    if (!interactive) {
      p <- p + mytheme
    }
  }

  # biomolecule fold change heatmap
  else if ("fcheatmap" %in% plot_type) {
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab))
      "Biomolecule" else
      x_lab
    the_y_label <- if (is.null(y_lab))
      "Comparison" else
      y_lab
    the_title_label <- if (is.null(title_lab))
      "Log Fold Change" else
      title_lab
    the_legend_label <- if (is.null(legend_lab))
      "Fold Change" else
      legend_lab

    # For now just consider biomolecules significant with respect to ANOVA
    volcano <- dplyr::filter(volcano, Type %in% c(
      "ANOVA", "EdgeR_F",
      "Voom_T", "DESeq_Wald",
      "DESeq_LRT"
    ))

    volcano_sigs <- dplyr::filter(volcano, P_value < attr(x, "pval_thresh"))
    if (!(nrow(volcano_sigs)) > 0)
      warning("No molecules significant at the provided p-value threshold")
    colnames(volcano_sigs)[1] <- "Biomolecule"

    if (cluster) {
      wide <- volcano_sigs %>% 
        tidyr::pivot_wider(
          id_cols = Biomolecule,
          names_from = Comparison,
          values_from = Fold_change
        ) %>% data.frame
      wide <- wide[!is.na(wide$Biomolecule), ]
      row.names(wide) <- wide$Biomolecule
      dist_mat <- dist(wide[-1])
      res_hclust <- hclust(dist_mat)
      order_biom <- rev(res_hclust$labels[res_hclust$order])

      volcano_sigs$Biomolecule <- as.character(volcano_sigs$Biomolecule)
      volcano_sigs$Biomolecule <- factor(
        volcano_sigs$Biomolecule,
        levels = unique(paste0(order_biom, volcano_sigs$Biomolecule))
      )
    } else {
      volcano_sigs$Biomolecule <- as.factor(volcano_sigs$Biomolecule)
    }

    p <- ggplot2::ggplot(
      volcano_sigs,
      ggplot2::aes(Biomolecule,
        Comparison,
        text = paste(
          "ID:",
          Biomolecule,
          "<br>",
          "Pval:",
          P_value
        )
      )
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = Fold_change), color = "white") +
      ggplot2::scale_fill_gradient(
        low = fc_colors[1],
        high = fc_colors[3],
        name = the_legend_label
      ) +
      ggplot2::xlab(the_x_label) +
      ggplot2::ylab(the_y_label) +
      ggplot2::ggtitle(the_title_label) +
      mytheme
    if (interactive && requirePlotly())
      return(plotly::ggplotly(p, tooltip = c("text"))) else
      return(p)
  } else if ("ma" %in% plot_type) {
    ## Color by significance
    comps <- strsplit(attr(x, "comparisons"), "_vs_")

    plotter <- purrr::map_dfr(1:length(comps), function(n_comp) {
      label <- attr(x, "comparisons")[n_comp]
      comp <- comps[[n_comp]]
      mean_df <- attr(x, "MA_means")
      pval <- grep(paste0("^P_value_", label), colnames(x), value = T)

      v1 <- mean_df[[comp[1]]]
      v2 <- mean_df[[comp[2]]]
      v3 <- x[[pval]]

      if (length(v1) == 0) {
        v1 <- NA
        v2 <- NA
        v3 <- NA
      }

      data.frame(var1 = v1, var2 = v2, pval = v3, comp = label)
    })

    p <- ggplot2::ggplot(
      plotter,
      ggplot2::aes(
        x = log2((var1 + var2) / 2),
        y = log2(var1 / var2), ## where mean is 0 or na in a group, goes to Inf
        color = pval < attr(x, "pval_thresh")
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~comp) +
      ggplot2::geom_segment(
        y = 0, yend = 0, linetype = "dashed", color = "red",
        x = min(log2((plotter$var1 + plotter$var2) / 2), na.rm = TRUE),
        xend = max(log2((plotter$var1 + plotter$var2) / 2), na.rm = TRUE)
      ) +
      ggplot2::labs(
        x = "A (Log2 Average Expression)",
        y = "M (Log2 Fold change)",
        color = paste("Significance < ", attr(x, "pval_thresh"))
      )

    if (bw_theme) p <- p +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))

    p <- p + mytheme

    if (interactive && requirePlotly())
      return(plotly::ggplotly(p, tooltip = c("text"))) else
      return(p)
  }

  return(p)
}

#' Extract flag columns from a statRes object
#'
#' Changes the flags columns from a statRes object into a format that the
#' statRes plot funcitons can handle. pmartR is an unruly beast that cannot be
#' tamed!!
#'
#' @param x A statRes object.
#' @param test character string indicating the type of test run.
#'
#' @return A data frame with the sample IDs and significance flags from a statistical test.
#' 
prep_flags <- function(x, test) {
  if (test == "anova") {
    # Assemble a data frame with the sample IDs and anova flags.
    da_flag <- data.frame(
      x[, 1, drop = FALSE],
      x[, grep("^Flag_A_", colnames(x))],
      check.names = FALSE
    )

    # Remove "Flag_A_" from column names. The first column name is removed
    # because it corresponds to the biomolecule ID column.
    colnames(da_flag)[-1] <- gsub(
      "^Flag_A_",
      "",
      colnames(x)[grep("^Flag_A_", colnames(x))]
    )
  } else if (test == "gtest") {
    # Throw an error at the user like Fezzik throwing rocks at the Dread Pirate
    # Roberts. It doesn't make sense to create a volcano plot with Gtest flags.
    # Coding would take so much less time if we didn't have to think of all the
    # ways a user can mess things up.
    stop("A volcano plot cannot be created with Gtest flags.")
  } else if (test %in% c("EdgeR_LRT", "EdgeR_F", "Voom_T", "DESeq_Wald", "DESeq_LRT")) {
    # Assemble a data frame with the sample IDs and flags.
    flag_cols <- grep("^Flag_", colnames(x))
    da_flag <- x[c(1, flag_cols)]

    # Remove "Flag_A_" from column names. The first column name is removed
    # because it corresponds to the biomolecule ID column.
    colnames(da_flag)[-1] <- gsub("^Flag_", "", colnames(da_flag)[-1])

    # Runs when the "combined" option was used.
  } else {
    # NOTE: It does not make sense to have Gtest flags present when creating a
    # volcano plot. For this reason we will only use the ANOVA flags when the
    # test is "combined".

    # Assemble a data frame with the sample IDs and anova flags.
    da_flag <- data.frame(
      x[, 1, drop = FALSE],
      x[, grep("^Flag_A_", colnames(x))],
      check.names = FALSE
    )

    # Remove "Flag_A_" from column names. The first column name is removed
    # because it corresponds to the biomolecule ID column.
    colnames(da_flag)[-1] <- gsub(
      "^Flag_A_",
      "",
      colnames(x)[grep("^Flag_A_", colnames(x))]
    )
  }

  return(da_flag)
}

#' Create a plotting dataframe for volcano plots and heatmaps
#'
#' A function internal to pmartR:::\link{plot.statRes} which creates the
#' dataframe necessary to construct volcano plots and heatmaps.
#'
#' @param x `statRes` object to be plotted, usually the result of `imd_anova`
#'
#' @returns `data.frame` object with plotting information about each biomolecule
#' such as missing counts per group, and p-values for t and g-tests.
#'
#' @keywords internal
#'
make_volcano_plot_df <- function(x) {
  # fold change values for volcano plot

  fc_data <- x[, c(1, grep("^Fold_change", colnames(x)))]

  colnames(fc_data) <-
    gsub(
      pattern = "^Fold_change_",
      replacement = "",
      x = colnames(fc_data)
    )

  fc_data <- fc_data %>%
    tidyr::pivot_longer(
      -1,
      names_to = "Comparison",
      values_to = "Fold_change",
      cols_vary = "slowest"
    ) %>% data.frame

  # Run the cmbn_flags function here.

  # fold change flags for coloring
  fc_flags <- prep_flags(
    x = x,
    test = attr(x, "statistical_test")
  )
  fc_flags <- fc_flags %>%
    tidyr::pivot_longer(
      -1,
      names_to = "Comparison",
      values_to = "Fold_change_flag",
      cols_vary = "slowest"
    ) %>%
    dplyr::mutate(
      Fold_change_flag = as.character(Fold_change_flag)
    ) %>%
    data.frame

  # p values for labeling and y axis in anova volcano plot
  p_data <- x[, c(1, grep("^P_value", colnames(x)))]
  pvals <- p_data %>%
    tidyr::pivot_longer(
      -1,
      names_to = "Comparison",
      values_to = "P_value",
      cols_vary = "slowest"
    ) %>% data.frame

  # grouping column based on test type
  if (attr(x, "statistical_test") == "combined") {
    pvals$Type <- "G-test"
    pvals$Type[grep(pattern = "^P_value_A_", x = pvals$Comparison)] <-
      "ANOVA"
    pvals$Comparison <-
      gsub(
        pattern = "^P_value_(G|A)_",
        replacement = "",
        pvals$Comparison,
        perl = T
      )
  } else if (attr(x, "statistical_test") == "gtest") {
    pvals$Type <- "G-test"
    pvals$Comparison <-
      gsub(
        pattern = "^P_value_G_",
        replacement = "",
        pvals$Comparison,
        perl = T
      )
  } else if (attr(x, "statistical_test") == "anova") {
    pvals$Type <- "ANOVA"
    pvals$Comparison <-
      gsub(
        pattern = "^P_value_A_",
        replacement = "",
        pvals$Comparison,
        perl = T
      )
  } else {
    pvals$Type <- attr(x, "statistical_test")
    pvals$Comparison <-
      gsub(
        pattern = "^P_value_",
        replacement = "",
        pvals$Comparison,
        perl = T
      )
  }

  volcano <-
    merge(merge(fc_data, pvals, all = TRUE), fc_flags, all = TRUE)

  # levels of comparison now of the form 'GROUPNAME_X vs GROUPNAME_Y'
  volcano$Comparison <-
    gsub(
      pattern = "_vs_",
      replacement = " vs ",
      volcano$Comparison
    )

  # create counts for gtest plot (number present in each group)
  if (attr(x, "statistical_test") %in% c("gtest", "combined")) {
    counts <-
      x[c(1, grep("^Count_", colnames(x)))]

    # trim column names so they are just group names
    colnames(counts) <-
      gsub("^Count_", replacement = "", colnames(counts))

    counts_df <- data.frame()
    for (comp in as.character(unique(volcano$Comparison))) {
      # create a vector of the two group names being compared
      groups = strsplit(comp, " vs ")[[1]]
      gsize_1 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[1]))
      gsize_2 = nrow(attr(x, "group_DF") %>% dplyr::filter(Group == groups[2]))

      # will contain ID column and count column corresponding to the two groups
      temp_df <-
        counts[c(
          1,
          which(colnames(counts) == groups[1]),
          which(colnames(counts) == groups[2])
        )]
      temp_df$Comparison <- comp

      # rename the columns to something static so they can be rbind-ed
      colnames(temp_df)[which(colnames(temp_df) == groups[1])] <- "Count_First_Group"
      colnames(temp_df)[which(colnames(temp_df) == groups[2])] <- "Count_Second_Group"

      # store proportion of nonmissing to color g-test values in volcano plot
      temp_df$Prop_First_Group <- temp_df$Count_First_Group / gsize_1
      temp_df$Prop_Second_Group <- temp_df$Count_Second_Group / gsize_2

      counts_df <- rbind(counts_df, temp_df)
    }

    # should automatically left join by ID AND Comparison
    suppressWarnings(volcano <-
      volcano %>% dplyr::left_join(counts_df))
  }

  return(volcano)
}

#' Fold change barplots for statres objects
#'
#' Plots a bar-chart with bar heights indicating the number of significant
#' biomolecules, grouped by test type and fold change direction.
#'
#' @param x,stacked,fc_colors passed from
#'   \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}
#' @param text_size An integer specifying the size of the text (number of
#'   non-missing values) within the plot. The default is 3.
#' @param display_count logical value. Indicates whether the non-missing counts will
#'   be displayed on the bar plot. The default is TRUE.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label.
#' @param title_lab character string specifying the plot title.
#' @param legend_lab character string specifying the legend title.
#'
#' @return `ggplot` object - barplot.
#'
#' @keywords internal
#'
statres_barplot <- function(x,
                            stacked,
                            fc_colors,
                            text_size,
                            display_count,
                            x_lab,
                            y_lab,
                            title_lab,
                            legend_lab,
                            free_y_axis) {
  # If free_y_axis is true then each bar plot can have its own y axis range.
  if (free_y_axis) {
    da_scales <- "free"
  } else {
    da_scales <- NULL
  }

  # Farm boy, do all the tedious label crap. As you wish.
  the_x_label <- if (is.null(x_lab))
    "Statistical test, by group comparison" else
    x_lab
  the_y_label <- if (is.null(y_lab))
    "Count of Biomolecules" else
    y_lab
  the_title_label <- if (is.null(title_lab))
    "Number of DE Biomolecules Between Groups" else
    title_lab
  the_legend_label <- if (is.null(legend_lab))
    "Fold Change Sign" else
    legend_lab

  comp_df <- attr(x, "number_significant")

  comp_df_melt <- comp_df %>%
    tidyr::pivot_longer(
      -Comparison,
      values_to = "Count",
      names_to = "Direction",
      cols_vary = "slowest"
    ) %>% data.frame
    
  levels(comp_df_melt$Comparison) <- gsub(
    pattern = "_",
    replacement = " ",
    levels(comp_df_melt$Comparison)
  )

  # Bar plots side-by-side, both going up
  # ggplot(data=comp_df_melt,aes(Comparison,Count,fill=Direction))+
  #   geom_bar(stat='identity',position='dodge')

  ## Up direction is positive, down direction is negative
  if (stacked)
    comp_df_melt[grep("Down", comp_df_melt$Direction), ]$Count <- (
      -comp_df_melt[grep("Down", comp_df_melt$Direction), ]$Count
    )

  # add whichtest, and posneg columns used for plot grouping and label
  # adjustment
  comp_df_melt <- comp_df_melt %>%
    dplyr::mutate(
      whichtest = ifelse(grepl("anova", Direction),
        "anova",
        ifelse(grepl("gtest", Direction),
          "gtest",
          "total"
        )
      ),
      posneg = ifelse(grepl("Up", Direction),
        "Positive",
        "Negative"
      )
    ) %>%
    dplyr::arrange(desc(posneg))

  # get only anova or only g-test rows if user did not specify combined
  if (attr(x, "statistical_test") %in% c("anova", "gtest")) {
    comp_df_melt <- comp_df_melt %>%
      dplyr::filter(whichtest == "total") %>%
      dplyr::mutate(whichtest = attr(x, "statistical_test"))
  }

  p <- ggplot2::ggplot(data = comp_df_melt, ggplot2::aes(Comparison, Count)) +
    ggplot2::geom_bar(
      ggplot2::aes(
        x = whichtest,
        fill = posneg
      ),
      stat = 'identity',
      position = if (stacked) {
        ggplot2::position_identity()
      } else {
        ggplot2::position_dodge(0.9)
      }
    ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour = 'gray50') +
    ggplot2::scale_fill_manual(
      values = c(fc_colors[1], fc_colors[3]),
      labels = c("Negative", "Positive"),
      name = the_legend_label
    ) +
    ggplot2::facet_wrap(~Comparison, scales = da_scales) +
    ggplot2::xlab(the_x_label) +
    ggplot2::ylab(the_y_label) +
    ggplot2::ggtitle(the_title_label)

  # Farm boy, display the counts on the plot. As you wish
  if (display_count) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          x = whichtest,
          group = posneg,
          label = abs(Count)
        ),
        size = text_size,
        position = if (stacked) {
          ggplot2::position_identity()
        } else {
          ggplot2::position_dodge(0.9)
        }
      )
  }

  return(p)
}

#' Plot a heatmap for the g-test results of imd-anova
#'
#' Plots a heatmap showing bins for combinations of # of biomolecules present
#' across groups.  Bins are colored by the number of biomolecules falling into
#' each bin, and have indicators for significance by g-test.
#'
#' @param volcano `data.frame` produced by pmartR:::\link{make_volcano_plot_df}
#' @param pval_thresh numeric value indicating the p-value threshold for
#'  significance.
#' @param show_sig Boolean whether to show the visual indicator that a certain
#'   bin combination is significant by the g-test
#' @param interactive passed from
#'   \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}. If T, will build
#'   a plotly version of the plot.  Defaults to FALSE.
#' @param color_low character string specifying the color of the gradient for
#'   low count values.
#' @param color_high character string specifying the color of the gradient for
#'   high count values.
#' @param plotly_layout A list of arguments, not including the plot, to be
#'   passed to plotly::layout if interactive = T.
#' @param text_size An integer specifying the size of the text (number of
#'   non-missing values) within the plot. The default is 3.
#' @param display_count logical value. Indicates whether the non-missing counts
#'   will be displayed on the bar plot. The default is TRUE.
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label.
#' @param title_lab character string specifying the plot title.
#' @param legend_lab character string specifying the legend title.
#'
#' @return `ggplot or plotly` object, depending on if \code{interactive} is set
#' to TRUE or FALSE respectively. A g-test heatmap
#'
#' @keywords internal
#'
gtest_heatmap <-
  function(volcano,
           pval_thresh,
           show_sig,
           interactive,
           color_low,
           color_high,
           plotly_layout,
           text_size,
           display_count,
           x_lab,
           y_lab,
           title_lab,
           legend_lab) {
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab))
      "Nonmissing (first group)" else
      x_lab
    the_y_label <- if (is.null(y_lab))
      "Nonmissing (second group)" else
      y_lab
    the_title_label <- if (is.null(title_lab))
      "Biomolecule counts (and # significant) for each combination of non-missing values" else
      title_lab
    the_legend_label <- if (is.null(legend_lab))
      "Number of biomolecules \nin this bin" else
      legend_lab

    temp_data_gtest <- volcano %>%
      dplyr::filter(Type == "G-test") %>%
      dplyr::mutate(
        Fold_change_flag = dplyr::case_when(
          Prop_First_Group > Prop_Second_Group &
            P_value <= pval_thresh ~ "2",
          Prop_First_Group < Prop_Second_Group &
            P_value <= pval_thresh ~ "-2",
          is.na(Fold_change) ~ "0",
          TRUE ~ Fold_change_flag
        )
      )

    # summarized data frame of # of biomolecules per count combination.
    gtest_counts <- temp_data_gtest %>%
      dplyr::group_by(Count_First_Group, Count_Second_Group, Comparison) %>%
      dplyr::summarise(
        n = dplyr::n(),
        n_sig = sum(na.omit(P_value < pval_thresh)),
        sig = any(na.omit(P_value < pval_thresh))
      ) %>%
      dplyr::mutate(
        Count_First_Group = as.character(Count_First_Group),
        Count_Second_Group = as.character(Count_Second_Group)
      )

    g1_counts <- unique(as.numeric(gtest_counts$Count_First_Group))
    g2_counts <- unique(as.numeric(gtest_counts$Count_Second_Group))

    all_counts <- expand.grid(
      as.character(g1_counts),
      as.character(g2_counts),
      unique(gtest_counts$Comparison),
      stringsAsFactors = F
    ) %>%
      `colnames<-`(c("Count_First_Group", "Count_Second_Group", "Comparison"))

    gtest_counts <- all_counts %>% dplyr::left_join(gtest_counts) %>%
      dplyr::mutate(
        Count_First_Group = factor(
            Count_First_Group,
            levels = sort(unique(as.numeric(Count_First_Group)))
          ),
        Count_Second_Group = factor(
            Count_Second_Group,
            levels = sort(unique(as.numeric(Count_Second_Group)))
          )
      )

    if (interactive && requirePlotly()) {
      comps <- unique(gtest_counts$Comparison)
      subplot_list <- list()
      limits = range(na.omit(gtest_counts$n))
      # legend entry will not appear for a single plot, so putting info in title
      subtext = if (length(comps) == 1)
        "\n(star indicates statistical significance)" else
        ""

      for (i in 1:length(comps)) {
        data <- gtest_counts %>% dplyr::filter(Comparison == comps[i])

        p <- plotly::plot_ly() %>%
          plotly::add_trace(
            data = data %>% dplyr::filter(sig),
            x = ~ as.numeric(as.character(Count_First_Group)),
            y = ~ as.numeric(as.character(Count_Second_Group)),
            type = 'scatter',
            mode = "markers",
            showlegend = i == length(comps),
            name = "Statistically significant",
            hoverinfo = "skip",
            marker = list(
              symbol = "star",
              color = "white",
              line = list(color = "black", width = 0.5)
            ),
            colors = grDevices::colorRamp(
              c(
                if (is.null(color_low))
                  "#132B43" else
                    color_low,
                if (is.null(color_high))
                  "#56B1F7" else
                    color_high
              )
            )
          ) %>%
          plotly::add_trace(
            data = data,
            x = ~ as.numeric(as.character(Count_First_Group)),
            y = ~ as.numeric(as.character(Count_Second_Group)),
            z = ~ n,
            type = "heatmap",
            hoverinfo = 'text',
            text = ~ sprintf("%s biomolecules in this bin. (%s) significant", n, n_sig),
            showscale = i == length(comps),
            xgap = 0.6, ygap = 0.6
          ) %>%
          plotly::add_annotations(
            text = comps[[i]],
            x = 0.5,
            y = 1,
            xref = "paper",
            yref = "paper",
            showarrow = FALSE,
            yanchor = "bottom",
            xanchor = "center",
            font = list(size = 14)
          )

        if (i == length(comps)) {
          p <- p %>% plotly::colorbar(
            title = "Count biomolecules \nin this bin.", limits = limits
          )
        }


        if (is.null(plotly_layout)) {
          p <- p %>%
            plotly::layout(
              plot_bgcolor = 'grey',
              title = list( 
                text = paste(the_title_label, subtext),
                font = list(size = 14),
                y = 0.95, yanchor = "top"
              ),
              margin = list(t = 65),
              xaxis = list(tickvals = g1_counts, zeroline = F),
              yaxis = list(tickvals = g2_counts, zeroline = F)
            )
        } else {
          p <- do.call(plotly::layout, c(list(p), plotly_layout))
        }

        subplot_list[[length(subplot_list) + 1]] <- p
      }

      p <- plotly::subplot(subplot_list, shareY = TRUE) %>%
        plotly::layout(
          xaxis = list(title = the_x_label),
          yaxis = list(title = the_y_label)
        )
    } else {
      # create the alternative column
      if (display_count) {
        gtest_counts <- gtest_counts %>% 
          dplyr::mutate(
            count_text = ifelse(
              is.na(n),
              NA,
              ifelse(n_sig > 0, sprintf("%s (%s)", n, n_sig), n)
            )
          ) 
      }
      
      p <- ggplot2::ggplot(gtest_counts) +
        ggplot2::theme_minimal() +
        ggplot2::geom_tile(
          ggplot2::aes(Count_First_Group,
            Count_Second_Group,
            fill = n
          ),
          color = "black"
        )

      if(show_sig) {
        sig_data <- gtest_counts %>% dplyr::filter(sig)
        
        if(nrow(sig_data) > 0) {
          p <- p + ggplot2::geom_point(
            data = sig_data,
            ggplot2::aes(Count_First_Group, Count_Second_Group, shape = "1"),
            fill = "white"
          ) +
          ggplot2::scale_shape_manual(name = "Statistically significant",
                                      labels = "",
                                      values = 21)
        }
      }

      if (display_count) {
        p <- p + ggplot2::geom_text(
          ggplot2::aes(
            Count_First_Group, 
            Count_Second_Group, 
            label = count_text
          ),
          nudge_x = -0.5, nudge_y = 0.5, hjust = -0.1, vjust = 1.5,
          color = "white", size = text_size
        )
      }

      p <- p +
        ggplot2::facet_wrap(~Comparison) +
        ggplot2::scale_fill_gradient(
          name = the_legend_label,
          low = if (is.null(color_low)) "#132B43" else color_low,
          high = if (is.null(color_high)) "#56B1F7" else color_high
        ) +
        ggplot2::xlab(the_x_label) +
        ggplot2::ylab(the_y_label) +
        ggplot2::ggtitle(the_title_label)
    }

    return(p)
  }

#' Volcano plot for the anova results of imd-anova
#'
#' Plots a volcano plot showing negative log10 p-values on the y axis and fold
#' change on the x axis. Each point is colored by fold change direction and
#' whether or not it was significant by ANOVA.
#'
#' @param volcano data frame produced by pmartR:::\link{make_volcano_plot_df}
#' @param pval_thresh numeric value between 0 and 1 for the alpha level to
#'   determine significance. Any values that are significant at this level will
#'   be colored based on fc_colors.
#' @param data_scale One of c("log2","log","log10"), for labeling purposes.
#' @param fc_colors,fc_threshold,interactive passed from
#'   \code{\link[pmartR:plot.statRes]{pmartR::plot.statRes()}}
#' @param x_lab character string specifying the x-axis label.
#' @param y_lab character string specifying the y-axis label.
#' @param title_lab character string specifying the plot title.
#' @param legend_lab character string specifying the legend title.
#'
#' @return `ggplot` object.  A volcano plot.
#'
#' @keywords internal
#'
statres_volcano_plot <-
  function(volcano,
           data_scale,
           pval_thresh,
           fc_colors,
           fc_threshold,
           interactive,
           x_lab,
           y_lab,
           title_lab,
           legend_lab) {
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab))
      paste("Fold-change (", data_scale, ")", sep = "") else
      x_lab
    the_y_label <- if (is.null(y_lab))
      "-log[10](p-value)" else
      y_lab
    the_title_label <- if (is.null(title_lab))
      "" else
      title_lab
    the_legend_label <- if (is.null(legend_lab))
      "Fold Change" else
      legend_lab

    # color vector which assigns black to gtest flag values (-2, 2)
    cols_anova <-
      c(
        "-2" = fc_colors[2],
        "-1" = fc_colors[1],
        "0" = fc_colors[2],
        "1" = fc_colors[3],
        "2" = fc_colors[2]
      )

    # temp data with rows only for ANOVA
    temp_data_anova <- volcano %>%
      dplyr::filter(Type %in% c(
        "ANOVA", "EdgeR_F", "Voom_T", "DESeq_Wald", "DESeq_LRT"
      )) %>%
      dplyr::mutate(
        Fold_change_flag = dplyr::case_when(
          is.na(Fold_change) |
            abs(Fold_change) < ifelse(length(fc_threshold) == 0,
              0,
              abs(fc_threshold)
            ) ~ "0",
          Fold_change > 0 &
            P_value <= pval_thresh ~ "1",
          Fold_change < 0 &
            P_value <= pval_thresh ~ "-1",
          TRUE ~ Fold_change_flag
        )
      )

    # interactive plots need manual text applied to prepare for ggplotly
    # conversion
    if (interactive && requirePlotly()) {
      p <-
        ggplot2::ggplot(temp_data_anova, ggplot2::aes(
          Fold_change,
          -log(P_value, base = 10),
          text = paste(
            "ID:",
            !!dplyr::sym(colnames(volcano)[1]),
            "<br>",
            "Pval:",
            round(P_value, 4)
          )
        ))
    } else {
      p <-
        ggplot2::ggplot(
          data = temp_data_anova,
          ggplot2::aes(
            Fold_change,
            -log(P_value, base = 10)
          )
        )
    }

    # draw vertical lines at +- fc threshold
    if (length(fc_threshold) > 0) {
      p <- p +
        ggplot2::geom_vline(ggplot2::aes(xintercept = abs(fc_threshold)),
          lty = 2
        ) +
        ggplot2::geom_vline(ggplot2::aes(xintercept = abs(fc_threshold) * (-1)),
          lty = 2
        )
    }

    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = Fold_change_flag), shape = 1) +
      ggplot2::facet_wrap(~Comparison) +
      ggplot2::ylab(the_y_label) +
      ggplot2::xlab(the_x_label) +
      ggplot2::ggtitle(the_title_label) +
      ggplot2::scale_color_manual(
        values = cols_anova,
        name = the_legend_label,
        labels = c(
          paste0("Neg(", unique(temp_data_anova$Type), ")"),
          "0", paste0("Pos(", unique(temp_data_anova$Type), ")")
        ),
        breaks = c("-1", "0", "1")
      )

    return(p)
  }

#' Require Plotly
#'
#' Loads plotly if installed, else prints a message telling the user to install
#' plotly.
#' 
#' @return whether plotly is installed
#' 
#' @keywords internal
#' 
#' @noRd
requirePlotly <- function() {
  if (requireNamespace("plotly", quietly = TRUE)) {
    return(TRUE)
  }
  
  warning("Package 'plotly' is not installed. Please install it to use interactive plots.")
  return(FALSE)
}
