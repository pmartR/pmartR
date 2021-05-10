#' Plots an object of class isobaricnormRes
#' 
#' For plotting an S3 object of type 'isobaricnormRes':
#' 
#' @param isobaricnormRes_object an object of type isobaricnormRes, created by \code{\link{normalize_isobaric}}  
#' @param x_lab character string to be used for x-axis label. Defaults to NULL  
#' @param ... further arguments 
#'
#' \tabular{ll}{
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to TRUE, in which case the ggplot2 default theme is used. \cr
#' \code{order} \tab logical indicates whether to order data by exp_cname'}
#' 
#' @return plots ggplot2 object
#' 
#' 
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(isobaric_object)
#' 
#' isobaric_object = edata_transform(isobaric_object, "log2")
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
#' @export
#'

plot.isobaricnormRes <- function (isobaricnormRes_object, x_lab = NULL, ...) {
  .plot.isobaricnormRes(isobaricnormRes_object, x_lab = NULL, ...)
}

.plot.isobaricnormRes <- function (isobaricnormRes_object, x_lab = NULL,
                                   y_lab = NULL, title_plot = NULL,
                                   title_size = 14, x_lab_size = 11,
                                   y_lab_size = 11, bw_theme = TRUE,
                                   order = FALSE) {
  
  # Check inputs ---------------------------------------------------------------
  
  if(!inherits(isobaricnormRes_object, "isobaricnormRes")) {
    
    stop("object must be of class 'isobaricnormRes'")
    
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
  
  # Combine info from e_data and f_data ----------------------------------------
  
  # extracting attributes from isobaricnormRes_object
  exp_cname = attr(isobaricnormRes_object, "isobaric_info")$exp_cname
  fdata_cname = attr(isobaricnormRes_object, "cnames")$fdata_cname
  edata_cname <- attr(isobaricnormRes_object, "cnames")$edata_cname
  
  # Transform the isobaricnormRes data frames into a format usable by ggplot2.
  tall_data <- prime_iso(isonormRes = isobaricnormRes_object,
                         exp_cname = exp_cname,
                         fdata_cname = fdata_cname,
                         edata_cname = edata_cname)
  
  # Create pretty plots --------------------------------------------------------
  
  # make labels
  xlabel <- ifelse(is.null(x_lab), "ref_sample", x_lab)
  ylabel <- ifelse(is.null(y_lab), "log Abundance", y_lab)
  plot_title <- ifelse(is.null(title_plot),
                       "Reference Sample Profile",
                       title_plot)
  
  
  # If order is TRUE order the box plots by experiment name/value.
  if (order == TRUE) {
    
    xlabel <- ifelse(is.null(x_lab),exp_cname, x_lab)
    
    p <- ggplot2::ggplot(data = tall_data,
                         ggplot2::aes(x = .data[[exp_cname]],
                                      y = values))
    
    # Otherwise separate the box plots by sample name.
  } else {
    
    p <- ggplot2::ggplot(data = tall_data,
                         ggplot2::aes(x = .data[[fdata_cname]],
                                      y = values))
    
  }
  
  p <- p + ggplot2::geom_boxplot() +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size),
                   axis.title.x = ggplot2::element_text(size = x_lab_size),
                   axis.title.y = ggplot2::element_text(size = y_lab_size)) 
  
  if (bw_theme == TRUE) {
    
    p <- p + ggplot2::theme_bw()
    
  }
  
  
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
