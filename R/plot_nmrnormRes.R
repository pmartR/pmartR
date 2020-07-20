#' Plots an object of class nmrnormRes
#' 
#' For plotting an S3 object of type 'nmrnormRes':
#' 
#' @param nmrnormRes_object an object of type nmrnormRes, created by \code{\link{normalize_nmr}}  
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
#' }
#' 
#' @return plots ggplot2 object
#' 
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object_identified, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, metabolite_name = "unkm1.53")
#' plot(nmr_norm)
#' 
#' # alternate specification: #
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, sample_property_cname = "Concentration")
#' plot(nmr_norm)
#'}
#' 
#' @rdname plot-nmrnormRes
#' @export
#'

plot.nmrnormRes <- function(nmrnormRes_object, x_lab = NULL, ...) {
  require(ggplot2)
  .plot.nmrnormRes(nmrnormRes_object, x_lab = NULL, ...)
}

.plot.nmrnormRes<- function(nmrnormRes_object, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = TRUE){
  
  #check for an nmrnormRes object #
  if(!inherits(nmrnormRes_object, "nmrnormRes")) stop("object must be of class 'nmrnormRes'")
  
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
  
  #extracting attributes from nmrnormRes_object
  sample_property_cname = attr(nmrnormRes_object, "nmr_info")$sample_property_cname
  metabolite_name <- attr(nmrnormRes_object, "nmr_info")$metabolite_name
  fdata_cname = attr(nmrnormRes_object, "cnames")$fdata_cname
  
  #organize nmrnormRes_object
  fdata_cname_ind = which(names(nmrnormRes_object) == fdata_cname)
  
  data = data.frame(Sample = nmrnormRes_object$Sample, value = nmrnormRes_object$value)
  
  #make labels
  xlabel <- ifelse(is.null(x_lab), "Sample ID", x_lab)
  if(is.null(sample_property_cname)){
    ylabel <- ifelse(is.null(y_lab), metabolite_name, y_lab)
    plot_title <- ifelse(is.null(title_plot), "Reference Metabolite Profile", title_plot)
  }else{
    if(is.null(metabolite_name)){
      ylabel <- ifelse(is.null(y_lab), sample_property_cname, y_lab)
      plot_title <- ifelse(is.null(title_plot), "Sample Property Profile", title_plot)
    }
  }
  
  p <- ggplot(data = data, aes(x = Sample, y = value))
  
  p <- p + geom_point() +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(plot_title) +
    theme(plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size)) 
  
  if(bw_theme == TRUE){
    p<- p + theme_bw()
  }
  
  
  return(p)  
}


