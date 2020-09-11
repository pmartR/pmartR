#' Plots an object of class isobaricnormRes
#' 
#' For plotting an S3 object of type 'isobaricnormRes':
#' 
#' @param an object of type isobaricnormRes, created by \code{\link{normalize_isobaric}}  
#' @param x_lab character string to be used for x-axis label. Defaults to NULL  
#' @param ... further arguments 
#'
#' \tabular{ll}{
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{order} \tab logical indicates whether to order data by exp_cname'}
#' 
#' @return plots ggplot2 object
#' 
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(isobaric_object)
#' 
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' result = normalize_isobaric(isobaric_object, apply_norm = F)
#' 
#' plot(result)
#'}
#' 
#' @rdname plot-isobaricnormRes
#' @export
#'

plot.isobaricnormRes <- function(isobaricnormRes_object, x_lab = NULL, ...) {
  require(ggplot2)
  .plot.isobaricnormRes(isobaricnormRes_object, x_lab = NULL, ...)
}

.plot.isobaricnormRes<- function(isobaricnormRes_object, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE, order = FALSE){
  
  #check for an isobaricnormRes object #
  if(!inherits(isobaricnormRes_object, "isobaricnormRes")) stop("object must be of class 'isobaricnormRes'")
  
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
  
  #extracting attributes from isobaricnormRes_object
  exp_cname = attr(isobaricnormRes_object, "isobaric_info")$exp_cname
  fdata_cname = attr(isobaricnormRes_object, "cnames")$fdata_cname
  
  #organize isobaricnormRes_object
  exp_cname_ind = which(names(isobaricnormRes_object) == exp_cname)
  fdata_cname_ind = which(names(isobaricnormRes_object) == fdata_cname)
  values_col_ind = which(names(isobaricnormRes_object) == "value")
  
  data = as.data.frame(isobaricnormRes_object[c(fdata_cname_ind, values_col_ind, exp_cname_ind)])
  data[[exp_cname]] = as.factor(data[[exp_cname]])
  
  #make labels
  xlabel <- ifelse(is.null(x_lab), "ref_sample", x_lab)
  ylabel <- ifelse(is.null(y_lab), "log2 Abundance", y_lab)
  plot_title <- ifelse(is.null(title_plot), "Reference Sample Profile", title_plot)
  
  p<- ggplot(data = data, aes(x = data[[fdata_cname]], y = value))
  
  if(order == TRUE){
    xlabel <- ifelse(is.null(x_lab),exp_cname, x_lab)
    p<- ggplot(data = data, aes(x = data[[exp_cname]], y = value)) 
    }
  
  p<- p + geom_boxplot() +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size)) 
  
  if(bw_theme == TRUE){
    p<- p + theme_bw()
  }
    

return(p)  
}



