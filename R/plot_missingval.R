#' plot.naRes
#' 
#' For plotting an S3 object of type 'naRes':
#' 
#'@rdname plot-naRes
#'@export
#'

plot.naRes <- function(naRes_object, x_lab = NULL, ...) {
  .plot.naRes(naRes_object, x_lab, ...)
}

.plot.naRes<- function(naRes_object, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE){
 
   # check for a naRes object #
  if(class(naRes_object)[1] != "naRes") stop("object must be of class 'naRes'")
  
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
  
  #extracting items from naRes_object
  na.by.sample<- naRes_object$na.by.sample
  na.by.molecule<- naRes_object$na.by.molecule
  edata_cname<- attr(naRes_object, "cnames")$edata_cname
  fdata_cname<- attr(naRes_object, "cnames")$fdata_cname
  
  
  # make labels #
  xlabel <- ifelse(is.null(x_lab), "Sample Name", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Missing values (count)", y_lab)
  plot_title <- ifelse(is.null(title_plot), "Missing Values by Sample", title_plot)
  
  #plots NA per sample
 if(bw_theme == FALSE){
   p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA)) + 
    ggplot2::geom_bar(stat = "identity", width = .8, color= "black", fill = "steelblue") +
    ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
   
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size))
 }
  else{p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA)) + 
    ggplot2::geom_bar(stat = "identity", width = .8, color= "black", fill = "steelblue") +
    ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
    bw_theme() +
    
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size))
  }
  
  

  
  return(p)

}
  #plots NA per molecule
 # g<- ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
 #   ggplot2::geom_histogram(binwidth = 1, aes(fill = ..count..)) + 
 #   ggplot2::xlab("Number of NA values per Molecule") +
 #   ggplot2::ylab("Molecules (count)") +
 #   ggplot2::ggtitle("Missing Values per Molecule") +
 #   ggplot2::theme(plot.title = ggplot2::element_text(size = 10), axis.title.x = ggplot2::element_text(size = 10), axis.title.y = ggplot2::element_text(size = 10))
    
  
