#' plot.naRes
#' 
#' For plotting an S3 object of type 'naRes':
#' 
#'
#'@param type is for specifying plot type, there are three options, "bySample", "byMolecule" and "Both"
#'
#'
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{legend_title} \tab character string to be used for legend_title label. Defaults to NULL \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bar_width} \tab integer value indicating the bar width in a barplot (when type = "bySample"). Defaults to .8. \cr
#' \code{binwidth} \tab integer value indicating the bin width in a histogram (when type = "byMolecule"). Defaults to 1. \cr
#' \code{palette} \tab character string indicating the name of the RColorBrewer palette to use. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
#' 
#' 
#'@rdname plot-naRes
#'@export
#'

plot.naRes <- function(naRes_object, type, x_lab = NULL, ...) {
  .plot.naRes(naRes_object, type, x_lab, ...)
}

.plot.naRes<- function(naRes_object, type, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bar_width = .8, binwidth = 1, bw_theme = FALSE, palette = "Spectral"){
 
   # check for a naRes object #
  if(class(naRes_object)[1] != "naRes") stop("object must be of class 'naRes'")
  
  #check that type is one of the three options
  if(!(type %in% c("bySample", "byMolecule", "Both"))) stop("type must be one of 'bySample', 'byMolecule', 'Both'")
  
  #check that if type "Both" is specified, then x_lab, y_lab and tile_plot are all NULL
  if(type == "Both" & (!is.null(x_lab)|!is.null(y_lab)|!is.null(title_plot))) stop("if type 'Both' is specified, x_lab, y_lab and title_plot cannot be specified")
  
  if(!is.null(title_plot)) {
    if(!is.character(title_plot)) stop("title_plot must be a character vector")
  }
  if(!is.null(x_lab)) {
    if(!is.character(x_lab)) stop("x_lab must be a character vector")
  }
  if(!is.null(y_lab)) {
    if(!is.character(y_lab)) stop("y_lab must be a character vector")
  }
  if(!is.null(legend_title)) {
    if(!is.character(legend_title)) stop("legend_title must be a character vector")
  }
  
  ## end of initial checks ##
  
  #extracting items from naRes_object
  na.by.sample<- naRes_object$na.by.sample
  na.by.molecule<- naRes_object$na.by.molecule
  edata_cname<- attr(naRes_object, "cnames")$edata_cname
  fdata_cname<- attr(naRes_object, "cnames")$fdata_cname
  
  if(type == "bySample"){
    # make labels #
    xlabel <- ifelse(is.null(x_lab), "Sample Name", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Missing values (count)", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Missing Values by Sample", title_plot)
    legendtitle<- ifelse(is.null(legend_title), "group", legend_title)
    
    
    #plots NA per sample
    if(bw_theme == FALSE){
      p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
        
          ggplot2::xlab(xlabel) +
          ggplot2::ylab(ylabel) +
          ggplot2::ggtitle(plot_title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) + scale_fill_brewer(palette = palette, name = legendtitle )

    }
    else{
      p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
          theme_bw() +
      
          ggplot2::xlab(xlabel) +
          ggplot2::ylab(ylabel) +
          ggplot2::ggtitle(plot_title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) + scale_fill_brewer(palette = palette, name = legendtitle)
    }
    return(p)
  }
 
  else if(type == "byMolecule"){
    # make labels #
    xlabel <- ifelse(is.null(x_lab), "Number of NA values per Molecule", x_lab)
    ylabel <- ifelse(is.null(y_lab), "Molecules (count)", y_lab)
    plot_title <- ifelse(is.null(title_plot), "Missing Values per Molecule", title_plot)
    legendtitle <- ifelse(is.null(legend_title), "Count", legend_title)
  
    #plots NA per molecule
    if(bw_theme == FALSE){
      p<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) + 
         
           ggplot2::xlab(xlabel) +
           ggplot2::ylab(ylabel) +
           ggplot2::ggtitle(plot_title) +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
           scale_fill_distiller(palette = palette, name = legendtitle)  
      }
    else{
      p<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) +      
           theme_bw() +
      
           ggplot2::xlab(xlabel) +
           ggplot2::ylab(ylabel) +
           ggplot2::ggtitle(plot_title) +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
           scale_fill_distiller(palette = palette, name = legendtitle)  
    }
    return(p)
  }
  
  #if type = "Both"
  else{
    
    #plots NA per sample
    if(bw_theme == FALSE){
      s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
        
          ggplot2::xlab("Sample Name") +
          ggplot2::ylab("Missing values (count)") +
          ggplot2::ggtitle("Missing Values by Sample") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) + scale_fill_brewer(palette = palette, name = "group")
    }
    else{
      s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "white", size = 3) +
          theme_bw() +
        
          ggplot2::xlab("Sample Name") +
          ggplot2::ylab("Missing values (count)") +
          ggplot2::ggtitle("Missing Values by Sample") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) + scale_fill_brewer(palette = palette, name = "group")
    }
    
    #plots NA per molecule
    if(bw_theme == FALSE){
      m<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) + 
        
           ggplot2::xlab("Number of NA values per Molecule") +
           ggplot2::ylab("Molecules (count)") +
           ggplot2::ggtitle("Missing Values per Molecule") +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
           scale_fill_distiller(palette = palette, name = "Count")
    }
    else{
      m<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) +      
           theme_bw() +
        
           ggplot2::xlab("Number of NA values per Molecule") +
           ggplot2::ylab("Molecules (count)") +
           ggplot2::ggtitle("Missing Values per Molecule") +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
           scale_fill_distiller(palette = palette, name = "Count")
       }
    
     grid.arrange(s, m, ncol = 2)
    
  }
}
