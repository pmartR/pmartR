#' Plots an object of class naRes
#' 
#' For plotting an S3 object of type 'naRes':
#' 
#' @param naRes_object a list of two data frames, one contains NA values by sample, the second contains NA values by molecule
#' @param type is for specifying plot type, there are three options, 'bySample': plots missing values per sample, 'byMolecule': plots missing values per molecule and 'Both': displays both 'bySample' and 'byMolecule' plots. 
#' @param x_lab character string to be used for x-axis label. Defaults to NULL  
#' @param ... further arguments 
#'
#' \tabular{ll}{
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{legend_title} \tab character string to be used for legend_title label. Defaults to NULL \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bar_width} \tab integer value indicating the bar width in a barplot (when type = "bySample"). Defaults to .8. \cr
#' \code{binwidth} \tab integer value indicating the bin width in a histogram (when type = "byMolecule"). Defaults to 1. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{palette} \tab palette is a character string indicating the name of the RColorBrewer palette to use; "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"\cr
#' \code{x_lab_angle} \tab integer value indicating the angle of x-axis labels \cr
#' \code{coordinate_flip} \tab logical indicates whether to flip cartesian coordinates so that horizontal becomes vertical and vise versa, defaults to false \cr
#'}
#' 
#' @return plots ggplot2 object
#' 
#' @details This function takes in an object of class naRes and creates either a bar chart or histogram of missing values. When parameter type = "bySample", a "sample name" by "missing values count" bar chart is returned. When parameter type = "byMolecule" a "number of missing values per molecule" by "Molecules" histogram is returned. If type is set to "Both", a bar chart and histogram will be returned. 
#' 
#' @examples
#' dontrun{
#' library(gridExtra)
#' library(pmartRdata)
#' data("lipid_object")
#' result<- missingval_result(lipid_object)
#' plot(result, type = "bySample", x_lab_angle = 50, palette = "Set1")
#' 
#' plot(result, type = "Both", palette = "Set2")
#'}
#' 
#' @rdname plot-naRes
#' @export
#'

plot.naRes <- function(naRes_object, type, x_lab = NULL, ...) {
  .plot.naRes(naRes_object, type, x_lab, ...)
}

.plot.naRes<- function(naRes_object, type, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bar_width = .8, binwidth = 1, bw_theme = FALSE, palette = "Spectral", x_lab_angle = 60, coordinate_flip = FALSE){
 
   # check for a naRes object #
  if(!inherits(naRes_object, "naRes")) stop("object must be of class 'naRes'")
  
  #check that type is one of the three options
  if(!(type %in% c("bySample", "byMolecule", "Both"))) stop("type must be one of 'bySample', 'byMolecule', 'Both'")
  
  #check that if type "Both" is specified, then x_lab, y_lab and tile_plot are all NULL
  if(type == "Both" & (!is.null(x_lab)|!is.null(y_lab)|!is.null(title_plot)| !is.null(legend_title))) stop("if type 'Both' is specified, x_lab, y_lab, legend_title and title_plot cannot be specified")
  
  if(!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", 
                      "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", 
                      "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) stop("palette must be one of RColorBrewer palettes")
  
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
  if(!(is.numeric(x_lab_angle))) stop("x_lab_angle must be numeric")
  
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
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "black", size = 3) +
          ggplot2::xlab(xlabel) +
          ggplot2::ylab(ylabel) +
          ggplot2::ggtitle(plot_title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = legendtitle) 
      
      if(coordinate_flip == TRUE){
        p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
            ggplot2::geom_bar(stat = "identity", width = bar_width) +
            ggplot2::geom_text(aes(label = num_NA), hjust = 2, color = "black", size = 3) +
            ggplot2::xlab(xlabel) +
            ggplot2::ylab(ylabel) +
            ggplot2::ggtitle(plot_title) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = legendtitle) + 
            ggplot2::coord_flip()
      }
      
    }
    else{
      p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "black", size = 3) +
          ggplot2::theme_bw() +
          ggplot2::xlab(xlabel) +
          ggplot2::ylab(ylabel) +
          ggplot2::ggtitle(plot_title) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = legendtitle)
         
          if(coordinate_flip == TRUE){
            p<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
                ggplot2::geom_bar(stat = "identity", width = bar_width) +
                ggplot2::geom_text(aes(label = num_NA), hjust = 2, color = "black", size = 3) +
                ggplot2::theme_bw() +
                ggplot2::xlab(xlabel) +
                ggplot2::ylab(ylabel) +
                ggplot2::ggtitle(plot_title) +
                ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = legendtitle) + 
                ggplot2::coord_flip()
          } 
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
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size),axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
           scale_fill_distiller(palette = palette, name = legendtitle)  
      
          if(coordinate_flip == TRUE){
            p = p + ggplot2::coord_flip()
          }
      }
    else{
      p<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) +      
           theme_bw() +
           ggplot2::xlab(xlabel) +
           ggplot2::ylab(ylabel) +
           ggplot2::ggtitle(plot_title) +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
           scale_fill_distiller(palette = palette, name = legendtitle)  
      
            if(coordinate_flip == TRUE){
              p = p + ggplot2::coord_flip()
            }
    }
    return(p)
  }
  
  #if type = "Both"
  else{
    
    #plots NA per sample
    if(bw_theme == FALSE){
      s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "black", size = 3) +
          ggplot2::xlab("Sample Name") +
          ggplot2::ylab("Missing values (count)") +
          ggplot2::ggtitle("Missing Values by Sample") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = "group")
      
          if(coordinate_flip == TRUE){
            s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
                ggplot2::geom_bar(stat = "identity", width = bar_width) +
                ggplot2::geom_text(aes(label = num_NA), hjust = 2, color = "black", size = 3) +
                ggplot2::xlab("Sample Name") +
                ggplot2::ylab("Missing values (count)") +
                ggplot2::ggtitle("Missing Values by Sample") +
                ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = "group") +
                ggplot2::coord_flip()
          }
    }
    else{
      s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
          ggplot2::geom_bar(stat = "identity", width = bar_width) +
          ggplot2::geom_text(aes(label = num_NA), vjust = 2, color = "black", size = 3) +
          theme_bw() +
          ggplot2::xlab("Sample Name") +
          ggplot2::ylab("Missing values (count)") +
          ggplot2::ggtitle("Missing Values by Sample") +
          ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = "group")
      
          if(coordinate_flip == TRUE){
            s<- ggplot2::ggplot(data=na.by.sample, aes(x=na.by.sample[[fdata_cname]], y = na.by.sample$num_NA, fill = group)) + 
                ggplot2::geom_bar(stat = "identity", width = bar_width) +
                ggplot2::geom_text(aes(label = num_NA), hjust = 2, color = "black", size = 3) +
                theme_bw() +
                ggplot2::xlab("Sample Name") +
                ggplot2::ylab("Missing values (count)") +
                ggplot2::ggtitle("Missing Values by Sample") +
                ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + scale_fill_brewer(palette = palette, name = "group") +
                ggplot2::coord_flip()
          }
    }
    
    #plots NA per molecule
    if(bw_theme == FALSE){
      m<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) + 
           ggplot2::xlab("Number of NA values per Molecule") +
           ggplot2::ylab("Molecules (count)") +
           ggplot2::ggtitle("Missing Values per Molecule") +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
           scale_fill_distiller(palette = palette, name = "Count")
      
          if(coordinate_flip == TRUE){
            m = m + ggplot2::coord_flip()
          }
    }
    else{
      m<-  ggplot2::ggplot(data= na.by.molecule, aes(x = na.by.molecule$num_NA)) +
           ggplot2::geom_histogram(binwidth = binwidth, aes(fill = ..count..)) +      
           theme_bw() +
           ggplot2::xlab("Number of NA values per Molecule") +
           ggplot2::ylab("Molecules (count)") +
           ggplot2::ggtitle("Missing Values per Molecule") +
           ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
           scale_fill_distiller(palette = palette, name = "Count")
     
           if(coordinate_flip == TRUE){
             m = m + ggplot2::coord_flip()
           }
       }
    
     grid.arrange(s, m, ncol = 2)
    
  }
}
