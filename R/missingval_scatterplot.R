#' missingval_scatterplot
#' 
#' mean intensity vs number of missing values (per molecule) 
#' 
#'@param palette is a character string indicating the name of the RColorBrewer palette to use; "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"
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
#' \code{point_size} \tab integer value indicating scatterplot point size, defaults to 3. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' \code{x_lab_angle} \tab integer value indicating the angle of x-axis labels \cr
#' \code{coordinate_flip} \tab logical indicates whether to flip cartesian coordinates so that horizontal becomes vertical and vise versa, defaults to false \cr
#' }
#' 
#' 
#'@rdname missingval_scatterplot
#'@export
#'

missingval_scatterplot <- function(omicsData, x_lab = NULL, y_lab = NULL, ...) {
  .missingval_scatterplot(omicsData, x_lab, y_lab, ...)
}

.missingval_scatterplot<- function(omicsData, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, point_size = 3, palette = "Spectral", bw_theme = FALSE, x_lab_angle = 0, coordinate_flip = FALSE){
  
  #check that omicsData is of correct class
  if(!(class(omicsData) %in% c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  #check that palette is in the list of RColorBrewer palettes
  if(!(palette %in% c("YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", 
                      "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", 
                      "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"))) stop("palette must be one of RColorBrewer palettes")
  
  #checking arguments are of correct class
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
  
  #extract from omicsData, e_data and attributes
  edata<- omicsData$e_data
  edata_cname<- attr(omicsData, "cnames")$edata_cname
  edata_cname_id<- which(names(edata) == edata_cname)
  
  #missing values results
  na_Res<- missingval_result(omicsData)
  bymolecule<- na_Res$na.by.molecule
  num_missing_vals<- bymolecule$num_NA
  
  # make labels #
  xlabel <- ifelse(is.null(x_lab), "Mean Intensity", x_lab)
  ylabel <- ifelse(is.null(y_lab), "Missing Values (per molecule)", y_lab)
  plot_title <- ifelse(is.null(title_plot), "Mean Intensity vs NA per Molecule", title_plot)
  legendtitle <- ifelse(is.null(legend_title), "Missing Values", legend_title)
  
  #check if group_designation has not been applied to omicsData
  if(is.null(attr(omicsData, "group_DF"))){
    
    mean_intensity<- apply(edata[, -edata_cname_id], 1, mean, na.rm = TRUE)
    
    plot_data<- as.data.frame(cbind(mean_intensity, num_missing_vals))
    
    if(bw_theme == FALSE){
    p<-ggplot2::ggplot(plot_data, aes(mean_intensity, num_missing_vals)) + geom_point(aes(color = num_missing_vals), size = point_size) + 
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) + 
      scale_color_distiller(palette = palette, name = legendtitle)
    
        if(coordinate_flip == TRUE){
          p = p + ggplot2::coord_flip()
        }
    
    }
    else{
      p<-ggplot2::ggplot(plot_data, aes(mean_intensity, num_missing_vals)) + geom_point(aes(color = num_missing_vals), size = point_size) + 
        ggplot2::theme_bw() +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab(ylabel) +
        ggplot2::ggtitle(plot_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
        scale_color_distiller(palette = palette, name = legendtitle)
      
        if(coordinate_flip == TRUE){
          p = p + ggplot2::coord_flip()
        }
      
    }
  }
  
  if(!is.null(attr(omicsData, "group_DF"))){
    
    group_df<-attr(omicsData, "group_DF")
    levels<- levels(group_df$Group)
    indices_list<- list()
    
    for(i in 1:length(levels)){
      inds<- which(group_df$Group == levels[i])
      indices_list[[i]]<- inds
    } 
    temp_edata<- edata[,-edata_cname_id]
    
    edata_by_group<- lapply(indices_list, function(x, temp_edata){temp_edata[,x]}, temp_edata = temp_edata)
    mean_by_group<- lapply(edata_by_group, rowMeans, na.rm = TRUE )
    names(mean_by_group)<- levels
    #NaN appear because the entire row to which rowMean is applied to is all NA values
    
    mean_intensity<- do.call(cbind, mean_by_group)
    
    plot_data<- cbind(num_missing_vals, mean_intensity)
    plot_data<- as.data.frame(plot_data)
    plot_data<- melt(plot_data, id.vars = "num_missing_vals")
    legendtitle <- ifelse(is.null(legend_title), "group", legend_title)
    
    if(bw_theme == FALSE){
      p<-ggplot2::ggplot(plot_data, aes(value, num_missing_vals)) + geom_point(aes(colour = variable), size = point_size) +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab(ylabel) +
        ggplot2::ggtitle(plot_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
        scale_color_brewer(palette = palette, name = legendtitle)
      
      if(coordinate_flip == TRUE){
        p = p + ggplot2::coord_flip()
      }
        
      
    }
    else{
      p<-ggplot2::ggplot(plot_data, aes(value, num_missing_vals)) + geom_point(aes(colour = variable), size = point_size) +
        ggplot2::theme_bw() +
        ggplot2::xlab(xlabel) +
        ggplot2::ylab(ylabel) +
        ggplot2::ggtitle(plot_title) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
        
        scale_color_brewer(palette = palette, name = legend_title)
      
        if(coordinate_flip == TRUE){
          p = p + ggplot2::coord_flip()
        }
    }
        
  }
  
   return(p) 
  
}

