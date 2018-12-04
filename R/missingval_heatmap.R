#' Produces a heatmap of missing data
#' 
#' This function takes in an omicsData object and returns a heatmap of omicsData$e_data
#'
#' @param omicsData an object of class "pepData", "proData", "metabData", or "lipidData", created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @param x_lab character string to be used for x-axis label. Defaults to NULL
#' @param y_lab character string to be used for y-axis label. Defaults to NULL 
#' @param ... further arguments
#' 
#' \tabular{ll}{
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{legend_title} \tab character string to be used for legend_title label. Defaults to NULL \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{palette} \tab is a character string indicating the name of the RColorBrewer palette to use; "YlOrRd", "YlOrBr", "YlGnBu", "YlGn", "Reds","RdPu", "Purples", "PuRd", "PuBuGn", "PuBu", "OrRd","Oranges", "Greys", "Greens", "GnBu", "BuPu","BuGn","Blues", "Set3", "Set2", "Set1", "Pastel2", "Pastel1", "Paired", "Dark2", "Accent", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr","PRGn", "PiYG", "BrBG"\cr
#' \code{x_lab_angle} \tab integer value indicating the angle of x-axis labels \cr
#' \code{coordinate_flip} \tab logical indicates whether to flip cartesian coordinates so that horizontal becomes vertical and vise versa, defaults to false \cr
#' }
#' 
#' @return plots ggplot2 object
#' 
#' @details This function takes in an omicsData object and returns a heatmap of omicsData$e_data, the colored tiles in the heatmap represent the intensity of the the e_data values, while the grey tiles represent NA values.  
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' library(ggplot2)
#' data("lipid_object")
#' 
#' lipid_object2 <- edata_transform(omicsData = lipid_object, data_scale="log2")
#' missingval_heatmap(lipid_object2, palette = "OrRd")
#' missingval_heatmap(lipid_object2, palette = "Reds", coordinate_flip = TRUE)
#'}
#' 
#' @rdname missingval_heatmap
#' @export
#'


missingval_heatmap <- function(omicsData, x_lab = NULL, y_lab = NULL, ...) {
  require(ggplot2)
  .missingval_heatmap(omicsData, x_lab, y_lab, ...)
}

.missingval_heatmap<- function(omicsData, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, palette = "YlOrRd", x_lab_angle = 60, coordinate_flip = FALSE){
  
  #check that omicsData is of correct class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  #check that omicsData data_scale attribute is log2
  if(attr(omicsData, "data_info")$data_scale != "log2") message("omicsData is not of log2 scale")
  
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

#here we call 'missingval_result' function that will give us an object of type naRes containing informaiton on the number of missing values per molecule   
na_Res<- missingval_result(omicsData)   
by_molecule<- na_Res$na.by.molecule
by_molecule<- by_molecule[order(by_molecule$num_NA),]

# make labels #
xlabel <- ifelse(is.null(x_lab), "Intensity", x_lab)
ylabel <- ifelse(is.null(y_lab), "Molecule", y_lab)
plot_title <- ifelse(is.null(title_plot), "Missing Values Heatmap", title_plot)
legendtitle<- ifelse(is.null(legend_title), "Intensity", legend_title)

#pull attr from omicsData
edata_cname<- attr(omicsData, "cnames")$edata_cname 
e_data<- omicsData$e_data

edata_melt<- melt(e_data, id.vars = edata_cname)
edata_melt[[edata_cname]]<- factor(edata_melt[[edata_cname]], levels = rev(by_molecule[[edata_cname]]))  
names(edata_melt)[1]<- "edata_cname"

  p <- ggplot(edata_melt, aes(x=variable, y = edata_cname)) + geom_tile(aes(fill = value)) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(plot_title) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size), axis.text.x = element_text(angle = x_lab_angle, hjust = 1)) +
    scale_fill_distiller(palette = palette, name = legendtitle)
  
  if(coordinate_flip == TRUE){
    p <- ggplot(edata_melt, aes(x=variable, y = edata_cname)) + geom_tile(aes(fill = value)) +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(plot_title) +
      coord_flip() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = title_size), axis.title.x = element_text(size = x_lab_size), axis.title.y = element_text(size = y_lab_size)) +
      scale_fill_distiller(palette = palette, name = legendtitle) 
  }

return(p)

}
