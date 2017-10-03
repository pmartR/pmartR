#' missingval_heatmapplot
#' 
#'takes in omicsData and returns a heatmap of omicsData$e_data
#' 
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
#' \code{palette} \tab character string indicating the name of the RColorBrewer palette to use. \cr
#' }
#' 
#' 
#'@rdname missingval_heatmapplot
#'@export
#'


missingval_heatmapplot <- function(omicsData, x_lab = NULL, y_lab = NULL, ...) {
  .missingval_heatmapplot(omicsData, x_lab, y_lab, ...)
}

.missingval_heatmapplot<- function(omicsData, x_lab = NULL, y_lab = NULL, title_plot = NULL, legend_title = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, palette = "YlOrRd"){
  
  #check that omicsData is of correct class
  if(!(class(omicsData) %in% c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
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

#here we call 'missingval_result' function that will give us an object of type naRes containing informaiton on the number of missing values per molecule   
na_Res<- missingval_result(omicsData)   
by_molecule<- na_Res$na.by.molecule
by_molecule<- by_molecule[order(by_molecule$num_NA),]

# make labels #
xlabel <- ifelse(is.null(x_lab), "Intensity", x_lab)
ylabel <- ifelse(is.null(y_lab), "Molecule", y_lab)
plot_title <- ifelse(is.null(title_plot), "Missing Values Heatmap", title_plot)
legendtitle<- ifelse(is.null(legend_title), "Value", legend_title)

#pull attr from omicsData
edata_cname<- attr(omicsData, "cnames")$edata_cname 
e_data<- omicsData$e_data

edata_melt<- melt(e_data, id.vars = edata_cname)
edata_melt[[edata_cname]]<- factor(edata_melt[[edata_cname]], levels = rev(by_molecule[[edata_cname]]))  
names(edata_melt)[1]<- "edata_cname"

  p <- ggplot2::ggplot(edata_melt, aes(x=variable, y = edata_cname)) + geom_tile(aes(fill = value), colour = "white") +
    ggplot2::xlab(xlabel) +
    ggplot2::ylab(ylabel) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
    scale_fill_distiller(palette = palette, name = legendtitle)


return(p)

}
