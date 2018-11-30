#' Plots an object of class dataRes
#' 
#' For plotting an S3 object of type 'dataRes'
#'
#' @param dataRes an object of the class 'dataRes, usually created by \code{\link{summarize}}.
#' @param metric a character string indicating which metric to use in plot, one of 'mean', 'median', 'sd, 'min' or 'max'
#' @param density logical default to FALSE, if set to true a density plot of the specified metric is returned
#' @param ncols defaults to NULL, is an integer the specifies the number of columns for the histogram facet_wrap
#' 
#' @details This function can only create plots of dataRes objects whose 'by' == 'molecule' and 'groupvar' attribute is non NULL
#' @return plots ggplot2 object
#'
#' @examples
#' dontrun{
#' data(lipid_object)
#' lipid_object = edata_transform(lipid_object, "log2")
#' result = edata_summary(omicsData = lipid_object, by = "molecule", groupvar = "Condition")
#' plot(result)
#'}
#'
#'@rdname plot-dataRes
#' @export

plot.dataRes<- function(dataRes, metric = NULL, density = FALSE, ncols = NULL){
  #check that attr(dataRes, "by") == "molecule"
  if(attr(dataRes, "by") != "molecule") stop("can only plot a dataRes object if its 'by' attribute is equal to 'molecule'")
  
  #check that attr(dataRes, "groupvar") is not NULL
  if(is.null(attr(dataRes, "groupvar"))) stop("can only plot a dataRes object if its 'groupvar' attribute is not NULL")
  
  #extract molecule name
  edata_cname = attr(dataRes, "cnames")$edata_cname
  
  if(is.null(metric)){
    #scatterplots
    #subseting dataRes object
    mean = dataRes$mean
    median = dataRes$median
    sd = dataRes$sd
    
    #melting data frames from dataRes object
    mean_melt = melt(mean, id.vars = edata_cname)
    names(mean_melt)[3]<- "mean"
    sd_melt = melt(sd, id.vars = edata_cname)
    names(sd_melt)[3]<- "sd"
    median_melt = melt(median, id.vars = edata_cname)
    names(median_melt)[3]<- "median"
    
    data_mean_sd<- merge(mean_melt, sd_melt, by = c(edata_cname, "variable"))
    data_mean_median<- merge(mean_melt, median_melt, by = c(edata_cname, "variable"))
    
    q<- ggplot2::ggplot(data_mean_sd, ggplot2::aes(x = mean, y = sd, color = variable)) + ggplot2::geom_point(shape = 1) + ggplot2::ggtitle("Mean x Standard Deviation") +
      theme_bw()
    
    p<- ggplot2::ggplot(data_mean_median, ggplot2::aes(x = mean, y = median, color = variable)) + ggplot2::geom_point(shape = 1) + ggplot2::ggtitle("Mean x Median") +
      theme_bw()
    
    gridExtra::grid.arrange(q, p, ncol = 2)
    
  }
  
  else if(!is.null(metric)){
    if(!(metric %in% c('mean', 'median','sd', 'pct_obs', 'min', 'max'))) stop("metric must be one of mean, median, sd, pct_obs, min or max")
    if(!is.logical(density)) stop("density must be either TRUE or FALSE")
    
    #if density == F, will plot faceted histograms
    if(density == FALSE){
      #subsetting dataRes object
      data = dataRes[[metric]]
      data_melt = melt(data, id.vars = edata_cname)
      
      r<- ggplot2::ggplot(data_melt, ggplot2::aes(x = value, fill = variable)) + ggplot2::geom_histogram(binwidth = .5, colour = "white") +
        ggplot2::facet_wrap(~variable, ncol = ncols) +
        ggplot2::ggtitle(paste("Histograms for ", metric, sep = "")) + ggplot2::theme_bw()
      
    }else{
      #if density == T, will plot geom_density
      data = dataRes[[metric]]
      data_melt = melt(data, id.vars = edata_cname)
      
      r<- ggplot2::ggplot(data_melt, ggplot2::aes(x = value, colour = variable)) + ggplot2::geom_density() +
        ggplot2::ggtitle(paste("Density plot for ", metric, sep = "")) + ggplot2::theme_bw()
    }

    return(r)
  }
}