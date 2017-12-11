#' Plots an object of class dataRes
#' 
#' For plotting an S3 object of type 'dataRes'
#'
#' @param dataRes an object of the class 'dataRes, usually created by \code{\link{summarize}}.
#' @param metric a character string indicating which metric to use in plot, one of 'mean', 'median', 'sd, 'min' or 'max'
#' 
#' @details This function can only create plots of dataRes objects whose 'by' == 'molecule' and 'groupvar' attribute is non NULL
#' @return plots ggplot2 object
#'
#' @examples
#' dontrun{
#' data(lipid_object)
#' lipid_object = edata_transform(lipid_object, "log2")
#' result = summarize(omicsData = lipid_object, by = "molecule", groupvar = "Condition")
#' plot(result)
#'}
#'
#'@rdname plot-dataRes
#' @export

plot.dataRes<- function(dataRes, metric = NULL){
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
    
    q<- ggplot(data_mean_sd, aes(x = mean, y = sd, color = variable)) + geom_point(shape = 1) + ggtitle("Mean x Standard Deviation") +
      theme_bw()
    
    p<- ggplot(data_mean_median, aes(x = mean, y = median, color = variable)) + geom_point(shape = 1) + ggtitle("Mean x Median") +
      theme_bw()
    
    grid.arrange(q, p, ncol = 2)
    
  }else{
    #histogram plots
    
    if(!(metric %in% c('mean', 'median','sd', 'pct_obs', 'min', 'max'))) stop("metric must be one of mean, median, sd, pct_obs, min or max")
    
    #subsetting dataRes object
    data = dataRes[[metric]]
    data_melt = melt(data, id.vars = edata_cname)
    
    #the number of groups, which is the number of histograms
    vars = levels(data_melt$variable)
    
    if(length(vars) == 1){
      r<- ggplot(data_melt, aes(x = value)) + geom_histogram(alpha = .5, fill = "red", binwidth = .5) + ggtitle(paste("Histogram for ", metric, sep = "")) + theme_bw()
    }
    else if(length(vars) == 2){
      r<- ggplot(data_melt, aes(value)) +
      geom_histogram(data = subset(data_melt, variable == vars[1]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[2]), aes(fill = variable), alpha = .2, binwidth = .5) +
      ggtitle(paste("Histogram for ", metric, sep = "")) + theme_bw()      
      }
    else if(length(vars) == 3){
      r<- ggplot(data_melt, aes(value)) + 
      geom_histogram(data = subset(data_melt, variable == vars[1]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[2]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[3]), aes(fill = variable), alpha = .2, binwidth = .5) +  
      ggtitle(paste("Histogram for ", metric, sep = "")) + theme_bw()
    }
    else if(length(vars) == 4){
      r<- ggplot(data_melt, aes(value)) +
      geom_histogram(data = subset(data_melt, variable == vars[1]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[2]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[3]), aes(fill = variable), alpha = .2, binwidth = .5) +
      geom_histogram(data = subset(data_melt, variable == vars[4]), aes(fill = variable), alpha = .2, binwidth = .5) +    
      ggtitle(paste("Histogram for ", metric, sep = "")) + theme_bw()
    }
    
    return(r)
  }
}