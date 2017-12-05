#' Creates a data frame displaying multiple metrics
#'
#' This function takes in an object of class 'dataRes' and returns a data frame displaying a combination of metrics. The six summarizing metrics include, mean, standard deviation, median, percent observed, minimum and maximum.   
#'
#' @param statRes an object of the class 'statRes', created by \code{\link{summarize}}.
#' @param minmax logical, for specifying whether or not to include minimum and maximum data in the returned data frame. Defaults to FALSE.
#' 
#' @details If the 'by' attribute of the input 'statRes' object is set to 'sample', then its groupvar attribute must be set to NULL
#' 
#' @return prints a data frame
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object2 <- edata_transform(omicsData = lipid_object, data_scale = "log2")
#' 
#' dataRes_sample = summarize(omicsData = lipid_object2, groupvar = NULL, by = "sample")
#' report_dataRes(dataRes_sample)
#'}
#'
#'
#' @export

report_dataRes<- function(dataRes, minmax = FALSE, digits = 2){
  #checks
  if(class(dataRes) != "dataRes") stop("dataRes must be of class dataRes")
  if(!is.logical(minmax)) stop("minmax must be of class 'logical'")
  if(digits <= 0) stop("digits must be an integer greater than zero")
  
  #extracting some dataRes attr
  edata_cname = attr(dataRes, "cnames")$edata_cname
  
  if(attr(dataRes, "by") == "molecule"){
      if(!is.null(attr(dataRes, "groupvar"))){
        
        #subsetting dataRes object
        avg = dataRes$mean
        avg[, 2:ncol(avg)]<- round(avg[,-1], digits = digits)
        
        std_div = dataRes$sd
        std_div[, 2:ncol(std_div)]<- round(std_div[, -1], digits = digits)
        
        minimum = dataRes$min
        minimum[, 2:ncol(minimum)]<- round(minimum[, -1], digits = digits)
        
        maximum = dataRes$max
        maximum[, 2:ncol(maximum)]<- round(maximum[, -1], digits = digits)
        
        #organizing data
        mean_melt = melt(avg, id.vars = edata_cname)
        names(mean_melt)[3]<- "mean"
        sd_melt = melt(std_div, id.vars = edata_cname)
        names(sd_melt)[3]<- "sd"
        min_melt = melt(minimum, id.vars = edata_cname)
        names(min_melt)[3]<- "min"
        max_melt = melt(maximum, id.vars = edata_cname)
        names(max_melt)[3]<- "max"
        
        data_mean_sd<- merge(mean_melt, sd_melt, by = c(edata_cname, "variable"))
        data_min_max<- merge(min_melt, max_melt, by = c(edata_cname, "variable"))
        all_data = merge(data_mean_sd, data_min_max, by = c(edata_cname, "variable"))
        mean_sd_cols = all_data[,c("mean", "sd")]
        first_cols = data_mean_sd[, -which(names(data_mean_sd) %in% names(mean_sd_cols))]
        new_col = apply(mean_sd_cols, 1, function(x){paste(x[1], " \u00B1 ", x[2], sep = "")})
        
        res = cbind(first_cols, new_col) 
        names(res)[3]<- "value"
        formula2 = paste(edata_cname, "~...", sep = "")
        res = dcast(res, formula = formula2)
        
        if(minmax == TRUE){
          min_max_cols = all_data[,c("min", "max")]
          new_col2 = apply(min_max_cols, 1, function(x){paste(x[1], x[2], sep = "/")})
          
          res2 = cbind(first_cols, paste(new_col, new_col2, sep = '\n')) 
          names(res2)[3]<- "value"
          res2 = dcast(res2, formula = formula2)
          
          res = res2
        }
      }else{
        #subsetting dataRes object
        avg = dataRes$mean
        avg$mean<- round(avg[,-which(names(avg) == edata_cname)], digits = digits)
        
        std_div = dataRes$sd
        std_div$sd<- round(std_div[,-which(names(std_div) == edata_cname)], digits = digits)
        
        minimum = dataRes$min
        minimum$min<- round(minimum[,-which(names(minimum) == edata_cname)], digits = digits)
        
        maximum = dataRes$max
        maximum$max<- round(maximum[,-which(names(maximum) == edata_cname)], digits = digits)
        
        data_mean_sd<- merge(avg, std_div, by = edata_cname)
        data_min_max = merge(minimum, maximum, by = edata_cname)
        all_data = merge(data_mean_sd, data_min_max, by = edata_cname)
        
        mean_sd_cols = all_data[,c("mean", "sd")]
        first_cols = as.data.frame(data_mean_sd[, -which(names(data_mean_sd) %in% names(mean_sd_cols))])
        names(first_cols)[1]<- edata_cname
        new_col = apply(mean_sd_cols, 1, function(x){paste(x[1], " \u00B1 ", x[2], sep = "")})
        
        res = cbind(first_cols, new_col) 
        names(res)[2]<- "value"
        
        if(minmax == TRUE){
          min_max_cols = all_data[,c("min", "max")]
          new_col2 = apply(min_max_cols, 1, function(x){paste(x[1], x[2], sep = "/")})
          res2 = cbind(first_cols, paste(new_col, new_col2, sep = "\n"))
          names(res2)[2]<- "value"
          res = res2
        }
      }
  }
  
  else if(attr(dataRes, "by") == "sample"){
    #check that groupvar is NULL
    if(!is.null(attr(dataRes, "groupvar"))) stop("if 'by' attribute is set to 'sample','groupvar' attribute must be NULL")
    
    #minmax defaults to FALSE
    
    #subsetting dataRes object and applying round function with digits argument
    avg = dataRes$mean
    avg$mean<- round(avg[,-1], digits = digits)
    
    std_div = dataRes$sd
    std_div$sd<- round(std_div[,-1], digits = digits)
  
    minimum = dataRes$min
    minimum$min = round(minimum[,-1], digits = digits)
    
    maximum = dataRes$max
    maximum$max = round(maximum[,-1], digits = digits)
    
    data_mean_sd<- merge(avg, std_div, by = "sample")
    data_min_max = merge(minimum, maximum, by = "sample")
    all_data = merge(data_mean_sd, data_min_max, by = "sample")
    
    mean_sd_cols = all_data[,c("mean", "sd")]
    first_cols = as.data.frame(data_mean_sd[, -which(names(data_mean_sd) %in% names(mean_sd_cols))])
    names(first_cols)<- "sample"
    new_col = apply(mean_sd_cols, 1, function(x){paste(x[1], " \u00B1 ", x[2], sep = "")})
    
    res = cbind(first_cols, new_col) 
    names(res)[2]<- paste("mean", " \u00B1 ", "sd", sep = "") 
    
    if(minmax == TRUE){
      min_max_cols = all_data[,c("min", "max")]
      new_col2 = apply(min_max_cols, 1, function(x){paste(x[1], x[2], sep = "/")})
      
      res2 = cbind(first_cols, paste(new_col, new_col2, sep = "\n")) 
      names(res2)[2]<- paste("mean", " \u00B1 ", "sd", " (min/max)", sep = "") 
      
      res = res2
    }
  }
  return(knitr::kable(res))
}