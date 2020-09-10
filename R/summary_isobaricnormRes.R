#' Summarizes an object of class isobaricnormRes
#' 
#' For creating a summary of an S3 object of type 'isobaricnormRes':
#' 
#' @param isobaricnormRes_object an object of type isobaricnormRes, created by \code{\link{normalize_isobaric}}  
#'
#' @return data.frame object
#' 
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(isobaric_object)
#' 
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' result = normalize_isobaric(isobaric_object, exp_cname = "Set", apply_norm = FALSE, channel_cname = "iTRAQ.Channel", refpool_channel = "116")
#' 
#' summary(result)
#'}
#' 
#'@export
#'@rdname summary-isobaricnormRes
#'@name summary-pmartR
#'
summary.isobaricnormRes<-function(isobaricnormRes_object){
  
  #check for an isobaricnormRes object #
  if(!inherits(isobaricnormRes_object, "isobaricnormRes")) stop("object must be of class 'isobaricnormRes'")
  
  #extract attributes from isobaricnomrRes_object
  exp_cname = attr(isobaricnormRes_object, "isobaric_info")$exp_cname
  
  exp_cname_ind = which(names(isobaricnormRes_object) %in% exp_cname)
  value_col_ind = which(names(isobaricnormRes_object) == "value")
  
  #subset data columns
  data = as.data.frame(isobaricnormRes_object[c(exp_cname_ind, value_col_ind)])
  
  split_data = split(data, data[[exp_cname]])
  
  res_median = lapply(split_data, function(item){median(item[["value"]], na.rm = T)})
  res_sd = lapply(split_data, function(item){sd(item[["value"]], na.rm = T)})
  
  final_res = as.data.frame(cbind(names(res_median), res_median, res_sd))
  row.names(final_res) = NULL
  names(final_res) = c(exp_cname, "Median", "SD")
  
  return(final_res)
}