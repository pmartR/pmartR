#' Summarizes an object of class nmrnormRes
#' 
#' For creating a summary of an S3 object of type 'nmrnormRes':
#' 
#' @param nmrnormRes_object an object of type nmrnormRes, created by \code{\link{normalize_nmr}}  
#'
#' @return data.frame object
#' 
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object_identified, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, metabolite_name = "unkm1.53")
#' summary(nmr_norm)
#' 
#' # alternate specification: #
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, sample_property_cname = "Concentration")
#' summary(nmr_norm)
#' 
#'}
#' 
#'@export
#'@rdname summary-nmrnormRes
#'@name summary-pmartR
#'
summary.nmrnormRes<-function(nmrnormRes_object){
  
  #check for an nmrnormRes object #
  if(!inherits(nmrnormRes_object, "nmrnormRes")) stop("object must be of class 'nmrnormRes'")
  
  #extract attributes from nmrnomrRes_object
  sample_property_cname = attr(nmrnormRes_object, "nmr_info")$sample_property_cname
  metabolite_name <- attr(nmrnormRes_object, "nmr_info")$metabolite_name
  if(is.null(sample_property_cname)){
    normalized_using <- metabolite_name
  }else{
    if(is.null(metabolite_name)){
      normalized_using <- sample_property_cname
    }
  }
  #value_col_ind = which(names(nmrnormRes_object) == "value")
  
  res_median = median(nmrnormRes_object$value, na.rm = TRUE)
  res_sd = sd(nmrnormRes_object$value, na.rm = TRUE)
  
  final_res = as.data.frame(cbind(normalized_using, res_median, res_sd))
  row.names(final_res) = NULL
  
  if(is.null(sample_property_cname)){
    #normalized_using <- metabolite_name
    names(final_res) = c("Metabolite", "Median", "SD")
  }else{
    if(is.null(metabolite_name)){
      #normalized_using <- sample_property_name
      names(final_res) = c("Sample Property", "Median", "SD")
    }
  }
  
  
  return(final_res)
}
