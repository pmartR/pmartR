#' Creates custom sample names for plots
#' 
#' This helper function creates custom sample names for plot data object function
#' 
#' @param omicsData omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{pepData}, \code{proData}, \code{metabData}, or \code{lipidData}, respectively.
#' @param firstn is an integer specifying the first n characters to keep as the sample name
#' @param from is an integer specifying the start of the range of characters to keep as the sample name
#' @param to is an integer specifying the end of the range of characters to keep as the sample name
#' @param delim is a delimiter to separate sample name components by
#' @param components an integer vector specifying which components separated by delim to keep as sample name
#' 
#' @details 
#' 
#' @examples
#' dontrun{
#' data(pep_object)
#' results = custom_sampnames(pep_object, firstn = 5)
#' } 
#' 
#' @export

custom_sampnames = function(omicsData, firstn = NULL, from = NULL, to = NULL, delim = NULL, components = NULL){
  
  #extract sample names from omicsData object
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  names = as.list(as.character(omicsData$f_data[[fdata_cname]]))
  
  if(!is.null(firstn)){
    if(!is.null(from) | !is.null(to) | !is.null(delim) | !is.null(components)) stop("only firstn argument is needed")
    if(!is.numeric(firstn)) stop("firstn must be a non-negative numeric value less than the number of characters in the smallest sample name")
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")[[1]]
                                       if(firstn > length(temp)) warning(paste("there are less than", firstn, "characters in a sample name, nothing was truncated", sep = " "))
                                       temp <- temp[1:firstn]
                                       temp2 = paste(temp[!is.na(temp)], collapse = "")
                                       return(temp2)})
    output = unlist(output)
  }
  
  else if((!is.null(from)) & (!is.null(to))){
    if(!is.null(firstn) | !is.null(delim) | !is.null(components)) stop("only from and to arguments are needed")
    if(!is.numeric(from) | !is.numeric(to)) stop("'from' and 'to' must be non-negative numeric values")
    if(to < from) stop("'to' must be less than 'from'")
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")
                                       if(!(from %in% 1:length(temp[[1]])) & !(to %in% 1:length(temp[[1]]))) stop(paste(from, to, "are not in the range of the length of at least one sample name", sep = " "))
                                       temp2 = paste(temp[[1]][from:min(to, length(temp[[1]]))], collapse = "")
                                       return(temp2)})
    
    output = unlist(output)
  }
 
  else if((!is.null(delim)) & (!is.null(components))){
    if(!is.null(firstn) | !is.null(from) | !is.null(to)) stop("only delim and components arguments are needed")
    
    output = lapply(names, function(x){temp = strsplit(x, split = delim)[[1]]
                                       if(length(components) > length(temp)) stop(paste("the length of 'components vector must be less than", length(temp), sep = " "))
                                       if(!any(components %in% 1:length(temp))) stop("none of the indices specified in 'components' match indices of the split sample name")
                                       temp <- temp[components]
                                       temp2 = paste(temp[!is.na(temp)], collapse = delim)
                                       return(temp2)})
    
    output = unlist(output)
    
  }
  
  omicsData$f_data[["VizSampNames"]] = output
  
  return(omicsData)
}