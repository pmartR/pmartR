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
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")
                                       if(firstn > length(temp[[1]])) stop(paste("there are less than", firstn, "characters in a sample name", sep = " "))
                                       temp2 = paste(temp[[1]][1:firstn], collapse = "")
                                       return(temp2)})
    output = unlist(output)
  }
  
  else if((!is.null(from)) & (!is.null(to))){
    if(!is.null(firstn) | !is.null(delim) | !is.null(components)) stop("only from and to arguments are needed")
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")
                                       if(!(from %in% 1:length(temp[[1]])) & !(to %in% 1:length(temp[[1]]))) stop(paste(from, to, "are not in the range of the length of at least one sample name", sep = " "))
                                       temp2 = paste(temp[[1]][from:to], collapse = "")
                                       return(temp2)})
    
    output = unlist(output)
  }
  
  else if((!is.null(delim)) & (!is.null(components))){
    if(!is.null(firstn) | !is.null(from) | !is.null(to)) stop("only delim and components arguments are needed")
    
    output = lapply(names, function(x){temp = strsplit(x, split = delim)
                                       if(length(components) > length(temp[[1]])) stop(paste("the length of 'components vector must be less than", lenght(temp[[1]]), sep = " "))
                                       temp2 = paste(temp[[1]][components], collapse = "")
                                       return(temp2)})
    
    output = unlist(output)
    
  }
  
  omicsData$f_data[["VizSampNames"]] = output
  
  return(omicsData)
}