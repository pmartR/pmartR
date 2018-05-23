#' Creates custom sample names for plots
#' 
#' This helper function creates custom sample names for plot data object function
#' 
#' @param omicsData omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{pepData}, \code{proData}, \code{metabData}, or \code{lipidData}, respectively.
#' @param sample_names is an optional argument that has three options for creating custom sample names, 'option1', 'option2','option3'. 'option1' takes the first n characters of each sample name to be the new sample names. Requires 'firstx' argument to be input. 'option2' takes a range of characters from each sample name to be the new sample names. Requires 'from' and 'to' arguments be input. 'option3' takes components separated by a delimiter in each sample name to be the new sample names. Requires 'delim' and 'components' arguments be input.
#' @param firstx is an integer specifying the first n characters to keep as the sample name
#' @param from is an integer specifying the start of the range of characters to keep as the sample name
#' @param to is an integer specifying the end of the range of characters to keep as the sample name
#' @param delim is a delimiter to separate sample name components by
#' @param components an integer vector specifying which components separated by delim to keep as sample name
#' 
#' @details 
#' 
#' @examples 
#' 
#' 

custom_sampnames = function(omicsData, sample_names = "option1", firstx = NULL, from = NULL, to = NULL, delim = NULL, components = NULL){
  #check that sample_names argument is one of "option1", "option2", "option3"
  if(!(sample_names %in% c("option1", "option2", "option3"))) stop("sample_names must be one of, 'option1', 'option2', or 'option3'")
  
  #extract sample names from omicsData object
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  names = as.list(as.character(omicsData$f_data[[fdata_cname]]))
  
  if(sample_names == "option1" & !is.null(firstx)){
    if(!is.null(from) | !is.null(to) | !is.null(delim) | !is.null(components)) stop("only firstx argument is needed with 'option1'")
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")
                                       if(firstx > length(temp[[1]])) stop(paste("there are less than", firstx, "characters in a sample name", sep = " "))
                                       temp2 = paste(temp[[1]][1:firstx], collapse = "")
                                       return(temp2)})
    output = unlist(output)
  }
  
  else if((sample_names == "option2") & (!is.null(from)) & (!is.null(to))){
    if(!is.null(firstx) | !is.null(delim) | !is.null(components)) stop("only from and to arguments are needed with 'option2'")
    
    output = lapply(names, function(x){temp = strsplit(x, split = "")
                                       if(!(from %in% 1:length(temp[[1]])) & !(to %in% 1:length(temp[[1]]))) stop(paste(from, to, "are not in the range of the length of at least one sample name", sep = " "))
                                       temp2 = paste(temp[[1]][from:to], collapse = "")
                                       return(temp2)})
    
    output = unlist(output)
  }
  
  else if(sample_names == "option3" & (!is.null(delim)) & (!is.null(components))){
    if(!is.null(firstx) | !is.null(from) | !is.null(to)) stop("only delim and components arguments are needed with 'option3'")
    
    output = lapply(names, function(x){temp = strsplit(x, split = delim)
                                       if(length(components) > length(temp[[1]])) stop(paste("the length of 'components vector must be less than", lenght(temp[[1]]), sep = " "))
                                       temp2 = paste(temp[[1]][components], collapse = "")
                                       return(temp2)})
    
    output = unlist(output)
    
  }
  
  omicsData$f_data[["plot_names"]] = output
  
  return(omicsData)
}