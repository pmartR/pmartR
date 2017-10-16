#' Return comparisons of statRes object
#' 
#' This function returns comparisons from statRes object
#' 
#' @param statRes is an object of class 'statRes'
#' 
#' @return returns a data frame with comparisons and their indices
#'
#' @examples 
#' dontrun{
#' get_comparisons(statRes_object)
#'
#'}
#' 
#' @rdname get_comparisons
#' @export
#'

get_comparisons<- function(statRes){

  #check that statRes object is of 'statRes' class
  if(class(statRes) != "statRes") stop("object must be of class 'statRes'")
  
  #pull comparisons attribute
  comp = attr(statRes, "comparisons")
  
  result = data.frame("comparisons" = as.character(comp), "index" = 1:length(comp), stringsAsFactors = FALSE)
  
  return(result)
  
}

