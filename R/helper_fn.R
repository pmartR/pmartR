#' get_comparisons
#' 
#' returns comparisons given statRes object
#' 
#' 
#'@rdname get_comparisons
#'@export
#'

get_comparisons<- function(statRes){

  #check that statRes object is of 'statRes' class
  if(class(statRes) != "statRes") stop("object must be of class 'statRes'")
  
  #pull comparisons attribute
  comp = attr(statRes, "comparisons")
  
  result = data.frame("comparisons" = as.character(comp), "index" = 1:length(comp), stringsAsFactors = FALSE)
  
  return(result)
  
}


