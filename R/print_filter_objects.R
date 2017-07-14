#' print.moleculeFilt
#' 
#' For printing an S3 object of type 'moleculeFilt':
#' 
#'@rdname print-moleculeFilt
#'@export
#'
print.mofleculeFilt<- function(filter_object){
  if(class(filter_object)[1] != "moleculeFilt") stop("filter object must be of the class 'moleculeFilt'")
  
  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))
  
  to_filter<- rbind(filter_object_head, blank_row, filter_object_tail)
 
  cat("to_filter\n")
  cat(capture.output(to_filter), sep = "\n")
  cat("\n")
}

#' print.proteomicsFilt
#' 
#' For printing an S3 object of type 'proteomicsFilt':
#' 
#'@rdname print-proteomicsFilt
#'@export
#'
print.proteomicsFilt<- function(filter_object){
  if(class(filter_object)[1] != "proteomicsFilt") stop("filter object must be of the class 'proteomicsFilt'")

  
  
  
  
}