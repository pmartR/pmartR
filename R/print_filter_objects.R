#' print.moleculeFilt
#' 
#' For printing an S3 object of type 'moleculeFilt':
#' 
#'@rdname print-moleculeFilt
#'@export
#'
print.moleculeFilt<- function(filter_object){
  if(class(filter_object)[1] != "moleculeFilt") stop("filter object must be of the class 'moleculeFilt'")
  filter_object<- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)
  
  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))
  
  result<- rbind(filter_object_head, blank_row, filter_object_tail)
  
  cat("moleculeFilt object\n")
  cat(capture.output(result), sep = "\n")
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
  
  counts_by_pep<- as.data.frame(lapply(filter_object$counts_by_pep, as.character), stringsAsFactors = FALSE)
  counts_by_pro<- as.data.frame(lapply(filter_object$counts_by_pro, as.character), stringsAsFactors = FALSE)
  
  counts_by_pep_head<- head(counts_by_pep, 4)[1:ncol(counts_by_pep)]
  counts_by_pep_tail<- tail(counts_by_pep, 4)[1:ncol(counts_by_pep)]
  blank_row = rep("---", ncol(counts_by_pep))
  
  counts_by_pro_head<- head(counts_by_pro, 4)[1:ncol(counts_by_pro)]
  counts_by_pro_tail<- tail(counts_by_pro, 4)[1:ncol(counts_by_pro)]
  blank_row2 = rep("---", ncol(counts_by_pro))
  
  bypep<- rbind(counts_by_pep_head, blank_row, counts_by_pep_tail)
  bypro<- rbind(counts_by_pro_head, blank_row2, counts_by_pro_tail)
  
  cat("proteomicsFilt object\n")
  cat("counts_by_pep\n")
  cat(capture.output(bypep), sep = "\n")
  cat("\n")
  
  cat("counts_by_pro\n")
  cat(capture.output(bypro), sep= "\n")
  cat("\n")
  
}

#' print.imdanovaFilt
#' 
#' For printing an S3 object of type 'imdanovaFilt':
#' 
#'@rdname print-imdanovaFilt
#'@export
#'
print.imdanovaFilt<- function(filter_object){
  if(class(filter_object)[1] != "imdanovaFilt") stop("filter object must be of the class 'imdanovaFilt'")
  filter_object<- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)
  
  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))
  
  result<- rbind(filter_object_head, blank_row, filter_object_tail)
  
  cat("imdanovaFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.rmdFilt
#' 
#' For printing an S3 object of type 'rmdFilt':
#' 
#'@rdname print-rmdFilt
#'@export
#'
print.rmdFilt<- function(filter_object){
  if(class(filter_object)[1] != "rmdFilt") stop("filter object must be of the class 'rmdFilt'")
  filter_object<- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)
  num_cols<- ncol(filter_object)
  
  filter_object_head = head(filter_object, 4)[, 1:min(num_cols, 5)]
  filter_object_tail = tail(filter_object, 4)[, 1:min(num_cols, 5)]
  blank_row = rep("---", ncol(filter_object))
  
  result<- rbind(filter_object_head, blank_row, filter_object_tail)
 
  if(num_cols > 5) message("only first 5 columns are shown")
  cat("rmdFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

#' print.cvFilt
#' 
#' For printing an S3 object of type 'cvFilt':
#' 
#'@rdname print-cvFilt
#'@export
#'
print.cvFilt<- function(filter_object){
  if(class(filter_object)[1] != "cvFilt") stop("filter object must be of the class 'cvFilt'")
  filter_object<- as.data.frame(lapply(filter_object, as.character), stringsAsFactors = FALSE)
  
  filter_object_head = head(filter_object, 4)[, 1:ncol(filter_object)]
  filter_object_tail = tail(filter_object, 4)[, 1:ncol(filter_object)]
  blank_row = rep("---", ncol(filter_object))
  
  result<- rbind(filter_object_head, blank_row, filter_object_tail)
  
  cat("cvFilt object\n")
  cat(capture.output(result), sep = "\n")
  cat("\n")
}

