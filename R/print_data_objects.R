#' print.pepData
#' 
#' For printing an S3 object of type 'pepData':
#' 
#'@rdname print-pepData
#'@export
#'
print.pepData<- function(pepData){
  if(class(pepData) != "pepData") stop("pep_object must be of the class 'pepData'")
  
  e_data<- as.data.frame(lapply(pepData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(pepData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(pepData$e_meta)){
    e_meta<- as.data.frame(lapply(pepData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    emeta<-rbind(e_meta_head, blank_row, e_meta_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.metabData
#' 
#' For printing an S3 object of type 'metabData':
#' 
#'@rdname print-metabData
#'@export
#'
print.metabData<- function(metabData){
  if(class(metabData) != "metabData") stop("metab_object must be of the class 'metabData'")
  
  e_data<- as.data.frame(lapply(metabData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(metabData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(metabData$e_meta)){
    e_meta<- as.data.frame(lapply(metabData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    emeta<-rbind(e_meta_head, blank_row, e_meta_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.proData
#' 
#' For printing an S3 object of type 'proData':
#' 
#'@rdname print-proData
#'@export
#'
print.proData<- function(proData){
  if(class(proData) != "proData") stop("pro_object must be of the class 'proData'")
  
  e_data<- as.data.frame(lapply(proData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(proData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(proData$e_meta)){
    e_meta<- as.data.frame(lapply(proData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    emeta<-rbind(e_meta_head, blank_row, e_meta_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.lipidData
#' 
#' For printing an S3 object of type 'lipidData':
#' 
#'@rdname print-lipidData
#'@export
#'
print.lipidData<- function(lipidData){
  if(class(lipidData) != "lipidData") stop("lipid_object must be of the class 'lipidData'")
  
  e_data<- as.data.frame(lapply(lipidData$e_data, as.character), stringsAsFactors = FALSE)
  f_data<- as.data.frame(lapply(lipidData$f_data, as.character), stringsAsFactors = FALSE)
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(lipidData$e_meta)){
    e_meta<- as.data.frame(lapply(lipidData$e_meta, as.character), stringsAsFactors = FALSE)
    emeta_ncols<- ncol(e_meta)
    
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    emeta<-rbind(e_meta_head, blank_row, e_meta_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
    
    if(emeta_ncols > 5) message("only first 5 columns are shown")
    cat("e_meta\n")
    cat(capture.output(emeta), sep = "\n")
    cat("\n")
    
  }else{
    e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)]
    e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)]
    f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)]
    f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)]
    blank_row = rep("----", 5)
    
    edata<-rbind(e_data_head, blank_row, e_data_tail)
    fdata<-rbind(f_data_head, blank_row, f_data_tail)
    
    if(edata_ncols > 5) message("only first 5 columns are shown")
    cat("e_data\n")
    cat(capture.output(edata), sep = "\n")
    cat("\n")
    
    if(fdata_ncols > 5) message("only first 5 columns are shown")
    cat("f_data\n")
    cat(capture.output(fdata), sep = "\n")
    cat("\n")
  }
}

#' print.normRes
#' 
#' For printing an S3 object of type 'normRes':
#' 
#'@rdname print-normRes
#'@export
#'
print.normRes<- function(normRes){
  if(class(normRes) != "normRes") stop("normRes object must be of the class 'normRes'") 
  
  attr(normRes, "class")<- NULL
  attr(normRes, "omicsData")<- NULL
  
  print(normRes)
}