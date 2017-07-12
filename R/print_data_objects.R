#' print.pepData
#' 
#' For printing an S3 object of type 'pepData':
#' 
#'@rdname print-pepData
#'@export
#'
print.pepData<- function(pepData){
  if(class(pepData) != "pepData") stop("pep_object must be of the class 'pepData'")
  
  e_data<- pepData$e_data
  f_data<- pepData$f_data
 
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(pepData$e_meta)){
    e_meta<- pepData$e_meta
    emeta_ncols<- ncol(e_meta)
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)], e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)], e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)])
  }else{
    e_meta<-NULL
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)])
    
  }
  
  if(max(edata_ncols, fdata_ncols, emeta_ncols) > 5) message("only first 5 columns are shown")

  print(out)
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
  
  e_data<- metabData$e_data
  f_data<- metabData$f_data
  
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(metabData$e_meta)){
    e_meta<- metabData$e_meta
    emeta_ncols<- ncol(e_meta)
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)], e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)], e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)])
  }else{
    e_meta<-NULL
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)])
  }
  
  if(max(edata_ncols, fdata_ncols, emeta_ncols) > 5) message("only first 5 columns are shown")
  
  print(out)
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
  
  e_data<- proData$e_data
  f_data<- proData$f_data
  
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(proData$e_meta)){
    e_meta<- proData$e_meta
    emeta_ncols<- ncol(e_meta)
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)], e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)], e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)])
  }else{
    e_meta<-NULL
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)])
  }
  
  if(max(edata_ncols, fdata_ncols, emeta_ncols) > 5) message("only first 5 columns are shown")
  
  print(out)
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
  
  e_data<- lipidData$e_data
  f_data<- lipidData$f_data
  
  edata_ncols<- ncol(e_data)
  fdata_ncols<- ncol(f_data)
  
  if(!is.null(lipidData$e_meta)){
    e_meta<- lipidData$e_meta
    emeta_ncols<- ncol(e_meta)
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)], e_meta_head = head(e_meta, 4)[, 1:min(emeta_ncols, 5)], e_meta_tail = tail(e_meta, 4)[, 1:min(emeta_ncols, 5)])
  }else{
    e_meta<-NULL
    out<- list(e_data_head = head(e_data, 4)[, 1:min(edata_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edata_ncols, 5)], f_data_head = head(f_data, 4)[, 1:min(fdata_ncols, 5)], f_data_tail = tail(f_data, 4)[, 1:min(fdata_ncols, 5)])
  }
  
  if(max(edata_ncols, fdata_ncols, emeta_ncols) > 5) message("only first 5 columns are shown")
  
  print(out)
}