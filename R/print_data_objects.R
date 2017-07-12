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
 
  edat_ncols<- ncol(e_data)
  fdat_ncols<- ncol(f_data)
  
  if(!is.null(pepData$e_meta)){
    e_meta<- pepData$e_meta
    emeta_ncols<- ncol(e_meta)
    
  }else{
    e_meta<-NULL
  }
  
  if(max(edat_ncols, fdat_ncols, emeta_ncols) > 5) message("only first 5 columns are shown")

  out<- list(e_data = head(e_data, 4)[, 1:min(edat_ncols, 5)], e_data_tail = tail(e_data, 4)[, 1:min(edat_ncols, 5)], f_data = head(f_data, 4)[, 1:min(ncol(f_data), 5)], e_meta = head(e_meta, 4)[, 1:min(ncol(e_meta), 5)])
  
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
  message("only first 6 rows were printed")
  
  if(!is.null(metabData$e_meta)){
    e_meta<- metabData$e_meta
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = head(e_meta)[, 1:min(ncol(e_meta), 6)])
  }else{
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = NULL) 
  }
  
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
  f_data<-proData$f_data
  
  ncols_edata<- ncol(e_data)
  ncols_fdata<- ncol(f_data)
  
  
  
  if(!is.null(proData$e_meta)){
    e_meta<- proData$e_meta
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = head(e_meta)[, 1:min(ncol(e_meta), 6)])
  }else{
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = NULL) 
  }
  
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
  f_data<-lipidData$f_data
  message("only first 6 rows were printed")
  
  if(!is.null(lipidData$e_meta)){
    e_meta<-lipidData$e_meta
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = head(e_meta)[, 1:min(ncol(e_meta), 6)])
  }else{
    out<- list(e_data = head(e_data)[, 1:min(ncol(e_data), 6)], f_data = head(f_data)[, 1:min(ncol(f_data), 6)], e_meta = NULL) 
  }
  
  print(out)
}