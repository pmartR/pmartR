#' prot_quant wrapper function
#' 
#' Takes in a pepData object, the argument 'method' determines the quantitation method to be applied to pepData.
#' The argument isoformRes defaults to NULL, is a list of data.frames which are results of applying bpquant to
#' the original pepData object for each protein. If isoformRes is provided then,   
#' 
#' @export

prot_quant<- function(pepData, method, isoformRes = NULL){
  
  #some checks
  if(class(pepData) != "pepData") stop("pepData must be an object of class pepData")
  if(!(method %in% c('mean', 'median', 'rrollup'))) stop("method must be one of, mean, median, rrollup")
  
  edata_cname<- attr(pepData, "cnames")$edata_cname
  fdata_cname<- attr(pepData, "cnames")$fdata_cname
  emeta_cname<- attr(pepData, "cnames")$emeta_cname
  
  f_data<- pepData$f_data
  e_meta<- pepData$e_meta
  
  edata_cname_id<- which(names(pepData$e_data) == edata_cname)

  if(is.null(isoformRes)){
    
    if(method == 'mean'){
      results<- pquant_mean(pepData)
    }
    if(method == 'median'){
     results<- pquant_median(pepData)
    }
    if(method == 'rrollup'){
      results<- rrollup(pepData)
    }
    
    
  }  
  
 if(!is.null(isoformRes)){
   
   #we will extract 'isoformRes_subset' attribute from isoformRes, which is all the proteins that mapped to a nonzero proteoformID
   isoformRes2<- attr(isoformRes, "isoformRes_subset")
   
   peptides<- which(pepData$e_data[,edata_cname] %in% isoformRes2[, edata_cname])
   temp_pepdata<- as.pepData(e_data = pepData$e_data[peptides,], f_data = f_data, e_meta = isoformRes2, edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = "Protein_Isoform" )
   
   if(method == 'mean'){
     results<- pquant_mean(temp_pepdata)
   }
   
   if(method == 'median'){
         results<- pquant_median(temp_pepdata)
   }
   
   if(method == 'rrollup'){
      results<- rrollup(temp_pepdata)
   }
  
 }
  
  return(results)
}

