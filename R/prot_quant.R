#Wrapper function prot-quant

prot_quant<- function(pepData, method, proteoformRes = NULL){
  
  #some checks
  if(class(pepData) != "pepData") stop("pepData must be an object of class pepData")
  if(!(method %in% c('mean', 'median', 'rrollup'))) stop("method must be one of, mean, median, rrollup")
  
  edata_cname<- attr(pepData, "cnames")$edata_cname
  fdata_cname<- attr(pepData, "cnames")$fdata_cname
  emeta_cname<- attr(pepData, "cnames")$emeta_cname
  
  f_data<- pepData$f_data
  e_meta<- pepData$e_meta
  
  edata_cname_id<- which(names(pepData$e_data) == edata_cname)

  if(is.null(proteoformRes)){
    
    if(method == 'mean'){
      results<- pquant_mean(pepData)
    }
    if(method == 'median'){
      #use pquant_median function
    }
    if(method == 'rrollup'){
      #use rrollup_new function
    }
    
    
  }  
  
 if(!is.null(proteoformRes)){
   
   if(method == 'rrollup'){
     #apply "myfunction" to proteoformRes to identify "Protein_Isoform"
     proteoformRes2<- lapply(proteoformRes, myfunction)
     proteoformRes2<- do.call(rbind, proteoformRes2)
    
     peptides<- which(pepData$e_data[,edata_cname] %in% proteoformRes2[, edata_cname])
     
     temp_pepdata<- as.pepData(e_data = pepData$e_data[peptides,], f_data = f_data, e_meta = proteoformRes2, edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = "Protein_Isoform" )
     
     results<- pmartRqc:::rrollup(temp_pepdata)
     
   }
    
 
  
 }
  return(results)
}