# prot-quant function

prot_quant<- function(pepData, method, proteoformRes, parallel = TRUE){
  
  #some checks
  if(class(pepData) != "pepData") stop("pepData must be an object of class pepData")
  if(!(method %in% c('mean', 'median', 'rrollup'))) stop("method must be one of, mean, median, rrollup")
  
  edata_cname<- attr(pepData, "cnames")$edata_cname
  fdata_cname<- attr(pepData, "cnames")$fdata_cname
  emeta_cname<- attr(pepData, "cnames")$emeta_cname
  
  f_data<- pepData$f_data
  e_meta<- pepData$e_meta
  
  edata_cname_id<- which(names(pepData$e_data) == edata_cname)
  
  if(method == 'rrollup'){

    library(doParallel)
    cl<- makeCluster(4)
    registerDoParallel(cl)
    
    final_res<- foreach(i=1:length(proteoformRes), .combine = 'rbind', .packages="foreach") %dopar%{
        pformRes_sub<- proteoformRes[[i]]
          
        foreach(j = 1:max(pformRes_sub$proteoformID), .combine = 'rbind') %do%{
           prot_name<- pformRes_sub$Protein[1]
           cur_subset<- pformRes_sub[which(pformRes_sub$proteoformID == j), ]
           inds<- which(pepData$e_data[[edata_cname]] %in% cur_subset[[edata_cname]])
           edata_subset<- pepData$e_data[inds, -edata_cname_id]
           prot<- paste(prot_name, j, sep=';')
           result<- rrollup_r(edata_subset)
           result<- result[,-1]
           result<- cbind(prot, result)
      }
    }
    stopCluster(cl)
    
    names(final_res)[1]<- emeta_cname
    proteoformRes_sub<-attr(proteoformRes, "proteoformRes_subset")
    
    #proData_object<- as.proData(e_data = final_res, f_data= f_data, e_meta= NULL, edata_cname = emeta_cname, fdata_cname = fdata_cname)
    #return(proData_object)
    return(final_res)
  }  
  
  
  
}