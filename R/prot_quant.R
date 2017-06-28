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
    
   res<- vector("list", length(proteoformRes))
    
    r<-foreach(i=1:length(proteoformRes),.export = c("rrollup_r")) %dopar%{
      
      rollup_res<- list()
      column<- c()
      pformRes_sub<- proteoformRes[[i]]
      pformRes_sub<- as.data.frame(pformRes_sub)
      
      prot_name<- pformRes_sub$Protein[1]
      
      max_proteo_id<- max(pformRes_sub$proteoformID)
      
      for(j in 1:max_proteo_id){
        cur_subset<- pformRes_sub[which(pformRes_sub$proteoformID == j), ]
        
        inds<- which(pepData$e_data[[edata_cname]] %in% cur_subset[[edata_cname]])
        
        edata_subset<- pepData$e_data[inds, ]
        edata_subset<- edata_subset[, -edata_cname_id]
        column[j]<- paste(prot_name, j, sep=';')
        
        result<- rrollup_r(edata_subset)
        rollup_res[[j]]<- result
        
      }
      
      rollup_res<-do.call(rbind, rollup_res)
      rollup_res<- rollup_res[, -1]
      rollup_res<- cbind(column, rollup_res)
      names(rollup_res)[1]<- emeta_cname
      
      
        res[[i]]<- rollup_res
      
    }
    
    stopCluster(cl)
    
    final_res<- do.call(rbind, r)
    
    #proData_object<- as.proData(e_data = final_res, f_data= f_data, e_meta= NULL, edata_cname = emeta_cname, fdata_cname = fdata_cname)
    #return(proData_object)
    return(final_res)
  }  
  
  
  
}