#' prot_quant wrapper function
#' 
#' Takes in a pepData object, the argument 'method' determines the quantitation method to be applied to pepData.
#' The argument proteoformRes defaults to NULL, is a list of data.frames which are results of applying bpquant to
#' the original pepData object for each protein. If proteoformRes is provided then,   
#' 
#' @export

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
     results<- pquant_median(pepData)
    }
    if(method == 'rrollup'){
      results<- rrollup(pepData)
    }
    
    
  }  
  
 if(!is.null(proteoformRes)){
   
   #apply "proteoformRes_func" to proteoformRes to identify "Protein_Isoform"
   proteoformRes2<- lapply(proteoformRes, proteoformRes_func)
   proteoformRes2<- do.call(rbind, proteoformRes2)
   
   peptides<- which(pepData$e_data[,edata_cname] %in% proteoformRes2[, edata_cname])
   temp_pepdata<- as.pepData(e_data = pepData$e_data[peptides,], f_data = f_data, e_meta = proteoformRes2, edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = "Protein_Isoform" )
   
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




# function to use with lapply on proteoformRes
proteoformRes_func<- function(df){
  temp<- vector("list", max(df$proteoformID))
  
  for(i in 1:max(df$proteoformID)){
    cur_subset<- df[which(df$proteoformID == i), ]
    new_df<- data.frame(cur_subset$Protein, paste(cur_subset$Protein, cur_subset$proteoformID, sep = ';'), as.numeric(cur_subset$Mass_Tag_ID))
    names(new_df)<- c(names(df)[1], "Protein_Isoform", names(df)[2])
    temp[[i]]<- new_df
  }
  
  do.call(rbind, temp)
}