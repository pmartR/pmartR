#work on normalize_isobaric function 

normalize_isobaric<- function(omicsData, apply_norm = FALSE){
  # initial checks #
  
  #check that omicsData is of correct class
  if(!inherits(omicsData, "isobaricpepData")) stop("omicsData must be of the class 'isobaricpepData'")
  
  #check that omicsData$e_data is log transformed
  #if(!(attr(omicsData, "data_info")$data_scale %in% c('log2', 'log10', 'log'))) stop("omicsData$e_data must be log transformed")
  
  #pull some attributes from omicsData
  edata = omicsData$e_data
  fdata = omicsData$f_data
  
  edata_cname = attr(omicsData, "cnames")$edata_cname
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  exp_cname = attr(omicsData, "isobaric_info")$exp_cname
  channel_cname = attr(omicsData, "isobaric_info")$channel_cname
  refpool_channel = attr(omicsData, "isobaric_info")$refpool_channel
  #looking at isobaric_info attr
  
  #case where refpool_channel is non-NULL and channel_cname is non-NULL
  if(!is.null(attr(omicsData, "isobaric_info")$refpool_channel) & !is.null(attr(omicsData, "isobaric_info")$channel_cname)){
    
    
  }
}



#rearranging data 
edata_melt = melt(edata, id.vars = edata_cname)
names(edata_melt)[2] = fdata_cname
edata_melt = merge.data.frame(edata_melt, fdata[,c(exp_cname, fdata_cname)], by = fdata_cname)


