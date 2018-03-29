#' Normalize an object of class isobaricpepData
#' 
#' Description of what this function does
#' 
#' @param omicsData an object of the class 'isobaricpepData'
#' @param apply_norm logical, indicates whether normalization should be applied to omicsData$e_data
#'
#' @examples  
#' dontrun{ 
#' 
#' 
#' }
#' 

normalize_isobaric<- function(omicsData, apply_norm = FALSE){
  # initial checks #
  
  #check that omicsData is of correct class
  if(!inherits(omicsData, "isobaricpepData")) stop("omicsData must be of the class 'isobaricpepData'")
  
  #check that omicsData$e_data is log transformed
  if(!(attr(omicsData, "data_info")$data_scale %in% c('log2', 'log10', 'log'))) stop("omicsData$e_data must be log transformed")

  #check that apply_norm is of class logical
  if(!is.logical(apply_norm)) stop("apply_norm must be of class 'logical'")
  
  #pull some attributes from omicsData
  edata = omicsData$e_data
  fdata = omicsData$f_data
  
  edata_cname = attr(omicsData, "cnames")$edata_cname
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  exp_cname = attr(omicsData, "isobaric_info")$exp_cname
  channel_cname = attr(omicsData, "isobaric_info")$channel_cname
  refpool_channel = attr(omicsData, "isobaric_info")$refpool_channel
  
  refpool_cname = attr(omicsData, "isobaric_info")$refpool_cname
  refpool_notation = attr(omicsData, "isobaric_info")$refpool_notation
  
  #rearranging edata 
  edata_melt = melt(edata, id.vars = edata_cname)
  names(edata_melt)[2] = fdata_cname
  edata_melt2 = merge.data.frame(edata_melt, fdata[,c(exp_cname, fdata_cname)], by = fdata_cname)
  split_data = split(edata_melt2, edata_melt2[[exp_cname]])
  
  #case where refpool_channel is non-NULL and channel_cname is non-NULL
  if(!is.null(refpool_channel) & !is.null(channel_cname)){
    #looking for ref_samples
    ref_samples = fdata[[fdata_cname]][which(fdata[[channel_cname]] == refpool_channel)] 
  }
  #case where refpool_cname is non-NULL and refpool_notation is non-NULL
  else if(!is.null(refpool_cname) & !is.null(refpool_notation)){
    #looking for ref_samples
    ref_samples = fdata[[fdata_cname]][which(fdata[[refpool_cname]] == refpool_notation)] 
  }
  
  #case where apply_norm is TRUE
  if(apply_norm == TRUE){
    #applying myfunc to split_data
    res = lapply(split_data, myfunc, ref_samples = ref_samples, edata_cname = edata_cname, fdata_cname = fdata_cname)
    
    #combining normalized split_data
    final = join_all(res, by =  edata_cname)
    
    #replace omicsData$e_data with normalized edata
    omicsData$e_data<- final
    
    #updata isobaric_norm attr
    attr(omicsData, "data_info")$isobaric_norm = TRUE
    
    result = omicsData
  }else{
    #subsetting edata_melt to just contain ref_samples
    edata_melt = edata_melt[which(edata_melt[[fdata_cname]] %in% ref_samples), ]
    
    #subsetting fdata to just contain ref_samples
    fdata = fdata[which(fdata[[fdata_cname]] %in% ref_samples), ]
    
    #subset specific columns of fdata
    if(!is.null(channel_cname)){
      fdata_cname_ind = which(names(fdata) %in% fdata_cname)
      exp_cname_ind = which(names(fdata) %in% exp_cname)
      channel_cname_ind = which(names(fdata) %in% channel_cname)
      
      fdata = fdata[, c(fdata_cname_ind, exp_cname_ind, channel_cname_ind)]
      
      #merging edata_melt and fdata
      result = merge(edata_melt, fdata, by = fdata_cname)
      class(result) = "isobaricnormRes"
      
      #adding some more attributes to result
      attr(result, "isobaric_info")$exp_cname = exp_cname
      attr(result, "isobaric_info")$refpool_channel = refpool_channel
      attr(result, "isobaric_info")$channel_cname = channel_cname
      attr(result, "isobaric_info")$refpool_cname = NULL
      attr(result, "isobaric_info")$refpool_notation = NULL
    }
    else if(!is.null(refpool_cname)){
      fdata_cname_ind = which(names(fdata) %in% fdata_cname)
      exp_cname_ind = which(names(fdata) %in% exp_cname)
      refpool_cname_ind = which(names(fdata) %in% refpool_cname)
      
      fdata = fdata[, c(fdata_cname_ind, exp_cname_ind, refpool_cname_ind)]
      
      #merging edata_melt and fdata
      result = merge(edata_melt, fdata, by = fdata_cname)
      class(result) = "isobaricnormRes"
      
      #adding some more attributes to result
      attr(result, "isobaric_info")$exp_cname = exp_cname
      attr(result, "isobaric_info")$refpool_channel = NULL
      attr(result, "isobaric_info")$channel_cname = NULL
      attr(result, "isobaric_info")$refpool_cname = refpool_cname
      attr(result, "isobaric_info")$refpool_notation = refpool_notation
    }

  }
  
  return(result)
}




#lapply function to work on list items in split_data
myfunc<- function(item, ref_samples, edata_cname, fdata_cname){
  #for formula argument in dcast function
  formula = paste(edata_cname, "~", fdata_cname, sep = "")
  
  #casting data associated with one exp_cname
  cast_data = dcast(item, formula)
  
  #colum number for reference sample
  ref_samp_ind = which(names(cast_data) %in% ref_samples)
  #column number for edata_cname
  edata_cname_ind = which(names(cast_data) %in% edata_cname)
  
  #subtracting refpool sample from non-refpool samples in cast_data
  cast_data[, -edata_cname_ind] = cast_data[, -edata_cname_ind] - cast_data[, ref_samp_ind]
  #removing refpool sample from cast_data
  cast_data = cast_data[, -which(names(cast_data) %in% ref_samples)]
  return(cast_data)
}



