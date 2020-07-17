#' Normalize an object of class isobaricpepData
#' 
#' The samples are normalized to their corresponding reference pool sample
#' 
#' @param omicsData an object of the class 'isobaricpepData'
#' @param apply_norm logical, indicates whether normalization should be applied to omicsData$e_data
#' #' @param exp_cname character string specifying the name of the column containing the experiment/plate information in \code{f_data}.
#' @param channel_cname optional character string specifying the name of the column containing the instrument channel a sample was run on in \code{f_data}. This argument is optional, see Details for how to specify information regarding reference pool samples. If using this argument, the 'refpool_channel' argument must also be specified; in this case, 'refpool_cname' and 'refpool_notation' should not be specified.
#' @param refpool_channel optional character string specifying which channel contained the reference pool sample, only used when this remains the same from experiment to experiment. This argument is optional, see Details for how to specify information regarding reference pool samples. If using this argument, the 'channel_cname' argument must also be specified; in this case, 'refpool_cname' and 'refpool_notation' should not be specified.
#' @param refpool_cname optional character string specifying the name of the column containing information about which samples are reference samples in \code{f_data}. This argument is optional, see Details for how to specify information regarding reference pool samples. If using this argument, the 'refpool_notation' argument must also be specified; in this case, 'channel_cname' and 'refpool_channel' should not be specified.
#' @param refpool_notation optional character string specifying the value in the refpool_channel column which denotes that a sample is a reference sample. This argument is optional, see Details for how to specify information regarding reference pool samples. If using this argument, the 'refpool_cname' argument must also be specified; in this case, 'channel_cname' and 'refpool_channel' should not be specified.
#' @details 
#' #' There are two ways to specify the information needed for identifying reference samples which should be used for normalization:
#' \enumerate{
#' \item specify \code{channel_cname} and \code{refpool_channel}. This should be used when the reference sample for each experiment/plate was always located in the same channel. Here \code{channel_cname} gives the column name for the column in \code{f_data} which gives information about which channel each sample was run on, and \code{refpool_channel} is a character string specifying the value in \code{channel_colname} that corresponds to the reference sample channel.
#' \item specify \code{refpool_cname} and \code{refpool_notation}. This should be used when the reference sample is not in a consistent channel across experiments/plates. Here, \code{refpool_cname} gives the name of the column in \code{f_data} which indicates whether a sample is a reference or not, and \code{refpool_notation} is a character string giving the value used to denote a reference sample in that column.
#' }
#' In both cases you must specify \code{exp_cname} which gives the column name for the column in \code{f_data} containing information about which experiment/plate a sample was run on.
#' 
#' See examples below.
#' @examples  
#' dontrun{ 
#' library(pmartRdata)
#' data(isobaric_object)
#' 
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set", apply_norm = TRUE, channel_cname = "iTRAQ.Channel", refpool_channel = "116")
#' 
#' # alternate specification: #
#' data(isobaric_object)
#' 
#' isobaric_object = edata_transform(isobaric_object, "log2")
#' isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set", apply_norm = TRUE, refpool_cname = "Reference", refpool_notation = "Yes")
#' 
#' }
#' 
#' @export
#'
  
normalize_isobaric<- function(omicsData, exp_cname = NULL, apply_norm = FALSE, channel_cname = NULL, refpool_channel = NULL, refpool_cname = NULL, refpool_notation = NULL){
  # initial checks #
  
  #check that omicsData is of correct class
  if(!inherits(omicsData, "isobaricpepData")) stop("omicsData must be of the class 'isobaricpepData'")
  
  #check that omicsData$e_data is log transformed
  if(!(attr(omicsData, "data_info")$data_scale %in% c('log2', 'log10', 'log'))) stop("omicsData$e_data must be log transformed")

  # check that exp_cname is in f_data #
  if(!(exp_cname %in% names(omicsData$f_data))) stop(paste("Experiment column", exp_cname, "is not found in f_data.", sep = " "))
  
  #check that apply_norm is of class logical
  if(!is.logical(apply_norm)) stop("apply_norm must be of class 'logical'")
  
  # check that channel_cname is in f_data, if not NULL #
  if(!is.null(channel_cname)){
    if(!(channel_cname %in% names(omicsData$f_data))) stop(paste("Channel column", channel_cname, "is not found in f_data. See details of as.isobaricpepData for specifying column names.", sep = " "))
  }
  
  # check that refpool_cname is in f_data, if not NULL #
  if(!is.null(refpool_cname)){
    if(!(refpool_cname %in% names(omicsData$f_data))) stop(paste("Reference pool column", refpool_cname, "is not found in f_data. See details of as.isobaricpepData for specifying column names.", sep = " "))
  }
  # make sure the reference pool info info is specified appropriately #
  # possibility 1: specify refpool_cname #
  poss1 = !is.null(refpool_cname) & !is.null(refpool_notation)
  # possibility 2: specify refpool_channel and channel_cname#
  poss2 = !is.null(refpool_channel) & !is.null(channel_cname) 
  
  # throw an error if neither or both of these are true #
  if((poss1 + poss2) != 1) stop("Reference samples information was not correctly specified. See Details and Examples for more information.")
  
  # get some values #

  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  
  # if possibility 1 is used, check that ref_cname column 
  if(poss1 == TRUE){
    if(!is.character(refpool_notation)) stop("refpool_notation must be of class 'character'")
    omicsData$f_data[ , refpool_cname] = as.character(omicsData$f_data[ , refpool_cname])
    omicsData$f_data[ , fdata_cname] = as.character(omicsData$f_data[ , fdata_cname])
    idx = split(omicsData$f_data[ , refpool_cname], omicsData$f_data[ , exp_cname])
    temp_check = lapply(idx, function(x) refpool_notation %in% x)
    if(sum(unlist(temp_check)) != length(temp_check)) stop(paste("'refpool_notation=", refpool_notation, " is not in every experiment. See Details and Examples for more information."))
  }
  
  # if possibility 2 is used, check that refpool_channel is a value seen in each experiment #
  if(poss2 == TRUE){
    if(!is.character(refpool_channel)) stop("refpool_channel must be of class 'character'")
    omicsData$f_data[ , channel_cname] = as.character(omicsData$f_data[ , channel_cname])
    omicsData$f_data[ , fdata_cname] = as.character(omicsData$f_data[ , fdata_cname])
    idx = split(omicsData$f_data[ , channel_cname], omicsData$f_data[ , exp_cname])
    temp_check = lapply(idx, function(x) refpool_channel %in% x)
    if(sum(unlist(temp_check)) != length(temp_check)) stop(paste("'refpool_channel=", refpool_channel, " is not in every experiment. See Details and Examples for more information."))
  }
  
  # set up attributes for reference pool information #
  attr(omicsData, "isobaric_info") = list(exp_cname = exp_cname, channel_cname = channel_cname, refpool_channel = refpool_channel, refpool_cname = refpool_cname, refpool_notation = refpool_notation, norm_info = list(is_normalized = FALSE))
  
  
  #pull some attributes from omicsData
  edata = omicsData$e_data
  fdata = omicsData$f_data
  
  edata_cname = attr(omicsData, "cnames")$edata_cname
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  
  
  channel_cname = attr(omicsData, "isobaric_info")$channel_cname
  refpool_channel = attr(omicsData, "isobaric_info")$refpool_channel
  
  refpool_cname = attr(omicsData, "isobaric_info")$refpool_cname
  refpool_notation = attr(omicsData, "isobaric_info")$refpool_notation
  
  #rearranging edata 
  edata_melt = reshape2::melt(edata, id.vars = edata_cname)
  names(edata_melt)[2] = fdata_cname
  edata_melt2 = merge.data.frame(edata_melt, fdata[,c(exp_cname, fdata_cname)], by = fdata_cname)
  split_data = split(edata_melt2, edata_melt2[[exp_cname]])
  
  #case where refpool_channel is non-NULL and channel_cname is non-NULL
  if(!is.null(refpool_channel) & !is.null(channel_cname)){
    #looking for ref_samples
    ref_samples = fdata[[fdata_cname]][which(fdata[[channel_cname]] == refpool_channel)] 
  }else{
    #case where refpool_cname is non-NULL and refpool_notation is non-NULL
    if(!is.null(refpool_cname) & !is.null(refpool_notation)){
      #looking for ref_samples
      ref_samples = fdata[[fdata_cname]][which(fdata[[refpool_cname]] == refpool_notation)] 
    }
  }
  
  
  #case where apply_norm is TRUE
  if(apply_norm == TRUE){
    #applying myfunc to split_data
    res = lapply(split_data, myfunc, ref_samples = ref_samples, edata_cname = edata_cname, fdata_cname = fdata_cname)
    
    #combining normalized split_data
    final = Reduce(function(...) merge(..., all=T), res)
    
    #replace omicsData$e_data with normalized edata
    omicsData$e_data<- final
    
    #remove refsamples from f_data 
    omicsData$f_data<- fdata[-which(fdata[[fdata_cname]] %in% ref_samples),]
    
    #update isobaric norm flag
    attr(omicsData, "isobaric_info")$norm_info$is_normalized = TRUE
    attr(omicsData, "isobaric_info")$exp_cname = exp_cname
    
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
      attr(result, "cnames")$edata_cname = edata_cname
      attr(result, "cnames")$fdata_cname = fdata_cname
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
      attr(result, "cnames")$edata_cname = edata_cname
      attr(result, "cnames")$fdata_cname = fdata_cname
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
  cast_data = reshape2::dcast(item, formula)
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



