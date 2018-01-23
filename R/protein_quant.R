#' protein_quant wrapper function
#' 
#' This function takes in a pepData object, method (quantitation method, mean, median or rrollup), and the optional argument isoformRes (defaults to NULL). An object of the class 'proData' is returned. 
#' 
#' @param pepData an omicsData object of the class 'pepData'.
#' @param method is one of four protein quantitation methods, 'rollup', 'rrollup', 'qrollup' and 'zrollup'. When 'rollup' is selected, combine_fn must also be provided and will determine whether pquant_mean or pquant_median function will be used.
#' @param isoformRes is a list of data frames, the result of applying the 'bpquant_loop' function to original pepData object. Defaults to NULL.
#' @param qrollup_thresh is a numeric value; is the peptide abundance cutoff value. Is an argument to qrollup function.
#' @param single_pep logical indicating whether or not to remove proteins that have just a single peptide mapping to them, defaults to FALSE.
#' @param single_observation logical indicating whether or not to remove peptides that have just a single observation, defaults to FALSE.
#' @param combine_fn can either be 'mean' or 'median'.
#' @param use_parallel logical indicating whether or not to use "doParallel" loop in applying rollup functions. Defaults to TRUE. Is an argument of rrollup, qrollup and zrollup functions.
#' 
#' @return an omicsData object of the class 'proData'
#' 
#' @details If isoformRes is provided then, a temporary pepData object is formed using the isoformRes information as the e_meta component and the original pepData object will be used for e_data and f_data components. The emeta_cname for the temporary pepData object will be the âprotein_isoformâ column of isoformRes. Then one of the three 'method' functions can be applied to the temporary pepData object to return a proData object. If isofromRes is left NULL, then depending on the input for 'method', the correct 'method' function is applied directly to the input pepData object and a proData object is returned.
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' 
#' case where isoformRes is NULL:
#' results<- prot_quant(pepData = pep_object, method = 'rollup', combine_fn = 'median', isoformRes = NULL)
#' 
#' case where isoformRes is provided:
#' results2 = prot_quant(pep_data = pep_object, method = 'rollup', combine_fn = 'mean', isoformRes = isoformRes_object)
#' }
#'
#' @rdname protein_quant 
#' @export

protein_quant<- function(pepData, method, isoformRes = NULL, qrollup_thresh = NULL, single_pep = FALSE, single_observation = FALSE, combine_fn = "median", use_parallel = TRUE){
  
  #some checks
  if(!inherits(pepData, "pepData")) stop("pepData must be an object of class pepData")
  if(!(method %in% c('rollup', 'rrollup', 'qrollup', 'zrollup'))) stop("method must be one of, rollup, rrollup, qrollup, zrollup")
  if(!(combine_fn %in% c('median', 'mean'))) stop("combine_fn must be on of mean or median")
  
  edata_cname<- attr(pepData, "cnames")$edata_cname
  fdata_cname<- attr(pepData, "cnames")$fdata_cname
  emeta_cname<- attr(pepData, "cnames")$emeta_cname
  
  f_data<- pepData$f_data
  e_meta<- pepData$e_meta
  
  edata_cname_id<- which(names(pepData$e_data) == edata_cname)
  
  #gives message if single_pep and single_observation are TRUE and method is not zrollup
  if(method != 'zrollup' & (single_pep == TRUE | single_observation == TRUE)) message("single_pep and single_observation will be ignored, as they are only applicable if method is zrollup")
  
  #gives message if qrollup_thresh is not NULL and method is not qrollup
  if(method != 'qrollup' & !is.null(qrollup_thresh)) message("qrollup_thresh argument will be ignored, as it is only applicable if method is qrollup")
  
  #gives message if use_parallel is TRUE and method is not a rollup
  #if(!(method %in% c('rrollup','qrollup', 'zrollup')) & use_parallel == TRUE) message("use_parallel argument will be ignored, as it is not applicable if method isn't rrollup, qrollup or zrollup")
  
  if(is.null(isoformRes)){
    if(method == 'rollup'){
      if(combine_fn == 'median'){
        results<- pquant_median(pepData)
      }else{results = pmartRqc:::pquant_mean(pepData)}
    }
    if(method == 'rrollup'){
      results<- rrollup(pepData, combine_fn = combine_fn, parallel = use_parallel)
    }
    if(method == 'qrollup'){
      results<- qrollup(pepData, qrollup_thresh = qrollup_thresh, combine_fn = combine_fn, parallel = use_parallel)
    }
    if(method == 'zrollup'){
      
      #check for single peps and single observations
      if(single_pep == FALSE & single_observation == FALSE){
        proteomicsfilt = proteomics_filter(pepData)
        moleculefilt = molecule_filter(pepData)
        
        if(any(moleculefilt$Num_Observations == 1) & any(proteomicsfilt$counts_by_pro == 1)) stop("There are peptides with less than 2 observations and proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides/proteins and then run zrollup, set both 'single_observation' and 'single_pep' arguments to TRUE")
        if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
        if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE")
      }
      
      if(single_pep == TRUE & single_observation == FALSE){
        proteomicsfilt = proteomics_filter(pepData)
        moleculefilt = molecule_filter(pepData)
        
        #since single_pep is TRUE we will remove proteins with single peptide mapped to them
        pepData = applyFilt(proteomicsfilt, pepData, min_num_peps = 2)
        
        if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
      }
      
      if(single_pep == FALSE & single_observation == TRUE){
        proteomicsfilt = proteomics_filter(pepData)
        moleculefilt = molecule_filter(pepData)
        
        #since single_observation is TRUE we will remove peptides with single observation
        pepData = applyFilt(moleculefilt, pepData)
        
        if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE.")
      }
      
      if(single_pep == TRUE & single_observation == TRUE){
        proteomicsfilt = proteomics_filter(pepData)
        moleculefilt = molecule_filter(pepData)
        
        pepData0 = applyFilt(moleculefilt, pepData)
        pepData = applyFilt(proteomicsfilt, pepData0, min_num_peps = 2)
      }
      
      results<- zrollup(pepData, combine_fn = combine_fn, parallel = use_parallel)
    }
  }  
  
  if(!is.null(isoformRes)){
    
    #we will extract 'isoformRes_subset' attribute from isoformRes, which is all the proteins that mapped to a nonzero proteoformID
    isoformRes2<- attr(isoformRes, "isoformRes_subset")
    
    #pulling more attributes from pepData, will be used as arguments in as.pepData
    data_scale = attr(pepData, "data_info")$data_scale
    data_norm = attr(pepData, "data_info")$data_norm
    
    peptides<- which(pepData$e_data[,edata_cname] %in% isoformRes2[, edata_cname])
    temp_pepdata<- as.pepData(e_data = pepData$e_data[peptides,], f_data = f_data, e_meta = isoformRes2, edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = "Protein_Isoform", data_scale = data_scale, data_norm = data_norm)
    
    #now we will update the attributes of temp_pepdata with attributes from original pepData
    attr(temp_pepdata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
    attr(temp_pepdata, "data_info")$data_types = attr(pepData, "data_info")$data_types
    attr(temp_pepdata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
    
    attr(temp_pepdata, "filters")<- attr(pepData, "filters")
    attr(temp_pepdata, "group_DF")<- attr(pepData, "group_DF")
    attr(temp_pepdata, "imdanova")<- attr(pepData, "imdanova")
    
    if(method == 'rollup'){
      if(combine_fn == 'median'){
        results<- pquant_median(temp_pepdata)
      }else{results<- pquant_mean(temp_pepdata)}
    }

    if(method == 'rrollup'){
      results<- rrollup(temp_pepdata, combine_fn = combine_fn, parallel = use_parallel)
    }
    
    if(method == 'qrollup'){
      results<- pmartRqc:::qrollup(temp_pepdata, qrollup_thresh = qrollup_thresh, combine_fn = combine_fn, parallel = use_parallel)
    }
    
    if(method == 'zrollup'){
      #check for single peps and single observations
      if(single_pep == FALSE & single_observation == FALSE){
        proteomicsfilt = proteomics_filter(temp_pepdata)
        moleculefilt = molecule_filter(temp_pepdata)
        
        if(any(moleculefilt$Num_Observations == 1) & any(proteomicsfilt$counts_by_pro == 1)) stop("There are peptides with less than 2 observations and proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides/proteins and then run zrollup, set both 'single_observation' and 'single_pep' arguments to TRUE.")
        if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. The method zrollup cannot be applied when this is the case. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
        if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE.")
      }
      
      if(single_pep == TRUE & single_observation == FALSE){
        proteomicsfilt = proteomics_filter(temp_pepdata)
        moleculefilt = molecule_filter(temp_pepdata)
        
        #since single_pep is TRUE we will remove proteins with single peptide mapped to them
        temp_pepdata = applyFilt(proteomicsfilt, temp_pepdata, min_num_peps = 2)
        
        if(any(moleculefilt$Num_Observations == 1)) stop("There are peptides with less than 2 observations. If you would like to filter out these peptides and then run zrollup, set the 'single_observation' argument to TRUE.")
      }
      
      if(single_pep == FALSE & single_observation == TRUE){
        proteomicsfilt = proteomics_filter(temp_pepdata)
        moleculefilt = molecule_filter(temp_pepdata)
        
        #since single_observation is TRUE we will remove peptides with single observation
        temp_pepdata = applyFilt(moleculefilt, temp_pepdata)
        
        if(any(proteomicsfilt$counts_by_pro == 1)) stop ("There are proteins with just a single peptide mapping to them. The method zrollup cannot be applied when this is the case. If you would like to filter out these single peptide to protein mappings and then run zrollup, set the 'single_pep' input argument to TRUE.")
      }
      
      if(single_pep == TRUE & single_observation == TRUE){
        proteomicsfilt = proteomics_filter(temp_pepdata)
        moleculefilt = molecule_filter(temp_pepdata)
        
        temp_pepdata0 = applyFilt(proteomicsfilt, temp_pepdata, min_num_peps = 2)
        temp_pepdata = applyFilt(moleculefilt, temp_pepdata0)
      }
      
      results<- zrollup(temp_pepdata, combine_fn = combine_fn, parallel = use_parallel)
    }
  }
  return(results)
}
