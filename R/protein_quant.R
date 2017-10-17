#' protein_quant wrapper function
#' 
#' This function takes in a pepData object, method (quantitation method, mean, median or rrollup), and the optional argument isoformRes (defaults to NULL). An object of the class 'proData' is returned. 
#' 
#' @param pepData an omicsData object of the class 'pepData'
#' @param method is one of the three protein quantitation methods, 'mean', 'median' and 'rrollup'
#' @param isoformRes is a list of data frames, the result of applying the 'bpquant_loop' function to original pepData object. Defaults to NULL
#' 
#' @return an omicsData object of the class 'proData'
#' 
#' @details If isoformRes is provided then, a temporary pepData object is formed using the isoformRes information as the e_meta component and the original pepData object will be used for e_data and f_data components. The emeta_cname for the temporary pepData object will be the “protein_isoform” column of isoformRes. Then one of the three 'method' functions can be applied to the temporary pepData object to return a proData object. If isofromRes is left NULL, then depending on the input for 'method', the correct 'method' function is applied directly to the input pepData object and a proData object is returned.
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' 
#' case where isoformRes is NULL:
#' results<- prot_quant(pepData = pep_object, method = 'median', isoformRes = NULL)
#' 
#' case where isoformRes is provided:
#' results2 = prot_quant(pep_data = pep_object, method = 'mean', isoformRes = isoformRes_object)
#' }
#'
#' @rdname protein_quant 
#' @export

protein_quant<- function(pepData, method, isoformRes = NULL){
  
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
