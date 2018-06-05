#' Protein Quantitation using Mean Peptide Abundances
#' 
#' This function takes in a pepData object and returns a proData object
#' 
#' @param pepData omicsData object of class 'pepData'
#' 
#' @return an omicsData object of class 'proData'
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = pquant_mean(pep_object) 
#'}
#' 
#' @rdname pquant_mean
#' @export


pquant_mean <- function(pepData){
  
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  
  pep_id = attr(pepData, "cnames")$edata_cname
  pro_id = attr(pepData, "cnames")$emeta_cname
  
  
  pep = data.table(pepData$e_data)
  pro = data.table(pepData$e_meta[,c(pep_id, pro_id)])
  temp = data.table:::merge.data.table(x = pro, y = pep, by = pep_id, all.x = F, all.y = T)
  temp = as.data.frame(temp, check.names=check_names)[,-which(names(temp)==pep_id)]
  DT = data.table(temp)
  res = as.data.frame(DT[,lapply(.SD, function(x){if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}}), by = pro_id], check.names=check_names)
  
  
  samp_id = attr(pepData, "cnames")$fdata_cname
  data_scale = attr(pepData, "data_info")$data_scale
  data_norm = attr(pepData, "data_info")$data_norm
  
  if(ncol(pepData$e_meta) == 2){
    e_meta = as.data.frame(pepData$e_meta[,-which(names(pepData$e_meta)==pep_id)])
    names(e_meta)<-pro_id
  }
  
  else {e_meta = pepData$e_meta[,-which(names(pepData$e_meta)==pep_id)]} 
  
  
  prodata = as.proData(e_data = data.frame(res, check.names=check_names), f_data = pepData$f_data,  e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, data_norm = data_norm, check.names = check_names)
  
  #check for isobaricpepData class
  if(inherits(pepData, "isobaricpepData")){
    #update attributes in prodata
    attr(prodata, "isobaric_info") = attr(pepData, "isobaric_info")
    attr(prodata, "data_info")$isobaric_norm = attr(pepData, "data_info")$isobaric_norm
  }
  
  #updating prodata attributes
  attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
  attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
  attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
  
  attr(prodata, "filters")<- attr(pepData, "filters")
  attr(prodata, "group_DF")<- attr(pepData, "group_DF")
  attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  
  
  return(prodata)
}
