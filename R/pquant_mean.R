#' Protein Quantitation using Mean Peptide Abundances
#' 
#' Description here
#' 
#' @export

#pquant_mean for pmartRqc pep_object
pquant_mean <- function(pepData){
  # check that pepData is of appropraite class #
  if(class(pepData) != "pepData") stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  
  pep_id = attr(pepData, "cnames")$edata_cname
  pro_id = attr(pepData, "cnames")$emeta_cname
  
  
  pep = data.table(pepData$e_data)
  pro = data.table(pepData$e_meta[,c(pep_id, pro_id)])
  temp = data.table:::merge.data.table(x = pro, y = pep, by = pep_id, all.x = F, all.y = T)
  temp = as.data.frame(temp, check.names=FALSE)[,-which(names(temp)==pep_id)]
  DT = data.table(temp)
  res = as.data.frame(DT[,lapply(.SD, mean, na.rm = T), by = pro_id], check.names=FALSE)
  
  samp_id = attr(pepData, "cnames")$fdata_cname
  data_scale = attr(pepData, "data_info")$data_scale
  data_norm = attr(pepData, "data_info")$data_norm
  
  if(ncol(pepData$e_meta) == 2){
    e_meta = as.data.frame(pepData$e_meta[,-which(names(pepData$e_meta)==pep_id)])
    names(e_meta)<-pro_id
  }
  
  else {e_meta = pepData$e_meta[,-which(names(pepData$e_meta)==pep_id)]} 
  
  
  prodata = as.proData(e_data = data.frame(res, check.names=FALSE), f_data = pepData$f_data, e_meta = NULL , fdata_cname = samp_id, edata_cname = pro_id, data_scale = data_scale, data_norm = data_norm)
  
  attr(prodata, "meta_info") = attr(pepData, "meta_info")
  
  
  return(prodata)
}