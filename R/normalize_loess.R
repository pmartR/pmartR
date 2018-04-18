# start on normalize_loess function

normalize_loess<- function(omicsData, method = "fast", span = .4){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # data should be on log scale #
  if(!(attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log"))) stop("omicsData must be on log scale")

  #check that norm attribute is False
  if(attr(omicsData, "data_info")$data_norm == TRUE) stop("data has already been normalized")
  
  #extract cname
  edata_cname = attr(omicsData, "cnames")$edata_cname
  
  #remove edata_cname col from e_data
  e_data = omicsData$e_data
  e_data = e_data[, -which(names(e_data) == edata_cname)]
  
  #apply normalizeCyclicLoess function to e_data
  result = limma::normalizeCyclicLoess(e_data, span = span, method = method)
  
  #update e_data with normalized result
  omicsData$e_data[,-which(names(omicsData$e_data) == edata_cname)] = result
  
  #update norm attribute in omicsData
  attr(omicsData, "data_info")$data_norm = TRUE
  
  return(omicsData)
}