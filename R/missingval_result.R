#' Creates an object of class naRes
#' 
#' This function takes in an omicsData object, outputs a list of two data frames, one contains NA by sample, the second contains NA by molecule
#'
#' @param omicsData an object of class "pepData", "proData", "metabData", "lipidData", or "nmrData", created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' 
#' @return outputs a list of two data frames, one contains NA by sample (missing values per sample), the second contains NA by molecule (missing values per molecule). The output is assigned class 'naRes'. 
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data("lipid_object")
#' data("metab_object")
#' 
#' result = missingval_result(lipid_object)
#' result2 = missingval_result(metab_object)
#'
#'}  
#' 
#' @rdname missingval_result
#' @export
#' 

missingval_result<- function(omicsData){
  
  #check for correct class
  if(!inherits(omicsData, c("pepData", "proData", "lipidData", "metabData", "nmrData"))) stop("omicsData must have class of the following, 'pepData', 'proData', 'lipidData', 'metabData', 'nmrData'")
  
  
  #pulling cname attr
  edata_cname<- attr(omicsData, "cnames")$edata_cname
  edata_cname_id<- which(names(omicsData$e_data)== edata_cname)
  fdata_cname<- attr(omicsData, "cnames")$fdata_cname
  
  #we'll count the number of NA values per column
  temp<- omicsData$e_data[, -edata_cname_id]
  na_per_col<- apply(is.na(temp), 2, sum)
  
  na_by_sample<- data.frame("sample_names" = names(temp), "num_NA" = as.numeric(na_per_col))
  names(na_by_sample)[1]<- fdata_cname
  na_by_sample<- merge(na_by_sample, omicsData$f_data, by = fdata_cname)
  names(na_by_sample)[3]<- "group"
  
  #now count the number of NA values per row
  na_per_row<- apply(is.na(temp), 1, sum)
  na_by_molecule<- data.frame("molecule"= omicsData$e_data[,edata_cname_id], "num_NA"= as.numeric(na_per_row))
  names(na_by_molecule)[1]<- edata_cname
  
  result<- list("na.by.sample" = na_by_sample, "na.by.molecule" = na_by_molecule)
  class(result)<- "naRes"
  
  attr(result, "cnames")<- list("edata_cname" = edata_cname, "fdata_cname" = fdata_cname)

  return(result)  
  
}