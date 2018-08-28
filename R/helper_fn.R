#' Return comparisons of statRes object
#' 
#' This function returns comparisons from statRes object
#' 
#' @param statRes is an object of class 'statRes'
#' @return returns a data frame with comparisons and their indices
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' my_prodata = group_designation(omicsData = pro_object, main_effects = c("Condition"))
#' 
#' imdanova_Filt = imdanova_filter(omicsData = my_prodata)
#' my_prodata = applyFilt(filter_object = imdanova_Filt, omicsData = my_prodata, min_nonmiss_anova=2) 
#' imd_anova_res = imd_anova(omicsData = my_prodata, test_method = 'comb', pval_adjust='bon')
#' 
#' result = get_comparisons(imd_anova_res)
#'}
#' @rdname get_comparisons
#' @export
get_comparisons<- function(statRes){

  #check that statRes object is of 'statRes' class
  if(!inherits(statRes, "statRes")) stop("object must be of class 'statRes'")
  
  #pull comparisons attribute
  comp = attr(statRes, "comparisons")
  
  result = data.frame("comparisons" = as.character(comp), "index" = 1:length(comp), stringsAsFactors = FALSE)
  
  return(result)
  
}


#' Return check.names attribute of omicsData object
#' 
#' This function returns check.names attribute from omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return returns logical, either TRUE or FALSE
#' @examples 
#' dontrun{
#' getchecknames(omicsData)
#'}
#' @rdname getchecknames
#' @export
getchecknames<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  result = attr(omicsData, "check.names")
  
  return(result)
}


#' Set check.names attribute of omicsData object
#' 
#' This function sets the check.names attribute of an omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @param set_to logical indicating what to set check.names attribute to. Defaults to TRUE.
#' @return updated omicsData object with check.names attribute
#' @examples 
#' dontrun{
#' setchecknames(omicsData, set_to = TRUE)
#'}
#' @rdname setchecknames
#' @export
setchecknames<- function(omicsData, set_to = TRUE){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  #check that set_to is logical
  if(!is.logical(set_to)) stop ("set_to must be of class 'logical' ")
  
  attr(omicsData, "check.names")<- set_to
  
  return(omicsData)
}


#' Get group info of omicsData object
#' 
#' This function returns the "group_DF" attribute of an omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a data.frame containing omicsData object group info
#' @examples 
#' dontrun{
#' get_group_info(omicsData)
#'}
#' @rdname get_group_info
#' @export
get_group_info<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  if(is.null(attr(omicsData, "group_DF"))) stop("group_designation has not been run on omicsData")
  
  res = attr(omicsData, "group_DF")
  
  return(res)
}


#' Get group table
#' 
#' This function returns a table with number of samples per group
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a table containing number of samples per group
#' @examples 
#' dontrun{
#' get_group_table(omicsData)
#'}
#' @rdname get_group_table
#' @export
get_group_table<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  if(is.null(attr(omicsData, "group_DF"))) stop("group_designation has not been run on omicsData")
  
  #should the result be constructed from omicsData$f_data or from the group_DF attr? if so...
  group = attr(omicsData, "group_DF")$Group
  
  return(table(group))
}


#' Get data scale
#' 
#' This function returns data scale attribute
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a character string describing data scale
#' @examples 
#' dontrun{
#' get_data_scale(omicsData)
#'}
#' @rdname get_data_scale
#' @export
get_data_scale<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  data_scale = attr(omicsData, "data_info")$data_scale
  return(data_scale)
}


#' Get data norm
#' 
#' This function returns data norm attribute
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a character string describing data norm
#' @examples 
#' dontrun{
#' get_data_norm(omicsData)
#'}
#' @rdname get_data_norm
#' @export
get_data_norm<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  data_norm = attr(omicsData, "data_info")$data_norm
  
  if(inherits(omicsData, "isobaricpepData")){
    data_norm =  data_norm = attr(omicsData, "data_info")$isobaric_norm
  }
  
  return(data_norm)
}

#' Get isobaric norm
#' 
#' This function returns the isobaric norm attribute
#' 
#' @param omicsData an object of the class 'pepData', 'isobaricpepData' or 'proData', usually created by \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}}.
#' @return a character string describing isobaric norm attribute
#' @examples 
#' dontrun{
#' get_isobaric_norm(omicsData)
#'}
#' @rdname get_isobaric_norm
#' @export
get_isobaric_norm<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "isobaricpepData"))) stop("omicsData must be of class 'pepData', 'proData' or 'isobaricpepData'")
  
  isobaric_norm = attr(omicsData, "data_info")$isobaric_norm
  
  return(isobaric_norm)
}

#' Get e_data cname
#' 
#' This function returns e_data cname
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a character string describing e_data cname
#' @examples 
#' dontrun{
#' get_edata_cname(omicsData)
#'}
#' @rdname get_edata_cname
#' @export
get_edata_cname<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  edata_cname = attr(omicsData, "cnames")$edata_cname
  return(edata_cname)
}

#' Get f_data cname
#' 
#' This function returns f_data cname
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a character string describing f_data cname
#' @examples 
#' dontrun{
#' get_fdata_cname(omicsData)
#'}
#' @rdname get_fdata_cname
#' @export
get_fdata_cname<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  fdata_cname = attr(omicsData, "cnames")$fdata_cname
  return(fdata_cname)
}

#' Get e_meta cname
#' 
#' This function returns e_meta cname
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a character string describing e_meta cname
#' @examples 
#' dontrun{
#' get_emeta_cname(omicsData)
#'}
#' @rdname get_emeta_cname
#' @export
get_emeta_cname<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  #check if emeta is null
  if(is.null(omicsData$e_meta)) stop("e_meta is NULL in omicsData, thus emeta_cname is also NULL")
  
  emeta_cname = attr(omicsData, "cnames")$emeta_cname
  return(emeta_cname)
}

#' Get Filters
#' 
#' This function returns filters attribute
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @return a vector of filter names
#' @examples 
#' dontrun{
#' get_filters(omicsData)
#'}
#' @rdname get_filters
#' @export
get_filters<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")
  
  filters = names(attr(omicsData, "filters"))
  
  if(is.null(filters)) stop("no filters have been applied")
  
  return(filters)
}

#rollup combine_fn functions
combine_fn_mean<- function(x){
  if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}
}

combine_fn_median<- function(x){median(x, na.rm = T)}

