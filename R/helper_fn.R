# Functions to get omicsData attributes ----------------------------------------

# Functions to set omicsData attributes ----------------------------------------

# Create a set function that will return the value for each attribute. For 
# example, set_data_info will perform all of the calculations to fill in the
# data_info attribute. These functions will be called in the as.xxx functions
# to create a xxData object but can be used individually to update any one of
# the attributes at a later time.

set_data_info <- function (e_data,
                           edata_cname,
                           data_scale,
                           data_types,
                           norm_info,
                           is_normalized) {
  
  # Identify the column number that contains the IDs.
  id_col <- which(names(e_data) == edata_cname)
  
  # Count the number of missing values in e_data.
  num_miss_obs <- sum(is.na(e_data[, -id_col]))
  
  # Calculate the proportion of missing data in e_data.
  prop_missing <- num_miss_obs / prod(dim(e_data[, -id_col]))
  
  # Extract the number of things (e.g., peptides, lipids, metabolites, ...).
  num_edata <- nrow(e_data)
  
  # Procure the number of samples.
  num_samps <- ncol(e_data) - 1
  
  # Set data normalization information.
  norm_info$is_normalized <- is_normalized
  
  # Return all of the information that belongs in the data_info attribute.
  return (list(data_scale = data_scale,
               norm_info = norm_info,
               num_edata = num_edata,
               num_miss_obs = num_miss_obs,
               prop_missing = prop_missing,
               num_samps = num_samps,
               data_types = data_types))
  
}

set_meta_info <- function (e_meta,
                           emeta_cname) {
  
  # Determine if meta data is present.
  meta_data <- ifelse(is.null(e_meta), FALSE, TRUE)
  
  # Test if emeta_cname is not null.
  if (!is.null(emeta_cname)) {
    
    # Enumerate the number of unique proteins that map to a peptide in e_data.
    num_emeta <- length(unique(e_meta[, emeta_cname]))
    
  } else {
    
    # If emeta_cname is null set the number of proteins to null.
    num_emeta <- NULL
    
  }
  
  # Return the list of meta attributes.
  return (list(meta_data = meta_data,
               num_emeta = num_emeta))
  
}

# This is kind of a stupid function. Don't know why I wrote it.
set_isobaric_info <- function (exp_cname, 
                               channel_cname, 
                               refpool_channel, 
                               refpool_cname, 
                               refpool_notation, 
                               norm_info,
                               isobaric_norm) {
  
  # Set the elements of the norm_info list.
  norm_info$is_normalized <- isobaric_norm
  
  # Return the list of isobaric_info attributes.
  return (list(exp_cname = exp_cname,
               channel_cname = channel_cname, 
               refpool_channel = refpool_channel, 
               refpool_cname = refpool_cname, 
               refpool_notation = refpool_notation, 
               norm_info = norm_info))
  
}

# Another stupid function but what are you going to do?
set_nmr_info <- function (metabolite_name,
                          sample_property_cname,
                          norm_info,
                          nmr_norm,
                          backtransform) {
  
  # Set the elements of the norm_info list.
  norm_info$is_normalized <- nmr_norm
  norm_info$backtransform <- backtransform
  
  # Return the list of objects belonging to this attribute.
  return(list(metabolite_name = metabolite_name,
              sample_property_cname = sample_property_cname,
              norm_info = norm_info))
  
}

# This function will create a filter class object. The output will always have
# the same attributes but not all of them will be used for every data type.
set_filter <- function (threshold,
                        filtered,
                        filter_method = NULL) {
  
  # Create an object that will have the filter elements and class added to it
  # later (in the next 10 lines or so).
  filta <- list()
  
  # Add threshold to filta.
  filta$threshold <- threshold
  
  # Append filtered to the filta list.
  filta$filtered <- filtered
  
  # Affix the filter method to filta. This is only used with imdanovaFilt.
  filta$method <- filter_method
  
  # Create the filter class.
  class(filta) = "filter"
  
  return (filta)
  
}

# To be sorted ---------------

#' Return comparisons of statRes object
#' 
#' This function returns comparisons from statRes or trellData object
#' 
#' @param compObj is an object with the comparison attribute; specifically objects of class 'statRes' and 'trellData' objects derived from 'statRes' objects in \code{\link{format_data}}
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
get_comparisons<- function(compObj){

  #check that compObj object is of 'statRes' or 'trellData' class
  if(!inherits(compObj, c("statRes", "trellData"))) stop("object must be of class 'statRes' or 'trellData'")
  
  #check that compObj object is of 'statRes' or 'trellData' class
  if(inherits(compObj, "trellData") && is.null(attr(compObj, "comparisons"))) stop("trellData object did not inherit 'comparisons' attribute; value is NULL")
  
  #pull comparisons attribute
  comp = attr(compObj, "comparisons")
  
  result = data.frame("comparisons" = as.character(comp), "index" = 1:length(comp), stringsAsFactors = FALSE)
  
  return(result)

}

#' Return data_class of statRes or trellData object
#' 
#' This function returns data_class attribute from statRes or trellData object, inherited from the omicsData used in \code{\link{imd_anova}} or \code{\link{format_data}}
#' 
#' @param dcObj an object of class 'statRes' or 'trellData'
#' @return returns the data_class attribute from a 'statRes' or 'trellData' object
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
#' result = get_data_class(imd_anova_res)
#'}
#' @rdname get_data_class
#' @export
get_data_class<- function(dcObj){
  
  #check that compObj object is of 'statRes' class
  if(!inherits(dcObj, c("statRes", "trellData"))) stop("dcObj object must be of class 'statRes' or 'trellData'")
  
  result = attr(dcObj, "data_class")
  
  return(result)
  
}


#' Return check.names attribute of omicsData object
#' 
#' This function returns check.names attribute from omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @return returns logical, either TRUE or FALSE
#' @examples 
#' dontrun{
#' getchecknames(omicsData)
#'}
#' @rdname getchecknames
#' @export
getchecknames<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  result = attr(omicsData, "check.names")
  
  return(result)
}


#' Set check.names attribute of omicsData object
#' 
#' This function sets the check.names attribute of an omicsData object
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
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
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  #check that set_to is logical
  if(!is.logical(set_to)) stop ("set_to must be of class 'logical' ")
  
  attr(omicsData, "check.names")<- set_to
  
  return(omicsData)
}


#' Get group info of omicsObject object
#' 
#' This function returns the "group_DF" attribute of an omicsObject object
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a data.frame containing omicsObject object group info
#' @examples 
#' dontrun{
#' get_group_info(omicsObject)
#'}
#' @rdname get_group_info
#' @export
get_group_info<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  if(is.null(attr(omicsObject, "group_DF"))) stop("group_designation has not been run for omicsObject")
  
  res = attr(omicsObject, "group_DF")
  
  return(res)
}


#' Get group table
#' 
#' This function returns a table with number of samples per group
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a table containing number of samples per group
#' @examples 
#' dontrun{
#' get_group_table(omicsObject)
#'}
#' @rdname get_group_table
#' @export
get_group_table<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  if(is.null(attr(omicsObject, "group_DF"))) stop("group_designation has not been run for omicsObject")
  
  #should the result be constructed from omicsObject$f_data or from the group_DF attr? if so...
  group = attr(omicsObject, "group_DF")$Group
  
  return(table(group))
}


#' Get data scale
#' 
#' This function returns data scale attribute
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing data scale
#' @examples 
#' dontrun{
#' get_data_scale(omicsObject)
#'}
#' @rdname get_data_scale
#' @export
get_data_scale<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  data_scale = attr(omicsObject, "data_info")$data_scale
  return(data_scale)
}


#' Get data norm
#' 
#' This function returns data norm attribute
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing data norm
#' @examples 
#' dontrun{
#' get_data_norm(omicsObject)
#'}
#' @rdname get_data_norm
#' @export
get_data_norm<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  data_norm = attr(omicsObject, "data_info")$norm_info$is_normalized
  
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
  
  isobaric_norm = attr(omicsData, "isobaric_info")$norm_info$is_normalized
  
  return(isobaric_norm)
}

#' Get e_data cname
#' 
#' This function returns e_data cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing e_data cname
#' @examples 
#' dontrun{
#' get_edata_cname(omicsObject)
#'}
#' @rdname get_edata_cname
#' @export
get_edata_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  edata_cname = attr(omicsObject, "cnames")$edata_cname
  return(edata_cname)
}

#' Get f_data cname
#' 
#' This function returns f_data cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing f_data cname
#' @examples 
#' dontrun{
#' get_fdata_cname(omicsObject)
#'}
#' @rdname get_fdata_cname
#' @export
get_fdata_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData'")
  
  fdata_cname = attr(omicsObject, "cnames")$fdata_cname
  return(fdata_cname)
}

#' Get e_meta cname
#' 
#' This function returns e_meta cname
#' 
#' @param omicsObject an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', 'statRes', or 'trellData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{imd_anova}}, or \code{\link{format_data}} respectively.
#' @return a character string describing e_meta cname
#' @examples 
#' dontrun{
#' get_emeta_cname(omicsObject)
#'}
#' @rdname get_emeta_cname
#' @export
get_emeta_cname<- function(omicsObject){
  # check that omicsObject is of appropriate class #
  if(!inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData", "statRes", "trellData"))) stop("omicsObject must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nrmData', 'statRes', or 'trellData'")
  
  #check if emeta is null
  if(is.null(omicsObject$e_meta) && 
     inherits(omicsObject, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop(
       "e_meta is NULL in omicsObject, thus emeta_cname is also NULL")
  
  if(is.null(omicsObject$e_meta) && 
     inherits(omicsObject, "statRes")) stop(
       "emeta_cname of input statsRes object is dependent on omicsData used in imd_anova; emeta_cname is NULL")
  
  emeta_cname = attr(omicsObject, "cnames")$emeta_cname
  return(emeta_cname)
}

#' Get Filters
#' 
#' This function returns filters attribute
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @return a vector of filter names
#' @examples 
#' dontrun{
#' get_filters(omicsData)
#'}
#' @rdname get_filters
#' @export
get_filters<- function(omicsData){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData", "nmrData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  filters = names(attr(omicsData, "filters"))
  
  if(is.null(filters)) stop("no filters have been applied")
  
  return(filters)
}

#rollup combine_fn functions
combine_fn_mean<- function(x){
  if(all(is.na(x))){mean(x)}else{mean(x, na.rm = T)}
}

combine_fn_median<- function(x){median(x, na.rm = T)}

