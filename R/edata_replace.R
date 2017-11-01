#' Replace Values Equal to x with y
#'
#' This function finds all values of x in the e_data element of omicsData and replaces them with y
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{pepData}, \code{proData}, \code{metabData}, or \code{lipidData}, respectively.
#' @param x value to be replaced, usually numeric or NA
#' @param y replacment value, usually numeric or NA
#'
#' @details This function is often used to replace any 0 values in peptide, protein, metabolite, or lipid data with NA's.
#'
#' @return data object of the same class as omicsData
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' metab_object2 <- edata_replace(omicsData = metab_object, x=0, y=NA)
#'}
#' @author Kelly Stratton
#'
#' @export
edata_replace <- function(omicsData, x, y){
  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("pepData", "proData", "metabData", "lipidData")) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")

  UseMethod("edata_replace")
}


# function for peptides #
#' @export
#' @name edata_replace
#' @rdname edata_replace
edata_replace.pepData <- function(omicsData, x, y){
  edata_id = attr(omicsData, "cnames")$edata_cname
  check_names = getchecknames(omicsData)

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  e_meta <- omicsData$e_meta

  # get count of the number of values replaced #
  if(is.na(x)){
    truefalse <- is.na(edata)
  }else{
    truefalse <- edata==x
  }
  num_replaced <- sum(truefalse, na.rm=TRUE)

  # apply vector_replace to the cols of edata (identifier column has been pulled off) #
  edata_new <- apply(edata, 2, vector_replace, x=x, y=y)

  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  # replace e_data in omicsData with edata_new #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated metabData object (look at as.metabData function to see how the attributes are accessed) #
  updated_data <- as.pepData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, data_scale = attr(omicsData, "data_info")$data_scale, data_norm = attr(omicsData, "data_info")$data_norm, norm_info=attr(omicsData, "data_info")$norm_info, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  message(paste(num_replaced, "instances of", x, "have been replaced with", y, sep=" "))
  return(updated_data)
}



# function for proteins #
#' @export
#' @name edata_replace
#' @rdname edata_replace
edata_replace.proData <- function(omicsData, x, y){
  edata_id = attr(omicsData, "cnames")$edata_cname
  check_names = getchecknames(omicsData)

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  e_meta <- omicsData$e_meta

  # get count of the number of values replaced #
  if(is.na(x)){
    truefalse <- is.na(edata)
  }else{
    truefalse <- edata==x
  }
  num_replaced <- sum(truefalse, na.rm=TRUE)

  # apply vector_replace to the cols of edata (identifier column has been pulled off) #
  edata_new <- apply(edata, 2, vector_replace, x=x, y=y)

  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  # replace e_data in omicsData with edata_new #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated metabData object (look at as.metabData function to see how the attributes are accessed) #
  updated_data <- as.proData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, data_scale = attr(omicsData, "data_info")$data_scale, data_norm = attr(omicsData, "data_info")$data_norm, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  message(paste(num_replaced, "instances of", x, "have been replaced with", y, sep=" "))
  return(updated_data)
}


# function for metabolites #
#' @export
#' @name edata_replace
#' @rdname edata_replace
edata_replace.metabData <- function(omicsData, x, y){

  edata_id = attr(omicsData, "cnames")$edata_cname
  check_names = getchecknames(omicsData)

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  e_meta <- omicsData$e_meta

  # get count of the number of values replaced #
  if(is.na(x)){
    truefalse <- is.na(edata)
  }else{
    truefalse <- edata==x
  }
  num_replaced <- sum(truefalse, na.rm=TRUE)

  # apply vector_replace to the cols of edata (identifier column has been pulled off) #
  edata_new <- apply(edata, 2, vector_replace, x=x, y=y)

  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  # replace e_data in omicsData with edata_new #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated metabData object (look at as.metabData function to see how the attributes are accessed) #
  updated_data <- as.metabData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, data_scale = attr(omicsData, "data_info")$data_scale, data_norm = attr(omicsData, "data_info")$data_norm, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  message(paste(num_replaced, "instances of", x, "have been replaced with", y, sep=" "))
  return(updated_data)
}



# function for lipids #
#' @export
#' @name edata_replace
#' @rdname edata_replace
edata_replace.lipidData <- function(omicsData, x, y){
  edata_id = attr(omicsData, "cnames")$edata_cname
  check_names = getchecknames(omicsData)

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  e_meta <- omicsData$e_meta

  # get count of the number of values replaced #
  if(is.na(x)){
    truefalse <- is.na(edata)
  }else{
    truefalse <- edata==x
  }
  num_replaced <- sum(truefalse, na.rm=TRUE)

  # apply vector_replace to the cols of edata (identifier column has been pulled off) #
  edata_new <- apply(edata, 2, vector_replace, x=x, y=y)

  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  # replace e_data in omicsData with edata_new #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated metabData object (look at as.metabData function to see how the attributes are accessed) #
  updated_data <- as.lipidData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, data_scale = attr(omicsData, "data_info")$data_scale, data_norm = attr(omicsData, "data_info")$data_norm, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  message(paste(num_replaced, "instances of", x, "have been replaced with", y, sep=" "))
  return(updated_data)
}




#' Replace x with y for a single vector
#'
#' @param one_vector numeric vector
#' @param x value to be replaced
#' @param y replacement value
#'
#' @return numeric vector
#'
#' @author Kelly Stratton
#'
vector_replace <- function(one_vector, x, y){
  # find indices where the value is x #
  if(is.na(x)){
    inds <- is.na(one_vector)
  }else{
    inds <- which(one_vector==x)
  }


  # initialize a new vector, which will be returned after we replace x with y #
  new_vector <- one_vector

  # replace x with y #
  new_vector[inds] <- y

  return(new_vector)
}




