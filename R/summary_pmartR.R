#' Produce a basic summary of a pmartR omicsData S3 Object
#'
#' This function will provide basic summary statistics for omicsData objects from the pmartR package.
#'
#' @param omicsData an object of the class 'lipidData', 'metabData', 'pepData', 'proData', or 'nmrData' usually created by \code{\link{as.lipidData}}, \code{\link{as.metabData}}, \code{\link{as.pepData}}, \code{\link{as.proData}}, or \code\{\link{as.nmrData}}, respectively.
#'
#' @return a summary table for the pmartR omicsData object. If assigned to a variable, the elements of the summary table are saved in a list format.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' summary(pep_object)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton, Thomas Johansen
#'
#' @export



#'@rdname summary-omicsData
#'@name summary-omicsData
summary.pepData <- function(omicsData) {
  
  verify_data_info(omicsData)
  
  # get values #
  res <- list(class = class(omicsData), num_samps = attr(omicsData, "data_info")$num_samps,
              num_edata = attr(omicsData, "data_info")$num_edata,
              num_emeta = attr(omicsData, "data_info")$num_emeta,
              num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
              prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3))

  # construct output #
  newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(newres, use.names = FALSE))

  # assemble text strings #
  edata_name <- paste("Unique ", attributes(omicsData)$cnames$edata_cname, "s (e_data)", sep="")
  fdata_name <- paste("Unique ", attributes(omicsData)$cnames$fdata_cname, "s (f_data)", sep="")
  if(!is.null(attributes(omicsData)$cnames$emeta_cname)){
    emeta_name <- paste("Unique ", attributes(omicsData)$cnames$emeta_cname, "s (e_meta)", sep="")
  }else{
    emeta_name <- "Rows (e_meta) "
  }
  
  colnames(catmat) <- NULL
  rownames(catmat) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ")
  
  #if group_DF attr is present 
  if(!is.null(attr(omicsData, "group_DF"))){
    group_vec<- attr(omicsData, "group_DF")$Group
    levels<- unique(attr(omicsData, "group_DF")$Group)
    counts <- vector(mode="numeric", length=length(levels))
    
    for(i in 1:length(levels)){
      counts[i]<- length(which(group_vec == levels[i]))
    }
    res2<- as.list(counts)
    names(res2)<- levels
    newres2 <- lapply(res2, function(x) ifelse(is.null(x), "NA", as.character(x)))
    
    newres2<- c(newres, newres2)
    catmat2 <- data.frame(unlist(newres2, use.names = FALSE))
    colnames(catmat2) <- NULL
    rownames(catmat2) <- c("Class",fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ", paste("Samples per group:", levels, sep = " "))
    
    catmat<- catmat2

  }
  return(catmat)
}

#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.proData <- function(omicsData) {

  verify_data_info(omicsData)
  
  # get values #
  res <- list(class = class(omicsData), num_samps = attr(omicsData, "data_info")$num_samps,
              num_edata = attr(omicsData, "data_info")$num_edata,
              num_emeta = attr(omicsData, "data_info")$num_emeta,
              num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
              prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3))

  # construct output #
  newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(newres, use.names = FALSE))

  # assemble text strings #
  edata_name <- paste("Unique ", attributes(omicsData)$cnames$edata_cname, "s (e_data)", sep="")
  fdata_name <- paste("Unique ", attributes(omicsData)$cnames$fdata_cname, "s (f_data)", sep="")
  if(!is.null(attributes(omicsData)$cnames$emeta_cname)){
    emeta_name <- paste("Unique ", attributes(omicsData)$cnames$emeta_cname, "s (e_meta)", sep="")
  }else{
    emeta_name <- "Rows (e_meta) "
  }

  colnames(catmat) <- NULL
  rownames(catmat) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ")

  #if group_DF attr is present 
  if(!is.null(attr(omicsData, "group_DF"))){
    group_vec<- attr(omicsData, "group_DF")$Group
    levels<- unique(attr(omicsData, "group_DF")$Group)
    counts <- vector(mode="numeric", length=length(levels))
    
    for(i in 1:length(levels)){
      counts[i]<- length(which(group_vec == levels[i]))
    }
    res2<- as.list(counts)
    names(res2)<- levels
    newres2 <- lapply(res2, function(x) ifelse(is.null(x), "NA", as.character(x)))
    
    newres2<- c(newres, newres2)
    catmat2 <- data.frame(unlist(newres2, use.names = FALSE))
    colnames(catmat2) <- NULL
    rownames(catmat2) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ", paste("Samples per group:", levels, sep = " "))
    
    catmat<- catmat2
  }

  return(catmat)
}


#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.lipidData <- function(omicsData) {
  
  verify_data_info(omicsData)
  
  # get values #
  res <- list(class = class(omicsData), num_samps = attr(omicsData, "data_info")$num_samps,
              num_edata = attr(omicsData, "data_info")$num_edata,
              num_emeta = attr(omicsData, "data_info")$num_emeta,
              num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
              prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3))

  # construct output #
  newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(newres, use.names = FALSE))

  # assemble text strings #
  edata_name <- paste("Unique ", attributes(omicsData)$cnames$edata_cname, "s (e_data)", sep="")
  fdata_name <- paste("Unique ", attributes(omicsData)$cnames$fdata_cname, "s (f_data)", sep="")
  if(!is.null(attributes(omicsData)$cnames$emeta_cname)){
    emeta_name <- paste("Unique ", attributes(omicsData)$cnames$emeta_cname, "s (e_meta)", sep="")
  }else{
    emeta_name <- "Rows (e_meta) "
  }

  colnames(catmat) <- NULL
  rownames(catmat) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ")

  #if group_DF attr is present 
  if(!is.null(attr(omicsData, "group_DF"))){
    group_vec<- attr(omicsData, "group_DF")$Group
    levels<- unique(attr(omicsData, "group_DF")$Group)
    counts <- vector(mode="numeric", length=length(levels))
    
    for(i in 1:length(levels)){
      counts[i]<- length(which(group_vec == levels[i]))
    }
    res2<- as.list(counts)
    names(res2)<- levels
    newres2 <- lapply(res2, function(x) ifelse(is.null(x), "NA", as.character(x)))
    
    newres2<- c(newres, newres2)
    catmat2 <- data.frame(unlist(newres2, use.names = FALSE))
    colnames(catmat2) <- NULL
    rownames(catmat2) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ", paste("Samples per group:", levels, sep = " "))
    
    catmat<- catmat2
  }

  return(catmat)
}



#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.metabData <- function(omicsData) {
  
  verify_data_info(omicsData)
  
  # get values #
  res <- list(class = class(omicsData), num_samps = attr(omicsData, "data_info")$num_samps,
              num_edata = attr(omicsData, "data_info")$num_edata,
              num_emeta = attr(omicsData, "data_info")$num_emeta,
              num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
              prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3))

  # construct output #
  newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(newres, use.names = FALSE))

  # assemble text strings #
  edata_name <- paste("Unique ", attributes(omicsData)$cnames$edata_cname, "s (e_data)", sep="")
  fdata_name <- paste("Unique ", attributes(omicsData)$cnames$fdata_cname, "s (f_data)", sep="")
  if(!is.null(attributes(omicsData)$cnames$emeta_cname)){
    emeta_name <- paste("Unique ", attributes(omicsData)$cnames$emeta_cname, "s (e_meta)", sep="")
  }else{
    emeta_name <- "Rows (e_meta) "
  }

  colnames(catmat) <- NULL
  rownames(catmat) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ")

  #if group_DF attr is present 
  if(!is.null(attr(omicsData, "group_DF"))){
    group_vec<- attr(omicsData, "group_DF")$Group
    levels<- unique(attr(omicsData, "group_DF")$Group)
    counts <- vector(mode="numeric", length=length(levels))
    
    for(i in 1:length(levels)){
      counts[i]<- length(which(group_vec == levels[i]))
    }
    res2<- as.list(counts)
    names(res2)<- levels
    newres2 <- lapply(res2, function(x) ifelse(is.null(x), "NA", as.character(x)))
    
    newres2<- c(newres, newres2)
    catmat2 <- data.frame(unlist(newres2, use.names = FALSE))
    colnames(catmat2) <- NULL
    rownames(catmat2) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ", paste("Samples per group:", levels, sep = " "))
    
    catmat<- catmat2
  }

  return(catmat)
}


#'@export
#'@rdname summary-omicsData
#'@name summary-omicsData
summary.nmrData <- function(omicsData) {
  
  verify_data_info(omicsData)
  
  # get values #
  res <- list(class = class(omicsData), num_samps = attr(omicsData, "data_info")$num_samps,
              num_edata = attr(omicsData, "data_info")$num_edata,
              num_emeta = attr(omicsData, "data_info")$num_emeta,
              num_miss_obs = attr(omicsData, "data_info")$num_miss_obs,
              prop_missing = round(attr(omicsData, "data_info")$prop_missing, 3))
  
  # construct output #
  newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(newres, use.names = FALSE))
  
  # assemble text strings #
  edata_name <- paste("Unique ", attributes(omicsData)$cnames$edata_cname, "s (e_data)", sep="")
  fdata_name <- paste("Unique ", attributes(omicsData)$cnames$fdata_cname, "s (f_data)", sep="")
  if(!is.null(attributes(omicsData)$cnames$emeta_cname)){
    emeta_name <- paste("Unique ", attributes(omicsData)$cnames$emeta_cname, "s (e_meta)", sep="")
  }else{
    emeta_name <- "Rows (e_meta) "
  }
  
  colnames(catmat) <- NULL
  rownames(catmat) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ")
  
  #if group_DF attr is present 
  if(!is.null(attr(omicsData, "group_DF"))){
    group_vec<- attr(omicsData, "group_DF")$Group
    levels<- unique(attr(omicsData, "group_DF")$Group)
    counts <- vector(mode="numeric", length=length(levels))
    
    for(i in 1:length(levels)){
      counts[i]<- length(which(group_vec == levels[i]))
    }
    res2<- as.list(counts)
    names(res2)<- levels
    newres2 <- lapply(res2, function(x) ifelse(is.null(x), "NA", as.character(x)))
    
    newres2<- c(newres, newres2)
    catmat2 <- data.frame(unlist(newres2, use.names = FALSE))
    colnames(catmat2) <- NULL
    rownames(catmat2) <- c("Class", fdata_name, edata_name, emeta_name, "Missing Observations ", "Proportion Missing ", paste("Samples per group:", levels, sep = " "))
    
    catmat<- catmat2
  }
  
  return(catmat)
}
