#' Normalize an object of class nmrData
#' 
#' The data is normalied either to a spiked in metabolite or to a property 
#' taking sample-specific values
#' 
#' @param omicsData an object of the class 'nmrData'
#' @param apply_norm logical, indicates whether normalization should be applied 
#' to omicsData$e_data. Defaults to FALSE. If TRUE, the normalization is applied 
#' to the data and an S3 object of the same class as \code{omicsData} (e.g. 'nmrData') 
#' with normalized values in \code{e_data} is returned. If FALSE, the normalization 
#' is not applied and an S3 object of the class \code{nmrnormRes} is returned, 
#' allowing some exploratory data analysis prior to subsequently applying the 
#' normalization.
#' @param backtransform logical argument indicating if parameters for back 
#' transforming the data, after normalization, should be calculated. Defaults to 
#' FALSE. If TRUE, the parameters for back transforming the data after normalization 
#' will be calculated, and subsequently included in the data normalization if 
#' \code{apply_norm} is TRUE or \code{\link{apply.normRes}} is used downstream. 
#' See details for an explanation of how these factors are calculated.
#' @param metabolite_name optional character string specifying the name of the 
#' (spiked in) metabolite in \code{e_data} to use for instrument normalization 
#' of the nmrData object. These values will be used to divide the raw abundance 
#' of the corresponding sample in e_data (if e_data is log transformed, this 
#' function accounts for that and returns normalized data on the same scale it 
#' was provided). If using this argument, the 'sample_property_cname' argument 
#' should not be specified.
#' @param sample_property_cname optional character string specifying the name of 
#' the column in f_data containing information to use for instrument normalization 
#' of the nmrData object, such as a concentration. These values will be used to 
#' divide the raw abundance of the corresponding sample in e_data (if e_data is 
#' log transformed, this function accounts for that by temporarily un-log transforming the data 
#' and then returning normalized data on the same scale it was provided). If 
#' using this argument, the 'metabolite_name' argument should not be specified.
#' @details 
#' #' There are two ways to specify the information needed for performing 
#' instrument normalization on an nmrData object:
#' \enumerate{
#' \item specify \code{metabolite_name}. This should be used when normalization 
#' to a spiked in standard is desired. Here \code{metabolite_name} gives the name 
#' of the metabolite in e_data (and e_meta, if present) corresponding to the 
#' spiked in standard. If any samples have a missing value for this metabolite, 
#' an error is returned.
#' \item specify \code{sample_property_cname}. This should be used when 
#' normalization to a sample property, such as concentration, is desired. Here, 
#' \code{sample_property_cname} gives the name of the column in \code{f_data} 
#' which contains the property to use for normalization. If any samples have a 
#' missing value for this column, and error is returned.
#' }
#' @section Backtransform:
#' The purpose of back transforming data is to ensure values are on a scale 
#' similar to their raw values before normalization. The following values are 
#' calculated and/or applied for backtransformation purposes:
#' \tabular{ll}{
#' If normalization using a metabolite in \code{e_data} is specified \tab location 
#' parameter is the median of the values for  \code{metabolite_name} \cr
#' \tab \cr
#' If normalization using a sample property in \code{f_data} is specified \tab 
#' location parameter is the median of the values in \code{sample_property} \cr
#' }
#' See examples below.
#' @examples  
#' dontrun{ 
#' library(pmartRdata)
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object_identified, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = TRUE, metabolite_name = "unkm1.53")
#' 
#' # alternate specification: #
#' data(nmr_object_identified)
#' 
#' nmr_object = edata_transform(nmr_object, "log2")
#' nmr_norm = normalize_nmr(nmr_object, apply_norm = TRUE, sample_property_cname = "Concentration")
#' 
#' }
#' 
#' @export
#'

normalize_nmr<- function(omicsData, apply_norm = FALSE, backtransform = FALSE, metabolite_name = NULL, sample_property_cname = NULL){
  
  ### --------------------------------------------- ###
  ### initial checks ###
  
  # check that omicsData is of correct class
  if(!inherits(omicsData, "nmrData")) stop("omicsData must be of the class 'nmrData'")
  
  # check whether omicsData$e_data is log transformed
  if(!(attr(omicsData, "data_info")$data_scale %in% c('log2', 'log10', 'log'))){
    is_log <- FALSE
  }else{
    is_log <- TRUE
  }
  
  # check that apply_norm is of class logical
  if(!is.logical(apply_norm)) stop("apply_norm must be of class 'logical'")
  
  #check that backtransform is T/F #
  if(!inherits(backtransform, "logical")) stop("backtransform must be of class 'logical'")
  if(backtransform == TRUE){
    backtransform_value <- NA # set it later
  }else{
    # backtransform is FALSE #
    message("backtransform is set to FALSE. Examine the distribution of your data to ensure this is reasonable.")
  }
  
  # if apply_norm is FALSE and backtransform is TRUE, throw error #
  if(apply_norm == FALSE & backtransform == TRUE){stop("apply_norm is set to FALSE and backtransform is set to TRUE; a backtransform cannot be applied if the normalization is not being applied.")}
  
  # pull some attributes from omicsData #
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  
  edata <- omicsData$e_data
  fdata <- omicsData$f_data
  
  # check that metabolite_cname is in e_data, if not NULL #
  if(!is.null(metabolite_name)){
    if(!(metabolite_name %in% omicsData$e_data[, edata_cname])) stop(paste("Metabolite", metabolite_cname, "is not found in e_data. See details of as.nmrData for specifying column names.", sep = " "))
  }
  
  # check that sample_property_cname is in f_data, if not NULL #
  if(!is.null(sample_property_cname)){
    if(!(sample_property_cname %in% names(omicsData$f_data))) stop(paste("Sample characteristic to use for normalization", sample_property_cname, "is not found in f_data. See details of as.nmrData for specifying column names.", sep = " "))
  }
  # make sure the reference pool info info is specified appropriately #
  # possibility 1: specify metabolite_cname #
  poss1 = !is.null(metabolite_name)
  # possibility 2: specify sample_property_cname #
  poss2 = !is.null(sample_property_cname) 
  
  # throw an error if neither or both of these are true #
  if((poss1 + poss2) != 1) stop("Reference metabolite or sample property information was not correctly specified. See Details and Examples for more information.")
  
  # if possibility 1 is used, check that metabolite_name corresponds to a metabolite seen in every sample #
  if(poss1 == TRUE){
    if(!is.character(metabolite_name)) stop("metabolite_name must be of class 'character'")
    reference_metabolite <- omicsData$e_data[grep(metabolite_name, omicsData$e_data) , 2:ncol(omicsData$e_data)]
    if(!all(unlist(lapply(reference_metabolite, class)) == "numeric")){
      # if TRUE, there's a problem with the values since they aren't all numeric #
      stop(paste("Values for", metabolite_name, "must be numeric"))
    }
    reference_metabolite <- as.numeric(omicsData$e_data[grep(metabolite_name, omicsData$e_data) , 2:ncol(omicsData$e_data)])
    backtransform_value <- median(reference_metabolite, na.rm = TRUE)
    
    # make sure none are NA #
    if(any(is.na(reference_metabolite))){stop(paste("Metabolite = ", metabolite_name, " is not observed in every sample. See Details and Examples for more information."))}
    
    # make sure none are 0 #
    if(any(reference_metabolite == 0)){stop(paste("Metabolite = ", metabolite_name, " is not non-zero in every sample. See Details and Examples for more information."))}
  }
  
  # if possibility 2 is used, check that sample_property_cname is present for every sample AND is not zero #
  if(poss2 == TRUE){
    if(!is.character(sample_property_cname)) stop("sample_property_cname must be of class 'character'")
    
    # make sure the values in the column are numeric! #
    sample_property <- omicsData$f_data[, sample_property_cname]
    if(!all(unlist(lapply(sample_property, class)) == "numeric")){
      # if TRUE, there's a problem with the values since they aren't all numeric #
      stop(paste("Values for", sample_property_cname, "must be numeric"))
    }
    
    sample_property <- as.numeric(sample_property)
    backtransform_value <- median(sample_property, na.rm = TRUE)
    
    # make sure none are NA #
    if(any(is.na(sample_property))){stop(paste("sample_property_cname=", sample_property_cname, " is not present for every sample. See Details and Examples for more information."))}
    
    # make sure none are 0 #
    if(any(sample_property == 0)){stop(paste("sample_property_cname=", sample_property_cname, " is not non-zero for every sample. See Details and Examples for more information."))}
  }
  
  ### --------------------------------------------- ###
  
  
  # set up attributes for nmr normalization information #
  attr(omicsData, "nmr_info") = list(metabolite_name = metabolite_name, sample_property_cname = sample_property_cname, norm_info = list(is_normalized = FALSE, backtransform = backtransform))
  
  
  ### --------------------------------------------- ###
  ### use a reference metabolite to normalize the other metabolites ###
  
  if(poss1 == TRUE){
    
    if(apply_norm == TRUE){
      # apply the normalization; remove the reference metabolite from edata as well as the edata_cname column #
      edata_new <- edata[-grep(metabolite_name, edata[, edata_cname]), -grep(edata_cname, names(edata))]
      
      if(is_log == TRUE){
        # we need to subtract reference_metabolite from all the others #
        # note that if we are on log scale, then the reference_metabolite AND the
        # backtransform_value are both also on the log scale #
        edata_new <- edata_new - reference_metabolite
        if(backtransform == TRUE){
          # which log scale is the data on? #
          data_scale <- attributes(omicsData)$data_info$data_scale
          if(data_scale == "log2"){
            edata_new <- edata_new + backtransform_value
          }else{
            if(data_scale == "log10"){
              edata_new <- edata_new + backtransform_value
            }else{
              if(data_scale == "ln"){
                edata_new <- edata_new + backtransform_value
              }else{
                stop("data_scale is not recognized as valid. Check attributes(omcisData)$data_info$data_scale")
              }
            }
          }
        }
      }else{
        # is_log == FALSE, so we need to divide all the others by reference_metabolite #
        edata_new <- edata_new / reference_metabolite
        if(backtransform == TRUE){
          edata_new <- edata_new * backtransform_value
        }
      }
      # reattach the edata_cname column to edata_new
      metabolites_new <- omicsData$e_data[, edata_cname]
      metabolites_new <- metabolites_new[-which(omicsData$e_data[, edata_cname] == metabolite_name)]
      edata_new <- cbind(metabolites_new, edata_new)
      names(edata_new)[1] <- edata_cname
      
      norm_string <- paste0("nmrObject was normalized using metabolite_name: ", metabolite_name)
      norm_params <- reference_metabolite
      
      # replace omicsData$e_data with normalized edata
      omicsData$e_data <- edata_new
      
      # need to update attributes containing summary info, since number of metabolites will change in the case that a row in edata was used to normalize ... probably need to re-create the data object and set the other attributes to their previous values #
      attributes(omicsData)$data_info$num_edata <- length(unique(omicsData$e_data[, edata_cname]))
      attributes(omicsData)$data_info$num_miss_obs <- sum(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
      attributes(omicsData)$data_info$prop_missing <- mean(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
      attributes(omicsData)$data_info$num_samps <- ncol(omicsData$e_data) - 1
      
      
      # update nmr norm flag
      attr(omicsData, "nmr_info")$norm_info$is_normalized = TRUE
      
      # update attribute stating how data was normalized #
      attr(omicsData, "nmr_info")$norm_info$how_normalized <- norm_string
      
      result <- omicsData
      
    }else{
      # apply_norm is FALSE in the reference metabolite case #
      # want a list: Sample (fdata_cname col from f_data), Metabolite (single metabolite name), value (values for that metabolite from e_data)
      rowind <- which(omicsData$e_data[, edata_cname] == metabolite_name)
      result <- list(Sample = omicsData$f_data[, fdata_cname], Metabolite = metabolite_name, value = as.numeric(omicsData$e_data[rowind, -which(names(omicsData$e_data) == edata_cname)]))
      
      # set class #
      class(result) <- c("nmrnormRes", "list")
      
      # set attributes #
      attr(result, "cnames")$edata_cname = edata_cname
      attr(result, "cnames")$fdata_cname = fdata_cname
      attr(result, "nmr_info")$sample_property_cname <- NULL
      attr(result, "nmr_info")$metabolite_name <- metabolite_name
      
    }
  }
  
  ### --------------------------------------------- ###
  ### use a sample property from fdata to normalize the samples ###
  
  if(poss2 == TRUE){
    
    if(apply_norm == TRUE){
      
      edata_new <- edata[ , -grep(edata_cname, names(edata))]
      if(is_log == TRUE){
        # NEW: we need to un-log the data and divide the samples by their 
        # corresponding sample property values; if backtransforming, multiply
        # by that value; then re-log the data
        
        # which log scale is the data on? #
        data_scale <- attributes(omicsData)$data_info$data_scale

        if(backtransform == TRUE){
          if(data_scale == "log2"){
            # un-log2 the data and divide by sample property
            edata_new <- 2^edata_new / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
            # do the backtransform
            edata_new <- edata_new * backtransform_value
            # re-log2 to the data
            # edata_new < log2(edata_new) # this doesn't work...spits out matrix of FALSE and the edata_new itself doesn't get changed from before
            edata_new <- apply(edata_new, 2, log2)
          }else{
            if(data_scale == "log10"){
              # un-log10 the data and divide by sample property
              edata_new <- 10^edata_new / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
              # do the backtransform
              edata_new <- edata_new * log10(backtransform_value)
              # re-log10 the data
              # edata_new <- log10(edata_new)
              edata_new <- apply(edata_new, 2, log10)
            }else{
              if(data_scale == "ln"){
                # un ln the data and divide by sample property
                edata_new <- exp(edata_new) / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
                # do the backtransform
                edata_new <- edata_new * ln(backtransform_value)
                # re-ln the data
                # edata_new <- log(edata_new)
                edata_new <- apply(edata_new, 2, log)
              }else{
                stop("data_scale is not recognized as valid. Check attributes(omcisData)$data_info$data_scale")
              }
            }
          }
        }else{ # backtransform == FALSE
          if(data_scale == "log2"){
            # un-log2 the data and divide by sample property
            edata_new <- 2^edata_new / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
            # re-log2 to the data
            # edata_new < log2(edata_new)
            edata_new <- apply(edata_new, 2, log2)
          }else{
            if(data_scale == "log10"){
              # un-log10 the data and divide by sample property
              edata_new <- 10^edata_new / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
              # re-log10 the data
              # edata_new <- log10(edata_new)
              edata_new <- apply(edata_new, 2, log10)
            }else{
              if(data_scale == "ln"){
                # un ln the data and divide by sample property
                edata_new <- exp(edata_new) / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
                # re-ln the data
                # edata_new <- log(edata_new)
                edata_new <- apply(edata_new, 2, log)
              }else{
                stop("data_scale is not recognized as valid. Check attributes(omcisData)$data_info$data_scale")
              }
            }
          }
        }
        
        edata_new <- data.frame(edata_new) # convert back to data.frame since the apply statements above switch it to matrix, array
        
        
      }else{
        # is_log == FALSE, so we need to divide samples by the sample_property
        # and do the backtransform, if specified to do so #
        edata_new <- edata_new / matrix(sample_property, nrow = nrow(edata_new), ncol = ncol(edata_new), byrow = TRUE)
        if(backtransform == TRUE){
          edata_new <- edata_new * backtransform_value
        }
      }
      edata_new <- cbind(edata[, grep(edata_cname, names(edata))], edata_new)
      names(edata_new)[1] <- edata_cname
      
      norm_string <- paste0("nmrObject was normalized using sample property: ", sample_property_cname)
      norm_params <- sample_property
      
      # replace omicsData$e_data with normalized edata
      omicsData$e_data <- edata_new
      
      # need to update attributes containing summary info, since number of metabolites will change in the case that a row in edata was used to normalize ... probably need to re-create the data object and set the other attributes to their previous values #
      attributes(omicsData)$data_info$num_edata <- length(unique(omicsData$e_data[, edata_cname]))
      attributes(omicsData)$data_info$num_miss_obs <- sum(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
      attributes(omicsData)$data_info$prop_missing <- mean(is.na(omicsData$e_data[,-which(names(omicsData$e_data)==edata_cname)]))
      attributes(omicsData)$data_info$num_samps <- ncol(omicsData$e_data) - 1
      
      
      # update nmr norm flag
      attr(omicsData, "nmr_info")$norm_info$is_normalized = TRUE
      
      # update attribute stating how data was normalized #
      attr(omicsData, "nmr_info")$norm_info$how_normalized <- norm_string
      
      result <- omicsData
      
    }else{
      # apply_norm is FALSE in the sample property case #  
      # want a list: Sample (fdata_cname col from f_data), sample property name (single column name from f_data), value (values for that metabolite from e_data)
      result <- list(Sample = omicsData$f_data[, fdata_cname], Property = sample_property_cname, value = omicsData$f_data[, sample_property_cname])
      
      # set class #
      class(result) <- c("nmrnormRes", "list")
      
      # set attributes #
      attr(result, "cnames")$edata_cname = edata_cname
      attr(result, "cnames")$fdata_cname = fdata_cname
      attr(result, "nmr_info")$sample_property_cname <- sample_property_cname
      attr(result, "nmr_info")$metabolite_name <- NULL
    }
  }
  
  return(result)
}

