#' Normalize an Object of Class nmrData
#'
#' The data is normalized either to a spiked-in metabolite or to a 
#' sample-specific property
#'
#' @param omicsData an object of the class 'nmrData'
#'
#' @param apply_norm logical, indicates whether normalization should be applied
#'   to omicsData$e_data. Defaults to FALSE. If TRUE, the normalization is
#'   applied to the data and an S3 object of the same class as \code{omicsData}
#'   (e.g. 'nmrData') with normalized values in \code{e_data} is returned. If
#'   FALSE, the normalization is not applied and an S3 object of the class
#'   \code{nmrnormRes} is returned, allowing some exploratory data analysis
#'   prior to subsequently applying the normalization.
#'
#' @param backtransform logical argument indicating if parameters for back
#'   transforming the data, after normalization, should be calculated. Defaults
#'   to FALSE. If TRUE, the parameters for back transforming the data after
#'   normalization will be calculated, and subsequently included in the data
#'   normalization if \code{apply_norm} is TRUE or \code{\link{apply.normRes}}
#'   is used downstream. See details for an explanation of how these factors are
#'   calculated.
#'
#' @param metabolite_name optional character string specifying the name of the
#'   (spiked in) metabolite in \code{e_data} to use for instrument normalization
#'   of the nmrData object. These values will be used to divide the raw
#'   abundance of the corresponding sample in e_data (if e_data is log
#'   transformed, this function accounts for that and returns normalized data on
#'   the same scale it was provided). If using this argument, the
#'   'sample_property_cname' argument should not be specified.
#' 
#' @param sample_property_cname optional character string specifying the name of
#'   the column in f_data containing information to use for instrument
#'   normalization of the nmrData object, such as a concentration. These values
#'   will be used to divide the raw abundance of the corresponding sample in
#'   e_data (if e_data is log transformed, this function accounts for that by
#'   temporarily un-log transforming the data and then returning normalized data
#'   on the same scale it was provided). If using this argument, the
#'   'metabolite_name' argument should not be specified.
#'
#' @details There are two ways to specify the information needed for
#' performing instrument normalization on an nmrData object: \enumerate{ \item
#' specify \code{metabolite_name}. This should be used when normalization to a
#' spiked in standard is desired. Here \code{metabolite_name} gives the name of
#' the metabolite in e_data (and e_meta, if present) corresponding to the spiked
#' in standard. If any samples have a missing value for this metabolite, an
#' error is returned. \item specify \code{sample_property_cname}. This should be
#' used when normalizing to a sample property, such as concentration, is
#' desired. Here, \code{sample_property_cname} gives the name of the column in
#' \code{f_data} which contains the property to use for normalization. If any
#' samples have a missing value for this column, and error is returned.}
#' 
#' @section Backtransform: The purpose of back transforming data is to ensure
#' values are on a scale similar to their raw values before normalization. The
#' following values are calculated and/or applied for backtransformation
#' purposes: \tabular{ll}{ If normalization using a metabolite in
#' \code{e_data} is specified \tab location parameter is the median of the
#' values for \code{metabolite_name} \cr \tab \cr If normalization using a
#' sample property in \code{f_data} is specified \tab location parameter is
#' the median of the values in \code{sample_property} \cr } See examples below.
#' 
#' @examples
#' library(pmartRdata)
#'
#' # Normalize using a metabolite (this is merely an example of how to use this specification; the metabolite used was not actually spiked-in for the purpose of normalization)
#' mynmr <- edata_transform(omicsData = nmr_identified_object, 
#'                          data_scale = "log2")
#' nmr_norm <- normalize_nmr(omicsData = mynmr, apply_norm = TRUE,
#'                           metabolite_name = "unkm1.53", 
#'                           backtransform = TRUE)
#'
#' # Normalization using a sample property
#' mynmr <- edata_transform(omicsData = nmr_identified_object, 
#'                          data_scale = "log2")
#' nmr_norm <- normalize_nmr(omicsData = mynmr, apply_norm = TRUE,
#'                          sample_property_cname = "Concentration",
#'                          backtransform = TRUE)
#'
#' @export
#' 
normalize_nmr<- function (omicsData,
                          apply_norm = FALSE,
                          backtransform = FALSE,
                          metabolite_name = NULL,
                          sample_property_cname = NULL) {
  
  ### --------------------------------------------- ###
  ### initial checks ###
  
  # check that omicsData is of correct class
  if (!inherits(omicsData, "nmrData")) {
    
    stop("omicsData must be of the class 'nmrData'")
    
  }
  
  # check that data has not already been nmr normalized
  if (attr(omicsData, "nmr_info")$norm_info$is_normalized == TRUE) {
    stop("omicsData has already undergone NMR normalization")
  }
  
  # Check the current and original data scale of omicsData.
  if (get_data_scale(omicsData) %in% c('log2', 'log10', 'log')) {
    is_log <- TRUE
  } else {
    is_log <- FALSE
  }
  if (get_data_scale_orig(omicsData) %in% c("log", "log2", "log10")) {
    is_log_orig <- TRUE
  } else {
    is_log_orig <- FALSE
  }
  
  # check that apply_norm is of class logical
  if(!is.logical(apply_norm)) stop("apply_norm must be of class 'logical'")
  
  #check that backtransform is T/F #
  if(!inherits(backtransform, "logical")) {
    
    stop("backtransform must be of class 'logical'")
    
  }
  
  # if apply_norm is FALSE and backtransform is TRUE, throw error #
  if(apply_norm == FALSE & backtransform == TRUE) {
    
    stop (paste("apply_norm is set to FALSE and backtransform is set to TRUE;",
                "a backtransform cannot be applied if the normalization is not",
                "being applied.",
                sep = " "))
    
  }
  
  # pull some attributes from omicsData #
  edata_cname <- get_edata_cname(omicsData)
  edata_idx <- which(names(omicsData$e_data) == edata_cname)
  emeta_idx <- which(names(omicsData$e_meta) == edata_cname)
  fdata_cname <- get_fdata_cname(omicsData)
  
  # check that metabolite_cname is in e_data, if not NULL #
  if(!is.null(metabolite_name)){
    if(!(metabolite_name %in% omicsData$e_data[, edata_idx])) {
      stop(paste("Metabolite", metabolite_cname,
                 "is not found in e_data. See details of as.nmrData for",
                 "specifying column names.", sep = " "))
    }
  }
  
  # check that sample_property_cname is in f_data, if not NULL #
  if(!is.null(sample_property_cname)){
    if(!(sample_property_cname %in% names(omicsData$f_data))) {
      stop(paste("Sample characteristic to use for normalization",
                 sample_property_cname,
                 "is not found in f_data. See details of as.nmrData for",
                 "specifying column names.", sep = " "))
    }
  }
  
  # make sure the reference pool info info is specified appropriately #
  # possibility 1: specify metabolite_cname #
  poss1 = !is.null(metabolite_name)
  # possibility 2: specify sample_property_cname #
  poss2 = !is.null(sample_property_cname) 
  
  # throw an error if neither or both of these are true #
  if ((poss1 + poss2) != 1) {
    
    stop (paste("Reference metabolite or sample property information was not",
                "correctly specified. See Details and Examples for more",
                "information.",
                sep = " "))
    
  }
  
  # Prepare possibility 1 info -------------------------------------------------
  
  # if possibility 1 is used, check that metabolite_name corresponds to a
  # metabolite seen in every sample #
  if (poss1 == TRUE) {
    
    # Make sure the name of the reference metabolite is a character string.
    if (!is.character(metabolite_name)) {
      
      stop ("metabolite_name must be of class 'character'")
      
    }
    
    # Fish out the row index of the reference metabolite.
    rowind <- which(omicsData$e_data[, edata_idx] == metabolite_name)
    
    # Extract the row vector (minus the metabolite ID column) corresponding to
    # the reference metabolite.
    reference_metabolite <- omicsData$e_data[rowind , -edata_idx]
    
    # Ensure all values in the row vector are numeric.
    if(!all(unlist(lapply(reference_metabolite, class)) == "numeric")){
      # if TRUE, there's a problem with the values since they aren't all numeric
      stop(paste("Values for", metabolite_name, "must be numeric"))
    }
    
    # make sure none of the reference metabolite values are NA #
    if (any(is.na(reference_metabolite))) {
      
      stop (paste("Metabolite =",
                  metabolite_name,
                  "is not observed in every sample. See Details and Examples",
                  "for more information.",
                  sep = " "))
      
    }
    
    # make sure none of the reference metabolite values are 0 #
    if (any(reference_metabolite == 0)) {
      
      stop (paste("Metabolite =",
                  metabolite_name,
                  "is not non-zero in every sample. See Details and Examples",
                  "for more information.",
                  sep = " "))
      
    }
    
    # Convert to numeric vector from data.frame with single row. This removes
    # column names, row indices, and matrix/data frame dimensions.
    reference_metabolite <- as.numeric(reference_metabolite)
    
    # Compute the median for backtransforming purposes.
    backtransform_value <- median(reference_metabolite, na.rm = TRUE)
    
    # Remove the reference metabolite row from e_data.
    omicsData$e_data <- omicsData$e_data[-rowind, ]
    
    # Find the row index that corresponds to the reference metabolite in e_meta.
    emeta_rowind <- which(omicsData$e_meta[, emeta_idx] == metabolite_name)
    
    # Remove the reference metabolite row from e_meta.
    omicsData$e_meta <- omicsData$e_meta[-emeta_rowind, ]
    
  }
  
  # Prepare possibility 2 info -------------------------------------------------
  
  # if possibility 2 is used, check that sample_property_cname is present for
  # every sample AND is not zero #
  if (poss2 == TRUE) {
    
    # Make sure the sample property column name input is a character string.
    if (!is.character(sample_property_cname)) {
      
      stop ("sample_property_cname must be of class 'character'")
      
    }
    
    # Extract the column in f_data containing the normalizing values.
    sample_property <- omicsData$f_data[, sample_property_cname]
    
    # make sure the values in the column are numeric! #
    if(!all(unlist(lapply(sample_property, class)) == "numeric")){
      # if TRUE, there's a problem with the values since they aren't all numeric
      stop(paste("Values for", sample_property_cname, "must be numeric"))
    }
    
    # make sure none of the normalizing values are NA #
    if (any(is.na(sample_property))) {
      
      stop (paste("sample_property_cname =",
                  sample_property_cname,
                  "is not present for every sample. See Details and Examples",
                  "for more information.",
                  sep = " "))
      
    }
    
    # make sure none of the normalizing values are 0 #
    if (any(sample_property == 0)) {
      
      stop (paste("sample_property_cname =",
                  sample_property_cname,
                  "is not non-zero for every sample. See Details and Examples",
                  "for more information.",
                  sep = " "))
      
    }
    
    # Convert the column vector to a numeric vector. This removes any column
    # names, row indices, and matrix/data frame dimensions. Change the scale if
    # the original and current data scales differ. If the two data scales are
    # the same no change will be made to the scale of the vector.
    sample_property <- mutate_fdata(
      ds = get_data_scale(omicsData),
      ds_orig = get_data_scale_orig(omicsData),
      is_log = is_log,
      is_log_orig = is_log_orig,
      sample_property = as.numeric(sample_property)
    )
    
    # Compute the median for backtransforming purposes.
    backtransform_value <- median(sample_property, na.rm = TRUE)
    
  }
  
  # Normalize using possibility 1 ----------------------------------------------
  
  if (poss1 == TRUE) {
    
    # Check if the normalization will be applied to the data.
    if (apply_norm) {
      
      # Check if the data is on the log scale.
      if (is_log) {
        
        # we need to subtract reference_metabolite from all the others #
        # note that if we are on log scale, then the reference_metabolite AND
        # the backtransform_value are both also on the log scale #
        omicsData$e_data[, -edata_idx] <- (omicsData$e_data[, -edata_idx] -
                                             rep(
                                               reference_metabolite,
                                               each = nrow(omicsData$e_data)
                                             ))
        
        # Check if we should backtransform the data when on the log scale.
        if (backtransform) {
          
          # Add the backtransform value to the normalized data.
          omicsData$e_data[, -edata_idx] <- (omicsData$e_data[, -edata_idx] +
                                               backtransform_value)
          
          # The backtransform argument is set to FALSE.
        } else {
          
          # Throw down a message to let the user know they better be thinking
          # long and hard about the decisions they are making. 
          message(paste("backtransform is set to FALSE. Examine the",
                        "distribution of your data to ensure this is",
                        "reasonable.",
                        sep = " "))
          
        }
        
        # Runs when is_log is FALSE
      } else if (!is_log) {
        
        # Divide all other metabolites by reference_metabolite.
        omicsData$e_data[, -edata_idx] <- (omicsData$e_data[, -edata_idx] /
                                             rep(
                                               reference_metabolite,
                                               each = nrow(omicsData$e_data)
                                             ))
        
        # Check if we need to backtransmute when not on the log scale.
        if (backtransform) {
          
          # Multiply the non-log normalized data by the backtransmogrify value.
          omicsData$e_data[, -edata_idx] <- (omicsData$e_data[, -edata_idx] *
                                               backtransform_value)
          
          # The backtransmogrify argument is set to FALSE.
        } else {
          
          # Throw down a message to let the user know they better be thinking
          # long and hard about the decisions they are making. 
          message(paste("backtransform is set to FALSE. Examine the",
                        "distribution of your data to ensure this is",
                        "reasonable.",
                        sep = " "))
          
        }
        
      }
      
      # Forge the nmr_info attribute.
      attr(omicsData, "nmr_info") <- list(
        metabolite_name = metabolite_name,
        sample_property_cname = sample_property_cname,
        norm_info = list(
          is_normalized = TRUE,
          backtransform = backtransform,
          norm_method = paste("nmrObject was normalized using",
                              "metabolite_name:",
                              metabolite_name,
                              sep = " "),
          norm_params = reference_metabolite
        )
      )
      
      # Extract data_info attribute from omicsData. Some of the elements will be
      # used to update this attribute.
      dInfo <- get_data_info(omicsData)
      # dInfo <- attr(omicsData, 'data_info')
      
      # Update the data_info attribute because we removed a row from e_data.
      attr(omicsData, "data_info") <- set_data_info(
        e_data = omicsData$e_data,
        edata_cname = edata_cname,
        data_scale_orig = get_data_scale_orig(omicsData),
        data_scale = get_data_scale(omicsData),
        data_types = dInfo$data_types,
        norm_info = dInfo$norm_info,
        is_normalized = dInfo$norm_info$is_normalized,
        batch_info = dInfo$batch_info,
        is_bc = dInfo$batch_info$is_bc
      )
      
      # Update the meta_info attribute because we removed a row from e_meta.
      attr(omicsData, 'meta_info') <- set_meta_info(
        e_meta = omicsData$e_meta,
        emeta_cname = get_emeta_cname(omicsData)
      )
      
      # Return the normalized nmrData object for possibility 1.
      return (omicsData)
      
      # Create a nmrnormRes object: apply_norm is FALSE.
    } else {
      
      # apply_norm is FALSE in the reference metabolite case # want a list:
      # Sample (fdata_cname col from f_data), Metabolite (single metabolite
      # name), value (values for that metabolite from e_data)
      result <- list(Sample = omicsData$f_data[, fdata_cname],
                     Metabolite = metabolite_name,
                     value = reference_metabolite)
      
      # Lay down the class.
      class(result) <- c("nmrnormRes", "list")
      
      # Set column name attributes.
      attr(result, "cnames") <- list(edata_cname = edata_cname,
                                     fdata_cname = fdata_cname)
      
      # Prepare the nmr_info attribute.
      attr(result, "nmr_info") <- list(metabolite_name = metabolite_name,
                                       sample_property_cname = NULL)
      
      # Return the un-normalized nmrnormRes object for possibility 1.
      return (result)
      
    }
    
  }
  
  # Normalize using possibility 2 ----------------------------------------------
  
  if (poss2 == TRUE) {
    
    # Check if the normalization will be applied to the data.
    if (apply_norm) {
      
      # NMR normalize the heck out of the data with possibility 2.
      omicsData$e_data <- nmrnorm(omicsData = omicsData,
                                  is_log = is_log,
                                  is_log_orig = is_log_orig,
                                  backtransform = backtransform,
                                  edata_idx = edata_idx,
                                  backtransform_value = backtransform_value,
                                  sample_property = sample_property)
      
      # Forge the nmr_info attribute.
      attr(omicsData, "nmr_info") <- list(
        metabolite_name = metabolite_name,
        sample_property_cname = sample_property_cname,
        norm_info = list(
          is_normalized = TRUE,
          backtransform = backtransform,
          norm_method = paste("nmrObject was normalized using",
                              "sample property:",
                              sample_property_cname,
                              sep = " "),
          norm_params = sample_property
        )
      )
      
      #!#!#!#!#!
      # NOTE: No need to recalculate data_info and meta_info attributes for
      # possibility 2 because no rows were removed from e_data or e_meta.
      #!#!#!#!#!
      
      # Return the normalized nmrData object for possibility 1.
      return (omicsData)
      
      # Create a nmrnormRes object: apply_norm is FALSE.
    } else {
      
      # apply_norm is FALSE in the sample property case # want a list: Sample
      # (fdata_cname col from f_data), sample property name (single column name
      # from f_data), value (values for that metabolite from e_data)
      result <- list(Sample = omicsData$f_data[, fdata_cname],
                     Property = sample_property_cname,
                     value = omicsData$f_data[, sample_property_cname])
      
      # set class #
      class(result) <- c("nmrnormRes", "list")
      
      # Set column name attributes.
      attr(result, "cnames") <- list(edata_cname = edata_cname,
                                     fdata_cname = fdata_cname)
      
      # Prepare the nmr_info attribute.
      attr(result, "nmr_info") <- list(
        metabolite_name = NULL,
        sample_property_cname = sample_property_cname
      )
      
      # Return the un-normalized nmrnormRes object for possibility 2.
      return (result)
      
    }
    
  }
  
}

nmrnorm <- function (omicsData,
                     is_log,
                     is_log_orig,
                     backtransform,
                     edata_idx,
                     backtransform_value,
                     sample_property) {
  
  # Scenario 1: Normalize on abundance scale -----------------------------------
  
  # In this scenario the current data scale is abundance and the original data
  # scale is either abundance or one of the three log scales. If the original
  # data scale is abundance we will normalize the data without any reservations
  # because both data scales are the same. However, if the original data scale
  # is one of the log scales the sample property column in f_data was previously
  # converted to the abundance scale (instead of converting the entire e_data
  # data frame to the log scale). The data will then be normalized on the
  # abundance scale.
  if ((!is_log_orig && !is_log) || (is_log_orig && !is_log)) {
    
    # Normalize the data on the abundance scale.
    omicsData$e_data[, -edata_idx] <- (omicsData$e_data[ , -edata_idx] /
                                         rep(sample_property,
                                             each = nrow(omicsData$e_data)))
    
    # Check if we should backtransmute the data when on the abundance scale.
    if (backtransform) {
      
      # Backtransfigure the data.
      omicsData$e_data[, -edata_idx] <- (omicsData$e_data[ , -edata_idx] *
                                           backtransform_value)
      
      # The backtransform argument is set to FALSE.
    } else {
      
      # Throw down a message to let the user know they better be thinking
      # long and hard about the decisions they are making. 
      message(paste("backtransform is set to FALSE. Examine the",
                    "distribution of your data to ensure this is",
                    "reasonable.",
                    sep = " "))
      
    }
    
  }
  
  # Scenario 2: Normalize on log scale -----------------------------------------
  
  # In this scenario the current data scale is one of the three log scales and
  # the original data scale is either abundance or one of the three log scales.
  # If the original data is on the abundance scale the sample property column in
  # f_data was previously converted to the log scale and the data are normalized
  # on the log scale. If the original data is also on the log scale no changes
  # were made to the sample property column in f_data and the normalization can
  # proceed on the log scale.
  if ((!is_log_orig && is_log) || (is_log_orig && is_log)) {
    
    # Normalize the data on the log scale.
    omicsData$e_data[, -edata_idx] <- (omicsData$e_data[ , -edata_idx] -
                                         rep(sample_property,
                                             each = nrow(omicsData$e_data)))
    
    # Check if we should backtransmute the data when on the log scale.
    if (backtransform) {
      
      # Backtransfigure the data.
      omicsData$e_data[, -edata_idx] <- (omicsData$e_data[ , -edata_idx] +
                                           backtransform_value)
      
      # The backtransform argument is set to FALSE.
    } else {
      
      # Throw down a message to let the user know they better be thinking
      # long and hard about the decisions they are making. 
      message(paste("backtransform is set to FALSE. Examine the",
                    "distribution of your data to ensure this is",
                    "reasonable.",
                    sep = " "))
      
    }
    
  }
  
  # Return the normalized data.
  return (omicsData$e_data)
  
}

# Convert the column in f_data used for normalizing to the correct scale. We are
# assuming that the data in this column of f_data is on the same sacle as the
# data in e_data when it was read in. For example, if the original data scale of
# e_data is abundance and the current data scale is log we will assume the data
# in f_data is on the abundance scale. The data in f_data will then need to be
# converted to the log scale.
# ds - A character vector indicating the current data scale.
# ds_orig - A character vector indicating the original data scale.
# is_log - Logical. Indicates whether the current data is on a log scale.
# is_log_orig - Logical. Indicates whether the original data is on a log scale.
# sample_property - A column of f_data used to normalize e_data.
mutate_fdata <- function (ds,
                          ds_orig,
                          is_log,
                          is_log_orig,
                          sample_property) {
  
  # Original data scale - abundance and current data scale - log.
  if (!is_log_orig && is_log) {
    
    # Convert the sample property to the correct log scale.
    switch (ds,
            
            "log" = {
              
              return (log(sample_property))
              
            },
            
            "log2" = {
              
              return (log2(sample_property))
              
            },
            
            "log10" = {
              
              return (log10(sample_property))
              
            })
    
    # Original data scale - log and current data scale - abundance.
  } else if (is_log_orig && !is_log) {
    
    # Convert the sample property to the abundance scale depending on which log
    # scale the data was originally on.
    switch (ds_orig,
            
            "log" = {
              
              return (exp(sample_property))
              
              },
            
            "log2" = {
              
              return (2^sample_property)
              
              },
            
            "log10" = {
              
              return (10^sample_property)
              
              })
    
    # Original data is on a log scale and the current data is on a different log
    # scale.
  } else if (is_log_orig && is_log) {
    
    # Check if both log scales are the same. If they are then the
    # sample_property column can be returned unaltered.
    if (ds_orig == ds) {
      
      return (sample_property)
      
      # The two log scales are different from each other.
    } else {
      
      # Convert sample_property from one log scale to another.
      return (logify(ds = ds,
                     ds_orig = ds_orig,
                     sample_property = sample_property))
      
    }
    
    # Both data_scale_orig and data_scale are on the abundance scale.
  } else {
    
    return (sample_property)
    
  }
  
}

# Convert the sample_property vector from one log scale to another.
# ds - A character vector indicating the current data scale.
# ds_orig - A character vector indicating the original data scale.
# sample_property - A column of f_data used to normalize e_data.
logify <- function (ds,
                    ds_orig,
                    sample_property) {
  
  # Convert the original data from the natural log scale to another log scale.
  if (ds_orig == "log") {
    
    switch (ds,
            
            "log2" = {
              
              return (log2(exp(sample_property)))
              
            },
            
            "log10" = {
              
              return (log10(exp(sample_property)))
              
            })
    
  }
  
  # Convert the original data from the log 2 scale to another log scale.
  if (ds_orig == "log2") {
    
    switch (ds,
            
            "log" = {
              
              return (log(2^sample_property))
              
            },
            
            "log10" = {
              
              return (log10(2^sample_property))
              
            })
    
  }
  
  # Convert the original data from the log 10 scale to another log scale.
  if (ds_orig == "log10") {
    
    switch (ds,
            
            "log" = {
              
              return (log(10^sample_property))
              
            },
            
            "log2" = {
              
              return (log2(10^sample_property))
              
            })
    
  }
  
}
