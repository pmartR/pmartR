#' Apply a Transformation to the Data
#'
#' This function applies a transformation to the e_data element of omicsData
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively.
#' @param data_scale a character string indicating the type of transformation to be applied to the data. Valid values are: 'log2', 'log', 'log10', or 'abundance'. A value of 'abundance' indicates the data has previously undergone one of the log transformations and should be transformed back to raw values with no transformation applied.
#' 
#' @details This function is intended to be used before analysis of the data begins. Data are typically analyzed on a log scale.
#'
#' @return data object of the same class as omicsData
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' metab_object2 <- edata_transform(omicsData = metab_object, data_scale="log2")
#' attr(metab_object2, "data_info")$data_scale
#'}
#' @author Kelly Stratton, Natalie Heller
#'
#' @export
edata_transform <- function(omicsData, data_scale){
  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData"))) stop("omicsData must be of class 'pepData', 'proData', 'metabData', or 'lipidData'")

  # check that data_scale is one of the acceptable options #
  if(!(data_scale %in% c('log2', 'log10', 'log', 'abundance'))) stop(paste(data_scale, " is not a valid option for 'data_scale'. See details of as.pepData for specifics.", sep=""))

  # if desired scale is log scale, check to make sure data is not already on log scale #
  if(attr(omicsData, "data_info")$data_scale == "log2" & data_scale == "log2"){
    stop("Data is already on log2 scale.")
  }
  if(attr(omicsData, "data_info")$data_scale == "log10" & data_scale == "log10"){
    stop("Data is already on log10 scale.")
  }
  if(attr(omicsData, "data_info")$data_scale == "log" & data_scale == "log"){
    stop("Data is already on (natural) log scale.")
  }

  # if desired scale is abundance, check to make sure data is not already on abundance scale #
  if(data_scale=="abundance" & attr(omicsData, "data_info")$data_scale == "abundance"){
    stop("Data is already on abundance scale.")
  }

  UseMethod("edata_transform")
}



# function for peptides #
#' @export
#' @name edata_transform
#' @rdname edata_transform
edata_transform.pepData <- function(omicsData, data_scale){
  
  check_names = getchecknames(omicsData)
  edata_id = attr(omicsData, "cnames")$edata_cname

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  # apply the transformation #
  if(attr(omicsData, "data_info")$data_scale == "abundance"){
    if(data_scale == "log"){
      edata_new <- log(edata)
    }else{
      if(data_scale == "log2"){
        edata_new <- log2(edata)
      }else{
        if(data_scale == "log10"){
          edata_new <- log10(edata)
        }
      }
    }
  }else{
    if(attr(omicsData, "data_info")$data_scale == "log"){
      if(data_scale == "abundance"){
        edata_new <- exp(edata)
      }else{
        if(data_scale == "log2"){
          edata_new <- log2(exp(edata))
        }else{
          if(data_scale == "log10"){
            edata_new <- log10(exp(edata))
          }
        }
      }
    }else{
      if(attr(omicsData, "data_info")$data_scale == "log2"){
        if(data_scale == "abundance"){
          edata_new <- 2^(edata)
        }else{
          if(data_scale == "log"){
            edata_new <- log(2^(edata))
          }else{
            if(data_scale == "log10"){
              edata_new <- log10(2^(edata))
            }
          }
        }
      }else{
        if(attr(omicsData, "data_info")$data_scale == "log10"){
          if(data_scale == "abundance"){
            edata_new <- 10^(edata)
          }else{
            if(data_scale == "log"){
              edata_new <- log(10^(edata))
            }else{
              if(data_scale == "log2"){
                edata_new <- log2(10^(edata))
              }
            }
          }
        }
      }
    }
  }

  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  if(is.null(omicsData$e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated pepData object #
  updated_data <- as.pepData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, techrep_cname = attr(omicsData, "cnames")$techrep_cname,
                             data_scale = data_scale, is_normalized = attr(omicsData, "data_info")$norm_info$is_normalized, norm_info=attr(omicsData, "data_info")$norm_info, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)
  
  #check for isobaricpepData class
  if(inherits(omicsData, "isobaricpepData")){
    updated_data <- as.isobaricpepData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, techrep_cname = attr(omicsData, "cnames")$techrep_cname, exp_cname = attr(omicsData, "isobaric_info")$exp_cname, channel_cname = attr(omicsData, "isobaric_info")$channel_cname, refpool_channel = attr(omicsData, "isobaric_info")$refpool_channel, refpool_cname = attr(omicsData, "isobaric_info")$refpool_cname, refpool_notation = attr(omicsData, "isobaric_info")$refpool_notation, data_scale = data_scale, is_normalized = attr(omicsData, "data_info")$norm_info$is_normalized, isobaric_norm = attr(omicsData, "isobaric_info")$norm_info$is_normalized, norm_info=attr(omicsData, "data_info")$norm_info, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)
  }
  
  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  return(updated_data)
}



# function for proteins #
#' @export
#' @name edata_transform
#' @rdname edata_transform
edata_transform.proData <- function(omicsData, data_scale){

  check_names = getchecknames(omicsData)
  edata_id = attr(omicsData, "cnames")$edata_cname

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  # apply the transformation #
  if(attr(omicsData, "data_info")$data_scale == "abundance"){
    if(data_scale == "log"){
      edata_new <- log(edata)
    }else{
      if(data_scale == "log2"){
        edata_new <- log2(edata)
      }else{
        if(data_scale == "log10"){
          edata_new <- log10(edata)
        }
      }
    }
  }else{
    if(attr(omicsData, "data_info")$data_scale == "log"){
      if(data_scale == "abundance"){
        edata_new <- exp(edata)
      }else{
        if(data_scale == "log2"){
          edata_new <- log2(exp(edata))
        }else{
          if(data_scale == "log10"){
            edata_new <- log10(exp(edata))
          }
        }
      }
    }else{
      if(attr(omicsData, "data_info")$data_scale == "log2"){
        if(data_scale == "abundance"){
          edata_new <- 2^(edata)
        }else{
          if(data_scale == "log"){
            edata_new <- log(2^(edata))
          }else{
            if(data_scale == "log10"){
              edata_new <- log10(2^(edata))
            }
          }
        }
      }else{
        if(attr(omicsData, "data_info")$data_scale == "log10"){
          if(data_scale == "abundance"){
            edata_new <- 10^(edata)
          }else{
            if(data_scale == "log"){
              edata_new <- log(10^(edata))
            }else{
              if(data_scale == "log2"){
                edata_new <- log2(10^(edata))
              }
            }
          }
        }
      }
    }
  }


  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  if(is.null(omicsData$e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated proData object #
  updated_data <- as.proData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, techrep_cname = attr(omicsData, "cnames")$techrep_cname, data_scale = data_scale, is_normalized = attr(omicsData, "data_info")$norm_info$is_normalized, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  return(updated_data)
}



# function for metabolites #
#' @export
#' @name edata_transform
#' @rdname edata_transform
edata_transform.metabData <- function(omicsData, data_scale){

  check_names = getchecknames(omicsData)
  edata_id = attr(omicsData, "cnames")$edata_cname

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  # apply the transformation #
  if(attr(omicsData, "data_info")$data_scale == "abundance"){
    if(data_scale == "log"){
      edata_new <- log(edata)
    }else{
      if(data_scale == "log2"){
        edata_new <- log2(edata)
      }else{
        if(data_scale == "log10"){
          edata_new <- log10(edata)
        }
      }
    }
  }else{
    if(attr(omicsData, "data_info")$data_scale == "log"){
      if(data_scale == "abundance"){
        edata_new <- exp(edata)
      }else{
        if(data_scale == "log2"){
          edata_new <- log2(exp(edata))
        }else{
          if(data_scale == "log10"){
            edata_new <- log10(exp(edata))
          }
        }
      }
    }else{
      if(attr(omicsData, "data_info")$data_scale == "log2"){
        if(data_scale == "abundance"){
          edata_new <- 2^(edata)
        }else{
          if(data_scale == "log"){
            edata_new <- log(2^(edata))
          }else{
            if(data_scale == "log10"){
              edata_new <- log10(2^(edata))
            }
          }
        }
      }else{
        if(attr(omicsData, "data_info")$data_scale == "log10"){
          if(data_scale == "abundance"){
            edata_new <- 10^(edata)
          }else{
            if(data_scale == "log"){
              edata_new <- log(10^(edata))
            }else{
              if(data_scale == "log2"){
                edata_new <- log2(10^(edata))
              }
            }
          }
        }
      }
    }
  }


  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  if(is.null(omicsData$e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated metabData object #
  updated_data <- as.metabData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, techrep_cname = attr(omicsData, "cnames")$techrep_cname, data_scale = data_scale, is_normalized = attr(omicsData, "data_info")$norm_info$is_normalized, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  return(updated_data)
}


# function for lipids #
#' @export
#' @name edata_transform
#' @rdname edata_transform
edata_transform.lipidData <- function(omicsData, data_scale){

  check_names = getchecknames(omicsData)
  edata_id = attr(omicsData, "cnames")$edata_cname

  edata <- omicsData$e_data
  feature_names <- edata[which(names(edata)==edata_id)]
  # pull off the identifier column #
  edata <- edata[, -which(colnames(edata)==edata_id)]

  # apply the transformation #
  if(attr(omicsData, "data_info")$data_scale == "abundance"){
    if(data_scale == "log"){
      edata_new <- log(edata)
    }else{
      if(data_scale == "log2"){
        edata_new <- log2(edata)
      }else{
        if(data_scale == "log10"){
          edata_new <- log10(edata)
        }
      }
    }
  }else{
    if(attr(omicsData, "data_info")$data_scale == "log"){
      if(data_scale == "abundance"){
        edata_new <- exp(edata)
      }else{
        if(data_scale == "log2"){
          edata_new <- log2(exp(edata))
        }else{
          if(data_scale == "log10"){
            edata_new <- log10(exp(edata))
          }
        }
      }
    }else{
      if(attr(omicsData, "data_info")$data_scale == "log2"){
        if(data_scale == "abundance"){
          edata_new <- 2^(edata)
        }else{
          if(data_scale == "log"){
            edata_new <- log(2^(edata))
          }else{
            if(data_scale == "log10"){
              edata_new <- log10(2^(edata))
            }
          }
        }
      }else{
        if(attr(omicsData, "data_info")$data_scale == "log10"){
          if(data_scale == "abundance"){
            edata_new <- 10^(edata)
          }else{
            if(data_scale == "log"){
              edata_new <- log(10^(edata))
            }else{
              if(data_scale == "log2"){
                edata_new <- log2(10^(edata))
              }
            }
          }
        }
      }
    }
  }


  # add the identifier column back #
  edata_new <- data.frame(edata_id=feature_names, edata_new, check.names=check_names)

  if(is.null(omicsData$e_meta)){
    emeta_cname = NULL
  }else{
    emeta_cname = attr(omicsData, "cnames")$emeta_cname
  }

  # create an updated lipidData object #
  updated_data <- as.lipidData(e_data = edata_new, f_data = omicsData$f_data, e_meta = omicsData$e_meta, edata_cname = edata_id, emeta_cname=emeta_cname, fdata_cname = attr(omicsData, "cnames")$fdata_cname, techrep_cname = attr(omicsData, "cnames")$techrep_cname, data_scale = data_scale, is_normalized = attr(omicsData, "data_info")$norm_info$is_normalized, data_types=attr(omicsData, "data_info")$data_types, check.names = check_names)

  attributes(updated_data)$group_DF <- attributes(omicsData)$group_DF
  attributes(updated_data)$filters <- attributes(omicsData)$filters
  attributes(updated_data)$meta_info <- attributes(omicsData)$meta_info

  return(updated_data)
}







