#'Create a `multiData` object from multiple omicsData objects
#'
#'@param ... At least two objects of type 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}
#'@param f_meta A data.frame containing sample and group information for all 
#'omicsData objects supplied to the function.
#'@param sample_intersect Should only the samples that are common across all
#'datasets be kept in f_meta?  See details for how samples will be dropped.
#'@param combine_lipids Whether to combine lipid-data objects if two are present
#'defaults to FALSE.
#'
#'@export
as.multiData <- function(..., f_meta = NULL, sample_intersect = F, combine_lipids = F) {
  omicsData_objects <- list(...)
  
  if(length(omicsData_objects) < 2) stop("Must provide at least two datasets.")
  
  # combine lipids, return if there is only one lipid object.
  if (combine_lipids & sum(obj_types %in% c("lipidData")) == 2) {
    obj_types <- sapply(omicsData_objects, class)
    lipid_objects = omicsData_objects[which(obj_types %in% c("lipidData"))]
    combined_lipids = combine_omicsData(lipid_objects[[1]], lipid_objects[[2]])
    
    omicsData_objects = omicsData_objects[-which(!(obj_types %in% c("lipidData")))]
    omicsData_objects[[length(omicsData_objects) + 1]] = combined_lipids
    
    if(length(omicsData_objects) < 2){
      res <- list("omicsData" = omicsData_objects, "f_meta" = NULL)
      class(res) <- "multiData"
      return(res)
    }
  }
  
  # validate object types
  for(obj in omicsData_objects){
    if(!inherits(obj, c('pepData', 'proData', 'metabData','lipidData', 'nmrData'))){
      stop(sprintf("Object was expected to have one of type 'pepData', 'proData',
                   'metabData','lipidData', or 'nmrData', but was of type %s",
                   paste(class(obj), collapse = ", ")))
    }
  }
  
  ## Check data scale and normalization status are identical across all objects.
  data_scales <- sapply(omicsData_objects, function(obj) {
    attr(obj, "data_info")$data_scale
  })
  
  if(length(unique(data_scales)) != 1) {
    stop(sprintf("Expected all data to be on the same scale, got data scales:
                 %s", paste(data_scales, collapse = ", ")))
  } 
  
  is_normed <- sapply(omicsData_objects, function(obj) {
    attr(obj, "data_info")$norm_info$is_normalized
  })
  
  if(length(unique(is_normed)) != 1) {
    stop(sprintf("Expected all data to be either normalized or unnormalized, 
                 got normalizations statuses: %s", paste(is_normed, collapse = ", ")))
  }
  ##
  
  # Check that there are an appropriate number of data types.
  obj_types <- sapply(omicsData_objects, class)
  
  if(sum(obj_types %in% c("pepData", "proData")) > 1) {
    stop("There must be no more than 1 object total from types 'pepData' or 'proData'")
  }
  if(sum(obj_types %in% c("lipidData")) > 2) {
    stop("There must be no more than 2 objects total of type 'lipidData'")
  }
  if(sum(obj_types %in% c("metabData")) > 1) {
    stop("There must be no more than 1 object of type 'metabData'")
  }
  if(sum(obj_types %in% c("nmrData")) > 1) {
    stop("There must be no more than 1 object of type 'nmrData'")
  }
  
  if(!is.null(f_meta)) {
    res <- fmeta_matches(omicsData_objects, f_meta)
    if(any(sapply(res, length) == 0)) {
      bad_object_classes = classes[sapply(res, length) == 0]
      stop(sprintf("Objects of the following types did not have a column in
                   f_meta that contained all samples: %s", 
                   paste(bad_object_classes, collapse = " | ")))
    }
  } else {
    message("Manually combining sample information to make f_meta, this assumes 
            that your sample information is row-aligned.")
    
    fmeta_cols <- lapply(omicsData_objects, function(obj) {
      obj$f_data[,get_fdata_cname(obj)]
    })
    
    # pad the length of each sample info vector the max length
    maxlen = max(sapply(fmeta_cols, length))
    
    f_meta_cols <- lapply(fmeta_cols, function(x) {
      length(x) <- maxlen
      return(x)
    })
    
    f_meta_cnames <- sapply(omicsData_objects, function(obj) {
      paste(get_fdata_cname(obj), class(obj), sep = "_")
    })
    
    if(sample_intersect) {
      allsamps <- unique(unlist(f_meta))
      allsamps <- allsamps[!is.na(allsamps)]
      
      for(col in f_meta_cols) {
        allsamps <- intersect(allsamps, col)
      }
      
      if(length(allsamps) < 3) stop("There were fewer than 3 samples that appear in all datasets.")
      
      # f_meta will simply be the concatenation of the same vector, which is the
      # intersection of all sample ids.
      f_meta <- cbind.data.frame(
        lapply(1:length(f_meta_cnames), function(x) allsamps)
      ) 
      
      # apply a custom filter to all datasets, keeping only the intersect of
      # all samples
      omicsData_objects <- lapply(omicsData_objects, function(obj) {
        filt_ <- custom_filter(obj, f_data_keep = allsamps)
        applyFilt(filt_, obj)
      })
      
    } else {
      f_meta <- cbind.data.frame(fmeta_cols) 
    }
  
    colnames(f_meta) <- f_meta_cnames
    
    # Do we want to left join the f_data information from all datasets?
  }
  
  res <- list("omicsData" = omicsData_objects, "f_meta" = f_meta)
  class(res) <- "multiData"
  
  return(res)
  
}

#' Check that the f_meta file contains a column aligned to each omicsData objects
#' 
#' @param omicsDataObjects A list of omicsdata objects containing sample 
#' information matching that in f_meta
#' @param f_meta passed from \code{as.multiData}
#' 
#' @keywords internal
fmeta_matches <- function(omicsData_objects, f_meta){
  res <- lapply(omicsData_objects, function(obj) {
    has_col = sapply(f_meta, function(col) {
      all(obj$f_data[,get_fdata_cname(obj)] %in% col)
    }) 
    
    if(sum(has_col) > 0) colnames(f_meta)[has_col] else NULL
    
  })
}

#' 
#'@export
print.multiData <- function(multiData, ...) {
  classes <- sapply(multiData$omicsData, class)
  
  cat(sprintf("multiData object containing %s omicsData objects\n", length(multiData$omicsData)))
  cat(sprintf("Object Types:  %s\n", paste(classes, collapse = ", ")))
  cat("Sample alignment:\n")
  cat(capture.output(multiData$f_meta), sep = "\n")
  
}

#'
#'@export
summary.multiData <- function(multiData, ...) {
  # Assume data scale and norm status will be consistent across all objects.
  data_scale = unique(sapply(multiData$omicsData, function(x) attr(x, "data_info")$data_scale))
  is_normed <- all(sapply(multiData$omicsData, function(x) attr(x, "data_info")$norm_info$is_normalized))
  
  cat(sprintf("multiData object containing %s %s omicsData objects on the %s scale\n", 
              length(multiData$omicsData),
              if(is_normed) "normalized" else "unnormalized",
              paste(data_scale, collapse = ", ")))
  cat(sprintf("Object Types:  %s\n", paste(classes, collapse = ", ")))
  
  cat("Sample alignment:\n")
  cat(capture.output(multiData$f_meta), sep = "\n")
}
