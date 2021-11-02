#' Create a `multiData` object from multiple omicsData objects
#'
#' @param ... At least two objects of type 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}
#' @param f_meta A data.frame containing sample and group information for all
#' omicsData objects supplied to the function.
#' @param sample_intersect Should only the samples that are common across all
#' datasets be kept in f_meta?  See details for how samples will be dropped.
#' @param keep_sample_info Whether to attempt to append sample information 
#' contained in the objects f_data to the final f_meta via a series of left 
#' joins.  Defaults to FALSE.
#' @param auto_fmeta Whether to attempt to automatically construct f_meta from
#' objects sample information.  Defaults to FALSE.
#' @param match_samples If auto_fmeta = T, whether to attempt to match the 
#' names in the sample columns in f_data across all objects in an attempt to 
#' align them in f_meta.  Defaults to TRUE.
#'
#' @return Object of class 'multiData' containing the omicsData objects, and the
#'sample alignment information f_meta.
#'
#' @details 
#' Object limits:  Currently, as.multiData accepts at most one object from each 
#' of classes 'pepData/proData', 'metabData', 'nmrData', and at most two objects 
#' of class 'lipidData'. 
#' 
#' \code{sample_intersect} will auto-align samples that occur in all datasets.
#' Specifically, it creates a vector of all samples that are common across all
#' datasets, and simply create an f_meta by copying this vector for each
#' dataset and column-binding them.
#' 
#' 
#' @seealso \code{\link{combine_lipidData}} If you want to combine lipidData objects
#' before providing them to as.multiData
#' 
#' @examples
#'
#' \dontrun{
#' library(pmartRdata)
#' library(pmartR)
#' 
#' # Combine lipid and protein object into multidata, both must be log2 + normalized.
#' mylipid_object <- edata_transform(lipid_object, "log2")
#' mylipid_object <- normalize_global(mylipid_object, "all", "median", apply_norm = T)
#' 
#' # Combine without specifically supplying f_meta, either directly, or as one
#' # of the f_datas in any object.
#' mymultidata <- as.multiData(pro_object, mylipid_object, auto_fmeta = T)
#' 
#' # Manually supply an f_meta
#' f_meta <- data.frame(
#' "Proteins" = c(paste0("Mock", 1:3), paste0("Infection", c(1:7)), NA,  "Infection9"),
#' "Lipids" = c(paste0("Mock", 1:3), paste0("Infection", c(1:4)), NA, paste0("Infection", c(6:9))),
#' "Metabolites" = c(paste0("Mock", 1:3), paste0("Infection", c(1:9))),
#' "Condition" = c(rep("A", 3), rep("B", 9))
#' )
#' 
#' mymetab_object <- edata_transform(metab_object, "log2")
#' mymetab_object <- normalize_global(mymetab_object, "all", "median", apply_norm = T)
#' 
#' as.multiData(mylipid_object, pro_object, mymetab_object, f_meta = f_meta)
#' # remove samples that are not common across all data.
#' as.multiData(mylipid_object, pro_object, mymetab_object, f_meta = f_meta, sample_intersect = T)
#' }
#' 
#' @export
#'
as.multiData <-
  function(...,
           f_meta = NULL,
           sample_intersect = F,
           match_samples = T,
           keep_sample_info = F,
           auto_fmeta = F) {
  omicsData_objects <- list(...)
    
  if(length(omicsData_objects) < 2) stop("Must provide at least two datasets.")
  
  # Objects must either be all ungrouped ...
  is_grouped <- sapply(omicsData_objects, 
                         function(x) !is.null(attr(x, "group_DF")))
  if(length(unique(is_grouped)) != 1) {
    stop("All objects must be grouped or ungrouped")
  }
  
  # ... or all grouped and have the same sample names.
  if(all(is_grouped)) {
    for(i in 1:(length(omicsData_objects) - 1)) {
      g1 = attr(omicsData_objects[[i]], "group_DF")$Group
      g2 = attr(omicsData_objects[[i + 1]], "group_DF")$Group
      grp_diff = setdiff(
        union(g1, g2),
        intersect(g1, g2)
      )
      
      if(length(grp_diff) > 0) {
        stop("If objects are grouped, they must have the same group assignments")
      }
    }
  }
  
  # validate object types
  for(obj in omicsData_objects){
    if(!inherits(obj, c('pepData', 'proData', 'metabData','lipidData', 'nmrData'))){
      stop(strwrap( 
        sprintf("Object was expected to have one of type 'pepData', 'proData', 
                'metabData','lipidData', or 'nmrData', but was of type %s",
                 paste(class(obj), collapse = ", ")),
        prefix = " ", initial = ""
        ))
    }
  }
  
  ## Check data scale and normalization status are identical across all objects.
  data_scales <- sapply(omicsData_objects, function(obj) {
    attr(obj, "data_info")$data_scale
  })
  
  if(length(unique(data_scales)) != 1) {
    stop(sprintf("Expected all data to be on the same scale, got data scales: %s", 
                 paste(data_scales, collapse = ", ")))
  } 
  
  is_normed <- sapply(omicsData_objects, function(obj) {
    attr(obj, "data_info")$norm_info$is_normalized
  })
  
  if(length(unique(is_normed)) != 1) {
    stop(strwrap(
          sprintf("Expected all data to be either normalized or unnormalized, 
                  got normalizations statuses: %s", 
                 paste(is_normed, collapse = ", ")),
          prefix = " ", initial = ""
          ))
  }
  ##
  
  obj_types <- sapply(omicsData_objects, class)
  
  # Check that there are an appropriate number of data types.
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
  
  # special check for isobaric data
  for(obj in omicsData_objects) {
    if(inherits(obj, "isobaricpepData") & 
       !isTRUE(attr(obj, "isobaric_info")$norm_info$is_normalized)){
      stop("Isobaric peptide data must be reference pool normalized first.")
    }
  }
  
  ## f_meta construction
  if(!is.null(f_meta)) {
    res <- fmeta_matches(omicsData_objects, f_meta)
    
    if(any(sapply(res, length) == 0)) {
      bad_object_classes = classes[sapply(res, length) == 0]
      stop(
        strwrap(sprintf(
          "Objects of the following types did not have a column in f_meta that 
          contained all samples: %s", 
          paste(bad_object_classes, collapse = " | ")
        )),
        prefix = " ", initial = ""
      )
    }
    
    fmeta_cnames <- find_fmeta_cnames(res)
    
  } else if (auto_fmeta){
    message("Manually combining sample information to make f_meta.")
    fmeta_cols <- lapply(omicsData_objects, function(obj) {
      obj$f_data[,get_fdata_cname(obj)]
    })
    
    # pad the length of each sample info vector to the max length
    maxlen = max(sapply(fmeta_cols, length))
    
    fmeta_cols <- lapply(fmeta_cols, function(x) {
      length(x) <- maxlen
      return(x)
    })
    
    # only match samples in auto_fmeta mode, trust that data frames with sample
    # information are properly aligned
    if(match_samples) {
      allsamps <- unique(unlist(fmeta_cols))
      allsamps <- allsamps[!is.na(allsamps)]
      
      shared_samps <- allsamps
      for(col in fmeta_cols) {
        shared_samps <- intersect(shared_samps, col)
      }
      
      extra_samps = setdiff(allsamps, shared_samps)
      
      fmeta_cols <- lapply(fmeta_cols, function(col) {
        append_samps = extra_samps
        append_samps[which(!(extra_samps %in% col))] <- NA
        c(shared_samps, append_samps)
      })
    } else {
      wrap_message(
        "You chose not to match samples across datasets when creating f_meta 
        from sample information.  This assumes your sample identifiers are 
        row-aligned.")
    }
    
    # 
    fmeta_cnames <- sapply(omicsData_objects, function(obj) {
      paste(get_fdata_cname(obj), class(obj), sep = "_")
    }) %>% make.unique()
    
    
    f_meta <- cbind.data.frame(fmeta_cols) 
    
    colnames(f_meta) <- fmeta_cnames
  } else {
    check_fdatas <- lapply(omicsData_objects, function(obj) {
      fmeta_matches(omicsData_objects, obj$f_data)
    })
    
    unique_cols <- lapply(check_fdatas, function(x) {
      if(all(sapply(x, length) > 0)) {
        unique(unlist(x))
      } else NULL
    })
    
    if(!any(!is.null(unlist(unique_cols)))) {
      stop(strwrap("No f_meta was provided, and none of the sample information
                were valid f_meta.  Either provide a valid f_meta, or specify
                auto_fmeta = T to try and have an f_meta constructed from 
                combined sample information.", prefix = " ", initial = ""))
    }
    
    max_vals = which.max(sapply(unique_cols, length))
    res <- check_fdatas[[max_vals]]
    
    fmeta_cnames <- find_fmeta_cnames(res)
    
    f_meta <- omicsData_objects[[max_vals]] %>% 
      dplyr::select(fmeta_cnames)
    
  }
  
  # 
  if(sample_intersect) {
    allsamps <- unique(unlist(f_meta[,fmeta_cnames]))
    allsamps <- allsamps[!is.na(allsamps)]
    
    shared_samps <- allsamps
    for(col in dplyr::select(f_meta, dplyr::one_of(fmeta_cnames))) {
      shared_samps <- intersect(shared_samps, col)
    }
    
    # apply a custom filter to all datasets, keeping only the intersect of
    # all samples
    omicsData_objects <- lapply(omicsData_objects, function(obj) {
      filt_ <- custom_filter(obj, f_data_keep = shared_samps)
      applyFilt(filt_, obj)
    })
    
    # f_meta will just be a data frame with identical columns
    f_meta <- data.frame(setNames(
      rep(list(shared_samps), length(fmeta_cnames)), fmeta_cnames))
    
  } else {
    if(length(unique(unlist(f_meta[,fmeta_cnames]))) != nrow(f_meta)) {
      wrap_message( 
        "Some samples are not present across all datasets, consider keeping 
          only the intersect with sample_intersect = TRUE")
    }
  }
  
  if(any(sapply(f_meta, function(x) sum(!is.na(x))) < 3)) {
    stop("There were fewer than 3 samples that appear in all datasets.")
  }
  
  # left join sample info across all objects
  if(keep_sample_info) {
    for(i in 1:length(omicsData_objects)) {
      f_meta <- f_meta %>% 
        dplyr::left_join(
          omicsData_objects[[i]]$f_data, 
          by = setNames(get_fdata_cname(omicsData_objects[[i]]), fmeta_cnames[i])
        )
    }
  }
  
  res <- list("omicsData" = omicsData_objects, "f_meta" = f_meta)
  attr(res, "fmeta_samp_cname") <- fmeta_cnames
  class(res) <- "multiData"
  
  return(res)
  
}

#' Check that the f_meta file contains a column aligned to each omicsData objects
#' 
#' @param omicsDataObjects A list of omicsdata objects containing sample 
#' information matching that in f_meta
#' @param f_meta passed from \code{as.multiData}
#' 
#' @return A list, each element of which contains a character vector of column
#' name matches in f_meta for each omicsData object.
#' 
#' @keywords internal
fmeta_matches <- function(omicsData_objects, f_meta){
  res <- lapply(omicsData_objects, function(obj) {
    has_col = sapply(f_meta, function(col) {
      all(obj$f_data[,get_fdata_cname(obj)] %in% col)
    }) 
    
    if(sum(has_col) > 0) colnames(f_meta)[has_col] else NULL
    
  })
  
  return(res)
}

#' Find column names in f_meta for each object.  May return in two objects
#' sharing the same column name.
#' 
#' @param res The output of \code{fmeta_matches}, A list, the i-th element of 
#' which contains a character vector of column name matches in f_meta for the 
#' i-th omicsData object.
#' 
#' @return Character vector containing a column name in f_meta that the i-th 
#' object matches to
#' 
#' @keywords internal
find_fmeta_cnames <- function(res) {
  fmeta_cnames <- character(length(res))
  
  for(i in order(sapply(res, length))){
    rem_cols <- setdiff(res[[i]], fmeta_cnames)
    fmeta_cnames[i] <- if(length(rem_cols) == 0) res[[i]][1] else rem_cols[1]
  }
  return(fmeta_cnames)
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
