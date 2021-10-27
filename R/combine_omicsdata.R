#' Combines two omicsdata objects with identical sample information.
#'
#' @param obj_1/obj_2 Two omicsData objects of the same supported type, 
#' currently "lipidData".  See details for more requirements.
#' @param retain_groups Whether to attempt to apply existing group information 
#' to the new object (defaults to FALSE).
#' @param retain_filters Whether to retain filter information in the new object 
#' (defaults to FALSE).
#'
#' @return An object of the same type as the two input objects, with their 
#' combined data.
#'
#' @details
#' General requirements:
#' 
#' * sample names:  These must be identical for both objects (column names of e_data, and sample identifiers in f_data).
#' * data attributes:  Objects must be on the same scale and both be either normalized or unnormalized.  
#' * group designation:  Objects must have the same grouping structure if retain_groups = T.
#' 
#' @md
#' @export
combine_omicsData <- function(obj_1, obj_2, retain_groups = FALSE, retain_filters = FALSE) {
  if (class(obj_1) != class(obj_2)) {
    stop(sprintf(
      "Objects must be of the same class, found %s and %s",
      class(obj_1),
      class(obj_2)
    ))
  }
  
  # Check that it is among supported objects
  if(!(class(obj_1) %in% c("lipidData"))) stop("Currently only support lipidData")
  
  if(attr(obj_1, "data_info")$norm_info$is_normalized != 
     attr(obj_2, "data_info")$norm_info$is_normalized) {
    stop("Both objects must have the same normalization status (normalized/unnormalized)") 
  }
  
  if(attr(obj_1, "data_info")$data_scale != 
     attr(obj_2, "data_info")$data_scale) {
    stop(sprintf(
      "Objects must be on the same scale, found %s and %s"),
      attr(obj_1, "data_info")$data_scale,
      attr(obj_2, "data_info")$data_scale
    ) 
  }
  
  same_nsamps = attr(obj_1, "data_info")$num_samps == attr(obj_2, "data_info")$num_samps
  samps_diff = setdiff(
    union(
      obj_1$f_data[, get_fdata_cname(obj_1)],
      obj_1$f_data[, get_fdata_cname(obj_2)]
    ),
    intersect(
      obj_1$f_data[, get_fdata_cname(obj_1)],
      obj_1$f_data[, get_fdata_cname(obj_2)]
    )
  )
  
  if(!same_nsamps) {
    stop(sprintf("Number of samples must be the same in both objects, found %s and %s", 
                 attr(obj_1, "data_info")$num_samps, attr(obj_2, "data_info")$num_samps))
  }
  
  if(length(samps_diff) != 0) {
    stop(sprintf("Your sample names did not match, samples found across both datasets: %s",
                 paste(samps_diff, collapse = ", ")))
  }
  
  ## Create the combined e_data
  
  #' we will use the e_data cname from the first object, we did not require that
  #' they both have the same e_data_cname.
  new_edata_cname = get_edata_cname(obj_1)
  
  # bind the two data frames
  new_edata <- dplyr::bind_rows(
    obj_1$e_data, 
    obj_2$e_data %>% 
      dplyr::rename(setNames(
        get_edata_cname(obj_2), 
        get_edata_cname(obj_1)
      ))
  )
  
  #' Combined f_data is simply a left join, since we require the sample names
  #' are the same.
  new_fdata <- obj_1$f_data %>% 
    dplyr::left_join(obj_2$f_data, 
                     by = setNames(get_fdata_cname(obj_2), get_fdata_cname(obj_1)))
  
  #' Combine e_meta in the same way as e_data if it exists in both datasets.
  if(!is.null(obj_1$e_meta) & !is.null(obj_2$e_meta)) {
    new_emeta_cname = get_emeta_cname(obj_1)
    
    new_emeta <- dplyr::bind_rows(
      obj_1$e_meta, 
      obj_2$e_meta %>% 
        rename(setNames(
          get_emeta_cname(obj_2), 
          get_emeta_cname(obj_1)
        ))
    )
    
    #' Check and warn about non-unique e_meta identifiers, this is pre-empting a
    #' situation where this function can take objects with pepData-like e_meta.
    new_emeta_ids = new_emeta[,new_emeta_cname]
    emeta_ids_1 = obj_1$e_meta[,get_emeta_cname(obj_1)]
    emeta_ids_2 = obj_2$e_meta[,get_emeta_cname(obj_2)]
    
    if (length(unique(new_emeta_ids)) != 
        length(unique(emeta_ids_1)) + length(unique(emeta_ids_2))) {
      wrap_warning("There were e_meta identifiers that occurred in both datasets,
              they have been duplicated in the new object's e_meta.")
    }
    
  } else{
    new_emeta_cname = new_emeta = NULL
  }
  
  # Construct the new object using the appropriate type.
  constructor_fn <- get(sprintf("as.%s", class(obj_1)))
  
  new_object <- constructor_fn(
    e_data = new_edata,
    edata_cname = new_edata_cname,
    f_data = new_fdata,
    fdata_cname = get_fdata_cname(obj_1),
    e_meta = new_emeta, 
    emeta_cname = new_emeta_cname,
    data_scale = attr(obj_1, "data_info")$data_scale,
    is_normalized = attr(obj_1, "data_info")$norm_info$is_normalized
  )
  
  # Retain filter information and store it in the new object
  if(retain_filters) {
    filters <- c(
      attr(obj_1, "filters"),
      attr(obj_2, "filters")
    )
    
    attr(new_object, "filters") = filters
  }
  
  #' Set the group designation of the new object, we assume that the grouping is
  #' consistent across both objects.
  if(retain_groups) {
    if(is.null(attr(obj_1, "group_DF")) & is.null(attr(obj_2, "group_DF"))) {
      stop("At least one object must have group information")
    }
    
    if(is.null(attr(obj_1, "group_DF"))) {
      group_df = attr(obj_2, "group_DF")
      tmp_fdata = obj_2$f_data
    } else {
      group_df = attr(obj_1, "group_DF")
      tmp_fdata = obj_1$f_data
    }
    
    samp_info = new_object$f_data[,-which(colnames(new_object$f_data) == 
                                            get_fdata_cname(new_object))]
    
    main_matches = lapply(attr(group_df, "main_effects"), function(x) {
      column_matches_exact(samp_info, tmp_fdata[,x])[1]
    }) %>% unlist()
    
    covariate_matches = if(!is.null(attr(group_df, "covariates"))) {
      lapply(attr(group_df, "covariates"), function(x) {
        column_matches_exact(samp_info, tmp_fdata[,x])[1]
      }) %>% unlist()
    } else NULL
    
    message(sprintf(
      "Grouping new object with main effects: %s.%s",
      paste(main_matches, collapse = ", "),
      if (is.null(covariate_matches))
        ""
      else
        sprintf("  Covariates: %s", paste(covariate_matches, collapse = ", "))
    ))
    
    new_object <- group_designation(
      new_object, main_effects = main_matches, covariates = covariate_matches)
  }
  
  return(new_object)
}