#' Combines two omicsData objects with identical sample information.
#'
#' @param obj_1 omicsData object of the same supported type as obj_2, currently
#'   "lipidData". See details for more requirements.
#' @param obj_2 omicsData object of the same supported type as obj_1, currently
#'   "lipidData". See details for more requirements.
#' @param retain_groups logical indicator of whether to attempt to apply
#'   existing group information to the new object. Defaults to FALSE.
#' @param retain_filters Whether to retain filter information in the new object 
#' (defaults to FALSE).
#' @param ... Extra arguments, not one of 'omicsData', 'main_effects', or
#'   'covariates' to be passed to `pmartR::group_designation`.
#'   
#' @return An object of the same type as the two input objects, with their 
#' combined data.
#'
#' @details
#' General requirements:
#'
#' * sample names: These must be identical for both objects (column names of
#' e_data, and sample identifiers in f_data)
#' * data attributes: Objects must be on the same scale and both be either
#' normalized or unnormalized
#' * group designation: Objects must have the same grouping structure if
#' retain_groups = T
#'
#' @examples
#' library(pmartRdata)
#' 
#' obj_1 <- lipid_neg_object
#' obj_2 <- lipid_pos_object
#' 
#' # de-duplicate any duplicate edata identifiers
#' all(obj_2$e_data[,get_edata_cname(obj_2)] == obj_2$e_meta[,get_edata_cname(obj_2)])
#' obj_2$e_data[,get_edata_cname(obj_2)] <- paste0("obj_2_", obj_2$e_data[,get_edata_cname(obj_2)])
#' obj_2$e_meta[,get_edata_cname(obj_2)] <- obj_2$e_data[,get_edata_cname(obj_2)]
#' 
#' combine_object <- combine_lipidData(obj_1 = obj_1, obj_2 = obj_2)
#' 
#' # preprocess and group the data and keep filters/grouping structure
#' 
#' obj_1 <- edata_transform(omicsData = obj_1, data_scale = "log2")
#' obj_1 <- normalize_global(omicsData = obj_1, subset_fn = "all", norm_fn = "median", apply_norm = TRUE)
#' obj_2 <- edata_transform(omicsData = obj_2, data_scale = "log2")
#' obj_2 <- normalize_global(omicsData = obj_2, subset_fn = "all", norm_fn = "median", apply_norm = TRUE)
#' 
#' obj_1 <- group_designation(omicsData = obj_1, main_effects = "Virus")
#' obj_2 <- group_designation(omicsData = obj_2, main_effects = "Virus")
#' 
#' obj_1 <- applyFilt(filter_object = molecule_filter(omicsData = obj_1), omicsData = obj_1, min_num = 2)
#' obj_2 <- applyFilt(filter_object = cv_filter(omicsData = obj_2),obj_2, cv_thresh = 60)
#' 
#' combine_object_later <- combine_lipidData(obj_1 = obj_1, obj_2 = obj_2, retain_groups = T, retain_filters = T)
#' 
#' @export
#' 
combine_lipidData <- function(obj_1, obj_2, retain_groups = FALSE, retain_filters = FALSE, drop_duplicate_emeta = TRUE, ...) {
  if (class(obj_1) != class(obj_2)) {
    stop(sprintf(
      "Objects must be of the same class, found %s and %s",
      class(obj_1),
      class(obj_2)
    ))
  }
  
  # Check that it is among supported objects
  if(!(class(obj_1) %in% c("lipidData"))) stop("Currently only support lipidData")
  
  if(get_data_norm(obj_1) != get_data_norm(obj_2))
    stop("Both objects must have the same normalization status (normalized/unnormalized)")
  # if(attr(obj_1, "data_info")$norm_info$is_normalized != 
  #    attr(obj_2, "data_info")$norm_info$is_normalized) {
    # stop("Both objects must have the same normalization status (normalized/unnormalized)") 
  # }
  
  if(attr(obj_1, "data_info")$data_scale != 
     attr(obj_2, "data_info")$data_scale) {
    stop(sprintf(
      "Objects must be on the same scale, found %s and %s",
      attr(obj_1, "data_info")$data_scale,
      attr(obj_2, "data_info")$data_scale
    ))
  }
  
  same_nsamps = attr(obj_1, "data_info")$num_samps == attr(obj_2, "data_info")$num_samps
  samps_diff = setdiff(
    union(
      obj_1$f_data[, get_fdata_cname(obj_1)],
      obj_2$f_data[, get_fdata_cname(obj_2)]
    ),
    intersect(
      obj_1$f_data[, get_fdata_cname(obj_1)],
      obj_2$f_data[, get_fdata_cname(obj_2)]
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
  
  # we will use the e_data cname from the first object, we did not require that
  # they both have the same e_data_cname.
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
  
  molnames <- new_edata[,get_edata_cname(obj_1)]
  
  if(length(molnames) != length(unique(molnames))) {
    warning("Duplicate molecule identifiers were found in your combined data.")
  }
  
  # Combined fdata will keep all columns from the first dataset in the case of
  # duplicates.
  obj_1_fdata_colnames <- obj_1$f_data %>% 
    dplyr::select(-dplyr::one_of(get_fdata_cname(obj_1))) %>% 
    colnames()
  
  new_fdata <- obj_1$f_data %>% 
    dplyr::left_join(
      dplyr::select(obj_2$f_data, -dplyr::one_of(obj_1_fdata_colnames)), 
      by = setNames(get_fdata_cname(obj_2), get_fdata_cname(obj_1))
    )
  
  # Combine e_meta in the same way as e_data if it exists in both datasets.
  if(!is.null(obj_1$e_meta) & !is.null(obj_2$e_meta)) {
    new_emeta_cname = get_emeta_cname(obj_1)
    
    new_emeta <- dplyr::bind_rows(
      obj_1$e_meta, 
      obj_2$e_meta %>% 
        dplyr::rename(setNames(
          get_emeta_cname(obj_2), 
          get_emeta_cname(obj_1)
        ))
    )
    
    # Check and warn about non-unique e_meta identifiers, this is pre-empting a
    # situation where this function can take objects with pepData-like e_meta.
    new_emeta_ids = new_emeta[,new_emeta_cname]
    emeta_ids_1 = obj_1$e_meta[,get_emeta_cname(obj_1)]
    emeta_ids_2 = obj_2$e_meta[,get_emeta_cname(obj_2)]
    
    if (length(unique(new_emeta_ids)) != 
        length(unique(emeta_ids_1)) + length(unique(emeta_ids_2))) {
      
      if(drop_duplicate_emeta) {
        warning(
          "There were non-unique molecule identifiers in e_meta, dropping these duplicates, some meta-data information may be lost."
        )
        new_emeta <-
          new_emeta %>% dplyr::distinct(!!rlang::sym(new_edata_cname), .keep_all = TRUE)
      } else {
        warning(
          "There were non-unique molecule identifiers in e_meta, this may cause the object construction to fail if edata_cname and emeta_cname do not specify unique rows in the combined e_meta"
        )
      }
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
    # data_scale = attr(obj_1, "data_info")$data_scale,
    data_scale = get_data_scale(obj_1),
    # is_normalized = attr(obj_1, "data_info")$norm_info$is_normalized
    is_normalized = get_data_norm(obj_1)
  )
  
  # Retain filter information and store it in the new object
  if(retain_filters) {
    filters <- c(
      attr(obj_1, "filters"),
      attr(obj_2, "filters")
    )
    
    attr(new_object, "filters") = filters
  }
  
  # Set the group designation of the new objects
  if(retain_groups) {
    if(is.null(attr(obj_1, "group_DF")) | is.null(attr(obj_2, "group_DF"))) {
      stop("Both objects must be grouped.")
    }
    
    # check that main effects are functionally the same
    n_orig_groups <- attr(obj_1, "group_DF") %>% 
      dplyr::group_by(Group) %>% 
      attributes() %>% 
      `[[`("groups") %>% 
      nrow()
    
    n_combined_groups <- attr(obj_1, "group_DF") %>% 
      dplyr::left_join(
        attr(obj_2, "group_DF"), 
        by = setNames(get_fdata_cname(obj_2), get_fdata_cname(obj_1))
      ) %>% 
      dplyr::group_by(Group.x, Group.y) %>% 
      attributes() %>% 
      `[[`("groups") %>% 
      nrow()
    
    if(n_orig_groups != n_combined_groups) {
      stop("The main effect structures of the two omicsData objects were not identical.")
    }
    
    ## check covariates ##
    # 1. Rename covariates in each fdata to a temp name
    # 2. Join the two f_datas into a dataframe with the rename covariates
    # 3. Group by the just the first objects covariates and then both the first and second,
    # the number of groups should be the same in both cases if the covariate structure
    # is the same.  If not, throw an error.
    # 4. Run group_designation on the combined object with the first object's main effects/covariates.
    covariates_1 <- attr(obj_1, "group_DF") %>% 
      attributes() %>% 
      `[[`("covariates") %>% 
      {`[`(., -which(colnames(.) == get_fdata_cname(obj_1)))} %>% 
      colnames()
    
    covariates_2 <- attr(obj_2, "group_DF") %>% 
      attributes() %>% 
      `[[`("covariates") %>% 
      {`[`(., -which(colnames(.) == get_fdata_cname(obj_2)))} %>% 
      colnames()
    
    if(all(!sapply(list(covariates_1, covariates_2), is.null))) {
      # renaming ...
      tmp_covar_names_1 <- paste0("_COVARS_1_", 1:length(covariates_1))
      tmp_covar_names_2 <- paste0("_COVARS_2_", 1:length(covariates_2))
      
      rename_map_1 <- setNames(covariates_1, tmp_covar_names_1)
      rename_map_2 <- setNames(covariates_2, tmp_covar_names_2)
      
      tmp_fdata1 <- obj_1$f_data %>% 
        dplyr::rename(!!!rename_map_1)
      tmp_fdata2 <- obj_2$f_data %>% 
        dplyr::rename(!!!rename_map_2)
      
      # ... to perform a join ...
      combined_fdatas <- tmp_fdata1 %>% 
        dplyr::left_join(
          tmp_fdata2,
          by = setNames(get_fdata_cname(obj_2), get_fdata_cname(obj_1))
       )
      
      # ... and check that both objects have the same covariate structure.
      n_orig_covariate_levels_1 <- combined_fdatas %>% 
        dplyr::group_by(
          dplyr::across(dplyr::one_of(tmp_covar_names_1))
        ) %>% 
        attributes() %>% 
        `[[`("groups") %>% 
        nrow()
      
      n_orig_covariate_levels_2 <- combined_fdatas %>% 
        dplyr::group_by(
          dplyr::across(dplyr::one_of(tmp_covar_names_2))
        ) %>% 
        attributes() %>% 
        `[[`("groups") %>% 
        nrow()
      
      n_comb_covariate_levels <- combined_fdatas %>% 
        dplyr::group_by(
          dplyr::across(dplyr::one_of(c(tmp_covar_names_1, tmp_covar_names_2)))
        ) %>% 
        attributes() %>% 
        `[[`("groups") %>% 
        nrow()
      
      if (n_orig_covariate_levels_1 != n_comb_covariate_levels | 
          n_orig_covariate_levels_2 != n_comb_covariate_levels) {
        stop("The covariate structure of both omicsData objects was not identical.")
      }
    }
    
    main_effects <-  attr(obj_1, "group_DF") %>% 
      attributes() %>% 
      `[[`("main_effects")

    message(sprintf(
      "Grouping new object with main effects: %s.%s",
      paste(main_effects, collapse = ", "),
      if (is.null(covariates_1))
        ""
      else
        sprintf("  Covariates: %s", paste(covariates_1, collapse = ", "))
    ))
    
    new_object <- group_designation(
      new_object, 
      main_effects = main_effects, 
      covariates = covariates_1,
      ...)
  }
  
  return(new_object)
}