context('normalize: global')

test_that('normalize_global produces the correct output',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Create a pepData object with the reduced data set.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")
  
  # Natural logate the data.
  pdata <- edata_transform(omicsData = pdata,
                           data_scale = "log")
  
  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = 'Condition')
  
  # Create subset vectors ------------------------------------------------------
  
  # Extract peptide IDs for "all" (keep everything) ---------------
  kp_all <- as.character(pdata$e_data[, 1])
  
  # Extract peptide IDs for "los" (keeping top 10%) ---------------
  
  # Construct an empty vector that will hold the peptide IDs for the top 10%
  # order statistics for each column in e_data.
  kp_los <- c()
  
  # Loop through each sample in e_data.
  for (e in 2:13) {
    
    # Extract the eth sample from pdata$e_data.
    e_sample <- pdata$e_data[, e]
    
    # Name each element with the peptide IDs from e_data.
    names(e_sample) <- as.character(pdata$e_data[, 1])
    
    # Order the absolute value of the data in a manner such that the largest
    # value is first, the next largest value is second, the third largest value
    # is third, and so on :)
    o_data <- sort(abs(e_sample), decreasing = TRUE)
    
    # Add the peptide IDs to the kp_los vector. These are the IDs that will be
    # used to calculate the normalizing parameters.
    kp_los <- c(kp_los, names(o_data)[1:15])
    
  }
  
  # Only keep the unique peptide IDs.
  kp_los <- unique(kp_los)
  
  # Extract peptide IDs for "ppp" (threshold = 0.5) ---------------
  
  # Determine which elements of e_data are not NA.
  extant <- !is.na(pdata$e_data[, -1])
  
  # Compute the proportion of peptides present by row.
  ratio <- rowSums(extant)/12
  
  # Divine which peptides have a proportion present above the threshold.
  kp_ppp <- as.character(pdata$e_data[which(ratio >= 0.5), 1])
  
  # Extract peptide IDs for "complete" ---------------
  
  # Deduce which rows have complete data.
  kp_complete <- as.character(
    pdata$e_data[which(complete.cases(pdata$e_data[, -1])), 1]
  )
  
  # Extract peptide IDs for "rip" (threshold = 0.2) ---------------
  
  # Pluck out the rows in e_data that have no missing values (complete rows).
  c_idx <- which(complete.cases(pdata_gdf$e_data[, -1]))
  
  # Generate a vector for the row-wise p-values. These will be calculated by
  # group membership.
  pvals <- rep(NA, length(c_idx))
  
  # Fish out the groups
  groupie <- attr(pdata_gdf, "group_DF")$Group
  
  # Loop through each row with complete data and Kruskal-Wallisate the data.
  for (e in 1:length(c_idx)) {
    
    # Extract the p-value from the Kruskal-Wallis test.
    pvals[[e]] <- kruskal.test(
      as.numeric(pdata_gdf$e_data[c_idx[[e]], -1]) ~ as.factor(groupie)
    )$p.value
    
  }
  
  # Determine which peptide IDs to keep based on the p-value from the KW test.
  kp_rip <- as.character(pdata_gdf$e_data[c_idx[which(pvals > 0.2)], 1])
  
  # Extract peptide IDs for "ppp_rip" (threshold = 0.5, 0.2) ---------------
  
  # Fish out the rows whose proportion of missing values is above the threshold.
  c_idx <- which(pdata_gdf$e_data[, 1] %in% kp_ppp)
  
  # Generate a vector for the row-wise p-values. These will be calculated by
  # group membership.
  pvals <- rep(NA, length(c_idx))
  
  # Loop through each row with complete data and Kruskal-Wallisate the data.
  for (e in 1:length(c_idx)) {
    
    # Extract the p-value from the Kruskal-Wallis test. If all values are
    # missing from a group a zero will be returned. This will essentially
    # prevent these rows from ever being selected.
    pvals[[e]] <- tryCatch (kruskal.test(
      as.numeric(pdata_gdf$e_data[c_idx[[e]], -1]) ~ as.factor(groupie)
    )$p.value,
    error = function(e) {0})
    
  }
  
  # Determine which peptide IDs to keep based on the p-value from the KW test.
  kp_pip <- as.character(pdata_gdf$e_data[c_idx[which(pvals > 0.2)], 1])
  
  # Test normalize_global: normRes ---------------------------------------------
  
  # Fashion a normRes object with "all" and "median".
  norm_all_med <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "median",
                                   apply_norm = FALSE,
                                   backtransform = FALSE)
  
  # Inspect the class of norm_all_med
  expect_s3_class(norm_all_med,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_all_med, "omicsData"),
                   pdata)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_all_med$subset_fn,
                   "all")
  expect_identical(norm_all_med$norm_fn,
                   "median")
  expect_equal(
    norm_all_med$parameters,
    list(normalization = list(scale = NULL,
                              location = apply(pdata$e_data[, -1], 2,
                                               median, na.rm = TRUE)),
         backtransform = list(scale = NULL,
                              location = NULL))
  )
  expect_equal(norm_all_med$n_features_calc,
               150)
  expect_equal(norm_all_med$feature_subset,
               kp_all)
  expect_identical(norm_all_med$prop_features_calc,
                   1)
  
  # Fashion a normRes object with "ppp" and "median".
  norm_ppp_med <- normalize_global(omicsData = pdata,
                                   subset_fn = "ppp",
                                   norm_fn = "median",
                                   apply_norm = FALSE,
                                   backtransform = FALSE)
  
  # Inspect the class of norm_ppp_med
  expect_s3_class(norm_ppp_med,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_ppp_med, "omicsData"),
                   pdata)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_ppp_med$subset_fn,
                   "ppp")
  expect_identical(norm_ppp_med$norm_fn,
                   "median")
  expect_equal(
    norm_ppp_med$parameters,
    list(normalization = list(scale = NULL,
                              location = apply(
                                pdata$e_data[which(kp_all %in% kp_ppp), -1],
                                2, median, na.rm = TRUE
                              )),
         backtransform = list(scale = NULL,
                              location = NULL))
  )
  expect_equal(norm_ppp_med$n_features_calc,
               length(kp_ppp))
  expect_equal(norm_ppp_med$feature_subset,
               kp_ppp)
  expect_identical(norm_ppp_med$prop_features_calc,
                   length(kp_ppp) / length(kp_all))

  # Construct a normRes object with "los" and "mean".
  norm_los_mean <- normalize_global(omicsData = pdata,
                                    subset_fn = "los",
                                    norm_fn = "mean",
                                    params = list(los = 0.1),
                                    apply_norm = FALSE,
                                    backtransform = FALSE)
  
  # Inspect the class of norm_los_mean
  expect_s3_class(norm_los_mean,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_los_mean, "omicsData"),
                   pdata)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_los_mean$subset_fn,
                   "los")
  expect_identical(norm_los_mean$norm_fn,
                   "mean")
  expect_equal(
    norm_los_mean$parameters,
    list(normalization = list(scale = NULL,
                              location = colMeans(
                                pdata$e_data[which(kp_all %in% kp_los), -1],
                                na.rm = TRUE
                              )),
         backtransform = list(scale = NULL,
                              location = NULL))
  )
  expect_equal(norm_los_mean$n_features_calc,
               length(kp_los))
  expect_equal(norm_los_mean$feature_subset,
               kp_los)
  expect_identical(norm_los_mean$prop_features_calc,
                   length(kp_los) / length(kp_all))
  
  # Produce a normRes object with "ppp" and "zscore".
  norm_ppp_z <- normalize_global(omicsData = pdata,
                                 subset_fn = "ppp",
                                 norm_fn = "zscore",
                                 params = list(ppp = 0.5),
                                 apply_norm = FALSE,
                                 backtransform = FALSE)
  
  # Inspect the class of norm_ppp_z
  expect_s3_class(norm_ppp_z,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_ppp_z, "omicsData"),
                   pdata)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_ppp_z$subset_fn,
                   "ppp")
  expect_identical(norm_ppp_z$norm_fn,
                   "zscore")
  expect_equal(
    norm_ppp_z$parameters,
    list(normalization = list(
      scale = apply(pdata$e_data[which(kp_all %in% kp_ppp), -1],
                    2,
                    sd,
                    na.rm = TRUE),
      location = colMeans(pdata$e_data[which(kp_all %in% kp_ppp), -1],
                          na.rm = TRUE)
    ),
         backtransform = list(scale = NULL,
                              location = NULL))
  )
  expect_equal(norm_ppp_z$n_features_calc,
               length(kp_ppp))
  expect_equal(norm_ppp_z$feature_subset,
               kp_ppp)
  expect_identical(norm_ppp_z$prop_features_calc,
                   length(kp_ppp) / length(kp_all))
  
  # Construct a normRes object with "complete" and "mean".
  norm_complete_mean <- normalize_global(omicsData = pdata,
                                         subset_fn = "complete",
                                         norm_fn = "mean",
                                         apply_norm = FALSE,
                                         backtransform = FALSE)
  
  # Inspect the class of norm_complete_mean
  expect_s3_class(norm_complete_mean,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_complete_mean, "omicsData"),
                   pdata)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_complete_mean$subset_fn,
                   "complete")
  expect_identical(norm_complete_mean$norm_fn,
                   "mean")
  expect_equal(
    norm_complete_mean$parameters,
    list(normalization = list(
      scale = NULL,
      location = colMeans(pdata$e_data[which(kp_all %in% kp_complete), -1],
                          na.rm = TRUE)
    ),
    backtransform = list(scale = NULL,
                         location = NULL))
  )
  expect_equal(norm_complete_mean$n_features_calc,
               length(kp_complete))
  expect_equal(norm_complete_mean$feature_subset,
               kp_complete)
  expect_identical(norm_complete_mean$prop_features_calc,
                   length(kp_complete) / length(kp_all))
  
  # Manufacture a normRes object with "rip" and "mad".
  norm_rip_mad <- normalize_global(omicsData = pdata_gdf,
                                   subset_fn = "rip",
                                   norm_fn = "mad",
                                   params = list(rip = 0.2),
                                   apply_norm = FALSE,
                                   backtransform = FALSE)
  
  # Inspect the class of norm_rip_mad
  expect_s3_class(norm_rip_mad,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_rip_mad, "omicsData"),
                   pdata_gdf)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_rip_mad$subset_fn,
                   "rip")
  expect_identical(norm_rip_mad$norm_fn,
                   "mad")
  expect_equal(
    norm_rip_mad$parameters,
    list(normalization = list(
      scale = apply(pdata_gdf$e_data[which(kp_all %in% kp_rip), -1],
                    2,
                    mad,
                    na.rm = TRUE),
      location = apply(pdata_gdf$e_data[which(kp_all %in% kp_rip), -1],
                       2,
                       median,
                       na.rm = TRUE)
    ),
    backtransform = list(scale = NULL,
                         location = NULL))
  )
  expect_equal(norm_rip_mad$n_features_calc,
               length(kp_rip))
  expect_equal(norm_rip_mad$feature_subset,
               kp_rip)
  expect_identical(norm_rip_mad$prop_features_calc,
                   length(kp_rip) / length(kp_all))
  
  # Manufacture a normRes object with "ppp_rip" and "mad".
  norm_pip_mad <- normalize_global(omicsData = pdata_gdf,
                                   subset_fn = "ppp_rip",
                                   norm_fn = "mad",
                                   params = list(ppp_rip = list(ppp = 0.5,
                                                                rip = 0.2)),
                                   apply_norm = FALSE,
                                   backtransform = FALSE)
  
  # Inspect the class of norm_pip_mad
  expect_s3_class(norm_pip_mad,
                  "normRes")
  
  # Survey the normRes attribute.
  expect_identical(attr(norm_pip_mad, "omicsData"),
                   pdata_gdf)
  
  # Examine the elements of the normRes object.
  expect_identical(norm_pip_mad$subset_fn,
                   "ppp_rip")
  expect_identical(norm_pip_mad$norm_fn,
                   "mad")
  expect_equal(
    norm_pip_mad$parameters,
    list(normalization = list(
      scale = apply(pdata_gdf$e_data[which(kp_all %in% kp_pip), -1],
                    2,
                    mad,
                    na.rm = TRUE),
      location = apply(pdata_gdf$e_data[which(kp_all %in% kp_pip), -1],
                       2,
                       median,
                       na.rm = TRUE)
    ),
    backtransform = list(scale = NULL,
                         location = NULL))
  )
  expect_equal(norm_pip_mad$n_features_calc,
               length(kp_pip))
  expect_equal(norm_pip_mad$feature_subset,
               kp_pip)
  expect_identical(norm_pip_mad$prop_features_calc,
                   length(kp_pip) / length(kp_all))
  
  # Test normalize_global: apply normalization ---------------------------------
  
  # Normalize the data using "all" and "median".
  norm_all_med <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "median",
                                   apply_norm = TRUE,
                                   backtransform = FALSE)
  
  # Verify the data was normalized.
  expect_identical(
    norm_all_med$e_data[, -1],
    pdata$e_data[, -1] - rep(apply(pdata$e_data[, -1], 2, median, na.rm = TRUE),
                             each = nrow(pdata$e_data))
  )
  
  # Confirm the attributes that shouldn't have changed didn't change.
  expect_identical(attr(norm_all_med, "cnames"),
                   attr(pdata, "cnames"))
  expect_identical(attr(norm_all_med, "meta_info"),
                   attr(pdata, "meta_info"))
  expect_identical(attr(norm_all_med, "filters"),
                   attr(pdata, "filters"))
  expect_identical(attr(norm_all_med, "data_info")[1:2],
                   attr(pdata, "data_info")[1:2])
  expect_identical(attr(norm_all_med, "data_info")[4:8],
                   attr(pdata, "data_info")[4:8])
  
  # Explore the data_info$norm_info attribute.
  expect_equal(
    attr(norm_all_med, "data_info")$norm_info,
    list(is_normalized = TRUE,
         norm_type = "global",
         subset_fn = "all",
         subset_params = NULL,
         norm_fn = "median",
         n_features_calc = 150,
         prop_features_calc = 1,
         params = list(
           norm_scale = NULL,
           norm_location = apply(pdata$e_data[, -1], 2, median, na.rm = TRUE),
           bt_scale = NULL,
           bt_location = NULL
         ))
  )
  
  # Normalize the data using "all" and "mean".
  norm_all_mea <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "mean",
                                   apply_norm = TRUE,
                                   backtransform = FALSE)
  
  # Verify the data was normalized.
  expect_identical(
    norm_all_mea$e_data[, -1],
    pdata$e_data[, -1] - rep(colMeans(pdata$e_data[, -1], na.rm = TRUE),
                             each = nrow(pdata$e_data))
  )
  
  # Confirm the attributes that shouldn't have changed didn't change.
  expect_identical(attr(norm_all_mea, "cnames"),
                   attr(pdata, "cnames"))
  expect_identical(attr(norm_all_mea, "meta_info"),
                   attr(pdata, "meta_info"))
  expect_identical(attr(norm_all_mea, "filters"),
                   attr(pdata, "filters"))
  expect_identical(attr(norm_all_mea, "data_info")[1:2],
                   attr(pdata, "data_info")[1:2])
  expect_identical(attr(norm_all_mea, "data_info")[4:8],
                   attr(pdata, "data_info")[4:8])
  
  # Explore the data_info$norm_info attribute.
  expect_equal(
    attr(norm_all_mea, "data_info")$norm_info,
    list(is_normalized = TRUE,
         norm_type = "global",
         subset_fn = "all",
         subset_params = NULL,
         norm_fn = "mean",
         n_features_calc = 150,
         prop_features_calc = 1,
         params = list(
           norm_scale = NULL,
           norm_location = colMeans(pdata$e_data[, -1], na.rm = TRUE),
           bt_scale = NULL,
           bt_location = NULL
         ))
  )
  
  # Normalize the data using "all" and "zscore".
  norm_all_z <- normalize_global(omicsData = pdata,
                                 subset_fn = "all",
                                 norm_fn = "zscore",
                                 apply_norm = TRUE,
                                 backtransform = FALSE)
  
  # Verify the data was normalized.
  expect_identical(
    norm_all_z$e_data[, -1],
    (pdata$e_data[, -1] - rep(colMeans(pdata$e_data[, -1], na.rm = TRUE),
                              each = nrow(pdata$e_data))) /
       rep(apply(pdata$e_data[, -1], 2, sd, na.rm = TRUE),
           each = nrow(pdata$e_data))
  )
  
  # Confirm the attributes that shouldn't have changed didn't change.
  expect_identical(attr(norm_all_z, "cnames"),
                   attr(pdata, "cnames"))
  expect_identical(attr(norm_all_z, "meta_info"),
                   attr(pdata, "meta_info"))
  expect_identical(attr(norm_all_z, "filters"),
                   attr(pdata, "filters"))
  expect_identical(attr(norm_all_z, "data_info")[1:2],
                   attr(pdata, "data_info")[1:2])
  expect_identical(attr(norm_all_z, "data_info")[4:8],
                   attr(pdata, "data_info")[4:8])
  
  # Explore the data_info$norm_info attribute.
  expect_equal(
    attr(norm_all_z, "data_info")$norm_info,
    list(is_normalized = TRUE,
         norm_type = "global",
         subset_fn = "all",
         subset_params = NULL,
         norm_fn = "zscore",
         n_features_calc = 150,
         prop_features_calc = 1,
         params = list(
           norm_scale = apply(pdata$e_data[, -1], 2, sd, na.rm = TRUE),
           norm_location = colMeans(pdata$e_data[, -1], na.rm = TRUE),
           bt_scale = NULL,
           bt_location = NULL
         ))
  )
  
  # Normalize the data using "all" and "mad".
  norm_all_mad <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "mad",
                                   apply_norm = TRUE,
                                   backtransform = FALSE)
  
  # Verify the data was normalized.
  expect_identical(
    norm_all_mad$e_data[, -1],
    (pdata$e_data[, -1] - rep(apply(pdata$e_data[, -1], 2, median, na.rm = TRUE),
                              each = nrow(pdata$e_data))) /
      rep(apply(pdata$e_data[, -1], 2, mad, na.rm = TRUE),
          each = nrow(pdata$e_data))
  )
  
  # Confirm the attributes that shouldn't have changed didn't change.
  expect_identical(attr(norm_all_mad, "cnames"),
                   attr(pdata, "cnames"))
  expect_identical(attr(norm_all_mad, "meta_info"),
                   attr(pdata, "meta_info"))
  expect_identical(attr(norm_all_mad, "filters"),
                   attr(pdata, "filters"))
  expect_identical(attr(norm_all_mad, "data_info")[1:2],
                   attr(pdata, "data_info")[1:2])
  expect_identical(attr(norm_all_mad, "data_info")[4:8],
                   attr(pdata, "data_info")[4:8])
  
  # Explore the data_info$norm_info attribute.
  expect_equal(
    attr(norm_all_mad, "data_info")$norm_info,
    list(is_normalized = TRUE,
         norm_type = "global",
         subset_fn = "all",
         subset_params = NULL,
         norm_fn = "mad",
         n_features_calc = 150,
         prop_features_calc = 1,
         params = list(
           norm_scale = apply(pdata$e_data[, -1], 2, mad, na.rm = TRUE),
           norm_location = apply(pdata$e_data[, -1], 2, median, na.rm = TRUE),
           bt_scale = NULL,
           bt_location = NULL
         ))
  )
  
  # Backtransformate -----------------------------------------------------------
  
  # Normalize the data using "all" and "median" with a backtransformation.
  back_all_med <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "median",
                                   apply_norm = TRUE,
                                   backtransform = TRUE)
  
  # Inspect the backtransformation.
  expect_equal(back_all_med$e_data[, -1],
               norm_all_med$e_data[, -1] +
                 median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Scrutinize the backtransform attributes.
  expect_null(attributes(back_all_med)$data_info$norm_info$params$bt_scale)
  expect_equal(attributes(back_all_med)$data_info$norm_info$params$bt_location,
               median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Verify the attributes are correct.
  expect_identical(attributes(back_all_med)$data_info[1:2],
                   attributes(norm_all_med)$data_info[1:2])
  expect_identical(attributes(back_all_med)$data_info[[3]][1:7],
                   attributes(norm_all_med)$data_info[[3]][1:7])
  expect_identical(attributes(back_all_med)$data_info[[3]][[8]][1:2],
                   attributes(norm_all_med)$data_info[[3]][[8]][1:2])
  expect_identical(attributes(back_all_med)$data_info[4:8],
                   attributes(norm_all_med)$data_info[4:8])
  expect_identical(attributes(back_all_med)$meta_info,
                   attributes(norm_all_med)$meta_info)
  expect_identical(attributes(back_all_med)$filters,
                   attributes(norm_all_med)$filters)
  
  # Normalize the data using "all" and "mean" with a backtransformation.
  back_all_mea <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "mean",
                                   apply_norm = TRUE,
                                   backtransform = TRUE)
  
  # Inspect the backtransformation.
  expect_equal(back_all_mea$e_data[, -1],
               norm_all_mea$e_data[, -1] +
                 median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Scrutinize the backtransform attributes.
  expect_null(attributes(back_all_mea)$data_info$norm_info$params$bt_scale)
  expect_equal(attributes(back_all_mea)$data_info$norm_info$params$bt_location,
               median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Verify the attributes are correct.
  expect_identical(attributes(back_all_mea)$data_info[1:2],
                   attributes(norm_all_mea)$data_info[1:2])
  expect_identical(attributes(back_all_mea)$data_info[[3]][1:7],
                   attributes(norm_all_mea)$data_info[[3]][1:7])
  expect_identical(attributes(back_all_mea)$data_info[[3]][[8]][1:2],
                   attributes(norm_all_mea)$data_info[[3]][[8]][1:2])
  expect_identical(attributes(back_all_mea)$data_info[4:8],
                   attributes(norm_all_mea)$data_info[4:8])
  expect_identical(attributes(back_all_mea)$meta_info,
                   attributes(norm_all_mea)$meta_info)
  expect_identical(attributes(back_all_mea)$filters,
                   attributes(norm_all_mea)$filters)
  
  # Normalize the data using "all" and "zscore" with a backtransformation.
  back_all_z <- normalize_global(omicsData = pdata,
                                 subset_fn = "all",
                                 norm_fn = "zscore",
                                 apply_norm = TRUE,
                                 backtransform = TRUE)
  
  # Calculate the pooled standard deviation.
  sample_var <- apply(pdata$e_data[, -1], 2, var, na.rm = TRUE)
  sample_nonmiss <- apply(pdata$e_data[, -1], 2, function(x) sum(!is.na(x)))
  pooled_var <- (sum((sample_nonmiss - 1) * sample_var)) /
    sum(sample_nonmiss - 1)
  pooled_sd <- sqrt(pooled_var)
  
  # Inspect the backtransformation.
  expect_equal(back_all_z$e_data[, -1],
               norm_all_z$e_data[, -1] * pooled_sd +
                 mean(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Scrutinize the backtransform attributes.
  expect_equal(attributes(back_all_z)$data_info$norm_info$params$bt_scale,
               pooled_sd)
  expect_equal(attributes(back_all_z)$data_info$norm_info$params$bt_location,
               mean(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Verify the attributes are correct.
  expect_identical(attributes(back_all_z)$data_info[1:2],
                   attributes(norm_all_z)$data_info[1:2])
  expect_identical(attributes(back_all_z)$data_info[[3]][1:7],
                   attributes(norm_all_z)$data_info[[3]][1:7])
  expect_identical(attributes(back_all_z)$data_info[[3]][[8]][1:2],
                   attributes(norm_all_z)$data_info[[3]][[8]][1:2])
  expect_identical(attributes(back_all_z)$data_info[4:8],
                   attributes(norm_all_z)$data_info[4:8])
  expect_identical(attributes(back_all_z)$meta_info,
                   attributes(norm_all_z)$meta_info)
  expect_identical(attributes(back_all_z)$filters,
                   attributes(norm_all_z)$filters)
  
  # Normalize the data using "all" and "mad" with a backtransformation.
  back_all_mad <- normalize_global(omicsData = pdata,
                                   subset_fn = "all",
                                   norm_fn = "mad",
                                   apply_norm = TRUE,
                                   backtransform = TRUE)
  
  # Compute the pooled mad.
  scale_param <- apply(pdata$e_data[, -1], 2, mad, na.rm = TRUE)
  sample_nonmiss <- apply(pdata$e_data[, -1], 2, function(x) sum(!is.na(x)))
  pooled_mad <- (sum(sample_nonmiss * scale_param)) / sum(sample_nonmiss)
  
  # Inspect the backtransformation.
  expect_equal(back_all_mad$e_data[, -1],
               norm_all_mad$e_data[, -1] * pooled_mad +
                 median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Scrutinize the backtransform attributes.
  expect_equal(attributes(back_all_mad)$data_info$norm_info$params$bt_scale,
               pooled_mad)
  expect_equal(attributes(back_all_mad)$data_info$norm_info$params$bt_location,
               median(as.matrix(pdata$e_data[, -1]), na.rm = TRUE))
  
  # Verify the attributes are correct.
  expect_identical(attributes(back_all_mad)$data_info[1:2],
                   attributes(norm_all_mad)$data_info[1:2])
  expect_identical(attributes(back_all_mad)$data_info[[3]][1:7],
                   attributes(norm_all_mad)$data_info[[3]][1:7])
  expect_identical(attributes(back_all_mad)$data_info[[3]][[8]][1:2],
                   attributes(norm_all_mad)$data_info[[3]][[8]][1:2])
  expect_identical(attributes(back_all_mad)$data_info[4:8],
                   attributes(norm_all_mad)$data_info[4:8])
  expect_identical(attributes(back_all_mad)$meta_info,
                   attributes(norm_all_mad)$meta_info)
  expect_identical(attributes(back_all_mad)$filters,
                   attributes(norm_all_mad)$filters)
  
})
