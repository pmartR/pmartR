context('normalize: zero one scaling')

test_that('normalize_zero_one_scaling produces the correct output', {
  # Load the nmr data frames ---------------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'
  ))
  
  # Create a pepData object with the reduced data set.
  pdata <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID",
    emeta_cname = "Protein"
  )
  
  expect_error(
    normalize_zero_one_scaling(omicsData = fdata),
    "omicsData must be of class "
  )
  
  expect_error(
    normalize_zero_one_scaling(omicsData = pdata),
    " should be log transformed"
  )
  
  # log2 the data.
  pdata <- edata_transform(
    omicsData = pdata,
    data_scale = "log2"
  )
  
  # Calculate the normalization standards --------------------------------------
  
  # Create a copy of edata. This will become the standard for the normalized
  # data.
  edata_stand <- edata
  
  # Normalize data: abundance scale --------------------------------------------

  # Normalize the nmr data on the abundance scale.
  a_norm <- normalize_zero_one_scaling(omicsData = pdata)
  
  # Put on my sleuth hat and start investigating the output.
  expect_identical(dim(a_norm$e_data), dim(edata_stand))
  expect_s3_class(a_norm, "pepData")
  
  # Check out the attributes.
  
  ## Norm info is different, so remove for comparison
  atr_norm <- attributes(a_norm)
  atr_norm$data_info$norm_info <- NULL
  atr_std <- attributes(pdata)
  atr_std$data_info$norm_info <- NULL
  
  expect_identical(atr_std, atr_norm)
  
  expect_true(get_data_norm(a_norm))
  
  ni <- get_data_info(a_norm)$norm_info
  
  expect_true(ni$norm_type == "zero_to_one")
  expect_true(ni$n_features_calc == 150)
  
  min_max <- rbind(apply(pdata$e_data[-1], 2, min, na.rm = T),
                   apply(pdata$e_data[-1], 2, max, na.rm = T))
  
  expect_identical(ni$params$norm_scale, min_max)
  
  ## Check NULLs
  null_atr <- purrr::map_lgl(
    c("subset_fn", "subset_params", "prop_features_calc"), function(str){
      is.null(ni[[str]])
  })
  expect_true(all(null_atr))
  null_atr <- purrr::map_lgl(
    c("norm_location", "bt_scale", "bt_location"), function(str){
      is.null(ni$params[[str]])
    })
  expect_true(all(null_atr))
  
  
  ## Basic norm expectations
  expect_true(all(a_norm$e_data[-1] <=1))
  expect_true(all(a_norm$e_data[-1] >=0))
  
  ## Check zeros are zeros
  min_val <- apply(edata_stand[-1], 2, min, na.rm = T)
  min_val <- purrr::map2_dfc(1:ncol(edata_stand[-1]), min_val, function(col, min){
    edata_stand[-1][col] != min
  })
  min_val[is.na(min_val)] <- FALSE
  min_val <- as.data.frame(min_val)
  a_norm_zeros <- as.matrix(a_norm$e_data[-1] != 0)
  row.names(a_norm_zeros) <- NULL
  
  expect_identical(a_norm_zeros, min_val == T)
  
  ## Check ones are ones
  max_val <- apply(edata_stand[-1], 2, max, na.rm = T)
  max_val <- purrr::map2_dfc(1:ncol(edata_stand[-1]), max_val, function(col, max){
    edata_stand[-1][col] != max
  })
  max_val[is.na(max_val)] <- TRUE
  max_val <- as.data.frame(max_val)
  a_norm_ones <- as.matrix(a_norm$e_data[-1] != 1)
  row.names(a_norm_ones) <- NULL
  
  expect_identical(a_norm_ones, max_val == T)
  
})
