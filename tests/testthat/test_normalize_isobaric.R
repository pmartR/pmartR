context('normalize: isobaric')

test_that('normalize_isobaric produces the correct output',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_isodata.RData',
                   package = 'pmartR'))
  
  # Create a isobaricpepData object with the reduced data set.
  isodata <- as.isobaricpepData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = "Peptide",
                                fdata_cname = "Sample",
                                emeta_cname = "Protein")
  
  # Natural logate the data.
  isodata <- edata_transform(omicsData = isodata,
                             data_scale = "log")
  
  # Calculate the normalization standards --------------------------------------
  
  # Subtract the reference values from all other samples within each set.
  set1 <- isodata$e_data[, 2:4] - isodata$e_data[, 5]
  set2 <- isodata$e_data[, 6:8] - isodata$e_data[, 9]
  set3 <- isodata$e_data[, 10:12] - isodata$e_data[, 13]
  
  # Combine each set by column.
  norm_standard <- cbind(isodata$e_data[, 1], set1, set2, set3)
  names(norm_standard)[1] <- names(isodata$e_data)[1]
  
  # Test normalize_isobaric: isobaricnormRes -----------------------------------
  
  # Use the first specification for identifying the reference samples.
  spec1 <- normalize_isobaric(isodata,
                              exp_cname = "Set",
                              apply_norm = FALSE,
                              channel_cname = "iTRAQ.Channel",
                              refpool_channel = "116")
  
  # Use the second specification for identifying the reference samples
  spec2 <- normalize_isobaric(isodata,
                              exp_cname = "Set",
                              apply_norm = FALSE,
                              refpool_cname = "Reference",
                              refpool_notation = "Yes")
  
  # Inspect the class of the two specifications.
  expect_s3_class(spec1, "isobaricnormRes")
  expect_s3_class(spec2, "isobaricnormRes")
  
  # Make sure the two specifications produce the same output.
  expect_identical(spec1$e_data, spec2$e_data)
  # Change all columns in f_data to character vectors because depending on which
  # specification is used some columns will/won't be converted to a character
  # vector within the normalize_isobaric function.
  expect_identical(apply(spec1$f_data, 2, as.character),
                   apply(spec2$f_data, 2, as.character))
  
  # Check e_data and f_data are correct.
  expect_identical(spec1$e_data,
                   isodata$e_data[, c(1, 5, 9, 13)])
  expect_identical(apply(spec1$f_data, 2, as.character),
                   apply(isodata$f_data[c(4, 8, 12), ], 2, as.character))
  expect_identical(spec2$e_data,
                   isodata$e_data[, c(1, 5, 9, 13)])
  expect_identical(apply(spec2$f_data, 2, as.character),
                   apply(isodata$f_data[c(4, 8, 12), ], 2, as.character))
  
  # Scan the attributes of the isobaricnormRes objects for deficiencies.
  expect_identical(attr(spec1, "cnames"),
                   list(edata_cname = "Peptide",
                        fdata_cname = "Sample"))
  expect_identical(attr(spec1, "isobaric_info"),
                   list(exp_cname = "Set",
                        channel_cname = "iTRAQ.Channel",
                        refpool_channel = "116",
                        refpool_cname = NULL,
                        refpool_notation = NULL,
                        norm_info = list(is_normalized = FALSE)))
  expect_identical(attr(spec2, "cnames"),
                   list(edata_cname = "Peptide",
                        fdata_cname = "Sample"))
  expect_identical(attr(spec2, "isobaric_info"),
                   list(exp_cname = "Set",
                        channel_cname = NULL,
                        refpool_channel = NULL,
                        refpool_cname = "Reference",
                        refpool_notation = "Yes",
                        norm_info = list(is_normalized = FALSE)))
  
  # Test normalize_isobaric: apply normalization -------------------------------
  
  # Use the first specification for identifying the reference samples.
  spec1 <- normalize_isobaric(isodata,
                              exp_cname = "Set",
                              apply_norm = TRUE,
                              channel_cname = "iTRAQ.Channel",
                              refpool_channel = "116")
  
  # Use the second specification for identifying the reference samples
  spec2 <- normalize_isobaric(isodata,
                              exp_cname = "Set",
                              apply_norm = TRUE,
                              refpool_cname = "Reference",
                              refpool_notation = "Yes")
  
  # Ensure the normalization is the same between the two specifications.
  expect_identical(spec1$e_data, spec2$e_data)
  # Change all columns in f_data to character vectors because depending on which
  # specification is used some columns will/won't be converted to a character
  # vector within the normalize_isobaric function.
  expect_identical(apply(spec1$f_data, 2, as.character),
                   apply(spec2$f_data, 2, as.character))
  expect_identical(spec1$e_meta, spec2$e_meta)
  
  # Compare the normalized data to the standard.
  expect_identical(spec1$e_data, norm_standard)
  expect_identical(spec2$e_data, norm_standard)
  
  # Change the "Group" column in f_data to "group" because the pre_flight
  # function does not allow any column in f_data to be named "Group".
  names(fdata)[5] <- "group"
  
  # Make sure f_data has the correct rows removed. Metamorphose all columns into
  # character vectors. This way the comparisons are apples to apples instead of
  # something like apples to zombies.
  expect_identical(apply(spec1$f_data, 2, as.character),
                   apply(fdata[-c(4, 8, 12), ], 2, as.character))
  expect_identical(apply(spec2$f_data, 2, as.character),
                   apply(fdata[-c(4, 8, 12), ], 2, as.character))
  
  # Have a looksie at the class.
  expect_s3_class(spec1, c("isobaricpepData", "pepData"))
  expect_s3_class(spec2, c("isobaricpepData", "pepData"))
  
  # See if the attributes that should not have changed did not change.
  expect_identical(attr(isodata, "cnames"),
                   attr(spec1, "cnames"))
  expect_identical(attr(isodata, "cnames"),
                   attr(spec2, "cnames"))
  expect_identical(attr(isodata, "meta_info"),
                   attr(spec1, "meta_info"))
  expect_identical(attr(isodata, "meta_info"),
                   attr(spec2, "meta_info"))
  expect_identical(attr(isodata, "filters"),
                   attr(spec1, "filters"))
  expect_identical(attr(isodata, "filters"),
                   attr(spec2, "filters"))
  
  # Sleuth around the data_info attribute.
  expect_equal(
    attr(spec1, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(spec1$e_data[, 1])),
         num_miss_obs = sum(is.na(spec1$e_data)),
         prop_missing = (sum(is.na(spec1$e_data)) /
                           prod(dim(spec1$e_data[, -1]))),
         num_samps = ncol(spec1$e_data[, -1]),
         data_types = NULL)
  )
  expect_equal(
    attr(spec2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(spec2$e_data[, 1])),
         num_miss_obs = sum(is.na(spec2$e_data)),
         prop_missing = (sum(is.na(spec2$e_data)) /
                           prod(dim(spec2$e_data[, -1]))),
         num_samps = ncol(spec2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Confirm the isobaric_info attributes are correct.
  expect_identical(attr(spec1, "isobaric_info"),
                   list(exp_cname = "Set",
                        channel_cname = "iTRAQ.Channel",
                        refpool_channel = "116",
                        refpool_cname = NULL,
                        refpool_notation = NULL,
                        norm_info = list(is_normalized = TRUE)))
  expect_identical(attr(spec2, "isobaric_info"),
                   list(exp_cname = "Set",
                        channel_cname = NULL,
                        refpool_channel = NULL,
                        refpool_cname = "Reference",
                        refpool_notation = "Yes",
                        norm_info = list(is_normalized = TRUE)))
  
})
