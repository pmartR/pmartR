context('populate the group_DF attribute')

test_that('the correct group data frame and attributes are created',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Fabricate objects and tests for 1 main effect ------------------------------
  
  # Construct a pepData object with the edata, fdata, and emeta data frames.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein',
                      data_scale_orig = "abundance")
  
  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = 'Condition')
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf$e_data),
               c(150, 13))
  expect_equal(dim(pdata_gdf$f_data),
               c(12, 2))
  expect_equal(dim(pdata_gdf$e_meta),
               c(150, 4))
  
  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf, 'group_DF')),
               data.frame(SampleID = as.character(fdata$SampleID),
                          Group = as.character(fdata$Condition),
                          stringsAsFactors = FALSE))
  
  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf, 'group_DF'))$main_effects,
               'Condition')
  expect_equal(attributes(attr(pdata_gdf, 'group_DF'))$nonsingleton_groups,
               c('Infection', 'Mock'))
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$time_course)
  
  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata, 'cnames'),
                   attr(pdata_gdf, 'cnames'))
  expect_identical(attr(pdata, 'data_info'),
                   attr(pdata_gdf, 'data_info'))
  expect_identical(attr(pdata, 'meta_info'),
                   attr(pdata_gdf, 'meta_info'))
  expect_identical(attr(pdata, 'filters'),
                   attr(pdata_gdf, 'filters'))
  
  # Generate objects and tests for 2 main effects ------------------------------
  
  # Add another column to fdata for a second main effect.
  fdata_2 <- data.frame(fdata,
                        Intensity = c('low', 'low', 'high', 'low', 'high',
                                      'high', 'high', 'high', 'low', 'none',
                                      'none', 'none'),
                        stringsAsFactors = FALSE)
  
  # Produce a pepData object with the new f_data data frame.
  pdata_2 <- as.pepData(e_data = edata,
                        f_data = fdata_2,
                        e_meta = emeta,
                        edata_cname = 'Mass_Tag_ID',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'Protein',
                        data_scale_orig = "abundance")
  
  # Run group_designation with two main effects.
  pdata_gdf_2 <- group_designation(omicsData = pdata_2,
                                   main_effects = c('Condition', 'Intensity'))
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf_2$e_data),
               c(150, 13))
  expect_equal(dim(pdata_gdf_2$f_data),
               c(12, 3))
  expect_equal(dim(pdata_gdf_2$e_meta),
               c(150, 4))
  
  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf_2, 'group_DF')),
               data.frame(SampleID = as.character(fdata_2$SampleID),
                          Group = as.character(paste(fdata_2$Condition,
                                                     fdata_2$Intensity,
                                                     sep = '_')),
                          Condition = as.character(fdata_2$Condition),
                          Intensity = as.character(fdata_2$Intensity),
                          stringsAsFactors = FALSE))
  
  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$main_effects,
               c('Condition', 'Intensity'))
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$nonsingleton_groups,
               c('Infection_high', 'Infection_low', 'Mock_none'))
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$time_course)
  
  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata_2, 'cnames'),
                   attr(pdata_gdf_2, 'cnames'))
  expect_identical(attr(pdata_2, 'data_info'),
                   attr(pdata_gdf_2, 'data_info'))
  expect_identical(attr(pdata_2, 'meta_info'),
                   attr(pdata_gdf_2, 'meta_info'))
  expect_identical(attr(pdata_2, 'filters'),
                   attr(pdata_gdf_2, 'filters'))
  
  # Add another column to fdata for a second main effect with some NA values.
  fdata_2 <- data.frame(fdata,
                        Intensity = c('low', 'low', 'high', 'low', 'high',
                                      'high', 'high', 'high', 'low', NA, NA,
                                      NA),
                        stringsAsFactors = FALSE)
  
  # Produce a pepData object with the new f_data data frame.
  pdata_2 <- as.pepData(e_data = edata,
                        f_data = fdata_2,
                        e_meta = emeta,
                        edata_cname = 'Mass_Tag_ID',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'Protein',
                        data_scale_orig = "abundance")
  
  # Run group_designation with two main effects.
  expect_warning(pdata_gdf_2 <- group_designation(omicsData = pdata_2,
                                                  main_effects = c('Condition',
                                                                   'Intensity')),
                 paste('The following 3 sample\\(s\\) has/have been removed',
                       'from the dataset due to missing group information:',
                       'Mock1, Mock2, Mock3',
                       sep = ' '))
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf_2$e_data),
               c(150, 10))
  expect_equal(dim(pdata_gdf_2$f_data),
               c(9, 3))
  expect_equal(dim(pdata_gdf_2$e_meta),
               c(150, 4))
  
  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf_2, 'group_DF')),
               data.frame(SampleID = as.character(fdata_2$SampleID[1:9]),
                          Group = as.character(paste(fdata_2$Condition[1:9],
                                                     fdata_2$Intensity[1:9],
                                                     sep = '_')),
                          Condition = as.character(fdata_2$Condition[1:9]),
                          Intensity = as.character(fdata_2$Intensity[1:9]),
                          stringsAsFactors = FALSE))
  
  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$main_effects,
               c('Condition', 'Intensity'))
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$nonsingleton_groups,
               c('Infection_high', 'Infection_low'))
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$time_course)
  
  # Checkerate the data_info attribute. The elements below should not be the
  # same as the original pepData object because the Mock samples were removed.
  expect_equal(get_data_info(pdata_gdf_2)$num_miss_obs,
               267)
  expect_equal(round(get_data_info(pdata_gdf_2)$prop_missing,
                     4),
               0.1978)
  expect_equal(get_data_info(pdata_gdf_2)$num_samps,
               9)
  
  # The following elements of data_info should be the same as the original
  # pepData object because they were not affected by removing the Mock samples.
  expect_identical(get_data_scale_orig(pdata),
                   get_data_scale_orig(pdata_gdf_2))
  expect_identical(get_data_scale(pdata),
                   get_data_scale(pdata_gdf_2))
  expect_identical(get_data_info(pdata)$norm_info$is_normalized,
                   get_data_info(pdata_gdf_2)$norm_info$is_normalized)
  expect_identical(get_data_info(pdata)$num_edata,
                   get_data_info(pdata_gdf_2)$num_edata)
  expect_identical(get_data_info(pdata)$data_types,
                   get_data_info(pdata_gdf_2)$data_types)
  
})
