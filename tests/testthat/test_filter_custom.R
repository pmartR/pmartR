context('filter by custom criteria')

test_that('custom_filter and applyFilt produce the correct output',{
  
  # Load data and standards and prepare omicsData objects ----------------------
  
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
  
  # Test custom_filter preliminary checks --------------------------------------
  
  expect_error(custom_filter(pdata),
               "No items have been identified for filtering.")
  
  expect_error(custom_filter(pdata,
                             f_data_remove = c("Infection3", "Infection5",
                                               "Mock2"),
                             e_meta_keep = c("GRP78_HUMAN", "GELS_HUMAN")),
               paste("Cannot have both remove arguments and keep arguments",
                     "be non-NULL. Create separate filter objects for the",
                     "remove arguments and keep arguments.",
                     sep = " "))
  
  # Test custom_filter remove --------------------------------------------------
  
  # e_data ---------------
  
  # Construct a remove filter object for e_data.
  edr <- custom_filter(pdata,
                       e_data_remove = c(6948849, 6948892, 6948891, 6948829,
                                         976139, 6948862, 6709059, 6948848,
                                         6948884, 6948913, 6706932, 6679059))
  
  # Examine the class of the filter object.
  expect_s3_class(edr,
                  c("customFilt", "list"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(edr$e_data_remove,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  
  # Make sure the other elements are NULL.
  expect_null(edr$f_data_remove)
  expect_null(edr$e_meta_remove)
  
  # Ensure custom filter attributes are correct.
  expect_equal(attr(edr, "num_samples"), 12)
  expect_equal(attr(edr, "num_edata"), 150)
  expect_equal(attr(edr, "num_emeta"), 83)
  expect_equal(attr(edr, "cnames")$edata_cname,
               "Mass_Tag_ID")
  expect_equal(attr(edr, "cnames")$fdata_cname,
               "SampleID")
  expect_equal(attr(edr, "cnames")$emeta_cname,
               "Protein")
  expect_equal(dim(attr(edr, "omicsData")$e_data),
               c(150, 13))
  expect_equal(dim(attr(edr, "omicsData")$f_data),
               c(12, 2))
  expect_equal(dim(attr(edr, "omicsData")$e_meta),
               c(150, 4))
  
  # f_data ---------------
  
  # Manufacture a remove filter object for f_data.
  fdr <- custom_filter(pdata,
                       f_data_remove = c("Infection3", "Infection5", "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(fdr$f_data_remove,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(fdr$e_data_remove)
  expect_null(fdr$e_meta_remove)
  
  # e_meta ---------------
  
  # Fashion a remove filter object for e_meta.
  emr <- custom_filter(pdata,
                       e_meta_remove = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(emr$e_meta_remove,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Make sure the other elements are NULL.
  expect_null(emr$e_data_remove)
  expect_null(emr$f_data_remove)
  
  # e_data and f_data ---------------
  
  # Produce a remove filter with e_data and f_data filter criteria.
  efdr <- custom_filter(pdata,
                        e_data_remove = c(6948849, 6948892, 6948891, 6948829,
                                          976139, 6948862, 6709059, 6948848,
                                          6948884, 6948913, 6706932, 6679059),
                        f_data_remove = c("Infection3", "Infection5", "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(efdr$e_data_remove,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(efdr$f_data_remove,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(efdr$e_meta_remove)
  
  # e_data and e_meta ---------------
  
  # Fabricate a remove filter with e_data and e_meta filter criteria.
  edemr <- custom_filter(pdata,
                         e_data_remove = c(6948849, 6948892, 6948891, 6948829,
                                           976139, 6948862, 6709059, 6948848,
                                           6948884, 6948913, 6706932, 6679059),
                         e_meta_remove = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(edemr$e_data_remove,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(edemr$e_meta_remove,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Make sure the other elements are NULL.
  expect_null(edemr$f_data_remove)
  
  # e_data, f_data, and e_meta ---------------
  
  # Design a remove filter with using all three filter criteria.
  efer <- custom_filter(pdata,
                        e_data_remove = c(6948849, 6948892, 6948891, 6948829,
                                          976139, 6948862, 6709059, 6948848,
                                          6948884, 6948913, 6706932, 6679059),
                        f_data_remove = c("Infection3", "Infection5", "Mock2"),
                        e_meta_remove = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(efer$e_data_remove,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(efer$f_data_remove,
                   c("Infection3", "Infection5", "Mock2"))
  expect_identical(efer$e_meta_remove,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Test custom_filter remove with groups --------------------------------------
  
  # Generate a pepData object with group information.
  pdata_grp <- group_designation(omicsData = pdata,
                                 main_effects = "Condition")
  
  # Manufacture a remove filter object for f_data.
  fdr_grp <- custom_filter(pdata_grp,
                           f_data_remove = c("Infection3", "Infection5",
                                             "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(fdr_grp$f_data_remove,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(fdr_grp$e_data_remove)
  expect_null(fdr_grp$e_meta_remove)
  
  # Test custom_filter keep ----------------------------------------------------
  
  # e_data ---------------
  
  # Construct a keep filter object for e_data.
  edk <- custom_filter(pdata,
                       e_data_keep = c(6948849, 6948892, 6948891, 6948829,
                                       976139, 6948862, 6709059, 6948848,
                                       6948884, 6948913, 6706932, 6679059))
  
  # Check the correct IDs are on the filter list.
  expect_identical(edk$e_data_keep,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  
  # Make sure the other elements are NULL.
  expect_null(edk$f_data_keep)
  expect_null(edk$e_meta_keep)
  
  # f_data ---------------
  
  # Manufacture a keep filter object for f_data.
  fdk <- custom_filter(pdata,
                       f_data_keep = c("Infection3", "Infection5", "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(fdk$f_data_keep,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(fdk$e_data_keep)
  expect_null(fdk$e_meta_keep)
  
  # e_meta ---------------
  
  # Fashion a keep filter object for e_meta.
  emk <- custom_filter(pdata,
                       e_meta_keep = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(emk$e_meta_keep,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Make sure the other elements are NULL.
  expect_null(emk$e_data_keep)
  expect_null(emk$f_data_keep)
  
  # e_data and f_data ---------------
  
  # Produce a keep filter with e_data and f_data filter criteria.
  efdk <- custom_filter(pdata,
                        e_data_keep = c(6948849, 6948892, 6948891, 6948829,
                                        976139, 6948862, 6709059, 6948848,
                                        6948884, 6948913, 6706932, 6679059),
                        f_data_keep = c("Infection3", "Infection5", "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(efdk$e_data_keep,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(efdk$f_data_keep,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(efdk$e_meta_keep)
  
  # e_data and e_meta ---------------
  
  # Fabricate a keep filter with e_data and e_meta filter criteria.
  edemk <- custom_filter(pdata,
                         e_data_keep = c(6948849, 6948892, 6948891, 6948829,
                                         976139, 6948862, 6709059, 6948848,
                                         6948884, 6948913, 6706932, 6679059),
                         e_meta_keep = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(edemk$e_data_keep,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(edemk$e_meta_keep,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Make sure the other elements are NULL.
  expect_null(edemk$f_data_keep)
  
  # e_data, f_data, and e_meta ---------------
  
  # Design a keep filter with using all three filter criteria.
  efek <- custom_filter(pdata,
                        e_data_keep = c(6948849, 6948892, 6948891, 6948829,
                                        976139, 6948862, 6709059, 6948848,
                                        6948884, 6948913, 6706932, 6679059),
                        f_data_keep = c("Infection3", "Infection5", "Mock2"),
                        e_meta_keep = c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(efek$e_data_keep,
                   c(6948849, 6948892, 6948891, 6948829, 976139, 6948862,
                     6709059, 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_identical(efek$f_data_keep,
                   c("Infection3", "Infection5", "Mock2"))
  expect_identical(efek$e_meta_keep,
                   c("GRP78_HUMAN", "GELS_HUMAN"))
  
  # Test custom_filter keep with groups ----------------------------------------
  
  # Manufacture a keep filter object for f_data.
  fdk_grp <- custom_filter(pdata_grp,
                           f_data_keep = c("Infection3", "Infection5", "Mock2"))
  
  # Check the correct IDs are on the filter list.
  expect_identical(fdk_grp$f_data_keep,
                   c("Infection3", "Infection5", "Mock2"))
  
  # Make sure the other elements are NULL.
  expect_null(fdk_grp$e_data_keep)
  expect_null(fdk_grp$e_meta_keep)
  
  # Test applyFilt remove no groups --------------------------------------------
  
  # e_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_edr <- applyFilt(filter_object = edr,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_edr, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_edr, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_edr))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_edr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edr, "filters")[[1]]$filtered$edata_filt,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_true(is.na(attr(filtered_edr, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_edr, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_edr, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_edr, "data_info")$num_edata,
               138)
  expect_equal(attr(filtered_edr, "data_info")$num_miss_obs,
               322)
  expect_equal(round(attr(filtered_edr, "data_info")$prop_missing, 4),
               0.1944)
  expect_equal(attr(filtered_edr, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_edr, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_edr, "meta_info")$meta_data)
  expect_equal(attr(filtered_edr, "meta_info")$num_emeta,
               81)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_edr$e_data),
               c(138, 13))
  expect_equal(dim(filtered_edr$f_data),
               c(12, 2))
  expect_equal(dim(filtered_edr$e_meta),
               c(138, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_edr, "group_DF"))
  
  # f_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_fdr <- applyFilt(filter_object = fdr,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_fdr, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_fdr, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_fdr))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdr, "filters")[[1]]$filtered$samples_filt,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdr, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_fdr, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_fdr, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_fdr, "data_info")$num_edata,
               150)
  expect_equal(attr(filtered_fdr, "data_info")$num_miss_obs,
               262)
  expect_equal(round(attr(filtered_fdr, "data_info")$prop_missing, 4),
               0.1941)
  expect_equal(attr(filtered_fdr, "data_info")$num_samps,
               9)
  expect_null(attr(filtered_fdr, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_fdr, "meta_info")$meta_data)
  expect_equal(attr(filtered_fdr, "meta_info")$num_emeta,
               83)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_fdr$e_data),
               c(150, 10))
  expect_equal(dim(filtered_fdr$f_data),
               c(9, 2))
  expect_equal(dim(filtered_fdr$e_meta),
               c(150, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_fdr, "group_DF"))
  
  # e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_emr <- applyFilt(filter_object = emr,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_emr, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_emr, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_emr))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_emr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_emr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_emr, "filters")[[1]]$filtered$emeta_filt,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_emr, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_emr, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_emr, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_emr, "data_info")$num_edata,
               141)
  expect_equal(attr(filtered_emr, "data_info")$num_miss_obs,
               323)
  expect_equal(round(attr(filtered_emr, "data_info")$prop_missing, 4),
               0.1909)
  expect_equal(attr(filtered_emr, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_emr, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_emr, "meta_info")$meta_data)
  expect_equal(attr(filtered_emr, "meta_info")$num_emeta,
               81)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_emr$e_data),
               c(141, 13))
  expect_equal(dim(filtered_emr$f_data),
               c(12, 2))
  expect_equal(dim(filtered_emr$e_meta),
               c(141, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_emr, "group_DF"))
  
  # e_data and f_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_efdr <- applyFilt(filter_object = efdr,
                             omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_efdr, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_efdr, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_efdr))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_efdr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efdr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efdr, "filters")[[1]]$filtered$edata_filt,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efdr, "filters")[[1]]$filtered$samples_filt,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_efdr, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_efdr, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_efdr, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_efdr, "data_info")$num_edata,
               138)
  expect_equal(attr(filtered_efdr, "data_info")$num_miss_obs,
               248)
  expect_equal(round(attr(filtered_efdr, "data_info")$prop_missing, 4),
               0.1997)
  expect_equal(attr(filtered_efdr, "data_info")$num_samps,
               9)
  expect_null(attr(filtered_efdr, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_efdr, "meta_info")$meta_data)
  expect_equal(attr(filtered_efdr, "meta_info")$num_emeta,
               81)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_efdr$e_data),
               c(138, 10))
  expect_equal(dim(filtered_efdr$f_data),
               c(9, 2))
  expect_equal(dim(filtered_efdr$e_meta),
               c(138, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_efdr, "group_DF"))
  
  # e_data and e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_edemr <- applyFilt(filter_object = edemr,
                              omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_edemr, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_edemr, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_edemr))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_edemr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edemr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edemr, "filters")[[1]]$filtered$edata_filt,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_edemr, "filters")[[1]]$filtered$emeta_filt,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_edemr, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_edemr, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_edemr, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_edemr, "data_info")$num_edata,
               130)
  expect_equal(attr(filtered_edemr, "data_info")$num_miss_obs,
               304)
  expect_equal(round(attr(filtered_edemr, "data_info")$prop_missing, 4),
               0.1949)
  expect_equal(attr(filtered_edemr, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_edemr, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_edemr, "meta_info")$meta_data)
  expect_equal(attr(filtered_edemr, "meta_info")$num_emeta,
               79)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_edemr$e_data),
               c(130, 13))
  expect_equal(dim(filtered_edemr$f_data),
               c(12, 2))
  expect_equal(dim(filtered_edemr$e_meta),
               c(130, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_edemr, "group_DF"))
  
  # e_data, f_data, and e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_efer <- applyFilt(filter_object = efer,
                             omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_efer, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_efer, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_efer))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_efer, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efer, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$edata_filt,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$samples_filt,
               c("Infection3", "Infection5", "Mock2"))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$emeta_filt,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_efer, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_efer, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_efer, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_efer, "data_info")$num_edata,
               130)
  expect_equal(attr(filtered_efer, "data_info")$num_miss_obs,
               234)
  expect_equal(round(attr(filtered_efer, "data_info")$prop_missing, 4),
               0.2)
  expect_equal(attr(filtered_efer, "data_info")$num_samps,
               9)
  expect_null(attr(filtered_efer, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_efer, "meta_info")$meta_data)
  expect_equal(attr(filtered_efer, "meta_info")$num_emeta,
               79)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_efer$e_data),
               c(130, 10))
  expect_equal(dim(filtered_efer$f_data),
               c(9, 2))
  expect_equal(dim(filtered_efer$e_meta),
               c(130, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_efer, "group_DF"))
  
  # Test applyFilt remove with groups ------------------------------------------
  
  # Apply the filter to the reduced pepData object.
  filtered_fdr_grp <- applyFilt(filter_object = fdr_grp,
                                omicsData = pdata_grp)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_grp, "cnames"),
                   attr(filtered_fdr_grp, "cnames"))
  expect_identical(attr(pdata_grp, "check.names"),
                   attr(filtered_fdr_grp, "check.names"))
  expect_identical(class(pdata_grp),
                   class(filtered_fdr_grp))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdr_grp, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdr_grp, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdr_grp, "filters")[[1]]$filtered$samples_filt,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdr_grp, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_fdr_grp, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_fdr_grp, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_fdr_grp, "data_info")$num_edata,
               150)
  expect_equal(attr(filtered_fdr_grp, "data_info")$num_miss_obs,
               262)
  expect_equal(round(attr(filtered_fdr_grp, "data_info")$prop_missing, 4),
               0.1941)
  expect_equal(attr(filtered_fdr_grp, "data_info")$num_samps,
               9)
  expect_null(attr(filtered_fdr_grp, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_fdr_grp, "meta_info")$meta_data)
  expect_equal(attr(filtered_fdr_grp, "meta_info")$num_emeta,
               83)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_fdr_grp$e_data),
               c(150, 10))
  expect_equal(dim(filtered_fdr_grp$f_data),
               c(9, 2))
  expect_equal(dim(filtered_fdr_grp$e_meta),
               c(150, 4))
  
  # Checkerate the group_DF attribute.
  expect_failure(expect_equal(attr(pdata_grp, "group_DF"),
                              attr(filtered_fdr_grp, "group_DF")))
  
  expect_equal(data.frame(attr(pdata_grp, "group_DF")[-c(3, 5, 11), ],
                          row.names = NULL),
               data.frame(attr(filtered_fdr_grp, "group_DF"),
                          row.names = NULL))
  
  # Test applyFilt keep no groups ----------------------------------------------
  
  # e_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_edk <- applyFilt(filter_object = edk,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_edk, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_edk, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_edk))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_edk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edk, "filters")[[1]]$filtered$edata_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_true(is.na(attr(filtered_edk, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_edk, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_edk, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_edk, "data_info")$num_edata,
               12)
  expect_equal(attr(filtered_edk, "data_info")$num_miss_obs,
               19)
  expect_equal(round(attr(filtered_edk, "data_info")$prop_missing, 4),
               0.1319)
  expect_equal(attr(filtered_edk, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_edk, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_edk, "meta_info")$meta_data)
  expect_equal(attr(filtered_edk, "meta_info")$num_emeta,
               10)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_edk$e_data),
               c(12, 13))
  expect_equal(dim(filtered_edk$f_data),
               c(12, 2))
  expect_equal(dim(filtered_edk$e_meta),
               c(12, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_edk, "group_DF"))
  
  # f_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_fdk <- applyFilt(filter_object = fdk,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_fdk, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_fdk, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_fdk))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdk, "filters")[[1]]$filtered$samples_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdk, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_fdk, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_fdk, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_fdk, "data_info")$num_edata,
               150)
  expect_equal(attr(filtered_fdk, "data_info")$num_miss_obs,
               79)
  expect_equal(round(attr(filtered_fdk, "data_info")$prop_missing, 4),
               0.1756)
  expect_equal(attr(filtered_fdk, "data_info")$num_samps,
               3)
  expect_null(attr(filtered_fdk, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_fdk, "meta_info")$meta_data)
  expect_equal(attr(filtered_fdk, "meta_info")$num_emeta,
               83)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_fdk$e_data),
               c(150, 4))
  expect_equal(dim(filtered_fdk$f_data),
               c(3, 2))
  expect_equal(dim(filtered_fdk$e_meta),
               c(150, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_fdk, "group_DF"))
  
  # e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_emk <- applyFilt(filter_object = emk,
                            omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_emk, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_emk, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_emk))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_emk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_emk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_emk, "filters")[[1]]$filtered$emeta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_emk, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_emk, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_emk, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_emk, "data_info")$num_edata,
               9)
  expect_equal(attr(filtered_emk, "data_info")$num_miss_obs,
               18)
  expect_equal(round(attr(filtered_emk, "data_info")$prop_missing, 4),
               0.1667)
  expect_equal(attr(filtered_emk, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_emk, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_emk, "meta_info")$meta_data)
  expect_equal(attr(filtered_emk, "meta_info")$num_emeta,
               2)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_emk$e_data),
               c(9, 13))
  expect_equal(dim(filtered_emk$f_data),
               c(12, 2))
  expect_equal(dim(filtered_emk$e_meta),
               c(9, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_emk, "group_DF"))
  
  # e_data and f_data ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_efdk <- applyFilt(filter_object = efdk,
                             omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_efdk, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_efdk, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_efdk))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_efdk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efdk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efdk, "filters")[[1]]$filtered$edata_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efdk, "filters")[[1]]$filtered$samples_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_efdk, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_efdk, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_efdk, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_efdk, "data_info")$num_edata,
               12)
  expect_equal(attr(filtered_efdk, "data_info")$num_miss_obs,
               5)
  expect_equal(round(attr(filtered_efdk, "data_info")$prop_missing, 4),
               0.1389)
  expect_equal(attr(filtered_efdk, "data_info")$num_samps,
               3)
  expect_null(attr(filtered_efdk, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_efdk, "meta_info")$meta_data)
  expect_equal(attr(filtered_efdk, "meta_info")$num_emeta,
               10)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_efdk$e_data),
               c(12, 4))
  expect_equal(dim(filtered_efdk$f_data),
               c(3, 2))
  expect_equal(dim(filtered_efdk$e_meta),
               c(12, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_efdk, "group_DF"))
  
  # e_data and e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_edemk <- applyFilt(filter_object = edemk,
                              omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_edemk, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_edemk, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_edemk))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_edemk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edemk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edemk, "filters")[[1]]$filtered$edata_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_edemk, "filters")[[1]]$filtered$emeta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_edemk, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_edemk, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_edemk, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_edemk, "data_info")$num_edata,
               12)
  expect_equal(attr(filtered_edemk, "data_info")$num_miss_obs,
               19)
  expect_equal(round(attr(filtered_edemk, "data_info")$prop_missing, 4),
               0.1319)
  expect_equal(attr(filtered_edemk, "data_info")$num_samps,
               12)
  expect_null(attr(filtered_edemk, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_edemk, "meta_info")$meta_data)
  expect_equal(attr(filtered_edemk, "meta_info")$num_emeta,
               2)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_edemk$e_data),
               c(12, 13))
  expect_equal(dim(filtered_edemk$f_data),
               c(12, 2))
  expect_equal(dim(filtered_edemk$e_meta),
               c(9, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_edemk, "group_DF"))
  
  # e_data, f_data, and e_meta ---------------
  
  # Apply the filter to the reduced pepData object.
  filtered_efek <- applyFilt(filter_object = efek,
                             omicsData = pdata)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(filtered_efek, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(filtered_efek, "check.names"))
  expect_identical(class(pdata),
                   class(filtered_efek))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_efek, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efek, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$edata_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$samples_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$emeta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_efek, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_efek, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_efek, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_efek, "data_info")$num_edata,
               12)
  expect_equal(attr(filtered_efek, "data_info")$num_miss_obs,
               5)
  expect_equal(round(attr(filtered_efek, "data_info")$prop_missing, 4),
               0.1389)
  expect_equal(attr(filtered_efek, "data_info")$num_samps,
               3)
  expect_null(attr(filtered_efek, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_efek, "meta_info")$meta_data)
  expect_equal(attr(filtered_efek, "meta_info")$num_emeta,
               2)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_efek$e_data),
               c(12, 4))
  expect_equal(dim(filtered_efek$f_data),
               c(3, 2))
  expect_equal(dim(filtered_efek$e_meta),
               c(9, 4))
  
  # Checkerate the group_DF attribute.
  expect_null(attr(filtered_efek, "group_DF"))
  
  # Test applyFilt keep with groups --------------------------------------------
  
  # Apply the filter to the reduced pepData object.
  filtered_fdk_grp <- applyFilt(filter_object = fdk_grp,
                                omicsData = pdata_grp)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_grp, "cnames"),
                   attr(filtered_fdk_grp, "cnames"))
  expect_identical(attr(pdata_grp, "check.names"),
                   attr(filtered_fdk_grp, "check.names"))
  expect_identical(class(pdata_grp),
                   class(filtered_fdk_grp))
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdk_grp, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdk_grp, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdk_grp, "filters")[[1]]$filtered$samples_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdk_grp, "filters")[[1]]$method))
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_fdk_grp, "data_info")$data_scale,
               "abundance")
  expect_false(attr(filtered_fdk_grp, "data_info")$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_fdk_grp, "data_info")$num_edata,
               150)
  expect_equal(attr(filtered_fdk_grp, "data_info")$num_miss_obs,
               79)
  expect_equal(round(attr(filtered_fdk_grp, "data_info")$prop_missing, 4),
               0.1756)
  expect_equal(attr(filtered_fdk_grp, "data_info")$num_samps,
               3)
  expect_null(attr(filtered_fdk_grp, "data_info")$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(filtered_fdk_grp, "meta_info")$meta_data)
  expect_equal(attr(filtered_fdk_grp, "meta_info")$num_emeta,
               83)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_fdk_grp$e_data),
               c(150, 4))
  expect_equal(dim(filtered_fdk_grp$f_data),
               c(3, 2))
  expect_equal(dim(filtered_fdk_grp$e_meta),
               c(150, 4))
  
  # Checkerate the group_DF attribute.
  expect_failure(expect_equal(attr(pdata_grp, "group_DF"),
                              attr(filtered_fdk_grp, "group_DF")))
  
  expect_equal(data.frame(attr(pdata_grp, "group_DF")[c(3, 5, 11), ],
                          row.names = NULL),
               data.frame(attr(filtered_fdk_grp, "group_DF"),
                          row.names = NULL))
  
})
