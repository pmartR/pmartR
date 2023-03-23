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
  expect_equal(attr(edr, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = "Protein",
                    fdata_cname = "SampleID"))
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
  expect_identical(class(pdata),
                   class(filtered_edr))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_edr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edr, "filters")[[1]]$filtered$e_data_remove,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_true(is.na(attr(filtered_edr, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_edr, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_edr$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_edr$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_edr$e_data[, -1])) /
                           prod(dim(filtered_edr$e_data[, -1]))),
         num_samps = ncol(filtered_edr$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_edr, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_edr$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_fdr))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdr, "filters")[[1]]$filtered$f_data_remove,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdr, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_fdr, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_fdr$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_fdr$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_fdr$e_data[, -1])) /
                           prod(dim(filtered_fdr$e_data[, -1]))),
         num_samps = ncol(filtered_fdr$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_fdr, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_fdr$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_emr))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_emr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_emr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_emr, "filters")[[1]]$filtered$e_meta_remove,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_emr, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_emr, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_emr$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_emr$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_emr$e_data[, -1])) /
                           prod(dim(filtered_emr$e_data[, -1]))),
         num_samps = ncol(filtered_emr$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_emr, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_emr$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_efdr))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_efdr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efdr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efdr, "filters")[[1]]$filtered$e_data_remove,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efdr, "filters")[[1]]$filtered$f_data_remove,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_efdr, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_efdr, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_efdr$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_efdr$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_efdr$e_data[, -1])) /
                           prod(dim(filtered_efdr$e_data[, -1]))),
         num_samps = ncol(filtered_efdr$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_efdr, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_efdr$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_edemr))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_edemr, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edemr, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edemr, "filters")[[1]]$filtered$e_data_remove,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_edemr, "filters")[[1]]$filtered$e_meta_remove,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_edemr, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_edemr, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_edemr$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_edemr$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_edemr$e_data[, -1])) /
                           prod(dim(filtered_edemr$e_data[, -1]))),
         num_samps = ncol(filtered_edemr$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_edemr, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_edemr$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_efer))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_efer, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efer, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$e_data_remove,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$f_data_remove,
               c("Infection3", "Infection5", "Mock2"))
  expect_equal(attr(filtered_efer, "filters")[[1]]$filtered$e_meta_remove,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_efer, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_efer, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_efer$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_efer$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_efer$e_data[, -1])) /
                           prod(dim(filtered_efer$e_data[, -1]))),
         num_samps = ncol(filtered_efer$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_efer, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_efer$e_meta$Protein)))
  )

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
  expect_identical(class(pdata_grp),
                   class(filtered_fdr_grp))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdr_grp, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdr_grp, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdr_grp, "filters")[[1]]$filtered$f_data_remove,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdr_grp, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_fdr_grp, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_fdr_grp$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_fdr_grp$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_fdr_grp$e_data[, -1])) /
                           prod(dim(filtered_fdr_grp$e_data[, -1]))),
         num_samps = ncol(filtered_fdr_grp$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_fdr_grp, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_fdr_grp$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_edk))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_edk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edk, "filters")[[1]]$filtered$e_data_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_true(is.na(attr(filtered_edk, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_edk, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_edk$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_edk$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_edk$e_data[, -1])) /
                           prod(dim(filtered_edk$e_data[, -1]))),
         num_samps = ncol(filtered_edk$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_edk, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_edk$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_fdk))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdk, "filters")[[1]]$filtered$f_data_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdk, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_fdk, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_fdk$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_fdk$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_fdk$e_data[, -1])) /
                           prod(dim(filtered_fdk$e_data[, -1]))),
         num_samps = ncol(filtered_fdk$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_fdk, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_fdk$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_emk))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_emk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_emk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_emk, "filters")[[1]]$filtered$e_meta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_emk, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_emk, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_emk$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_emk$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_emk$e_data[, -1])) /
                           prod(dim(filtered_emk$e_data[, -1]))),
         num_samps = ncol(filtered_emk$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_emk, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_emk$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_efdk))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_efdk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efdk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efdk, "filters")[[1]]$filtered$e_data_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efdk, "filters")[[1]]$filtered$f_data_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_efdk, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_efdk, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_efdk$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_efdk$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_efdk$e_data[, -1])) /
                           prod(dim(filtered_efdk$e_data[, -1]))),
         num_samps = ncol(filtered_efdk$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_efdk, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_efdk$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_edemk))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_edemk, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_edemk, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_edemk, "filters")[[1]]$filtered$e_data_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_edemk, "filters")[[1]]$filtered$e_meta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_edemk, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_edemk, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_edemk$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_edemk$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_edemk$e_data[, -1])) /
                           prod(dim(filtered_edemk$e_data[, -1]))),
         num_samps = ncol(filtered_edemk$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_edemk, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_edemk$e_meta$Protein)))
  )

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
  expect_identical(class(pdata),
                   class(filtered_efek))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_efek, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_efek, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$e_data_keep,
               c(6948849, 6948892, 6948891, 6948829, 976139, 6948862, 6709059,
                 6948848, 6948884, 6948913, 6706932, 6679059))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$f_data_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_equal(attr(filtered_efek, "filters")[[1]]$filtered$e_meta_keep,
               c("GRP78_HUMAN", "GELS_HUMAN"))
  expect_true(is.na(attr(filtered_efek, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_efek, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_efek$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_efek$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_efek$e_data[, -1])) /
                           prod(dim(filtered_efek$e_data[, -1]))),
         num_samps = ncol(filtered_efek$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_efek, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_efek$e_meta$Protein)))
  )

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
  expect_identical(class(pdata_grp),
                   class(filtered_fdk_grp))

  # Investigate the filters attribute.
  expect_equal(attr(filtered_fdk_grp, "filters")[[1]]$type,
               "customFilt")
  expect_true(is.na(attr(filtered_fdk_grp, "filters")[[1]]$threshold))
  expect_equal(attr(filtered_fdk_grp, "filters")[[1]]$filtered$f_data_keep,
               c("Infection3", "Infection5", "Mock2"))
  expect_true(is.na(attr(filtered_fdk_grp, "filters")[[1]]$method))

  # Examinate the data_info attribute.
  expect_equal(
    attr(filtered_fdk_grp, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_fdk_grp$e_data$Mass_Tag_ID)),
         num_miss_obs = sum(is.na(filtered_fdk_grp$e_data[, -1])),
         prop_missing = (sum(is.na(filtered_fdk_grp$e_data[, -1])) /
                           prod(dim(filtered_fdk_grp$e_data[, -1]))),
         num_samps = ncol(filtered_fdk_grp$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_fdk_grp, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_fdk_grp$e_meta$Protein)))
  )

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

  # Create paired objects ------------------------------------------------------

  load(system.file('testdata',
                   'little_pairdata.RData',
                   package = 'pmartR'))

  # Create a pepData object with the original main effect and pairing variable.
  pairdata <- as.pepData(e_data = edata,
                         f_data = fdata,
                         e_meta = emeta,
                         edata_cname = 'Mass_Tag_ID',
                         fdata_cname = 'Name',
                         emeta_cname = 'Protein')
  pairdata <- edata_transform(pairdata,
                              data_scale = "log")
  pairdata <- group_designation(pairdata,
                                main_effects = "Virus",
                                pair_id = "PairID",
                                pair_group = "Time",
                                pair_denom = "0")

  split_remove <- custom_filter(
    pairdata,
    f_data_remove = c("FM_0hr_3", "AM_18hr_5")
  )

  split_keep <- custom_filter(
    pairdata,
    f_data_keep = fdata$Name[c(1:10, 11:19, 21:30)]
  )

  # Test filter on paired data -------------------------------------------------

  expect_error(
    applyFilt(filter_object = split_remove, omicsData = pairdata),
    paste("The following samples should also be removed based on the input:",
    "FM_18hr_3 and AM_0hr_5. Samples in a pair must be removed together.",
    sep = " ")
  )

  expect_error(
    applyFilt(filter_object = split_keep, omicsData = pairdata),
    paste("The following samples should also be kept based on the input:",
          "FM_18hr_5. Samples in a pair must be kept together.",
          sep = " ")
  )
  
  # Expect warning if columns/rows have already been filtered ------------------
  
  # Columns
  
  fdr2_1 <- custom_filter(pdata,
                          f_data_remove = c("Infection3"))
  fdr2_2 <- custom_filter(pdata,
                          f_data_remove = c("Infection3"))
  fdr2_3 <- custom_filter(pdata,
                          f_data_remove = c("Infection3", "Mock2"))
  
  # Data should be identical if filtered columns don't exist
  filtered_fdr2_1 <- applyFilt(fdr2_1, pdata) 
  expect_warning(filtered_fdr2_2 <- applyFilt(fdr2_2, filtered_fdr2_1),
                 "Specified samples Infection3 were not found in the data")
  expect_identical(filtered_fdr2_1$e_data,
                   filtered_fdr2_2$e_data)
  
  # Columns that do exist should be removed, even if there are columns that
  # don't exist
  expect_warning(filtered_fdr2_3 <- applyFilt(fdr2_3, filtered_fdr2_1))
  expect_null(filtered_fdr2_3$e_data$Mock2)
  
  # Rows
  
  edr2_1 <- custom_filter(pdata,
                          e_data_remove = 6948849)
  edr2_2 <- custom_filter(pdata,
                          e_data_remove = 6948849)
  edr2_3 <- custom_filter(pdata,
                          e_data_remove = c(6948849, 6679059))
  
  # Data should be identical if filtered rows don't exist
  filtered_edr2_1 <- applyFilt(edr2_1, pdata)
  expect_warning(filtered_edr2_2 <- applyFilt(edr2_2, filtered_edr2_1),
                 "Specified biomolecules 6948849 were not found in e_data")
  expect_identical(filtered_edr2_1$e_data,
                   filtered_edr2_2$e_data)
  
  # Rows that do exist should be removed, even if there are rows that
  # don't exist
  expect_warning(filtered_edr2_3 <- applyFilt(edr2_3, filtered_edr2_1),
                 "Specified biomolecules 6948849 were not found in e_data")
  expect_false(any(filtered_edr2_3$e_data$Mass_Tag_ID == 6679059))
})
