context("filter by robust Mahalanobis distance")

test_that("rmd_filter and applyFilt produce the correct output",{

  # Load and prepare omicsData objects -----------------------------------------

  # Load peptide data.
  load(system.file("testdata",
                   "little_pdata.RData",
                   package = "pmartR"))

  # Fabricate a pepData object with the edata, fdata, and emeta data frames.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")

  # Log transfigure the peptide data.
  pdata <- edata_transform(omicsData = pdata,
                           data_scale = "log")

  # Group designate the peptide data.
  pdata <- group_designation(omicsData = pdata,
                             main_effects = "Condition")

  # Forge a pepData object with only one sample from Mock. This will create a
  # singleton group when group designating the data. sg: singleton group.
  pdata_sg <- as.pepData(e_data = edata[, 1:11],
                         f_data = fdata[1:10, ],
                         e_meta = emeta,
                         edata_cname = "Mass_Tag_ID",
                         fdata_cname = "SampleID",
                         emeta_cname = "Protein")

  # Log transmogrify the peptide data with a singleton group.
  pdata_sg <- edata_transform(omicsData = pdata_sg,
                              data_scale = "log")

  # Group designate the peptide data with a singleton group.
  pdata_sg <- group_designation(omicsData = pdata_sg,
                                main_effects = "Condition")

  # Load nmr data.
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))

  # Produce a nmrData object with the edata, fdata, and emeta data frames.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'nmrClass')

  # Log transmute the nmr data.
  nmrdata <- edata_transform(omicsData = nmrdata,
                             data_scale = "log")

  # Group designate the nmr data.
  nmrdata <- group_designation(omicsData = nmrdata,
                               main_effects = "Condition")

  # Test rmd_filter no singleton groups ----------------------------------------

  # Forge a filter object for the peptide data.
  pfilter_rmd <- rmd_filter(omicsData = pdata,
                            ignore_singleton_groups = FALSE)

  # Check the class of the filter object.
  expect_s3_class(pfilter_rmd,
                  c("rmdFilt", "data.frame"))

  # Checkify the peptide filter object data frame.
  expect_equal(dim(pfilter_rmd),
               c(12, 9))
  expect_equal(round(pfilter_rmd$Log2.md, 3),
               c(1.754, 1.983, 1.340, 2.158, 0.635, 1.697,
                 4.748, 2.210, 2.131, 0.787, 0.374, 1.477))
  expect_equal(round(pfilter_rmd$pvalue, 3),
               c(0.643, 0.556, 0.772, 0.485, 0.907, 0.663,
                 0.000, 0.463, 0.496, 0.886, 0.935, 0.733))
  expect_equal(round(pfilter_rmd$MAD, 3),
               c(1.331, 1.185, 1.241, 1.228, 1.228, 1.324,
                 1.203, 1.081, 1.169, 1.142, 1.274, 1.245))
  expect_equal(round(pfilter_rmd$Kurtosis, 3),
               c(1.250, 0.962, 1.085, -0.016, 0.523, -0.274,
                 -0.108, 0.174, -0.084, 0.286, 0.300, 0.605))
  expect_equal(round(pfilter_rmd$Skewness, 3),
               c(-0.429, -0.269, -0.296, 0.170, -0.045, 0.282,
                 -0.006, 0.286, 0.355, -0.143, -0.064, -0.148))
  expect_equal(round(pfilter_rmd$Corr, 3),
               c(0.971, 0.975, 0.971, 0.975, 0.974, 0.966,
                 0.960, 0.941, 0.966, 0.974, 0.982, 0.985))
  expect_equal(round(pfilter_rmd$Proportion_Missing, 3),
               c(0.173, 0.180, 0.180, 0.207, 0.187, 0.220,
                 0.153, 0.260, 0.220, 0.167, 0.160, 0.167))

  # Inspectify the filter object attributes.
  expect_equal(attr(pfilter_rmd, "sample_names"),
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection7", "Infection8",
                 "Infection9", "Mock1", "Mock2", "Mock3"))
  expect_equal(attr(pfilter_rmd, "group_DF"),
               attr(pdata, "group_DF"))
  expect_equal(attr(pfilter_rmd, "df"),
               5)
  expect_equal(attr(pfilter_rmd, "metrics"),
               c("MAD", "Kurtosis", "Skewness", "Corr", "Proportion_Missing"))

  # Test for an error with too few arguments in metrics (because there is no
  # missing data).
  expect_error(suppressWarnings(rmd_filter(omicsData = nmrdata,
                                           metrics = c("MAD",
                                                       "Proportion_Missing"))),
               paste("Vector of metrics must contain at least two valid",
                     "entries. Try including a metric other than",
                     "Proportion_Missing.",
                     sep = " "))

  # Test for a warning because there are few observations.
  expect_warning(rmd_filter(omicsData = nmrdata,
                            metrics = c("MAD",
                                        "Kurtosis",
                                        "Proportion_Missing")),
                 paste("Use the results of the RMD filter with caution due",
                       "to a small number of biomolecules \\(<50\\).",
                       sep = " "))

  # Test for a warning because there is no missing data.
  expect_warning(rmd_filter(omicsData = nmrdata,
                            metrics = c("MAD",
                                        "Kurtosis",
                                        "Proportion_Missing")),
                 paste("There are no missing values in e_data, therefore",
                       "Proportion_Missing will not be used as one of the",
                       "metrics for RMD-Runs.",
                       sep = " "))

  # Generate a filter object for the nmr data.
  nmrfilter_rmd <- suppressWarnings(rmd_filter(omicsData = nmrdata,
                                               metrics = c("MAD", "Kurtosis",
                                                           "Skewness",
                                                           "Correlation")))

  # Check the class of the filter object.
  expect_s3_class(nmrfilter_rmd,
                  c("rmdFilt", "data.frame"))

  # Checkify the nmr filter object data frame.
  expect_equal(dim(nmrfilter_rmd),
               c(41, 8))
  expect_equal(round(nmrfilter_rmd$Log2.md, 3),

               c(3.913, 2.919, 0.785, 3.937, 5.23, 0.897, 1.553, 1.644, 1.447,
                 2.126, 1.36, 1.995, 1.256, 0.095, 1.066, 0.999, 2.728, 1.966,
                 1.471, 3.022, 2.122, -1.172, 4.629, 2.262, 5.708, 3.163, 4.236,
                 1.666, 2.165, 1.107, 2.833, 2.18, 0.1, 0.131, 1.718, 1.712,
                 1.757, 3.113, 1.832, 2.577, 2.01))
  expect_equal(round(nmrfilter_rmd$pvalue, 3),
               c(0.005, 0.109, 0.786, 0.004, 0, 0.761, 0.569, 0.537, 0.605,
                 0.359, 0.633, 0.408, 0.665, 0.899, 0.719, 0.736, 0.157, 0.419,
                 0.597, 0.087, 0.36, 0.979, 0, 0.309, 0, 0.062, 0.001, 0.529,
                 0.344, 0.708, 0.13, 0.339, 0.899, 0.895, 0.511, 0.513, 0.496,
                 0.07, 0.469, 0.202, 0.402))
  expect_equal(round(nmrfilter_rmd$MAD, 3),
               c(0.922, 0.978, 1.411, 0.951, 1.047, 1.173, 1.245, 1.271, 1.176,
                 1.365, 1.234, 1.26, 1.27, 1.463, 1.461, 1.39, 1.328, 1.435,
                 1.632, 1.264, 1.65, 1.286, 1.018, 1.04, 0.788, 0.997, 0.967,
                 1.12, 1.197, 1.404, 1.619, 1.486, 1.409, 1.409, 1.413, 1.403,
                 1.654, 1.503, 1.529, 1.634, 1.502))
  expect_equal(round(nmrfilter_rmd$Kurtosis, 3),
               c(0.004, -0.022, -0.268, -0.114, -0.297, -0.094, -0.168, -0.371,
                 -0.082, -0.09, -0.105, -0.216, -0.372, -0.369, -0.529, -0.418,
                 -0.619, -0.376, -0.642, -0.752, -0.774, -0.378, -0.206, -0.374,
                 -0.173, -0.324, -0.283, -0.087, -0.282, -0.38, -0.58, -0.591,
                 -0.332, -0.43, -0.509, -0.194, -0.506, -0.578, -0.478, -0.652,
                 -0.374))
  expect_equal(round(nmrfilter_rmd$Skewness, 3),
               c(0.828, 0.874, 0.638, 0.73, 0.517, 0.75, 0.713, 0.575, 0.808,
                 0.788, 0.772, 0.695, 0.538, 0.578, 0.43, 0.551, 0.576, 0.704,
                 0.44, 0.49, 0.385, 0.578, 0.778, 0.62, 0.494, 0.618, 0.645,
                 0.868, 0.75, 0.534, 0.47, 0.371, 0.575, 0.517, 0.448, 0.794,
                 0.544, 0.602, 0.632, 0.508, 0.721))
  expect_equal(round(nmrfilter_rmd$Corr, 3),
               c(0.917, 0.922, 0.949, 0.914, 0.896, 0.954, 0.957, 0.951, 0.951,
                 0.955, 0.954, 0.957, 0.951, 0.94, 0.941, 0.946, 0.921, 0.932,
                 0.915, 0.907, 0.902, 0.939, 0.885, 0.918, 0.91, 0.913, 0.899,
                 0.931, 0.94, 0.951, 0.945, 0.947, 0.946, 0.943, 0.948, 0.932,
                 0.925, 0.93, 0.926, 0.926, 0.922))

  # Run standard diagnostics on the filter object attributes.
  expect_equal(attr(nmrfilter_rmd, "sample_names"),
               c("F3-049", "F3-097", "F3-002", "F3-050", "F3-098", "F3-013",
                 "F3-061", "F3-109", "F3-062", "F3-110", "F3-025", "F3-073",
                 "F3-121", "F3-026", "F3-074", "F3-122", "F3-085", "F3-133",
                 "F3-038", "F3-086", "F3-134", "F4-001", "F4-009", "F4-065",
                 "F4-005", "F4-069", "F4-037", "F4-017", "F4-045", "F4-021",
                 "F4-049", "F4-081", "F4-025", "F4-053", "F4-085", "F4-029",
                 "F4-057", "F4-089", "F4-033", "F4-061", "F4-093"))
  expect_equal(attr(nmrfilter_rmd, "group_DF"),
               attr(nmrdata, "group_DF"))
  expect_equal(attr(nmrfilter_rmd, "df"),
               4)
  expect_equal(attr(nmrfilter_rmd, "metrics"),
               c("MAD", "Kurtosis", "Skewness", "Corr"))

  # Test rmd_filter with singleton groups --------------------------------------

  # Forge a filter object for the peptide data.
  expect_warning(pfilter_rmd_sg <- rmd_filter(omicsData = pdata_sg,
                                              ignore_singleton_groups = TRUE),
                 paste("Use the results of the RMD filter with caution due to",
                       "a small number of samples relative to the number of",
                       "metrics being used. Consider reducing the number of",
                       "metrics being used.",
                       sep = " "))

  # Check the class of the filter object.
  expect_s3_class(pfilter_rmd_sg,
                  c("rmdFilt", "data.frame"))

  # Checkify the peptide filter object data frame.
  expect_equal(dim(pfilter_rmd_sg),
               c(9, 9))
  expect_equal(round(pfilter_rmd_sg$Log2.md, 3),
               c(6.321, 6.015, 6.225, 6.321, 5.608,
                 5.976, 8.923, 6.132, 6.380))
  expect_equal(round(pfilter_rmd_sg$pvalue, 16),
               c(9.000000e-16, 1.299400e-12, 1.030000e-14, 9.000000e-16,
                 2.459952e-09, 2.959300e-12, 0.000000e+00, 9.600000e-14,
                 2.000000e-16))
  expect_equal(round(pfilter_rmd_sg$MAD, 3),
               c(1.331, 1.185, 1.241, 1.228, 1.228,
                 1.324, 1.203, 1.081, 1.169))
  expect_equal(round(pfilter_rmd_sg$Kurtosis, 3),
               c(1.250, 0.962, 1.085, -0.016, 0.523,
                 -0.274, -0.108, 0.174, -0.084))
  expect_equal(round(pfilter_rmd_sg$Skewness, 3),
               c(-0.429, -0.269, -0.296, 0.170, -0.045,
                 0.282, -0.006, 0.286, 0.355))
  expect_equal(round(pfilter_rmd_sg$Corr, 3),
               c(0.971, 0.975, 0.971, 0.975, 0.974,
                 0.966, 0.960, 0.941, 0.966))
  expect_equal(round(pfilter_rmd_sg$Proportion_Missing, 3),
               c(0.173, 0.180, 0.180, 0.207, 0.187,
                 0.220, 0.153, 0.260, 0.220))

  # Inspectify the filter object attributes.
  expect_equal(attr(pfilter_rmd_sg, "sample_names"),
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection7", "Infection8",
                 "Infection9"))
  temp_grp_df <- attr(pdata_sg, "group_DF")[-10, ]
  expect_equal(attr(pfilter_rmd_sg, "group_DF"),
               temp_grp_df)
  expect_equal(attr(pfilter_rmd_sg, "df"),
               5)
  expect_equal(attr(pfilter_rmd_sg, "metrics"),
               c("MAD", "Kurtosis", "Skewness", "Corr", "Proportion_Missing"))

  # Forge a filter object for the peptide data. There is one singleton group but
  # it is not ignored for the following tests. _f: ignore singletons is FALSE
  pfilter_rmd_sg_f <- rmd_filter(omicsData = pdata_sg,
                               ignore_singleton_groups = FALSE)

  # Check the class of the filter object.
  expect_s3_class(pfilter_rmd_sg_f,
                  c("rmdFilt", "data.frame"))

  # Checkify the peptide filter object data frame.
  expect_equal(dim(pfilter_rmd_sg_f),
               c(10, 9))
  expect_equal(round(pfilter_rmd_sg_f$Log2.md, 3),
               c(2.292, 0.680, 1.661, 0.385, 1.255,
                 2.047, 1.655, 2.474, 2.772, 5.614))
  expect_equal(round(pfilter_rmd_sg_f$pvalue, 10),
               c(0.4285077890, 0.9010201653, 0.6750442797, 0.9343177198,
                 0.7934655431, 0.5305775013, 0.6770948538, 0.3520673347,
                 0.2333855866, 0.0000000022))
  expect_equal(round(pfilter_rmd_sg_f$MAD, 3),
               c(1.331, 1.185, 1.241, 1.228, 1.228,
                 1.324, 1.203, 1.081, 1.169, 1.142))
  expect_equal(round(pfilter_rmd_sg_f$Kurtosis, 3),
               c(1.250, 0.962, 1.085, -0.016, 0.523,
                 -0.274, -0.108, 0.174, -0.084, 0.286))
  expect_equal(round(pfilter_rmd_sg_f$Skewness, 3),
               c(-0.429, -0.269, -0.296, 0.170, -0.045,
                 0.282, -0.006, 0.286, 0.355, -0.143))
  expect_equal(round(pfilter_rmd_sg_f$Corr, 3),
               c(0.971, 0.975, 0.971, 0.975, 0.974,
                 0.966, 0.960, 0.941, 0.966, 0.903))
  expect_equal(round(pfilter_rmd_sg_f$Proportion_Missing, 3),
               c(0.173, 0.180, 0.180, 0.207, 0.187,
                 0.220, 0.153, 0.260, 0.220, 0.167))

  # Inspectify the filter object attributes.
  expect_equal(attr(pfilter_rmd_sg_f, "sample_names"),
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection7", "Infection8",
                 "Infection9", "Mock1"))
  expect_equal(attr(pfilter_rmd_sg_f, "group_DF"),
               attr(pdata_sg, "group_DF"))
  expect_equal(attr(pfilter_rmd_sg_f, "df"),
               5)
  expect_equal(attr(pfilter_rmd_sg_f, "metrics"),
               c("MAD", "Kurtosis", "Skewness", "Corr", "Proportion_Missing"))

  # Test applyFilt no singletons -----------------------------------------------

  # Peptide data ---------------

  # Apply the filter to the peptide data without singleton groups.
  pfiltered <- applyFilt(filter_object = pfilter_rmd,
                         omicsData = pdata,
                         pvalue_threshold = 0.0001,
                         min_num_biomolecules = 50)

  # Ensurify the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, "cnames"),
                   attr(pfiltered, "cnames"))
  expect_identical(attr(pdata, "check.names"),
                   attr(pfiltered, "check.names"))
  expect_identical(class(pdata),
                   class(pfiltered))

  # Examine the filters attribute.
  expect_equal(attr(pfiltered, "filters")[[1]]$type,
               "rmdFilt")
  expect_identical(attr(pfiltered, "filters")[[1]]$threshold,
                   0.0001)
  expect_equal(attr(pfiltered, "filters")[[1]]$filtered,
               "Infection7")
  expect_true(is.na(attr(pfiltered, "filters")[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(pfiltered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pfiltered$e_data[, 1])),
         num_miss_obs = sum(is.na(pfiltered$e_data)),
         prop_missing = (sum(is.na(pfiltered$e_data)) /
                           prod(dim(pfiltered$e_data[, -1]))),
         num_samps = ncol(pfiltered$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(pfiltered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pfiltered$e_meta$Protein)))
  )

  # Dissect the group_DF attribute.
  expect_equal(dim(attr(pfiltered, "group_DF")),
               c(11, 2))
  expect_equal(attr(pfiltered, "group_DF")$SampleID,
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection8", "Infection9",
                 "Mock1", "Mock2", "Mock3"))

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(pfiltered$e_data),
               c(150, 12))
  expect_equal(dim(pfiltered$f_data),
               c(11, 2))
  expect_equal(dim(pfiltered$e_meta),
               c(150, 4))

  # nmr data ---------------

  # Apply the filter to the nmr data with default min_num_biomolecules.
  expect_error(applyFilt(filter_object = nmrfilter_rmd,
                         omicsData = nmrdata,
                         pvalue_threshold = 0.0001,
                         min_num_biomolecules = 50),
               paste("There are fewer biomolecules in omicsData than",
                     "min_num_biomolecules \\(50\\).",
                     "See applyFilt for details.",
                     sep = " "))

  # Apply the filter to the nmr data with a lower biomolecule threshold.
  nmrfiltered <- applyFilt(filter_object = nmrfilter_rmd,
                           omicsData = nmrdata,
                           pvalue_threshold = 0.0001,
                           min_num_biomolecules = 30)

  # Ensurify the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(nmrdata, "cnames"),
                   attr(nmrfiltered, "cnames"))
  expect_identical(attr(nmrdata, "check.names"),
                   attr(nmrfiltered, "check.names"))
  expect_identical(class(nmrdata),
                   class(nmrfiltered))

  # Examine the filters attribute.
  expect_equal(attr(nmrfiltered, "filters")[[1]]$type,
               "rmdFilt")
  expect_identical(attr(nmrfiltered, "filters")[[1]]$threshold,
                   0.0001)
  expect_equal(attr(nmrfiltered, "filters")[[1]]$filtered,
               c("F3-098", "F4-009", "F4-005"))
  expect_true(is.na(attr(nmrfiltered, "filters")[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(nmrfiltered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrfiltered$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrfiltered$e_data)),
         prop_missing = (sum(is.na(nmrfiltered$e_data)) /
                           prod(dim(nmrfiltered$e_data[, -1]))),
         num_samps = ncol(nmrfiltered$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(nmrfiltered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(nmrfiltered$e_meta$nmrClass)))
  )

  # Dissect the group_DF attribute.
  expect_equal(dim(attr(nmrfiltered, "group_DF")),
               c(38, 2))
  expect_equal(attr(nmrfiltered, "group_DF")$SampleID,
               c("F3-049", "F3-097", "F3-002", "F3-050", "F3-013", "F3-061",
                 "F3-109", "F3-062", "F3-110", "F3-025", "F3-073", "F3-121",
                 "F3-026", "F3-074", "F3-122", "F3-085", "F3-133", "F3-038",
                 "F3-086", "F3-134", "F4-001", "F4-065", "F4-069", "F4-037",
                 "F4-017", "F4-045", "F4-021", "F4-049", "F4-081", "F4-025",
                 "F4-053", "F4-085", "F4-029", "F4-057", "F4-089", "F4-033",
                 "F4-061", "F4-093"))

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(nmrfiltered$e_data),
               c(38, 39))
  expect_equal(dim(nmrfiltered$f_data),
               c(38, 5))
  expect_equal(dim(nmrfiltered$e_meta),
               c(38, 3))

  # Test applyFilt with singletons ---------------------------------------------

  # Apply the filter to the peptide data without singleton groups.
  pfiltered_sg <- applyFilt(filter_object = pfilter_rmd_sg,
                         omicsData = pdata_sg,
                         pvalue_threshold = 0.0000000000000008,
                         min_num_biomolecules = 50)

  # Ensurify the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, "cnames"),
                   attr(pfiltered_sg, "cnames"))
  expect_identical(attr(pdata_sg, "check.names"),
                   attr(pfiltered_sg, "check.names"))
  expect_identical(class(pdata_sg),
                   class(pfiltered_sg))

  # Examine the filters attribute.
  expect_equal(attr(pfiltered_sg, "filters")[[1]]$type,
               "rmdFilt")
  expect_identical(attr(pfiltered_sg, "filters")[[1]]$threshold,
                   0.0000000000000008)
  expect_equal(attr(pfiltered_sg, "filters")[[1]]$filtered,
               c("Infection7", "Infection9"))
  expect_true(is.na(attr(pfiltered_sg, "filters")[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(pfiltered_sg, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pfiltered_sg$e_data[, 1])),
         num_miss_obs = sum(is.na(pfiltered_sg$e_data)),
         prop_missing = (sum(is.na(pfiltered_sg$e_data)) /
                           prod(dim(pfiltered_sg$e_data[, -1]))),
         num_samps = ncol(pfiltered_sg$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(pfiltered_sg, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pfiltered_sg$e_meta$Protein)))
  )

  # Dissect the group_DF attribute.
  expect_equal(dim(attr(pfiltered_sg, "group_DF")),
               c(8, 2))
  expect_equal(attr(pfiltered_sg, "group_DF")$SampleID,
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection8", "Mock1"))

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(pfiltered_sg$e_data),
               c(150, 9))
  expect_equal(dim(pfiltered_sg$f_data),
               c(8, 2))
  expect_equal(dim(pfiltered_sg$e_meta),
               c(150, 4))

  # Apply the filter to the peptide data without singleton groups.
  pfiltered_sg_f <- applyFilt(filter_object = pfilter_rmd_sg_f,
                              omicsData = pdata_sg,
                              pvalue_threshold = 0.0001,
                              min_num_biomolecules = 50)

  # Ensurify the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, "cnames"),
                   attr(pfiltered_sg_f, "cnames"))
  expect_identical(attr(pdata_sg, "check.names"),
                   attr(pfiltered_sg_f, "check.names"))
  expect_identical(class(pdata_sg),
                   class(pfiltered_sg_f))

  # Examine the filters attribute.
  expect_equal(attr(pfiltered_sg_f, "filters")[[1]]$type,
               "rmdFilt")
  expect_identical(attr(pfiltered_sg_f, "filters")[[1]]$threshold,
                   0.0001)
  expect_equal(attr(pfiltered_sg_f, "filters")[[1]]$filtered,
               "Mock1")
  expect_true(is.na(attr(pfiltered_sg_f, "filters")[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(pfiltered_sg_f, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pfiltered_sg_f$e_data[, 1])),
         num_miss_obs = sum(is.na(pfiltered_sg_f$e_data)),
         prop_missing = (sum(is.na(pfiltered_sg_f$e_data)) /
                           prod(dim(pfiltered_sg_f$e_data[, -1]))),
         num_samps = ncol(pfiltered_sg_f$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(pfiltered_sg_f, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pfiltered_sg_f$e_meta$Protein)))
  )

  # Dissect the group_DF attribute.
  expect_equal(dim(attr(pfiltered_sg_f, "group_DF")),
               c(9, 2))
  expect_equal(attr(pfiltered_sg_f, "group_DF")$SampleID,
               c("Infection1", "Infection2", "Infection3", "Infection4",
                 "Infection5", "Infection6", "Infection7", "Infection8",
                 "Infection9"))

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(pfiltered_sg_f$e_data),
               c(150, 10))
  expect_equal(dim(pfiltered_sg_f$f_data),
               c(9, 2))
  expect_equal(dim(pfiltered_sg_f$e_meta),
               c(150, 4))


  # Test VizSampNames attribute ------------------------------------------------

  # Assemble long sample names to make sure the VizSampNames attribute is
  # created correctly.
  names(pdata$e_data) <- c("Mass_Tag_ID",
                           "qwerty_one_infection1_asdf",
                           "qwerty_two_infection2_asdf",
                           "qwerty_three_infection3_asdf",
                           "qwerty_four_infection4_asdf",
                           "qwerty_five_infection5_asdf",
                           "qwerty_six_infection6_asdf",
                           "qwerty_seven_infection7_asdf",
                           "qwerty_eight_infection8_asdf",
                           "qwerty_nine_infection9_asdf",
                           "qwerty_one_mock1_asdf",
                           "qwerty_two_mock2_asdf",
                           "qwerty_three_mock3_asdf")

  # Change to the new long names in f_data
  pdata$f_data$SampleID <- c("qwerty_one_infection1_asdf",
                             "qwerty_two_infection2_asdf",
                             "qwerty_three_infection3_asdf",
                             "qwerty_four_infection4_asdf",
                             "qwerty_five_infection5_asdf",
                             "qwerty_six_infection6_asdf",
                             "qwerty_seven_infection7_asdf",
                             "qwerty_eight_infection8_asdf",
                             "qwerty_nine_infection9_asdf",
                             "qwerty_one_mock1_asdf",
                             "qwerty_two_mock2_asdf",
                             "qwerty_three_mock3_asdf")

  # Use delim and components to make tiny sample names.
  delim_u <- custom_sampnames(pdata, delim = "_", components = 3)

  # Don't forget to rerun the group_designation function!!
  delim_u <- group_designation(omicsData = delim_u,
                               main_effects = "Condition")

  # Forge a filter object for the peptide data with tiny sample names.
  delim_rmd <- rmd_filter(omicsData = delim_u,
                          ignore_singleton_groups = FALSE)

  # Sleuth around the sample names attributes.
  expect_equal(attr(delim_rmd, "sample_names"),
               c("qwerty_one_infection1_asdf", "qwerty_two_infection2_asdf",
                 "qwerty_three_infection3_asdf", "qwerty_four_infection4_asdf",
                 "qwerty_five_infection5_asdf", "qwerty_six_infection6_asdf",
                 "qwerty_seven_infection7_asdf", "qwerty_eight_infection8_asdf",
                 "qwerty_nine_infection9_asdf", "qwerty_one_mock1_asdf",
                 "qwerty_two_mock2_asdf", "qwerty_three_mock3_asdf"))
  expect_equal(attr(delim_rmd, "VizSampNames"),
               c("infection1", "infection2", "infection3", "infection4",
                 "infection5", "infection6", "infection7", "infection8",
                 "infection9", "mock1", "mock2", "mock3"))

  # Make sure the output is in the same order even though the names changed.
  expect_identical(delim_rmd$Group, pfilter_rmd$Group)
  expect_identical(delim_rmd$Log2.md, pfilter_rmd$Log2.md)
  expect_identical(delim_rmd$pvalue, pfilter_rmd$pvalue)
  expect_identical(delim_rmd$MAD, pfilter_rmd$MAD)
  expect_identical(delim_rmd$Kurtosis, pfilter_rmd$Kurtosis)
  expect_identical(delim_rmd$Skewness, pfilter_rmd$Skewness)
  expect_identical(delim_rmd$Corr, pfilter_rmd$Corr)
  expect_identical(delim_rmd$Proportion_Missing, pfilter_rmd$Proportion_Missing)

  # Test scenario when nothing is filtered -------------------------------------

  # Apply the filter with a value for pvalue_threshold that will not filter any
  # samples.
  expect_message(noFilta <- applyFilt(filter_object = pfilter_rmd,
                                      omicsData = pdata,
                                      pvalue_threshold = 0.0000000000000000001,
                                      min_num_biomolecules = 50),
                 paste("No samples were filtered with the value specified",
                       "for the pvalue_threshold argument.",
                       sep = " "))

  # The output of applyFilt should be the same as the omicsData object used as
  # the input because the filter was not applied. Therefore, the filters
  # attribute should remain how it was before running applyFilt.
  expect_identical(noFilta, pdata)

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
                                pair_id = "PairID",
                                pair_group = "Time",
                                pair_denom = "18")

  pair_filter <- rmd_filter(omicsData = pairdata)

  pair_filtered <- applyFilt(filter_object = pair_filter,
                             omicsData = pairdata,
                             pvalue_threshold = 0.001,
                             min_num_biomolecules = 50)

  # Holy paired rmd filter tests with no main effects, Batman! -----------------

  # Test all aspects of the rmd filter object.
  expect_equal(dim(pair_filter), c(30, 9))
  expect_equal(pair_filter$Name, fdata$Name)
  expect_equal(pair_filter$Group, rep("paired_diff", 30))
  expect_equal(round(pair_filter$Log2.md, 3),
               c(1.201, 0.455, 4.439, 2.086, 1.153, 0.794, 1.472, 0.316, 2.123,
                 1.651, 2.345, 3.44, 2.028, 1.374, 3.235, 2.474, 2.488, 2.533,
                 2.363, 2.58, 1.186, 1.056, 1.451, 1.504, -0.005, 2.21, 4.319,
                 1.654, 2.1, 2.47))
  expect_equal(round(pair_filter$pvalue, 3),
               c(0.807, 0.927, 0.001, 0.515, 0.817, 0.885, 0.735, 0.941, 0.499,
                 0.678, 0.406, 0.054, 0.538, 0.763, 0.094, 0.352, 0.346, 0.327,
                 0.399, 0.308, 0.81, 0.838, 0.741, 0.725, 0.963, 0.463, 0.001,
                 0.677, 0.509, 0.354))
  expect_equal(round(pair_filter$MAD, 3),
               c(0.1, 0.087, 0.09, 0.076, 0.094, 0.076, 0.095, 0.091, 0.091,
                 0.081, 0.075, 0.077, 0.07, 0.075, 0.091, 0.097, 0.089, 0.089,
                 0.084, 0.091, 0.088, 0.077, 0.09, 0.082, 0.088, 0.092, 0.099,
                 0.102, 0.105, 0.082))
  expect_equal(round(pair_filter$Kurtosis, 3),
               c(-0.212, -0.775, 0.196, -0.974, -0.511, -0.689, -0.784, -0.175,
                 -0.45, -0.248, -0.742, 0.279, -0.692, -0.763, 0.491, -0.939,
                 -0.962, 0.258, -0.767, -0.908, -0.434, -0.459, -0.43, -0.637,
                 -0.377, -0.646, 0.729, -0.238, -0.447, -0.799))
  expect_equal(round(pair_filter$Skewness, 3),
               c(-0.479, -0.123, -0.661, 0.104, -0.358, -0.094, -0.128, -0.423,
                 -0.331, -0.262, -0.128, -0.597, -0.05, 0.058, -0.685, -0.173,
                 0.109, -0.705, 0.057, -0.118, -0.25, -0.187, -0.433, -0.156,
                 -0.303, -0.322, -0.774, -0.435, -0.44, -0.215))
  expect_equal(round(pair_filter$Corr, 3),
               c(0.881, 0.867, 0.808, 0.896, 0.915, 0.88, 0.903, 0.889, 0.921,
                 0.86, 0.897, 0.86, 0.857, 0.873, 0.884, 0.862, 0.884, 0.898,
                 0.836, 0.877, 0.91, 0.899, 0.884, 0.892, 0.896, 0.915, 0.87,
                 0.901, 0.872, 0.877))
  expect_equal(round(pair_filter$Proportion_Missing, 3),
               c(0.22, 0.247, 0.353, 0.287, 0.267, 0.253, 0.273, 0.3, 0.187,
                 0.213, 0.207, 0.293, 0.267, 0.273, 0.233, 0.253, 0.18, 0.3,
                 0.26, 0.387, 0.207, 0.273, 0.287, 0.18, 0.213, 0.34, 0.333,
                 0.32, 0.233, 0.253))

  # Test the particulars of the filtered object.
  expect_equal(attr(pair_filtered, "filters")[[1]]$type,
               "rmdFilt")
  expect_identical(attr(pair_filtered, "filters")[[1]]$threshold,
                   0.001)
  expect_equal(attr(pair_filtered, "filters")[[1]]$filtered,
               c("Mock_0hr_3", "Mock_18hr_3"))
  expect_true(is.na(attr(pair_filtered, "filters")[[1]]$method))
  expect_equal(
    attr(pair_filtered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pair_filtered$e_data[, 1])),
         num_miss_obs = sum(is.na(pair_filtered$e_data)),
         prop_missing = (sum(is.na(pair_filtered$e_data)) /
                           prod(dim(pair_filtered$e_data[, -1]))),
         num_samps = ncol(pair_filtered$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  expect_equal(
    attr(pair_filtered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pair_filtered$e_meta$Protein)))
  )
  expect_equal(dim(attr(pair_filtered, "group_DF")),
               c(28, 2))
  expect_equal(attr(attr(pair_filtered, "group_DF"), "main_effects"),
               "no_main_effect")
  expect_equal(attr(attr(pair_filtered, "group_DF"), "nonsingleton_groups"),
               "paired_diff")
  expect_equal(dim(pair_filtered$e_data),
               c(150, 29))
  expect_equal(dim(pair_filtered$f_data),
               c(28, 6))
  expect_equal(dim(pair_filtered$e_meta),
               c(175, 6))

  # Expect warning if data has already been filtered ---------------------------
  
  # The original peptide data gets overwritten at some point. Load it as before.
  load(system.file("testdata",
                   "little_pdata.RData",
                   package = "pmartR"))
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")
  pdata <- edata_transform(omicsData = pdata,
                           data_scale = "log")
  pdata <- group_designation(omicsData = pdata,
                             main_effects = "Condition")
  
  # Create three identical filters using the same omicsData
  pfilter1 <- rmd_filter(omicsData = pdata,
                         ignore_singleton_groups = FALSE)
  pfilter2 <- rmd_filter(omicsData = pdata,
                         ignore_singleton_groups = FALSE)
  pfilter3 <- rmd_filter(omicsData = pdata,
                         ignore_singleton_groups = FALSE)
  
  # Apply two filters
  pfiltered1 <- applyFilt(filter_object = pfilter1,
                          omicsData = pdata,
                          pvalue_threshold = 0.0001,
                          min_num_biomolecules = 50)
  warnings <- capture_warnings(
    pfiltered2 <- applyFilt(filter_object = pfilter2,
                            omicsData = pfiltered1,
                            pvalue_threshold = 0.0001,
                            min_num_biomolecules = 50)
  )
  
  # The second applyFilt should generate warnings
  expect_match(
    warnings,
    "An RMD filter has already been applied to this data set.",
    all = FALSE
  )
  expect_match(
    warnings,
    "Specified samples Infection7 were not found in the e_data\\."
    all = FALSE
  )
  
  # Samples that do exist should be filtered even if there are samples that
  # don't exist
  warnings <- capture_warnings(
    pfiltered3 <- applyFilt(filter_object = pfilter3,
                            omicsData = pfiltered1,
                            pvalue_threshold = 0.48,
                            min_num_biomolecules = 50)
  )
  
  expect_match(
    warnings,
    "An RMD filter has already been applied to this data set.",
    all = FALSE
  )
  expect_match(
    warnings,
    "Specified samples Infection7 were not found in the e_data\\.",
    all = FALSE
  )
  
  expect_null(pfiltered3$e_data$Infection8)
})
