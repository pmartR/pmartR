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
               c(0.785,  0.897,  1.360,  0.095,  1.471,  3.913,  3.937,  1.553,
                 1.447,  1.995,  1.066,  2.728,  3.022,  2.919,  5.230,  1.644,
                 2.126,  1.256,  0.999,  1.966,  2.122, -1.172,  5.708,  4.629,
                 1.666,  1.107,  0.100,  1.712,  1.832,  4.236,  2.165,  2.833,
                 0.131,  1.757,  2.577,  2.262,  3.163,  2.180,  1.718,  3.113,
                 2.010))
  expect_equal(round(nmrfilter_rmd$pvalue, 3),
               c(0.786, 0.761, 0.633, 0.899, 0.597, 0.005, 0.004, 0.569, 0.605,
                 0.408, 0.719, 0.157, 0.087, 0.109, 0.000, 0.537, 0.359, 0.665,
                 0.736, 0.419, 0.360, 0.979, 0.000, 0.000, 0.529, 0.708, 0.899,
                 0.513, 0.469, 0.001, 0.344, 0.130, 0.895, 0.496, 0.202, 0.309,
                 0.062, 0.339, 0.511, 0.070, 0.402))
  expect_equal(round(nmrfilter_rmd$MAD, 3),
               c(1.411, 1.173, 1.234, 1.463, 1.632, 0.922, 0.951, 1.245, 1.176,
                 1.260, 1.461, 1.328, 1.264, 0.978, 1.047, 1.271, 1.365, 1.270,
                 1.390, 1.435, 1.650, 1.286, 0.788, 1.018, 1.120, 1.404, 1.409,
                 1.403, 1.529, 0.967, 1.197, 1.619, 1.409, 1.654, 1.634, 1.040,
                 0.997, 1.486, 1.413, 1.503, 1.502))
  expect_equal(round(nmrfilter_rmd$Kurtosis, 3),
               c(-0.268, -0.094, -0.105, -0.369, -0.642,  0.004, -0.114, -0.168,
                 -0.082, -0.216, -0.529, -0.619, -0.752, -0.022, -0.297, -0.371,
                 -0.090, -0.372, -0.418, -0.376, -0.774, -0.378, -0.173, -0.206,
                 -0.087, -0.380, -0.332, -0.194, -0.478, -0.283, -0.282, -0.580,
                 -0.430, -0.506, -0.652, -0.374, -0.324, -0.591, -0.509, -0.578,
                 -0.374))
  expect_equal(round(nmrfilter_rmd$Skewness, 3),
               c(0.638, 0.750, 0.772, 0.578, 0.440, 0.828, 0.730, 0.713, 0.808,
                 0.695, 0.430, 0.576, 0.490, 0.874, 0.517, 0.575, 0.788, 0.538,
                 0.551, 0.704, 0.385, 0.578, 0.494, 0.778, 0.868, 0.534, 0.575,
                 0.794, 0.632, 0.645, 0.750, 0.470, 0.517, 0.544, 0.508, 0.620,
                 0.618, 0.371, 0.448, 0.602, 0.721))
  expect_equal(round(nmrfilter_rmd$Corr, 3),
               c(0.949, 0.954, 0.954, 0.940, 0.915, 0.917, 0.914, 0.957, 0.951,
                 0.957, 0.941, 0.921, 0.907, 0.922, 0.896, 0.951, 0.955, 0.951,
                 0.946, 0.932, 0.902, 0.939, 0.910, 0.885, 0.931, 0.951, 0.946,
                 0.932, 0.926, 0.899, 0.940, 0.945, 0.943, 0.925, 0.926, 0.918,
                 0.913, 0.947, 0.948, 0.930, 0.922))

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
         data_types = NULL)
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
               c("F3-098", "F4-005", "F4-009"))
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
         data_types = NULL)
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
         data_types = NULL)
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
         data_types = NULL)
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

  # Test scenario when nothing is filtered -------------------------------------

  # Apply the filter with a value for pvalue_threshold that will not filter any
  # rows.
  expect_message(noFilta <- applyFilt(filter_object = pfilter_rmd,
                                      omicsData = pdata,
                                      pvalue_threshold = 0.00000001,
                                      min_num_biomolecules = 50),
                 paste("No samples were filtered with the value specified",
                       "for the pvalue_threshold argument.",
                       sep = " "))

  # The output of applyFilt should be the same as the omicsData object used as
  # the input because the filter was not applied. Therefore, the filters
  # attribute should remain how it was before running applyFilt.
  expect_identical(noFilta, pdata)

})
