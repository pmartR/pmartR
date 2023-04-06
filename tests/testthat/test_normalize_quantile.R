context('normalize: quantile')

test_that('normalize_quantile produces the correct output', {
  # Load the nmr data frames ---------------------------------------------------

  # Use nmr data because it does not contain any missing values and the quantile
  # normalization method cannot be used when there is missing data.
  load(system.file('testdata',
    'nmrData.RData',
    package = 'pmartR'
  ))

  # Create an nmrData object.
  nmrdata <- as.nmrData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = "Metabolite",
    fdata_cname = "SampleID",
    emeta_cname = "Bin"
  )

  # Copy nmrdata to create an nmrData object where the ID column is not the
  # first column in e_data.
  mixed_nmr <- nmrdata

  # Rearrange the columns of e_data.
  mixed_nmr$e_data <- nmrdata$e_data[, c(2:25, 1, 26:42)]

  # Natural logate the data.
  lognmrdata <- edata_transform(
    omicsData = nmrdata,
    data_scale = "log"
  )

  # Calculate the normalization standards --------------------------------------

  # Create a copy of edata. This will become the standard for the normalized
  # data.
  edata_stand <- edata

  # Sort the each value in ascending order by molecule (by row).
  x_sorted <- apply(
    t(edata_stand[, -1]), 2,
    function(c) sort(c, na.last = FALSE)
  )

  # Fish out the indices of the sorted data.
  x_ordered <- apply(
    t(edata_stand[, -1]), 2,
    function(c) order(c, na.last = FALSE)
  )

  # Compute the mean of each column. This will produce a vector of means that
  # are increasing.
  mean_sorted <- rowMeans(x_sorted)

  # Replace abundance values with their corresponding mean by metabolite and
  # reorder the means to the original abundance value order.
  for (i in 1:ncol(x_sorted)) {
    edata_stand[i, (x_ordered[, i] + 1)] <- mean_sorted
  }

  # Copy the normalized data. This will become the standard for the normalized
  # data on the log scale.
  edata_stand_log <- edata_stand

  # Natural logate the normalized standard.
  edata_stand_log[, -1] <- log(edata_stand[, -1])

  # test error for inccorect class
  expect_error(
    normalize_quantile(omicsData = fdata),
    "omicsData must be of class 'pepData', 'proData', "
  )


  # Normalize data: abundance scale --------------------------------------------

  # Normalize the nmr data on the abundance scale.
  a_norm <- normalize_quantile(omicsData = nmrdata)

  # Put on my sleuth hat and start investigating the output.
  expect_identical(a_norm$e_data, edata_stand)
  expect_s3_class(a_norm, "nmrData")

  # Check out the attributes.
  expect_identical(
    attr(a_norm, "cnames"),
    attr(nmrdata, "cname")
  )
  expect_identical(
    attr(a_norm, "nmr_info"),
    attr(nmrdata, "nmr_info")
  )
  expect_identical(
    attr(a_norm, "meta_info"),
    attr(nmrdata, "meta_info")
  )
  expect_identical(
    attr(a_norm, "filters"),
    attr(nmrdata, "filters")
  )
  expect_equal(
    attributes(a_norm)$data_info,
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(
        is_normalized = TRUE,
        norm_type = "quantile"
      ),
      num_edata = 38,
      num_miss_obs = 0,
      prop_missing = 0,
      num_samps = 41,
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Normalize the data with the mixed up columns.
  mixed_norm <- normalize_quantile(omicsData = mixed_nmr)

  # Look into the index of the ID column. It should now be the first column.
  expect_equal(which(names(mixed_norm$e_data) == "Metabolite"), 1)

  # Ensure the normalization is the same as the original nmr data.
  expect_identical(
    mixed_norm,
    a_norm
  )

  # Normalize data: log scale --------------------------------------------------

  # Normalize the nmr data on the abundance scale.
  l_norm <- normalize_quantile(omicsData = lognmrdata)

  # Put on my sleuth hat and start investigating the output.
  expect_equal(l_norm$e_data, edata_stand_log)
  expect_s3_class(l_norm, "nmrData")

  # Check out the attributes.
  expect_identical(
    attr(l_norm, "cnames"),
    attr(lognmrdata, "cname")
  )
  expect_identical(
    attr(l_norm, "nmr_info"),
    attr(lognmrdata, "nmr_info")
  )
  expect_identical(
    attr(l_norm, "meta_info"),
    attr(lognmrdata, "meta_info")
  )
  expect_identical(
    attr(l_norm, "filters"),
    attr(lognmrdata, "filters")
  )
  expect_equal(
    attributes(l_norm)$data_info,
    list(
      data_scale_orig = "abundance",
      data_scale = "log",
      norm_info = list(
        is_normalized = TRUE,
        norm_type = "quantile"
      ),
      num_edata = 38,
      num_miss_obs = 0,
      prop_missing = 0,
      num_samps = 41,
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Normailze data with missing values -----------------------------------------

  # Load the lipid data.
  load(system.file('testdata',
    'lipidData.RData',
    package = 'pmartR'
  ))

  # Create an lipidData object.
  ldata <- as.lipidData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = "LipidCommonName",
    fdata_cname = "Sample_Name",
    emeta_cname = "LipidClass"
  )

  # Try (and fail) to quantile normalize data with missing values.
  expect_error(
    normalize_quantile(omicsData = ldata),
    paste("The proportion of missing data is 0.55. Quantile",
      "normalization only works with complete data. Consider",
      "using SPANS to choose an appropriate normalization",
      "method for a dataset that includes missing values.",
      sep = " "
    )
  )
})
