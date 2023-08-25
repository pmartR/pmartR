context('class: isobaricpepData')

test_that('as.isobaricpepData returns the correct data frame and attributes', {
  # Load the reduced isobaric peptide data frames ------------------------------

  load(system.file('testdata',
    'little_isodata.RData',
    package = 'pmartR'
  ))

  # Run as.isobaricpepData with agreeable data frames --------------------------

  # Construct an isobaricpepData object with the edata, fdata, and emeta data
  # frames.
  isodata <- as.isobaricpepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Peptide',
    fdata_cname = 'Sample',
    emeta_cname = 'Protein'
  )

  # Check high level structure
  expect_equal(names(isodata), c("e_data", "f_data", "e_meta"))

  # Ensure the returned data frames are the correct dimension.
  expect_equal(
    dim(isodata$e_data),
    c(150, 13)
  )
  expect_equal(
    dim(isodata$f_data),
    c(12, 5)
  )
  expect_equal(
    dim(isodata$e_meta),
    c(150, 2)
  )

  # Confirm the correct attributes are present in the isobaricpepData object.
  expect_equal(
    names(attributes(isodata)),
    c(
      "names", "cnames", "data_info", "isobaric_info",
      "meta_info", "filters", "class"
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(isodata, "cnames"),
    list(
      edata_cname = "Peptide",
      emeta_cname = "Protein",
      fdata_cname = "Sample",
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(isodata, "data_info"),
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(isodata$e_data[, 1])),
      num_miss_obs = sum(is.na(isodata$e_data)),
      prop_missing = (sum(is.na(isodata$e_data)) /
        prod(dim(isodata$e_data[, -1]))),
      num_samps = ncol(isodata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Sleuth around the isobaric attributes.
  expect_equal(
    attr(isodata, "isobaric_info"),
    list(
      exp_cname = NA,
      channel_cname = NA,
      refpool_channel = NA,
      refpool_cname = NA,
      refpool_notation = NA,
      norm_info = list(is_normalized = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(isodata, "meta_info"),
    list(
      meta_data = TRUE,
      num_emeta = length(unique(isodata$e_meta$Protein))
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(isodata, "filters"), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(isodata, c("isobaricpepData", "pepData"))

  # Produce a isobaricpepData object only with the edata and fdata data frames.
  expect_message(
    isodata <- as.isobaricpepData(
      e_data = edata,
      f_data = fdata,
      e_meta = NULL,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Test'
    ),
    "emeta_cname set to NULL, no e_meta object was provided."
  )

  # Check high level structure
  expect_equal(names(isodata), c("e_data", "f_data", "e_meta"))

  # Ensure the returned data frames are the correct dimension.
  expect_equal(
    dim(isodata$e_data),
    c(150, 13)
  )
  expect_equal(
    dim(isodata$f_data),
    c(12, 5)
  )
  expect_null(isodata$e_meta)

  # Confirm the correct attributes are present in the isobaricpepData object.
  expect_equal(
    names(attributes(isodata)),
    c(
      "names", "cnames", "data_info", "isobaric_info",
      "meta_info", "filters", "class"
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(isodata, "cnames"),
    list(
      edata_cname = "Peptide",
      emeta_cname = NULL,
      fdata_cname = "Sample",
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(isodata, "data_info"),
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(isodata$e_data[, 1])),
      num_miss_obs = sum(is.na(isodata$e_data)),
      prop_missing = (sum(is.na(isodata$e_data)) /
        prod(dim(isodata$e_data[, -1]))),
      num_samps = ncol(isodata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Sleuth around the isobaric attributes.
  expect_equal(
    attr(isodata, "isobaric_info"),
    list(
      exp_cname = NA,
      channel_cname = NA,
      refpool_channel = NA,
      refpool_cname = NA,
      refpool_notation = NA,
      norm_info = list(is_normalized = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(isodata, "meta_info"),
    list(
      meta_data = FALSE,
      num_emeta = NULL
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(isodata, "filters"), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(isodata, c("isobaricpepData", "pepData"))

  # Run as.isobaricpepData with disagreeable data frames -----------------------

  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(
    as.isobaricpepData(
      e_data = data.frame(edata,
        Sample13 = edata[, 10]
      ),
      f_data = fdata,
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein'
    ),
    "1 samples from e_data not found in f_data"
  )

  # Check for an error when e_data has more rows than e_meta.
  expect_error(
    as.isobaricpepData(
      e_data = edata,
      f_data = fdata,
      e_meta = emeta[1:138, ],
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein'
    ),
    "Not all peptides in e_data are present in e_meta"
  )

  # Remove factors from the Sample column in fdata.
  fdata$Sample <- as.character(fdata$Sample)

  # Create an f_data object with an extra row.
  fdata_1 <- rbind(
    fdata,
    c('Sample13', 1, 113, 'Yes', 'A')
  )

  # Create a isobaricpepData object and check for a warning when the f_data
  # object has an extra row.
  expect_warning(
    isodata <- as.isobaricpepData(
      e_data = edata,
      f_data = fdata_1,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample'
    ),
    paste("Extra samples were found in f_data that were not in",
      "e_data. These have been removed from f_data.",
      sep = ' '
    )
  )

  # Check high level structure
  expect_equal(names(isodata), c("e_data", "f_data", "e_meta"))

  # Check the dimensions of the e_data and f_data data frames.
  expect_equal(
    dim(isodata$e_data),
    c(150, 13)
  )
  expect_equal(
    dim(isodata$f_data),
    c(12, 5)
  )
  expect_null(isodata$emeta)

  # Confirm the correct attributes are present in the isobaricpepData object.
  expect_equal(
    names(attributes(isodata)),
    c(
      "names", "cnames", "data_info", "isobaric_info",
      "meta_info", "filters", "class"
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(isodata, "cnames"),
    list(
      edata_cname = "Peptide",
      emeta_cname = NULL,
      fdata_cname = "Sample",
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(isodata, "data_info"),
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(isodata$e_data[, 1])),
      num_miss_obs = sum(is.na(isodata$e_data)),
      prop_missing = (sum(is.na(isodata$e_data)) /
        prod(dim(isodata$e_data[, -1]))),
      num_samps = ncol(isodata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Sleuth around the isobaric attributes.
  expect_equal(
    attr(isodata, "isobaric_info"),
    list(
      exp_cname = NA,
      channel_cname = NA,
      refpool_channel = NA,
      refpool_cname = NA,
      refpool_notation = NA,
      norm_info = list(is_normalized = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(isodata, "meta_info"),
    list(
      meta_data = FALSE,
      num_emeta = NULL
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(isodata, "filters"), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(isodata, c("isobaricpepData", "pepData"))

  # Generate a isobaricpepData object and check for a warning when the e_meta
  # object has more rows than e_data.
  expect_warning(
    isodata <- as.isobaricpepData(
      e_data = edata[1:117, ],
      f_data = fdata,
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein'
    ),
    paste('Extra peptides were found in e_meta that were not in',
      'e_data. These have been removed from e_meta.',
      sep = ' '
    )
  )

  # Check high level structure
  expect_equal(names(isodata), c("e_data", "f_data", "e_meta"))

  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(
    dim(isodata$e_data),
    c(117, 13)
  )
  expect_equal(
    dim(isodata$f_data),
    c(12, 5)
  )
  expect_equal(
    dim(isodata$e_meta),
    c(117, 2)
  )

  # Confirm the correct attributes are present in the isobaricpepData object.
  expect_equal(
    names(attributes(isodata)),
    c(
      "names", "cnames", "data_info", "isobaric_info",
      "meta_info", "filters", "class"
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(isodata, "cnames"),
    list(
      edata_cname = "Peptide",
      emeta_cname = "Protein",
      fdata_cname = "Sample",
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(isodata, "data_info"),
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(isodata$e_data[, 1])),
      num_miss_obs = sum(is.na(isodata$e_data)),
      prop_missing = (sum(is.na(isodata$e_data)) /
        prod(dim(isodata$e_data[, -1]))),
      num_samps = ncol(isodata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Sleuth around the isobaric attributes.
  expect_equal(
    attr(isodata, "isobaric_info"),
    list(
      exp_cname = NA,
      channel_cname = NA,
      refpool_channel = NA,
      refpool_cname = NA,
      refpool_notation = NA,
      norm_info = list(is_normalized = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(isodata, "meta_info"),
    list(
      meta_data = TRUE,
      num_emeta = length(unique(isodata$e_meta$Protein))
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(isodata, "filters"), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(isodata, c("isobaricpepData", "pepData"))

  # Check the technical replicates column for correct structure.
  expect_error(
    as.isobaricpepData(
      e_data = edata,
      f_data = data.frame(fdata,
        tReps = fdata[, 1]
      ),
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein',
      techrep_cname = 'tReps'
    ),
    paste('Specified technical replicate column had a unique value',
      'for each row.  Values should specify groups of technical',
      'replicates belonging to a biological sample.',
      sep = ' '
    )
  )

  set.seed(5)

  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(
    edata,
    edata[sample(1:150, 8), ]
  )

  # Create a isobaricpepData object with some of the rows of e_data repeated.
  isodata <- as.isobaricpepData(
    e_data = edata_1,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Peptide',
    fdata_cname = 'Sample',
    emeta_cname = 'Protein'
  )

  # Check high level structure
  expect_equal(names(isodata), c("e_data", "f_data", "e_meta"))

  # Verify that the returned data frames are the correct dimension.
  expect_equal(
    dim(isodata$e_data),
    c(150, 13)
  )
  expect_equal(
    dim(isodata$f_data),
    c(12, 5)
  )
  expect_equal(
    dim(isodata$e_meta),
    c(150, 2)
  )

  # Confirm the correct attributes are present in the isobaricpepData object.
  expect_equal(
    names(attributes(isodata)),
    c(
      "names", "cnames", "data_info", "isobaric_info",
      "meta_info", "filters", "class"
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(isodata, "cnames"),
    list(
      edata_cname = "Peptide",
      emeta_cname = "Protein",
      fdata_cname = "Sample",
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(isodata, "data_info"),
    list(
      data_scale_orig = "abundance",
      data_scale = "abundance",
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(isodata$e_data[, 1])),
      num_miss_obs = sum(is.na(isodata$e_data)),
      prop_missing = (sum(is.na(isodata$e_data)) /
        prod(dim(isodata$e_data[, -1]))),
      num_samps = ncol(isodata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Sleuth around the isobaric attributes.
  expect_equal(
    attr(isodata, "isobaric_info"),
    list(
      exp_cname = NA,
      channel_cname = NA,
      refpool_channel = NA,
      refpool_cname = NA,
      refpool_notation = NA,
      norm_info = list(is_normalized = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(isodata, "meta_info"),
    list(
      meta_data = TRUE,
      num_emeta = length(unique(isodata$e_meta$Protein))
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(isodata, "filters"), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(isodata, c("isobaricpepData", "pepData"))

  # Change the values in each of the samples for the repeated IDs.
  edata_1[151:158, 2:13] <- edata_1[151:158, 2:13] * 1.1

  # Forge a isobaricpepData object with some of the peptide IDs repeated but the
  # values for each of the repeated peptide IDs is different from the original
  # data.
  expect_error(
    as.isobaricpepData(
      e_data = edata_1,
      f_data = fdata,
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein'
    ),
    "The 'edata_cname' identifier is non-unique."
  )

  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(
    as.isobaricpepData(
      e_data = edata_1,
      f_data = fdata,
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = NULL
    ),
    'Since e_meta is non-NULL, emeta_cname must also be non-NULL.'
  )

  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(
    as.isobaricpepData(
      e_data = edata_1,
      f_data = fdata,
      e_meta = emeta,
      edata_cname = 'Peptide',
      fdata_cname = 'Sample',
      emeta_cname = 'Protein2'
    ),
    paste('Mapping variable column',
      'Protein2',
      'not found in e_meta. See details of as.isobaricpepData',
      'for specifying column names.',
      sep = ' '
    )
  )
})
