context('class: seqData')


### TO do: add emeta

test_that('as.seqData returns the correct data frame and attributes', {
  # Load the reduced peptide data frames ---------------------------------------

  load(system.file('testdata',
    'little_seqdata.RData',
    package = 'pmartR'
  ))

  # Run as.seqData with agreeable data frames ----------------------------------

  # Construct a seqData object with the edata, fdata, and emeta data frames.
  seqdata <- as.seqData(
    e_data = edata,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  )

  # Check high level structure
  expect_equal(names(seqdata), c("e_data", "f_data", "e_meta"))

  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(seqdata$e_data), c(1200, 41))
  expect_equal(dim(seqdata$f_data), c(40, 4))

  # Confirm the correct attributes are present in the seqData object.
  expect_equal(
    names(attributes(seqdata)),
    c(
      'names', 'cnames', 'data_info', 'meta_info',
      'filters', 'class'
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(seqdata, 'cnames'),
    list(
      edata_cname = 'ID_REF',
      emeta_cname = NULL,
      fdata_cname = 'Samples',
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(seqdata, 'data_info'),
    list(
      data_scale_orig = 'counts',
      data_scale = 'counts',
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(seqdata$e_data[, 1])),
      num_zero_obs = sum(seqdata$e_data == 0),
      prop_zeros = (sum(seqdata$e_data == 0) /
        prod(dim(seqdata$e_data[, -1]))),
      num_samps = ncol(seqdata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(seqdata, 'filters'), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(seqdata, 'seqData')

  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(seqdata$e_data), c(1200, 41))
  expect_equal(dim(seqdata$f_data), c(40, 4))

  # Confirm the correct attributes are present in the seqData object.
  expect_equal(
    names(attributes(seqdata)),
    c(
      'names', 'cnames', 'data_info', 'meta_info',
      'filters', 'class'
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(seqdata, 'cnames'),
    list(
      edata_cname = 'ID_REF',
      emeta_cname = NULL,
      fdata_cname = 'Samples',
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(seqdata, 'data_info'),
    list(
      data_scale_orig = 'counts',
      data_scale = 'counts',
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(seqdata$e_data[, 1])),
      num_zero_obs = sum(seqdata$e_data == 0),
      prop_zeros = (sum(seqdata$e_data == 0) /
        prod(dim(seqdata$e_data[, -1]))),
      num_samps = ncol(seqdata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(seqdata, 'meta_info'),
    list(
      meta_data = FALSE,
      num_emeta = NULL
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(seqdata, 'filters'), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(seqdata, 'seqData')

  # Run as.seqData with disagreeable data frames -------------------------------

  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(
    as.seqData(
      e_data = data.frame(edata,
        Mock4 = edata[, 10]
      ),
      f_data = fdata,
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples'
    ),
    '1 samples from e_data not found in f_data'
  )


  fdata

  # Create an f_data object with an extra row.
  fdata_1 <- rbind(fdata, data.frame(
    Samples = 'uterus_PBS_R5',
    Tissue = 'uterus',
    Treatment = 'PBS',
    Label = 'R5'
  ))

  # Create a seqData object and check for a warning when the f_data object has
  # an extra row.
  testthat::expect_warning(
    seqdata <- as.seqData(
      e_data = edata,
      f_data = fdata_1,
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples'
    ),
    paste('Extra samples were found in f_data that were not in',
      'e_data. These have been removed from f_data.',
      sep = ' '
    )
  )

  # Check high level structure
  expect_equal(names(seqdata), c("e_data", "f_data", "e_meta"))

  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(seqdata$e_data), c(1200, 41))
  expect_equal(dim(seqdata$f_data), c(40, 4))

  # Confirm the correct attributes are present in the seqData object.
  expect_equal(
    names(attributes(seqdata)),
    c(
      'names', 'cnames', 'data_info', 'meta_info',
      'filters', 'class'
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(seqdata, 'cnames'),
    list(
      edata_cname = 'ID_REF',
      emeta_cname = NULL,
      fdata_cname = 'Samples',
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(seqdata, 'data_info'),
    list(
      data_scale_orig = 'counts',
      data_scale = 'counts',
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(seqdata$e_data[, 1])),
      num_zero_obs = sum(seqdata$e_data == 0),
      prop_zeros = (sum(seqdata$e_data == 0) /
        prod(dim(seqdata$e_data[, -1]))),
      num_samps = ncol(seqdata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(seqdata, 'meta_info'),
    list(
      meta_data = FALSE,
      num_emeta = NULL
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(seqdata, 'filters'), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(seqdata, 'seqData')

  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(
    dim(seqdata$e_data),
    c(1200, 41)
  )
  expect_equal(dim(seqdata$f_data), c(40, 4))

  # Confirm the correct attributes are present in the seqData object.
  expect_equal(
    names(attributes(seqdata)),
    c(
      'names', 'cnames', 'data_info', 'meta_info',
      'filters', 'class'
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(seqdata, 'cnames'),
    list(
      edata_cname = 'ID_REF',
      emeta_cname = NULL,
      fdata_cname = 'Samples',
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(seqdata, 'data_info'),
    list(
      data_scale_orig = 'counts',
      data_scale = 'counts',
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(seqdata$e_data[, 1])),
      num_zero_obs = sum(seqdata$e_data == 0),
      prop_zeros = (sum(seqdata$e_data == 0) /
        prod(dim(seqdata$e_data[, -1]))),
      num_samps = ncol(seqdata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(seqdata, 'filters'), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(seqdata, 'seqData')

  # Check the technical replicates column for correct structure.
  expect_error(
    as.seqData(
      e_data = edata,
      f_data = data.frame(fdata,
        tReps = fdata[, 1]
      ),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      techrep_cname = 'tReps'
    ),
    paste('Specified technical replicate column had a unique value',
      'for each row.  Values should specify groups of technical',
      'replicates belonging to a biological sample.',
      sep = ' '
    )
  )

  set.seed(5)

  # Fabricate an e_data object with some of the IDs repeated.
  edata_1 <- rbind(
    edata,
    edata[sample(1:150, 8), ]
  )

  # Create a seqData object with some of the rows of e_data repeated.
  seqdata <- as.seqData(
    e_data = edata_1,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  )

  # Check high level structure
  expect_equal(names(seqdata), c("e_data", "f_data", "e_meta"))

  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(seqdata$e_data), c(1200, 41))
  expect_equal(dim(seqdata$f_data), c(40, 4))

  # Confirm the correct attributes are present in the seqData object.
  expect_equal(
    names(attributes(seqdata)),
    c(
      'names', 'cnames', 'data_info', 'meta_info',
      'filters', 'class'
    )
  )

  # Scrutinize the column names attribute.
  expect_equal(
    attr(seqdata, 'cnames'),
    list(
      edata_cname = 'ID_REF',
      emeta_cname = NULL,
      fdata_cname = 'Samples',
      techrep_cname = NULL
    )
  )

  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(seqdata, 'data_info'),
    list(
      data_scale_orig = 'counts',
      data_scale = 'counts',
      norm_info = list(is_normalized = FALSE),
      num_edata = length(unique(seqdata$e_data[, 1])),
      num_zero_obs = sum(seqdata$e_data == 0),
      prop_zeros = (sum(seqdata$e_data == 0) /
        prod(dim(seqdata$e_data[, -1]))),
      num_samps = ncol(seqdata$e_data[, -1]),
      data_types = NULL,
      batch_info = list(is_bc = FALSE)
    )
  )

  # Take a looksie at the filters attribute.
  expect_identical(attr(seqdata, 'filters'), list())

  # Ensure the omicsData object is classy.
  expect_s3_class(seqdata, 'seqData')

  # Change the values in each of the samples for the repeated IDs.
  # edata_1[151:158, 2:13] <- edata_1[151:158, 2:13] * 1.1
  edata_1 <- edata
  edata_1[2, 1] <- edata_1[1, 1]
  # Forge a seqData object with some of the peptide IDs repeated but the values
  # for each of the repeated peptide IDs is different from the original data.
  expect_error(
    as.seqData(
      e_data = edata_1,
      f_data = fdata,
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples'
    ),
    "The 'edata_cname' identifier is non-unique."
  )

  # Check for an error when e_meta is non-null but emeta_cname is null.
  # expect_error(as.seqData(e_data = edata_1,
  #                         f_data = fdata,
  #                         edata_cname = 'ID_REF',
  #                         fdata_cname = 'Samples'),
  #              'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')

  # Verify there is an error when emeta_cname is not a column name in e_meta.
  # expect_error(as.seqData(e_data = edata_1,
  #                         f_data = fdata,
  #                         e_meta = emeta,
  #                         edata_cname = 'ID_REF',
  #                         fdata_cname = 'Samples',
  #                         emeta_cname = 'Protein2'),
  #              paste('Mapping variable column',
  #                    'Protein2',
  #                    'not found in e_meta. See details of as.seqData for',
  #                    'specifying column names.',
  #                    sep = ' '))

  ## Warning for non-int (allowed for visualizing)
  testthat::expect_warning(
    as.seqData(
      e_data = cbind(edata[1], edata[-1] / 2),
      edata_cname = 'ID_REF',
      f_data = fdata,
      fdata_cname = "Samples"
    ),
    "Non-integers detected"
  )


  ## Error for non "count"
  testthat::expect_error(
    as.seqData(
      e_data = edata,
      edata_cname = 'ID_REF',
      f_data = fdata,
      fdata_cname = "Samples",
      data_scale = "log2"
    ),
    "data_scale must be 'counts' for as.seqData"
  )

  ## NA replaced by zeros
  e_data_temp <- edata
  e_data_temp[e_data_temp == 0] <- NA

  x <- testthat::expect_message(
    as.seqData(
      e_data = e_data_temp,
      edata_cname = 'ID_REF',
      f_data = fdata,
      fdata_cname = "Samples"
    )$e_data,
    "instances of NA have been replaced with 0"
  )
  testthat::expect_false(any(is.na(x)))

  ###### Make the class checks happen

  ## Supports Data table conversion
  expect_error(as.seqData(
    e_data = 5,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  ), "e_data must be of class 'data.frame'")

  expect_error(as.seqData(
    e_data = edata,
    f_data = 5,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  ), "f_data must be of class 'data.frame'")

  expect_error(as.seqData(
    e_data = edata,
    f_data = fdata,
    e_meta = 5,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  ), "e_meta must be of class 'data.frame'")

  emeta <- data.frame(edata[1], class = 1:nrow(edata))

  ## Supports Data table conversion also
  seqdata <- as.seqData(
    e_data = data.table::as.data.table(edata),
    f_data = data.table::as.data.table(fdata),
    e_meta = data.table::as.data.table(emeta),
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples',
    emeta_cname = "class"
  )

  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 5,
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "must be of the class 'character'"
  )

  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 5,
      emeta_cname = "class"
    ),
    "must be of the class 'character'"
  )


  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class",
      is_normalized = 5
    ),
    "must be of the class 'logical'"
  )

  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class",
      is_bc = 5
    ),
    "must be of the class 'logical'"
  )

  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class",
      batch_info = 5
    ),
    "must be of the class 'list'"
  )

  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class",
      data_types = 5
    ),
    "must be of the class 'character'"
  )

  ## Valid data_scale

  expect_error(
    as.nmrData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class",
      data_scale = 5
    ),
    "data_scale must be one of the"
  )

  expect_error(
    as.nmrData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      techrep_cname = 5,
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "techrep_cname must be a character string "
  )


  ## cnames correct

  # Bad edata
  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_RF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "not found in"
  )
  # Bad techrep
  expect_error(
    as.nmrData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      techrep_cname = 'Samples',
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "not found in"
  )

  # Bad fdata
  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Sampes',
      emeta_cname = "class"
    ),
    "not found in"
  )

  # Bad emeta
  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "cass"
    ),
    "not found in"
  )

  emeta2 <- emeta
  colnames(emeta2) <- c("t", "class")

  # Bad edata/emeta
  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata),
      e_meta = data.table::as.data.table(emeta2),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "not found in"
  )

  fdata2 <- fdata[1]

  # Bad fdata
  expect_error(
    as.seqData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata2),
      e_meta = data.table::as.data.table(emeta),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "must contain at least 2"
  )

  fdata3 <- fdata
  fdata3$techrep <- 1:nrow(fdata)

  ## bad techrep
  expect_error(
    as.nmrData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata3),
      e_meta = data.table::as.data.table(emeta),
      techrep_cname = 'techrep',
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "Specified technical replicate column"
  )

  # Bad emeta
  emeta2 <- rbind(emeta, emeta[1, ])

  expect_error(
    as.nmrData(
      e_data = data.table::as.data.table(edata),
      f_data = data.table::as.data.table(fdata3),
      e_meta = data.table::as.data.table(emeta2),
      edata_cname = 'ID_REF',
      fdata_cname = 'Samples',
      emeta_cname = "class"
    ),
    "Not all e_data cname and "
  )
})
