context('edata transform')

test_that('edata_transform correctly transforms the data',{
  
  # Load the data and create a pepData object ----------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Construct a pepData object.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein',
                      data_scale_orig = "abundance")
  
  # Create test standards ------------------------------------------------------
  
  # Natural log transfiguration ---------------
  
  # Apply a natrual log transmogrification to e_data.
  ledata <- edata
  ledata[, -1] <- log(ledata[, -1])
  lpdata <- as.pepData(e_data = ledata,
                       f_data = fdata,
                       e_meta = emeta,
                       edata_cname = 'Mass_Tag_ID',
                       fdata_cname = 'SampleID',
                       emeta_cname = 'Protein',
                       data_scale_orig = "log")
  
  # Log base 2 transfiguration ---------------
  
  # Apply a log 2 transmogrification to e_data.
  l2edata <- edata
  l2edata[, -1] <- log2(l2edata[, -1])
  l2pdata <- as.pepData(e_data = l2edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Mass_Tag_ID',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'Protein',
                        data_scale_orig = "log2")
  
  # Log base 10 transfiguration ---------------
  
  # Apply a log 10 transmogrification to e_data.
  l10edata <- edata
  l10edata[, -1] <- log10(l10edata[, -1])
  l10pdata <- as.pepData(e_data = l10edata,
                         f_data = fdata,
                         e_meta = emeta,
                         edata_cname = 'Mass_Tag_ID',
                         fdata_cname = 'SampleID',
                         emeta_cname = 'Protein',
                         data_scale_orig = "log10")
  
  # Fabricate an attribute list that will be used for testing the data_info
  # attributes (except for the data_scale attribute).
  data_info_standard <- list(norm_info = list(is_normalized = FALSE),
                             num_edata = 150,
                             num_miss_obs = 341,
                             prop_missing = attr(pdata,
                                                 'data_info')$prop_missing,
                             num_samps = 12,
                             data_types = NULL)
  
  # Run tests on the transformed data ------------------------------------------
  
  # Check for an error when the omicsData argument is not the correct class.
  expect_error(edata_transform(omicsData = pdata$e_data,
                               data_scale = 'abundance'),
               paste("omicsData must be of class 'pepData', 'proData',",
                     "'metabData', 'lipidData', or 'nmrData'",
                     sep = ' '))
  
  # Test for an error when the data_scale argument is not one of the specified
  # options.
  expect_error(edata_transform(omicsData = pdata,
                               data_scale = 'aaabundance'),
               paste("aaabundance is not a valid option for 'data_scale'.",
                     "See details of as.pepData for specifics.",
                     sep = ' '))
  
  # Test for errors when the selected data_scale is also the current data scale.
  expect_error(edata_transform(omicsData = pdata,
                               data_scale = 'abundance'),
               "Data is already on abundance scale.")
  expect_error(edata_transform(omicsData = lpdata,
                               data_scale = 'log'),
               "Data is already on log scale.")
  expect_error(edata_transform(omicsData = l2pdata,
                               data_scale = 'log2'),
               "Data is already on log2 scale.")
  expect_error(edata_transform(omicsData = l10pdata,
                               data_scale = 'log10'),
               "Data is already on log10 scale.")
  
  # Transmogrify the data and check results ------------------------------------
  
  # Test from abundance to log, log2, and log10 ---------------
  
  # Generate new object on the log scale.
  transmuted <- edata_transform(omicsData = pdata,
                                data_scale = 'log')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "abundance")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_identical(transmuted$e_data,
                   lpdata$e_data)
  
  # Forge a new object on the log base 2 scale.
  transmuted <- edata_transform(omicsData = pdata,
                                data_scale = 'log2')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "abundance")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log2")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_identical(transmuted$e_data,
                   l2pdata$e_data)
  
  # Manufacture a new object on the log base 10 scale.
  transmuted <- edata_transform(omicsData = pdata,
                                data_scale = 'log10')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "abundance")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log10")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_identical(transmuted$e_data,
                   l10pdata$e_data)
  
  # Ensure the remaining attributes have not changed.
  expect_equal(attr(transmuted, 'cnames'),
               attr(pdata, 'cnames'))
  expect_equal(attr(transmuted, 'check.names'),
               attr(pdata, 'check.names'))
  expect_equal(attr(transmuted, 'meta_info'),
               attr(pdata, 'meta_info'))
  expect_equal(attr(transmuted, 'group_DF'),
               attr(pdata, 'group_DF'))
  expect_equal(attr(transmuted, 'filters'),
               attr(pdata, 'filters'))
  
  # Test from log to abundance, log2, and log10 ---------------
  
  # Generate new object on the abundance scale.
  transmuted <- edata_transform(omicsData = lpdata,
                                data_scale = 'abundance')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "abundance")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               pdata$e_data)
  
  # Forge a new object on the log base 2 scale.
  transmuted <- edata_transform(omicsData = lpdata,
                                data_scale = 'log2')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log2")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               l2pdata$e_data)
  
  # Manufacture a new object on the log base 10 scale.
  transmuted <- edata_transform(omicsData = lpdata,
                                data_scale = 'log10')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log10")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               l10pdata$e_data)
  
  # Test from log2 to abundance, log, and log10 ---------------
  
  # Generate new object on the abundance scale.
  transmuted <- edata_transform(omicsData = l2pdata,
                                data_scale = 'abundance')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log2")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "abundance")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               pdata$e_data)
  
  # Forge a new object on the log scale.
  transmuted <- edata_transform(omicsData = l2pdata,
                                data_scale = 'log')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log2")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               lpdata$e_data)
  
  # Manufacture a new object on the log base 10 scale.
  transmuted <- edata_transform(omicsData = l2pdata,
                                data_scale = 'log10')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log2")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log10")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               l10pdata$e_data)
  
  # Test from log10 to abundance, log, and log2 ---------------
  
  # Generate new object on the abundance scale.
  transmuted <- edata_transform(omicsData = l10pdata,
                                data_scale = 'abundance')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log10")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "abundance")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               pdata$e_data)
  
  # Forge a new object on the log scale.
  transmuted <- edata_transform(omicsData = l10pdata,
                                data_scale = 'log')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log10")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               lpdata$e_data)
  
  # Manufacture a new object on the log base 2 scale.
  transmuted <- edata_transform(omicsData = l10pdata,
                                data_scale = 'log2')
  
  # Ensure the data_scale attributes are correct.
  expect_equal(attr(transmuted, "data_info")$data_scale_orig,
               "log10")
  expect_equal(attr(transmuted, 'data_info')$data_scale,
               "log2")
  
  # Check that the remaining data_info elements have not changed.
  expect_equal(attr(transmuted, 'data_info')[-c(1, 2)],
               data_info_standard)
  
  # Make sure the values of the e_data data frame match the standard.
  expect_equal(transmuted$e_data,
               l2pdata$e_data)
  
})
