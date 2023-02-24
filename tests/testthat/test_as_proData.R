context('class: proData')

test_that('as.proData returns the correct data frame and attributes',{
  
  # Load the reduced protein data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_prdata.RData',
                   package = 'pmartR'))
  
  # Run as.proData with agreeable data frames ----------------------------------
  
  # Produce a proData object with the edata, fdata, and emeta data frames.
  prdata <- as.proData(e_data = edata,
                       f_data = fdata,
                       e_meta = emeta,
                       edata_cname = 'Reference',
                       fdata_cname = 'SampleID',
                       emeta_cname = 'PClass')
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(prdata$e_data),
               c(150, 12))
  expect_equal(dim(prdata$f_data),
               c(11, 2))
  expect_equal(dim(prdata$e_meta),
               c(150, 2))
  
  # Confirm the correct attributes are present in the proData object.
  expect_equal(names(attributes(prdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "pro_quant_info", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(prdata, "cnames"),
               list(edata_cname = "Reference",
                    emeta_cname = "PClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(prdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata$e_data)),
         prop_missing = (sum(is.na(prdata$e_data)) /
                           prod(dim(prdata$e_data[, -1]))),
         num_samps = ncol(prdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(prdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(prdata$e_meta$PClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(prdata, "filters"), list())
  
  # Make sure the protein quantitation attribute is NA.
  expect_true(is.na(attr(prdata, "pro_quant_info")$method))
  
  # Ensure the omicsData object is classy.
  expect_s3_class(prdata, "proData")
  
  # Construct a proData object only with the edata and fdata data frames.
  expect_message(prdata <- as.proData(e_data = edata,
                                      f_data = fdata,
                                      edata_cname = 'Reference',
                                      fdata_cname = 'SampleID',
                                      emeta_cname = 'Test'),
                 "emeta_cname set to NULL, no e_meta object was provided.")
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(prdata$e_data),
               c(150, 12))
  expect_equal(dim(prdata$f_data),
               c(11, 2))
  expect_null(prdata$e_meta)
  
  # Confirm the correct attributes are present in the proData object.
  expect_equal(names(attributes(prdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "pro_quant_info", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(prdata, "cnames"),
               list(edata_cname = "Reference",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(prdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata$e_data)),
         prop_missing = (sum(is.na(prdata$e_data)) /
                           prod(dim(prdata$e_data[, -1]))),
         num_samps = ncol(prdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(prdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(prdata, "filters"), list())
  
  # Make sure the protein quantitation attribute is NA.
  expect_true(is.na(attr(prdata, "pro_quant_info")$method))
  
  # Ensure the omicsData object is classy.
  expect_s3_class(prdata, "proData")
  
  # Run as.proData with disagreeable data frames -------------------------------
  
  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(as.proData(e_data = data.frame(edata,
                                              Mock4 = edata[, 10]),
                          f_data = fdata,
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID'),
               "1 samples from e_data not found in f_data")
  
  # Verify an error is thrown when edata has more rows than emeta.
  expect_error(as.proData(e_data = edata,
                          f_data = fdata,
                          e_meta = emeta[1:130, ],
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'PClass'))
  
  # Forge an f_data object with an extra row.
  fdata_1 <- data.frame(SampleID = c(paste0('Mock', 1:3),
                                     paste0('Infection', 1:9)),
                        Condition = c(rep('Mock', 3),
                                      rep('Infection', 9)))
  
  # Create a proData object and check for a warning when the f_data object has
  # an extra row.
  expect_warning(prdata <- as.proData(e_data = edata,
                                      f_data = fdata_1,
                                      edata_cname = 'Reference',
                                      fdata_cname = 'SampleID'),
                 paste("Extra samples were found in f_data that were not in",
                       "e_data. These have been removed from f_data.",
                       sep = ' '))
  
  # Confirm the dimensions of the e_data and f_data data frames.
  expect_equal(dim(prdata$e_data),
               c(150, 12))
  expect_equal(dim(prdata$f_data),
               c(11, 2))
  expect_null(prdata$e_meta)
  
  # Confirm the correct attributes are present in the proData object.
  expect_equal(names(attributes(prdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "pro_quant_info", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(prdata, "cnames"),
               list(edata_cname = "Reference",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(prdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata$e_data)),
         prop_missing = (sum(is.na(prdata$e_data)) /
                           prod(dim(prdata$e_data[, -1]))),
         num_samps = ncol(prdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(prdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(prdata, "filters"), list())
  
  # Make sure the protein quantitation attribute is NA.
  expect_true(is.na(attr(prdata, "pro_quant_info")$method))
  
  # Ensure the omicsData object is classy.
  expect_s3_class(prdata, "proData")
  
  # Generate a proData object and check for a warning when the e_meta object
  # has more rows than e_data.
  expect_warning(prdata <- as.proData(e_data = edata[1:131, ],
                                      f_data = fdata,
                                      e_meta = emeta,
                                      edata_cname = 'Reference',
                                      fdata_cname = 'SampleID',
                                      emeta_cname = 'PClass'),
                 paste('Extra proteins were found in e_meta that were not in',
                       'e_data. These have been removed from e_meta.',
                       sep = ' '))
  
  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(dim(prdata$e_data),
               c(131, 12))
  expect_equal(dim(prdata$f_data),
               c(11, 2))
  expect_equal(dim(prdata$e_meta),
               c(131, 2))
  
  # Confirm the correct attributes are present in the proData object.
  expect_equal(names(attributes(prdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "pro_quant_info", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(prdata, "cnames"),
               list(edata_cname = "Reference",
                    emeta_cname = "PClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(prdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata$e_data)),
         prop_missing = (sum(is.na(prdata$e_data)) /
                           prod(dim(prdata$e_data[, -1]))),
         num_samps = ncol(prdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(prdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(prdata$e_meta$PClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(prdata, "filters"), list())
  
  # Make sure the protein quantitation attribute is NA.
  expect_true(is.na(attr(prdata, "pro_quant_info")$method))
  
  # Ensure the omicsData object is classy.
  expect_s3_class(prdata, "proData")
  
  # Check the technical replicates column for correct structure.
  expect_error(as.proData(e_data = edata,
                          f_data = data.frame(fdata,
                                              tReps = fdata[, 1]),
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID',
                          techrep_cname = 'tReps'),
               paste('Specified technical replicate column had a unique value',
                     'for each row.  Values should specify groups of technical',
                     'replicates belonging to a biological sample.',
                     sep = ' '))
  
  set.seed(3)
  
  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(edata,
                   edata[sample(1:150, 8), ])
  
  # Create a proData object with some of the rows of e_data repeated.
  prdata <- as.proData(e_data = edata_1,
                       f_data = fdata,
                       edata_cname = 'Reference',
                       fdata_cname = 'SampleID')
  
  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(prdata$e_data),
               c(150, 12))
  expect_equal(dim(prdata$f_data),
               c(11, 2))
  expect_null(prdata$e_meta)
  
  # Confirm the correct attributes are present in the proData object.
  expect_equal(names(attributes(prdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "pro_quant_info", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(prdata, "cnames"),
               list(edata_cname = "Reference",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(prdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata$e_data)),
         prop_missing = (sum(is.na(prdata$e_data)) /
                           prod(dim(prdata$e_data[, -1]))),
         num_samps = ncol(prdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(prdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(prdata, "filters"), list())
  
  # Make sure the protein quantitation attribute is NA.
  expect_true(is.na(attr(prdata, "pro_quant_info")$method))
  
  # Ensure the omicsData object is classy.
  expect_s3_class(prdata, "proData")
  
  # Change the values in each of the samples for the repeated IDs.
  edata_1[151:158, 2:12] <- edata_1[151:158, 2:12] * 1.1
  
  # Forge a proData object with some of the peptide IDs but the values for
  # each of the samples is different from the original data frame.
  expect_error(as.proData(e_data = edata_1,
                          f_data = fdata,
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID'),
               "The 'edata_cname' identifier is non-unique.")
  
  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(as.proData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID',
                          emeta_cname = NULL),
               'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')
  
  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(as.proData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Reference',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'PClass2'),
               paste('Mapping variable column',
                     'PClass2',
                     'not found in e_meta. See details of as.proData for',
                     'specifying column names.',
                     sep = ' '))
  
})
