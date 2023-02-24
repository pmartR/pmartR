context('class: pepData')

test_that('as.pepData returns the correct data frame and attributes',{

  # Load the reduced peptide data frames ---------------------------------------

  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))

  # Run as.pepData with agreeable data frames ----------------------------------

  # Construct a pepData object with the edata, fdata, and emeta data frames.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')

  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(pdata$e_data),
               c(150, 13))
  expect_equal(dim(pdata$f_data),
               c(12, 2))
  expect_equal(dim(pdata$e_meta),
               c(150, 4))

  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(pdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(pdata, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = "Protein",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(pdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata$e_data)),
         prop_missing = (sum(is.na(pdata$e_data)) /
                           prod(dim(pdata$e_data[, -1]))),
         num_samps = ncol(pdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(pdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pdata$e_meta$Protein)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(pdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(pdata, "pepData")
  
  # Produce a pepData object only with the edata and fdata data frames.
  expect_message(pdata <- as.pepData(e_data = edata,
                                     f_data = fdata,
                                     e_meta = NULL,
                                     edata_cname = 'Mass_Tag_ID',
                                     fdata_cname = 'SampleID',
                                     emeta_cname = 'Test'),
                 "emeta_cname set to NULL, no e_meta object was provided.")
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(pdata$e_data),
               c(150, 13))
  expect_equal(dim(pdata$f_data),
               c(12, 2))
  expect_null(pdata$e_meta)
  
  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(pdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(pdata, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(pdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata$e_data)),
         prop_missing = (sum(is.na(pdata$e_data)) /
                           prod(dim(pdata$e_data[, -1]))),
         num_samps = ncol(pdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(pdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(pdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(pdata, "pepData")

  # Run as.pepData with disagreeable data frames -------------------------------

  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(as.pepData(e_data = data.frame(edata,
                                              Mock4 = edata[, 10]),
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein'),
               "1 samples from e_data not found in f_data")

  # Check for an error when e_data has more rows than e_meta.
  expect_error(as.pepData(e_data = edata,
                          f_data = fdata,
                          e_meta = emeta[1:138, ],
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein'),
               "Not all peptides in e_data are present in e_meta")

  # Create an f_data object with an extra row.
  fdata_1 <- data.frame(SampleID = c(paste0('Infection', 1:9),
                                     paste0('Mock', 1:4)),
                        Condition = c(rep('Infection', 9),
                                      rep('Mock', 4)))

  # Create a pepData object and check for a warning when the f_data object has
  # an extra row.
  expect_warning(pdata <- as.pepData(e_data = edata,
                                     f_data = fdata_1,
                                     edata_cname = 'Mass_Tag_ID',
                                     fdata_cname = 'SampleID'),
                 paste("Extra samples were found in f_data that were not in",
                       "e_data. These have been removed from f_data.",
                       sep = ' '))

  # Check the dimensions of the e_data and f_data data frames.
  expect_equal(dim(pdata$e_data),
               c(150, 13))
  expect_equal(dim(pdata$f_data),
               c(12, 2))
  expect_null(pdata$e_meta)
  
  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(pdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(pdata, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(pdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata$e_data)),
         prop_missing = (sum(is.na(pdata$e_data)) /
                           prod(dim(pdata$e_data[, -1]))),
         num_samps = ncol(pdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(pdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(pdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(pdata, "pepData")

  # Generate a pepData object and check for a warning when the e_meta object has
  # more rows than e_data.
  expect_warning(pdata <- as.pepData(e_data = edata[1:117, ],
                                     f_data = fdata,
                                     e_meta = emeta,
                                     edata_cname = 'Mass_Tag_ID',
                                     fdata_cname = 'SampleID',
                                     emeta_cname = 'Protein'),
                 paste('Extra peptides were found in e_meta that were not in',
                       'e_data. These have been removed from e_meta.',
                       sep = ' '))

  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(dim(pdata$e_data),
               c(117, 13))
  expect_equal(dim(pdata$f_data),
               c(12, 2))
  expect_equal(dim(pdata$e_meta),
               c(117, 4))
  
  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(pdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(pdata, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = "Protein",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(pdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata$e_data)),
         prop_missing = (sum(is.na(pdata$e_data)) /
                           prod(dim(pdata$e_data[, -1]))),
         num_samps = ncol(pdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(pdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pdata$e_meta$Protein)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(pdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(pdata, "pepData")

  # Check the technical replicates column for correct structure.
  expect_error(as.pepData(e_data = edata,
                          f_data = data.frame(fdata,
                                              tReps = fdata[, 1]),
                          e_meta = emeta,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein',
                          techrep_cname = 'tReps'),
               paste('Specified technical replicate column had a unique value',
                     'for each row.  Values should specify groups of technical',
                     'replicates belonging to a biological sample.',
                     sep = ' '))

  set.seed(5)

  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(edata,
                   edata[sample(1:150, 8), ])

  # Create a pepData object with some of the rows of e_data repeated.
  pdata <- as.pepData(e_data = edata_1,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')

  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(pdata$e_data),
               c(150, 13))
  expect_equal(dim(pdata$f_data),
               c(12, 2))
  expect_equal(dim(pdata$e_meta),
               c(150, 4))

  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(pdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(pdata, "cnames"),
               list(edata_cname = "Mass_Tag_ID",
                    emeta_cname = "Protein",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(pdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata$e_data)),
         prop_missing = (sum(is.na(pdata$e_data)) /
                           prod(dim(pdata$e_data[, -1]))),
         num_samps = ncol(pdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(pdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(pdata$e_meta$Protein)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(pdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(pdata, "pepData")

  # Change the values in each of the samples for the repeated IDs.
  edata_1[151:158, 2:13] <- edata_1[151:158, 2:13] * 1.1

  # Forge a pepData object with some of the peptide IDs repeated but the values
  # for each of the repeated peptide IDs is different from the original data.
  expect_error(as.pepData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein'),
               "The 'edata_cname' identifier is non-unique.")
  
  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(as.pepData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = NULL),
               'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')
  
  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(as.pepData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein2'),
               paste('Mapping variable column',
                     'Protein2',
                     'not found in e_meta. See details of as.pepData for',
                     'specifying column names.',
                     sep = ' '))

})
