context('class: nmrData')

test_that('as.nmrData returns the correct data frame and attributes',{
  
  # Load the nmr data frames ---------------------------------------------------
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Run as.nmrData with agreeable data frames ----------------------------------
  
  # Produce a nmrData object with the edata, fdata, and emeta data frames.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'nmrClass')
  
  # Check high level structure
  expect_equal(names(nmrdata), c("e_data", "f_data", "e_meta"))
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(nmrdata$e_data),
               c(38, 42))
  expect_equal(dim(nmrdata$f_data),
               c(41, 5))
  expect_equal(dim(nmrdata$e_meta),
               c(38, 3))
  
  # Confirm the correct attributes are present in the nmrData object.
  expect_equal(names(attributes(nmrdata)),
               c("names", "cnames", "data_info", "nmr_info",
                 "meta_info", "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(nmrdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = "nmrClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(nmrdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata$e_data)),
         prop_missing = (sum(is.na(nmrdata$e_data)) /
                           prod(dim(nmrdata$e_data[, -1]))),
         num_samps = ncol(nmrdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Sleuth around the nmr attributes.
  expect_equal(attr(nmrdata, "nmr_info"),
               list(metabolite_name = NA,
                    sample_property_cname = NA,
                    norm_info = list(is_normalized = FALSE,
                                     backtransform = NA)))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(nmrdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(nmrdata$e_meta$nmrClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(nmrdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(nmrdata, "nmrData")
  
  # Construct a nmrData object only with the edata and fdata data frames.
  expect_message(nmrdata <- as.nmrData(e_data = edata,
                                       f_data = fdata,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'Test'),
                 "emeta_cname set to NULL, no e_meta object was provided.")
  
  # Check high level structure
  expect_equal(names(nmrdata), c("e_data", "f_data", "e_meta"))
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(nmrdata$e_data),
               c(38, 42))
  expect_equal(dim(nmrdata$f_data),
               c(41, 5))
  expect_null(nmrdata$e_meta)
  
  # Confirm the correct attributes are present in the nmrData object.
  expect_equal(names(attributes(nmrdata)),
               c("names", "cnames", "data_info", "nmr_info",
                 "meta_info", "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(nmrdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(nmrdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata$e_data)),
         prop_missing = (sum(is.na(nmrdata$e_data)) /
                           prod(dim(nmrdata$e_data[, -1]))),
         num_samps = ncol(nmrdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Sleuth around the nmr attributes.
  expect_equal(attr(nmrdata, "nmr_info"),
               list(metabolite_name = NA,
                    sample_property_cname = NA,
                    norm_info = list(is_normalized = FALSE,
                                     backtransform = NA)))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(nmrdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(nmrdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(nmrdata, "nmrData")
  
  # Run as.nmrData with disagreeable data frames -------------------------------
  
  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(as.nmrData(e_data = data.frame(edata,
                                              `F5-007` = edata[, 10], check.names = FALSE),
                          f_data = fdata,
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID'),
               "1 samples from e_data not found in f_data")
  
  # Verify an error is thrown when edata has more rows than emeta.
  expect_error(as.nmrData(e_data = edata,
                          f_data = fdata,
                          e_meta = emeta[1:15, ],
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'nmrClass'))
  
  # Forge an f_data object with an extra row.
  fdata_1 <- rbind(fdata,
                   c('F5-007', 'Shadow', 'J55', 53))
  
  # Create a nmrData object and check for a warning when the f_data object has
  # an extra row.
  expect_warning(nmrdata <- as.nmrData(e_data = edata,
                                       f_data = fdata_1,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID'),
                 paste("Extra samples were found in f_data that were not in",
                       "e_data. These have been removed from f_data.",
                       sep = ' '))
  
  # Check high level structure
  expect_equal(names(nmrdata), c("e_data", "f_data", "e_meta"))
  
  # Confirm the dimensions of the e_data and f_data data frames.
  expect_equal(dim(nmrdata$e_data),
               c(38, 42))
  expect_equal(dim(nmrdata$f_data),
               c(41, 5))
  expect_null(nmrdata$e_meta)
  
  # Confirm the correct attributes are present in the nmrData object.
  expect_equal(names(attributes(nmrdata)),
               c("names", "cnames", "data_info", "nmr_info",
                 "meta_info", "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(nmrdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(nmrdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata$e_data)),
         prop_missing = (sum(is.na(nmrdata$e_data)) /
                           prod(dim(nmrdata$e_data[, -1]))),
         num_samps = ncol(nmrdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Sleuth around the nmr attributes.
  expect_equal(attr(nmrdata, "nmr_info"),
               list(metabolite_name = NA,
                    sample_property_cname = NA,
                    norm_info = list(is_normalized = FALSE,
                                     backtransform = NA)))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(nmrdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(nmrdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(nmrdata, "nmrData")
  
  # Generate a nmrData object and check for a warning when the e_meta object
  # has more rows than e_data.
  expect_warning(nmrdata <- as.nmrData(e_data = edata[1:26, ],
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'nmrClass'),
                 paste('Extra metabolites were found in e_meta that were not',
                       'in e_data. These have been removed from e_meta.',
                       sep = ' '))
  
  # Check high level structure
  expect_equal(names(nmrdata), c("e_data", "f_data", "e_meta"))
  
  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(dim(nmrdata$e_data),
               c(26, 42))
  expect_equal(dim(nmrdata$f_data),
               c(41, 5))
  expect_equal(dim(nmrdata$e_meta),
               c(26, 3))
  
  # Confirm the correct attributes are present in the nmrData object.
  expect_equal(names(attributes(nmrdata)),
               c("names", "cnames", "data_info", "nmr_info",
                 "meta_info", "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(nmrdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = "nmrClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(nmrdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata$e_data)),
         prop_missing = (sum(is.na(nmrdata$e_data)) /
                           prod(dim(nmrdata$e_data[, -1]))),
         num_samps = ncol(nmrdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Sleuth around the nmr attributes.
  expect_equal(attr(nmrdata, "nmr_info"),
               list(metabolite_name = NA,
                    sample_property_cname = NA,
                    norm_info = list(is_normalized = FALSE,
                                     backtransform = NA)))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(nmrdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(nmrdata$e_meta$nmrClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(nmrdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(nmrdata, "nmrData")
  
  # Check the technical replicates column for correct structure.
  expect_error(as.nmrData(e_data = edata,
                          f_data = data.frame(fdata,
                                              tReps = fdata[, 1]),
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID',
                          techrep_cname = 'tReps'),
               paste('Specified technical replicate column had a unique value',
                     'for each row.  Values should specify groups of technical',
                     'replicates belonging to a biological sample.',
                     sep = ' '))
  
  set.seed(8)
  
  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(edata,
                   edata[sample(1:38, 8), ])
  
  # Create a nmrData object with some of the rows of e_data repeated.
  nmrdata <- as.nmrData(e_data = edata_1,
                        f_data = fdata,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID')
  
  # Check high level structure
  expect_equal(names(nmrdata), c("e_data", "f_data", "e_meta"))
  
  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(nmrdata$e_data),
               c(38, 42))
  expect_equal(dim(nmrdata$f_data),
               c(41, 5))
  expect_null(nmrdata$e_meta)
  
  # Confirm the correct attributes are present in the nmrData object.
  expect_equal(names(attributes(nmrdata)),
               c("names", "cnames", "data_info", "nmr_info",
                 "meta_info", "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(nmrdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(nmrdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata$e_data)),
         prop_missing = (sum(is.na(nmrdata$e_data)) /
                           prod(dim(nmrdata$e_data[, -1]))),
         num_samps = ncol(nmrdata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Sleuth around the nmr attributes.
  expect_equal(attr(nmrdata, "nmr_info"),
               list(metabolite_name = NA,
                    sample_property_cname = NA,
                    norm_info = list(is_normalized = FALSE,
                                     backtransform = NA)))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(nmrdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(nmrdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(nmrdata, "nmrData")
  
  # Change the values in each of the samples for the repeated IDs.
  edata_1[39:46, 2:42] <- edata_1[39:46, 2:42] * 1.1
  
  # Forge a nmrData object with some of the peptide IDs but the values for
  # each of the samples is different from the original data frame.
  expect_error(as.nmrData(e_data = edata_1,
                          f_data = fdata,
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID'),
               "The 'edata_cname' identifier is non-unique.")
  
  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(as.nmrData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID',
                          emeta_cname = NULL),
               'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')
  
  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(as.nmrData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'Metabolite',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'MetaboliteClass'),
               paste('Mapping variable column',
                     'MetaboliteClass',
                     'not found in e_meta. See details of as.nmrData for',
                     'specifying column names.',
                     sep = ' '))
  
})
