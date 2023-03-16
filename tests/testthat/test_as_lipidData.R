context('class: lipidData')

test_that('as.lipidData returns the correct data frame and attributes',{
  
  # Load the lipid data frames -------------------------------------------------
  
  load(system.file('testdata',
                   'lipidData.RData',
                   package = 'pmartR'))
  
  # Run as.lipidData with agreeable data frames --------------------------------
  
  # Produce a lipidData object with the edata, fdata, and emeta data frames.
  ldata <- as.lipidData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'LipidCommonName',
                        fdata_cname = 'Sample_Name',
                        emeta_cname = 'LipidClass')
  
  # Check high level structure
  expect_equal(names(ldata), c("e_data", "f_data", "e_meta"))
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(ldata$e_data),
               c(146, 12))
  expect_equal(dim(ldata$f_data),
               c(11, 2))
  expect_equal(dim(ldata$e_meta),
               c(146, 2))
  
  # Confirm the correct attributes are present in the lipidData object.
  expect_equal(names(attributes(ldata)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(ldata, "cnames"),
               list(edata_cname = "LipidCommonName",
                    emeta_cname = "LipidClass",
                    fdata_cname = "Sample_Name",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(ldata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata$e_data)),
         prop_missing = (sum(is.na(ldata$e_data)) /
                           prod(dim(ldata$e_data[, -1]))),
         num_samps = ncol(ldata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(ldata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(ldata$e_meta$LipidClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(ldata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(ldata, "lipidData")
  
  # Construct a lipidData object only with the edata and fdata data frames.
  expect_message(ldata <- as.lipidData(e_data = edata,
                                       f_data = fdata,
                                       edata_cname = 'LipidCommonName',
                                       fdata_cname = 'Sample_Name',
                                       emeta_cname = 'Test'),
                 "emeta_cname set to NULL, no e_meta object was provided.")
  
  # Check high level structure
  expect_equal(names(ldata), c("e_data", "f_data", "e_meta"))
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(ldata$e_data),
               c(146, 12))
  expect_equal(dim(ldata$f_data),
               c(11, 2))
  expect_null(ldata$e_meta)
  
  # Confirm the correct attributes are present in the lipidData object.
  expect_equal(names(attributes(ldata)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(ldata, "cnames"),
               list(edata_cname = "LipidCommonName",
                    emeta_cname = NULL,
                    fdata_cname = "Sample_Name",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(ldata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata$e_data)),
         prop_missing = (sum(is.na(ldata$e_data)) /
                           prod(dim(ldata$e_data[, -1]))),
         num_samps = ncol(ldata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(ldata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(ldata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(ldata, "lipidData")
  
  # Run as.lipidData with disagreeable data frames -----------------------------
  
  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(as.lipidData(e_data = data.frame(edata,
                                              Mock4 = edata[, 10]),
                            f_data = fdata,
                            edata_cname = 'LipidCommonName',
                            fdata_cname = 'Sample_Name'),
               "1 samples from e_data not found in f_data")
  
  # Verify an error is thrown when edata has more rows than emeta.
  expect_error(as.lipidData(e_data = edata,
                            f_data = fdata,
                            e_meta = emeta[1:130, ],
                            edata_cname = 'LipidCommonName',
                            fdata_cname = 'Sample_Name',
                            emeta_cname = 'LipidClass'))
  
  # Forge an f_data object with an extra row.
  fdata_1 <- data.frame(Sample_Name = c(paste0('Mock', 1:3),
                                        paste0('Infection', 1:10)),
                        Condition = c(rep('Mock', 3),
                                      rep('Infection', 10)))
  
  # Create a lipidData object and check for a warning when the f_data object has
  # an extra row.
  expect_warning(ldata <- as.lipidData(e_data = edata,
                                       f_data = fdata_1,
                                       edata_cname = 'LipidCommonName',
                                       fdata_cname = 'Sample_Name'),
                 paste("Extra samples were found in f_data that were not in",
                       "e_data. These have been removed from f_data.",
                       sep = ' '))
  
  # Check high level structure
  expect_equal(names(ldata), c("e_data", "f_data", "e_meta"))
  
  # Confirm the dimensions of the e_data and f_data data frames.
  expect_equal(dim(ldata$e_data),
               c(146, 12))
  expect_equal(dim(ldata$f_data),
               c(11, 2))
  expect_null(ldata$e_meta)
  
  # Confirm the correct attributes are present in the lipidData object.
  expect_equal(names(attributes(ldata)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(ldata, "cnames"),
               list(edata_cname = "LipidCommonName",
                    emeta_cname = NULL,
                    fdata_cname = "Sample_Name",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(ldata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata$e_data)),
         prop_missing = (sum(is.na(ldata$e_data)) /
                           prod(dim(ldata$e_data[, -1]))),
         num_samps = ncol(ldata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(ldata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(ldata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(ldata, "lipidData")
  
  # Generate a lipidData object and check for a warning when the e_meta object
  # has more rows than e_data.
  expect_warning(ldata <- as.lipidData(e_data = edata[1:117, ],
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'LipidCommonName',
                                       fdata_cname = 'Sample_Name',
                                       emeta_cname = 'LipidClass'),
                 paste('Extra lipids were found in e_meta that were not in',
                       'e_data. These have been removed from e_meta.',
                       sep = ' '))
  
  # Check high level structure
  expect_equal(names(ldata), c("e_data", "f_data", "e_meta"))
  
  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(dim(ldata$e_data),
               c(117, 12))
  expect_equal(dim(ldata$f_data),
               c(11, 2))
  expect_equal(dim(ldata$e_meta),
               c(117, 2))
  
  # Confirm the correct attributes are present in the lipidData object.
  expect_equal(names(attributes(ldata)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(ldata, "cnames"),
               list(edata_cname = "LipidCommonName",
                    emeta_cname = "LipidClass",
                    fdata_cname = "Sample_Name",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(ldata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata$e_data)),
         prop_missing = (sum(is.na(ldata$e_data)) /
                           prod(dim(ldata$e_data[, -1]))),
         num_samps = ncol(ldata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(ldata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(ldata$e_meta$LipidClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(ldata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(ldata, "lipidData")
  
  # Check the technical replicates column for correct structure.
  expect_error(as.lipidData(e_data = edata,
                          f_data = data.frame(fdata,
                                              tReps = fdata[, 1]),
                          edata_cname = 'LipidCommonName',
                          fdata_cname = 'Sample_Name',
                          techrep_cname = 'tReps'),
               paste('Specified technical replicate column had a unique value',
                     'for each row.  Values should specify groups of technical',
                     'replicates belonging to a biological sample.',
                     sep = ' '))
  
  set.seed(3)
  
  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(edata,
                   edata[sample(1:146, 8), ])
  
  # Create a lipidData object with some of the rows of e_data repeated.
  ldata <- as.lipidData(e_data = edata_1,
                      f_data = fdata,
                      edata_cname = 'LipidCommonName',
                      fdata_cname = 'Sample_Name')
  
  # Check high level structure
  expect_equal(names(ldata), c("e_data", "f_data", "e_meta"))
  
  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(ldata$e_data),
               c(146, 12))
  expect_equal(dim(ldata$f_data),
               c(11, 2))
  expect_null(ldata$e_meta)
  
  # Confirm the correct attributes are present in the lipidData object.
  expect_equal(names(attributes(ldata)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(ldata, "cnames"),
               list(edata_cname = "LipidCommonName",
                    emeta_cname = NULL,
                    fdata_cname = "Sample_Name",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(ldata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata$e_data)),
         prop_missing = (sum(is.na(ldata$e_data)) /
                           prod(dim(ldata$e_data[, -1]))),
         num_samps = ncol(ldata$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(ldata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(ldata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(ldata, "lipidData")
  
  # Change the values in each of the samples for the repeated IDs.
  edata_1[147:154, 2:12] <- edata_1[147:154, 2:12] * 1.1
  
  # Forge a lipidData object with some of the peptide IDs but the values for
  # each of the samples is different from the original data frame.
  expect_error(as.lipidData(e_data = edata_1,
                            f_data = fdata,
                            edata_cname = 'LipidCommonName',
                            fdata_cname = 'Sample_Name'),
               "The 'edata_cname' identifier is non-unique.")
  
  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(as.lipidData(e_data = edata_1,
                            f_data = fdata,
                            e_meta = emeta,
                            edata_cname = 'LipidCommonName',
                            fdata_cname = 'Sample_Name',
                            emeta_cname = NULL),
               'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')
  
  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(as.lipidData(e_data = edata_1,
                          f_data = fdata,
                          e_meta = emeta,
                          edata_cname = 'LipidCommonName',
                          fdata_cname = 'Sample_Name',
                          emeta_cname = 'Lipid'),
               paste('Mapping variable column',
                     'Lipid',
                     'not found in e_meta. See details of as.lipidData for',
                     'specifying column names.',
                     sep = ' '))
  
})
