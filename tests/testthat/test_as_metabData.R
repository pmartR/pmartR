context('class: metabData')

test_that('as.metabData returns the correct data frame and attributes',{
  
  # Load the metabolite data frames --------------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
  # Run as.metabData with agreeable data frames --------------------------------
  
  # Produce a metabData object with the edata, fdata, and emeta data frames.
  mdata <- as.metabData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'MClass',
                        data_scale_orig = "abundance")
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(mdata$e_data),
               c(80, 13))
  expect_equal(dim(mdata$f_data),
               c(12, 2))
  expect_equal(dim(mdata$e_meta),
               c(80, 2))
  
  # Confirm the correct attributes are present in the metabData object.
  expect_equal(names(attributes(mdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(mdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = "MClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(mdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata$e_data)),
         prop_missing = (sum(is.na(mdata$e_data)) /
                           prod(dim(mdata$e_data[, -1]))),
         num_samps = ncol(mdata$e_data[, -1]),
         data_types = NULL)
  )
  
  # Check the checkers.
  expect_true(attr(mdata, "check.names"))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(mdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(mdata$e_meta$MClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(mdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(mdata, "metabData")
  
  # Construct a metabData object only with the edata and fdata data frames.
  expect_message(mdata <- as.metabData(e_data = edata,
                                       f_data = fdata,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'Test',
                                       data_scale_orig = "abundance"),
                 "emeta_cname set to NULL, no e_meta object was provided.")
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(mdata$e_data),
               c(80, 13))
  expect_equal(dim(mdata$f_data),
               c(12, 2))
  expect_null(mdata$e_meta)
  
  # Confirm the correct attributes are present in the metabData object.
  expect_equal(names(attributes(mdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(mdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(mdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata$e_data)),
         prop_missing = (sum(is.na(mdata$e_data)) /
                           prod(dim(mdata$e_data[, -1]))),
         num_samps = ncol(mdata$e_data[, -1]),
         data_types = NULL)
  )
  
  # Check the checkers.
  expect_true(attr(mdata, "check.names"))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(mdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(mdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(mdata, "metabData")
  
  # Run as.metabData with disagreeable data frames -----------------------------
  
  # Check for an error when e_data has more columns than f_data has rows.
  expect_error(as.metabData(e_data = data.frame(edata,
                                                Mock4 = edata[, 10]),
                            f_data = fdata,
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            data_scale_orig = "abundance"),
               "1 samples from e_data not found in f_data")
  
  # Verify an error is thrown when edata has more rows than emeta.
  expect_error(as.metabData(e_data = edata,
                            f_data = fdata,
                            e_meta = emeta[1:72, ],
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            emeta_cname = 'MClass',
                            data_scale_orig = "abundance"))
  
  # Forge an f_data object with an extra row.
  fdata_1 <- data.frame(SampleID = c(paste0('Mock', 1:3),
                                     paste0('Infection', 1:10)),
                        Condition = c(rep('Mock', 3),
                                      rep('Infection', 10)))
  
  # Create a metabData object and check for a warning when the f_data object has
  # an extra row.
  expect_warning(mdata <- as.metabData(e_data = edata,
                                       f_data = fdata_1,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       data_scale_orig = "abundance"),
                 paste("Extra samples were found in f_data that were not in",
                       "e_data. These have been removed from f_data.",
                       sep = ' '))
  
  # Confirm the dimensions of the e_data and f_data data frames.
  expect_equal(dim(mdata$e_data),
               c(80, 13))
  expect_equal(dim(mdata$f_data),
               c(12, 2))
  expect_null(mdata$e_meta)
  
  # Confirm the correct attributes are present in the metabData object.
  expect_equal(names(attributes(mdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(mdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(mdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata$e_data)),
         prop_missing = (sum(is.na(mdata$e_data)) /
                           prod(dim(mdata$e_data[, -1]))),
         num_samps = ncol(mdata$e_data[, -1]),
         data_types = NULL)
  )
  
  # Check the checkers.
  expect_true(attr(mdata, "check.names"))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(mdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(mdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(mdata, "metabData")
  
  # Generate a metabData object and check for a warning when the e_meta object
  # has more rows than e_data.
  expect_warning(mdata <- as.metabData(e_data = edata[1:72, ],
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'MClass',
                                       data_scale_orig = "abundance"),
                 paste('Extra metabolites were found in e_meta that were not in',
                       'e_data. These have been removed from e_meta.',
                       sep = ' '))
  
  # Confirm the dimensions of the e_data and e_meta data frames.
  expect_equal(dim(mdata$e_data),
               c(72, 13))
  expect_equal(dim(mdata$f_data),
               c(12, 2))
  expect_equal(dim(mdata$e_meta),
               c(72, 2))
  
  # Confirm the correct attributes are present in the metabData object.
  expect_equal(names(attributes(mdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(mdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = "MClass",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(mdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata$e_data)),
         prop_missing = (sum(is.na(mdata$e_data)) /
                           prod(dim(mdata$e_data[, -1]))),
         num_samps = ncol(mdata$e_data[, -1]),
         data_types = NULL)
  )
  
  # Check the checkers.
  expect_true(attr(mdata, "check.names"))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(mdata, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(mdata$e_meta$MClass)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(mdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(mdata, "metabData")
  
  # Check the technical replicates column for correct structure.
  expect_error(as.metabData(e_data = edata,
                            f_data = data.frame(fdata,
                                                tReps = fdata[, 1]),
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            techrep_cname = 'tReps',
                            data_scale_orig = "abundance"),
               paste('Specified technical replicate column had a unique value',
                     'for each row.  Values should specify groups of technical',
                     'replicates belonging to a biological sample.',
                     sep = ' '))
  
  set.seed(3)
  
  # Fabricate an e_data object with some of the peptide IDs repeated.
  edata_1 <- rbind(edata,
                   edata[sample(1:80, 8), ])
  
  # Create a metabData object with some of the rows of e_data repeated.
  mdata <- as.metabData(e_data = edata_1,
                        f_data = fdata,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        data_scale_orig = "abundance")
  
  # Verify that the returned data frames are the correct dimension.
  expect_equal(dim(mdata$e_data),
               c(80, 13))
  expect_equal(dim(mdata$f_data),
               c(12, 2))
  expect_null(mdata$e_meta)
  
  # Confirm the correct attributes are present in the metabData object.
  expect_equal(names(attributes(mdata)),
               c("names", "cnames", "data_info", "check.names", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(mdata, "cnames"),
               list(edata_cname = "Metabolite",
                    emeta_cname = NULL,
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(mdata, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata$e_data)),
         prop_missing = (sum(is.na(mdata$e_data)) /
                           prod(dim(mdata$e_data[, -1]))),
         num_samps = ncol(mdata$e_data[, -1]),
         data_types = NULL)
  )
  
  # Check the checkers.
  expect_true(attr(mdata, "check.names"))
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(mdata, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(mdata, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(mdata, "metabData")
  
  # Change the values in each of the samples for the repeated IDs.
  edata_1[81:88, 2:13] <- edata_1[81:88, 2:13] * 1.1
  
  # Forge a metabData object with some of the peptide IDs but the values for
  # each of the samples is different from the original data frame.
  expect_error(as.metabData(e_data = edata_1,
                            f_data = fdata,
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            data_scale_orig = "abundance"),
               "The 'edata_cname' identifier is non-unique.")
  
  # Check for an error when e_meta is non-null but emeta_cname is null.
  expect_error(as.metabData(e_data = edata_1,
                            f_data = fdata,
                            e_meta = emeta,
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            emeta_cname = NULL,
                            data_scale_orig = "abundance"),
               'Since e_meta is non-NULL, emeta_cname must also be non-NULL.')
  
  # Verify there is an error when emeta_cname is not a column name in e_meta.
  expect_error(as.metabData(e_data = edata_1,
                            f_data = fdata,
                            e_meta = emeta,
                            edata_cname = 'Metabolite',
                            fdata_cname = 'SampleID',
                            emeta_cname = 'MetaboliteClass',
                            data_scale_orig = "abundance"),
               paste('Mapping variable column',
                     'MetaboliteClass',
                     'not found in e_meta. See details of as.metabData for',
                     'specifying column names.',
                     sep = ' '))
  
})
