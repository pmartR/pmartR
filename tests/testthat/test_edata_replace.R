context('edata replace')

test_that('edata_replace correctly replaces one value with another',{
  
  # Load: isobaricpepData ------------------------------------------------------
  
  load(system.file('testdata',
                   'little_isodata.RData',
                   package = 'pmartR'))
  
  # Test: isobaricpepData ------------------------------------------------------
  
  # Fabricate an isobaricpepData object.
  isodata <- as.isobaricpepData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Peptide',
                                fdata_cname = 'Sample',
                                emeta_cname = 'Protein',
                                data_scale_orig = "abundance")
  
  # Replace NAs with 0s in isodata.
  expect_message(isodata2 <- edata_replace(omicsData = isodata,
                                           x = NA,
                                           y = 0),
                 "400 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(isodata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(isodata2$e_data[, 1])),
         num_miss_obs = sum(is.na(isodata2$e_data)),
         prop_missing = (sum(is.na(isodata2$e_data)) /
                           prod(dim(isodata2$e_data[, -1]))),
         num_samps = ncol(isodata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(isodata, 'cnames'),
               attr(isodata2, 'cnames'))
  expect_equal(attr(isodata, 'isobaric_info'),
               attr(isodata2, 'isobaric_info'))
  expect_equal(attr(isodata, 'check.names'),
               attr(isodata2, 'check.names'))
  expect_equal(attr(isodata, 'meta_info'),
               attr(isodata2, 'meta_info'))
  expect_equal(attr(isodata, 'group_DF'),
               attr(isodata2, 'group_DF'))
  expect_equal(attr(isodata, 'filters'),
               attr(isodata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(isodata2$e_data == 0),
               400)
  
  # Examine the class of the new isodata object.
  expect_s3_class(isodata2,
                  c('isobaricpepData', 'pepData'))
  
  # Load: lipidData ------------------------------------------------------------
  
  load(system.file('testdata',
                   'lipidData.RData',
                   package = 'pmartR'))
  
  # Test: lipidData ------------------------------------------------------------
  
  # Produce a lipidData object.
  ldata <- as.lipidData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'LipidCommonName',
                        fdata_cname = 'Sample_Name',
                        emeta_cname = 'LipidClass',
                        data_scale_orig = "abundance")
  
  # Replace NAs with 0s in ldata.
  expect_message(ldata2 <- edata_replace(omicsData = ldata,
                                         x = NA,
                                         y = 0),
                 "884 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(ldata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(ldata2$e_data[, 1])),
         num_miss_obs = sum(is.na(ldata2$e_data)),
         prop_missing = (sum(is.na(ldata2$e_data)) /
                           prod(dim(ldata2$e_data[, -1]))),
         num_samps = ncol(ldata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(ldata, 'cnames'),
               attr(ldata2, 'cnames'))
  expect_equal(attr(ldata, 'check.names'),
               attr(ldata2, 'check.names'))
  expect_equal(attr(ldata, 'meta_info'),
               attr(ldata2, 'meta_info'))
  expect_equal(attr(ldata, 'group_DF'),
               attr(ldata2, 'group_DF'))
  expect_equal(attr(ldata, 'filters'),
               attr(ldata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(ldata2$e_data == 0),
               884)

  # Examine the class of the new ldata object.
  expect_s3_class(ldata2,
                  'lipidData')
  
  # Load: metabData ------------------------------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
  # Test: metabData ------------------------------------------------------------
  
  # Forge a metabData object.
  mdata <- as.metabData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'MClass',
                        data_scale_orig = "abundance")
  
  # Replace NAs with 0s in mdata.
  expect_message(mdata2 <- edata_replace(omicsData = mdata,
                                         x = NA,
                                         y = 0),
                 "148 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(mdata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(mdata2$e_data[, 1])),
         num_miss_obs = sum(is.na(mdata2$e_data)),
         prop_missing = (sum(is.na(mdata2$e_data)) /
                           prod(dim(mdata2$e_data[, -1]))),
         num_samps = ncol(mdata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(mdata, 'cnames'),
               attr(mdata2, 'cnames'))
  expect_equal(attr(mdata, 'check.names'),
               attr(mdata2, 'check.names'))
  expect_equal(attr(mdata, 'meta_info'),
               attr(mdata2, 'meta_info'))
  expect_equal(attr(mdata, 'group_DF'),
               attr(mdata2, 'group_DF'))
  expect_equal(attr(mdata, 'filters'),
               attr(mdata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(mdata2$e_data == 0),
               148)
  
  # Examine the class of the new mdata object.
  expect_s3_class(mdata2,
                  'metabData')
  
  # Load: nmrData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Alter the original nmr data frame to contain NAs ---------------------------
  
  set.seed(22)
  
  # Select row indices.
  rows <- sample(1:38, 338, replace = TRUE)
  
  # Select column indices.
  columns <- sample(2:42, 338, replace = TRUE)
  
  # Loop through each row column pair.
  for (e in 1:338) {
    
    # Replace current value in edata with NA.
    edata[rows[[e]], columns[[e]]] <- NA
    
  }
  
  # Test: nmrData --------------------------------------------------------------
  
  # Fabricate an nmrData object.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'nmrClass',
                        check.names = FALSE,
                        data_scale_orig = "abundance")
  
  # Replace NAs with 0s in nmrdata.
  expect_message(nmrdata2 <- edata_replace(omicsData = nmrdata,
                                           x = NA,
                                           y = 0),
                 "302 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(nmrdata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(nmrdata2$e_data[, 1])),
         num_miss_obs = sum(is.na(nmrdata2$e_data)),
         prop_missing = (sum(is.na(nmrdata2$e_data)) /
                           prod(dim(nmrdata2$e_data[, -1]))),
         num_samps = ncol(nmrdata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(nmrdata, 'cnames'),
               attr(nmrdata2, 'cnames'))
  expect_equal(attr(nmrdata, 'nmr_info'),
               attr(nmrdata2, 'nmr_info'))
  expect_equal(attr(nmrdata, 'check.names'),
               attr(nmrdata2, 'check.names'))
  expect_equal(attr(nmrdata, 'meta_info'),
               attr(nmrdata2, 'meta_info'))
  expect_equal(attr(nmrdata, 'group_DF'),
               attr(nmrdata2, 'group_DF'))
  expect_equal(attr(nmrdata, 'filters'),
               attr(nmrdata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(nmrdata2$e_data == 0),
               302)
  
  # Examine the class of the new nmrdata object.
  expect_s3_class(nmrdata,
                  'nmrData')
  
  # Load: pepData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Test: pepData --------------------------------------------------------------
  
  # Construct a pepData object.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein',
                      data_scale_orig = "abundance")
  
  # Replace NAs with 0s in pdata.
  expect_message(pdata2 <- edata_replace(omicsData = pdata,
                                         x = NA,
                                         y = 0),
                 "341 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(pdata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(pdata2$e_data[, 1])),
         num_miss_obs = sum(is.na(pdata2$e_data)),
         prop_missing = (sum(is.na(pdata2$e_data)) /
                           prod(dim(pdata2$e_data[, -1]))),
         num_samps = ncol(pdata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(pdata, 'cnames'),
               attr(pdata2, 'cnames'))
  expect_equal(attr(pdata, 'check.names'),
               attr(pdata2, 'check.names'))
  expect_equal(attr(pdata, 'meta_info'),
               attr(pdata2, 'meta_info'))
  expect_equal(attr(pdata, 'group_DF'),
               attr(pdata2, 'group_DF'))
  expect_equal(attr(pdata, 'filters'),
               attr(pdata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(pdata2$e_data == 0),
               341)
  
  # Examine the class of the new pdata object.
  expect_s3_class(pdata,
                  'pepData')
  
  # Load: proData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'little_prdata.RData',
                   package = 'pmartR'))
  
  # Test: proData --------------------------------------------------------------
  
  # Generate a proData object.
  prdata <- as.proData(e_data = edata,
                       f_data = fdata,
                       e_meta = emeta,
                       edata_cname = 'Reference',
                       fdata_cname = 'SampleID',
                       emeta_cname = 'PClass',
                       data_scale_orig = "abundance")
  
  # Replace NAs with 0s in prdata.
  expect_message(prdata2 <- edata_replace(omicsData = prdata,
                                          x = NA,
                                          y = 0),
                 "234 instances of NA have been replaced with 0")
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(
    attr(prdata2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(prdata2$e_data[, 1])),
         num_miss_obs = sum(is.na(prdata2$e_data)),
         prop_missing = (sum(is.na(prdata2$e_data)) /
                           prod(dim(prdata2$e_data[, -1]))),
         num_samps = ncol(prdata2$e_data[, -1]),
         data_types = NULL)
  )
  
  # Make sure all other attributes have not changed.
  expect_equal(attr(prdata, 'cnames'),
               attr(prdata2, 'cnames'))
  expect_equal(attr(prdata, 'check.names'),
               attr(prdata2, 'check.names'))
  expect_equal(attr(prdata, 'meta_info'),
               attr(prdata2, 'meta_info'))
  expect_equal(attr(prdata, 'group_DF'),
               attr(prdata2, 'group_DF'))
  expect_equal(attr(prdata, 'filters'),
               attr(prdata2, 'filters'))
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(prdata2$e_data == 0),
               234)
  
  # Examine the class of the new prdata object.
  expect_s3_class(prdata,
                  'proData')
  
})
