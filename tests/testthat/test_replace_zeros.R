context('replace missing values')

test_that('as.xxx correctly replaces missing values (0s) with NAs',{
  
  # Load: isobaricpepData ------------------------------------------------------
  
  load(system.file('testdata',
                   'little_isodata.RData',
                   package = 'pmartR'))
  
  # Mutate all NAs in edata.
  edata[is.na(edata)] <- 0
  
  # Test: isobaricpepData ------------------------------------------------------
  
  # Fabricate an isobaricpepData object.
  expect_message(isodata <- as.isobaricpepData(e_data = edata,
                                               f_data = fdata,
                                               e_meta = emeta,
                                               edata_cname = 'Peptide',
                                               fdata_cname = 'Sample',
                                               emeta_cname = 'Protein'),
                 '400 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(isodata)$data_info$num_edata,
               150)
  expect_equal(attributes(isodata)$data_info$num_miss_obs,
               400)
  expect_equal(round(attributes(isodata)$data_info$prop_missing, 4),
               0.2222)
  expect_equal(attributes(isodata)$data_info$num_samps,
               12)
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(is.na(isodata$e_data)),
               400)
  
  # Examine the class of the new isodata object.
  expect_s3_class(isodata,
                  c('isobaricpepData', 'pepData'))
  
  # Load: lipidData ------------------------------------------------------------
  
  load(system.file('testdata',
                   'lipidData.RData',
                   package = 'pmartR'))
  
  # Transmogrify all NAs in edata.
  edata[is.na(edata)] <- 0
  
  # Test: lipidData ------------------------------------------------------------
  
  # Produce a lipidData object.
  expect_message(ldata <- as.lipidData(e_data = edata,
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'LipidCommonName',
                                       fdata_cname = 'Sample_Name',
                                       emeta_cname = 'LipidClass'),
                 '884 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(ldata)$data_info$num_edata,
               146)
  expect_equal(attributes(ldata)$data_info$num_miss_obs,
               884)
  expect_equal(round(attributes(ldata)$data_info$prop_missing, 4),
               0.5504)
  expect_equal(attributes(ldata)$data_info$num_samps,
               11)
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(is.na(ldata$e_data)),
               884)
  
  # Examine the class of the new ldata object.
  expect_s3_class(ldata,
                  'lipidData')
  
  # Load: metabData ------------------------------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
  # Transmute all NAs in edata.
  edata[is.na(edata)] <- 0
  
  # Test: metabData ------------------------------------------------------------
  
  # Forge a metabData object.
  expect_message(mdata <- as.metabData(e_data = edata,
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'MClass'),
                 '148 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(mdata)$data_info$num_edata,
               80)
  expect_equal(attributes(mdata)$data_info$num_miss_obs,
               148)
  expect_equal(round(attributes(mdata)$data_info$prop_missing, 4),
               0.1542)
  expect_equal(attributes(mdata)$data_info$num_samps,
               12)
  
  # Ensure the number of replaced elements is accurate.
  expect_equal(sum(is.na(mdata$e_data)),
               148)
  
  # Examine the class of the new mdata object.
  expect_s3_class(mdata,
                  'metabData')
  
  # Load: nmrData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Alter the original nmr data frame to contain zeros
  set.seed(22)
  
  # Select row indices.
  rows <- sample(1:38, 338, replace = TRUE)
  
  # Select column indices.
  columns <- sample(2:42, 338, replace = TRUE)
  
  # Loop through each row column pair.
  for (e in 1:338) {
    
    # Replace current value in edata with zero.
    edata[rows[[e]], columns[[e]]] <- 0
    
  }
  
  # Test: nmrData --------------------------------------------------------------
  
  # Fabricate an nmrData object.
  expect_message(nmrdata <- as.nmrData(e_data = edata,
                                       f_data = fdata,
                                       e_meta = emeta,
                                       edata_cname = 'Metabolite',
                                       fdata_cname = 'SampleID',
                                       emeta_cname = 'nmrClass'),
                 '302 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(nmrdata)$data_info$num_edata,
               38)
  expect_equal(attributes(nmrdata)$data_info$num_miss_obs,
               302)
  expect_equal(round(attributes(nmrdata)$data_info$prop_missing, 4),
               0.1938)
  expect_equal(attributes(nmrdata)$data_info$num_samps,
               41)
  
  # Ensure the number of replaced elements is correct.
  expect_equal(sum(is.na(nmrdata$e_data)),
               302)
  
  # Examine the class of the new nmrdata object.
  expect_s3_class(nmrdata,
                  'nmrData')
  
  # Load: pepData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Transfigure all NAs in edata.
  edata[is.na(edata)] <- 0
  
  # Test: pepData --------------------------------------------------------------
  
  # Construct a pepData object.
  expect_message(pdata <- as.pepData(e_data = edata,
                                     f_data = fdata,
                                     e_meta = emeta,
                                     edata_cname = 'Mass_Tag_ID',
                                     fdata_cname = 'SampleID',
                                     emeta_cname = 'Protein'),
                 '341 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(pdata)$data_info$num_edata,
               150)
  expect_equal(attributes(pdata)$data_info$num_miss_obs,
               341)
  expect_equal(round(attributes(pdata)$data_info$prop_missing, 4),
               0.1894)
  expect_equal(attributes(pdata)$data_info$num_samps,
               12)
  
  # Ensure the number of replaced elements is accurate.
  expect_equal(sum(is.na(pdata$e_data)),
               341)
  
  # Examine the class of the new pdata object.
  expect_s3_class(pdata,
                  'pepData')
  
  # Load: proData --------------------------------------------------------------
  
  load(system.file('testdata',
                   'little_prdata.RData',
                   package = 'pmartR'))
  
  # Transmute all NAs in edata.
  edata[is.na(edata)] <- 0
  
  # Test: proData --------------------------------------------------------------
  
  # Generate a proData object.
  expect_message(prdata <- as.proData(e_data = edata,
                                      f_data = fdata,
                                      e_meta = emeta,
                                      edata_cname = 'Reference',
                                      fdata_cname = 'SampleID',
                                      emeta_cname = 'PClass'),
                 '234 instances of 0 have been replaced with NA')
  
  # Check that the elements of the data_info attribute are all correct.
  expect_equal(attributes(prdata)$data_info$num_edata,
               150)
  expect_equal(attributes(prdata)$data_info$num_miss_obs,
               234)
  expect_equal(round(attributes(prdata)$data_info$prop_missing, 4),
               0.1418)
  expect_equal(attributes(prdata)$data_info$num_samps,
               11)
  
  # Ensure the number of replaced elements is accurate.
  expect_equal(sum(is.na(prdata$e_data)),
               234)
  
  # Examine the class of the new prdata object.
  expect_s3_class(prdata,
                  'proData')

})
