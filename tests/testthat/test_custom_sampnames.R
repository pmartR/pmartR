context("custom_sampnames")

test_that("sample names are properly abridged",{
  
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
                      emeta_cname = 'Protein')
  
  # Holy short sample names, Batman! -------------------------------------------
  
  # firstn ---------------
  
  # Create short names with the first four characters of the sample names using
  # the firstn argument.
  firstn_4 <- custom_sampnames(pdata, firstn = 4)
  
  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(firstn_4))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(firstn_4$f_data),
               c(12, 3))
  expect_equal(firstn_4$f_data$VizSampNames,
               c(rep("Infe", 9), rep("Mock", 3)))
  
  # Create not so short sample names (the input to firstn is larger than the 
  # longest sample name).
  expect_warning(firstn_11 <- custom_sampnames(pdata, firstn = 11))
  
  # Make sure nothing changes if firstn is larger than largest sample name.
  expect_equal(as.character(firstn_11$f_data$SampleID),
               firstn_11$f_data$VizSampNames)
  
  # from to ---------------
  
  # Range completely excludes some or all sample names.
  expect_error(custom_sampnames(pdata, from = 6, to = 10))
  expect_error(custom_sampnames(pdata, from = 11, to = 12))
  
  # Non-numeric inputs for from and to.
  expect_error(custom_sampnames(pdata, from = "2", to = "10"))
  
  # The value for to is less than from.
  expect_error(custom_sampnames(pdata, from = 3, to = 2))
  
  # Create short names using the from/to combo.
  from_to <- custom_sampnames(pdata, from = 3, to = 5)
  
  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(from_to))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(from_to$f_data),
               c(12, 3))
  expect_equal(from_to$f_data$VizSampNames,
               c(rep("fec", 9), paste0("ck", 1:3)))
  
  # delim ---------------
  
  # Use delim and components to make tiny sample names.
  delim_c <- custom_sampnames(pdata, delim = "c", components = 1)
  delim_empty <- custom_sampnames(pdata, delim = "", components = 5)
  
  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(delim_c))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(delim_c$f_data),
               c(12, 3))
  expect_equal(delim_c$f_data$VizSampNames,
               c(rep("Infe", 9), rep("Mo", 3)))
  expect_equal(dim(delim_empty$f_data),
               c(12, 3))
  expect_equal(delim_empty$f_data$VizSampNames,
               c(rep("c", 9), as.character(1:3)))
  
  # Components does not specify an index created by delim that exists.
  expect_error(custom_sampnames(pdata, delim = "c", components = 4))
  
})
