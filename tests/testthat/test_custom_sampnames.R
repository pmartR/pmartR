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

  # Make long sample names that can be used to test each method for reducing the
  # length of the sample names.
  names(pdata$e_data) <- c("Mass_Tag_ID",
                           "qwerty_one_infection1_asdf",
                           "qwerty_two_infection2_asdf",
                           "qwerty_three_infection3_asdf",
                           "qwerty_four_infection4_asdf",
                           "qwerty_five_infection5_asdf",
                           "qwerty_six_infection6_asdf",
                           "qwerty_seven_infection7_asdf",
                           "qwerty_eight_infection8_asdf",
                           "qwerty_nine_infection9_asdf",
                           "qwerty_one_mock1_asdf",
                           "qwerty_two_mock2_asdf",
                           "qwerty_three_mock3_asdf")

  # Change to the new long names in f_data
  pdata$f_data$SampleID <- c("qwerty_one_infection1_asdf",
                             "qwerty_two_infection2_asdf",
                             "qwerty_three_infection3_asdf",
                             "qwerty_four_infection4_asdf",
                             "qwerty_five_infection5_asdf",
                             "qwerty_six_infection6_asdf",
                             "qwerty_seven_infection7_asdf",
                             "qwerty_eight_infection8_asdf",
                             "qwerty_nine_infection9_asdf",
                             "qwerty_one_mock1_asdf",
                             "qwerty_two_mock2_asdf",
                             "qwerty_three_mock3_asdf")

  # Holy short sample names, Batman! -------------------------------------------

  # firstn ---------------

  # Create short sample names with firstn that are not unique
  expect_error(custom_sampnames(pdata, firstn = 6),
               paste("The input used does not produce a unique sample name",
                     "for each sample.",
                     sep = " "))

  # Create short names with the first four characters of the sample names using
  # the firstn argument.
  firstn_16 <- custom_sampnames(pdata, firstn = 16)

  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(firstn_16))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(firstn_16$f_data),
               c(12, 3))
  expect_equal(firstn_16$f_data$VizSampNames,
               c("qwerty_one_infec",
                 "qwerty_two_infec",
                 "qwerty_three_inf",
                 "qwerty_four_infe",
                 "qwerty_five_infe",
                 "qwerty_six_infec",
                 "qwerty_seven_inf",
                 "qwerty_eight_inf",
                 "qwerty_nine_infe",
                 "qwerty_one_mock1",
                 "qwerty_two_mock2",
                 "qwerty_three_moc"))

  # Create not so short sample names (the input to firstn is larger than the
  # longest sample name).
  expect_warning(firstn_35 <- custom_sampnames(pdata, firstn = 35))

  # Make sure nothing changes if firstn is larger than largest sample name.
  expect_identical(as.character(firstn_35$f_data$SampleID),
                   firstn_35$f_data$VizSampNames)

  # from to ---------------

  # Range completely excludes some or all sample names.
  expect_error(custom_sampnames(pdata, from = 16, to = 22))
  expect_error(custom_sampnames(pdata, from = 30, to = 33))

  # Non-numeric inputs for from and to.
  expect_error(custom_sampnames(pdata, from = "2", to = "10"))

  # The value for to is less than from.
  expect_error(custom_sampnames(pdata, from = 3, to = 2))

  # Create short names using the from/to combo.
  from_to <- custom_sampnames(pdata, from = 8, to = 16)

  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(from_to))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(from_to$f_data),
               c(12, 3))
  expect_equal(from_to$f_data$VizSampNames,
               c("one_infec",
                 "two_infec",
                 "three_inf",
                 "four_infe",
                 "five_infe",
                 "six_infec",
                 "seven_inf",
                 "eight_inf",
                 "nine_infe",
                 "one_mock1",
                 "two_mock2",
                 "three_moc"))

  # delim ---------------

  # Use delim and components to make tiny sample names.
  delim_u <- custom_sampnames(pdata, delim = "_", components = 3)
  delim_empty <- custom_sampnames(pdata, delim = "", components = 8:16)

  # Make sure attributes have not been changed after adding short names.
  expect_identical(attributes(pdata), attributes(delim_u))
  # Make sure the short sample names are present and correct.
  expect_equal(dim(delim_u$f_data),
               c(12, 3))
  expect_equal(delim_u$f_data$VizSampNames,
               c(paste0("infection", 1:9), paste0("mock", 1:3)))
  expect_equal(dim(delim_empty$f_data),
               c(12, 3))
  expect_equal(delim_empty$f_data$VizSampNames,
               c("one_infec",
                 "two_infec",
                 "three_inf",
                 "four_infe",
                 "five_infe",
                 "six_infec",
                 "seven_inf",
                 "eight_inf",
                 "nine_infe",
                 "one_mock1",
                 "two_mock2",
                 "three_moc"))

  # Components does not specify an index created by delim that exists.
  expect_error(custom_sampnames(pdata, delim = "_", components = 20),
               paste("none of the indices specified in 'components' match",
                     "indices of the split sample name",
                     sep = " "))
  
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  myseqData <- as.seqData(e_data = edata,
                          f_data = fdata,
                          edata_cname = 'ID_REF',
                          fdata_cname = 'Samples'
  )
  
  err_1 <- "Please specify both 'from' and 'to' if either argument is used."
  err_2 <- "Please specify both 'delim' and 'components' if either argument is used."
  
  expect_identical(
    strtrim(myseqData$f_data$Samples, 13),
    suppressWarnings(custom_sampnames(myseqData, firstn = 13)$f_data$VizSampNames) 
    ## Might wanna think about a lastn as well
  )
  
  testthat::expect_error(custom_sampnames(myseqData, from = 1), err_1)
  testthat::expect_error(custom_sampnames(myseqData, to = 13), err_1)
  
  expect_identical(
    substring(myseqData$f_data$Samples, first = 3, last = 13),
    custom_sampnames(myseqData, from = 3, to = 13)$f_data$VizSampNames
  )
  
  testthat::expect_error(custom_sampnames(myseqData, delim = "_"), err_2)
  testthat::expect_error(custom_sampnames(myseqData, components = c(1,2,3)), err_2)
  
  ## No changes actually happening here, demo of arg use
  expect_identical(
    myseqData$f_data$Samples,
    custom_sampnames(myseqData, delim = "_", components = c(1,2,3))$f_data$VizSampNames,
    custom_sampnames(myseqData, pattern = ".+")$f_data$VizSampNames
  )

})
