context('plotting errors')

test_that('plotting errors are triggered correctly', {
  load(system.file('testdata',
    'little_seqdata.RData',
    package = 'pmartR'
  ))

  myseqData <- as.seqData(
    e_data = edata,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  )

  load(system.file('testdata',
    'little_isodata.RData',
    package = 'pmartR'
  ))

  isodata <- as.isobaricpepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Peptide',
    fdata_cname = 'Sample',
    emeta_cname = 'Protein'
  )

  err_1 <- "lkpm is not a valid option for 'transformation'"
  err_2 <- "log is not a valid option for 'transformation'. "

  testthat::expect_error(plot(myseqData, transformation = "lkpm"), err_1)
  testthat::expect_error(plot(myseqData, transformation = "log"), err_2)
  
  # "unused argument(s): lcpm"
  testthat::expect_warning(plot(isodata, transformation = "lcpm")) 
})
