context('plotting errors')

test_that('plotting errors are triggered correctly',{

  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  
  myseqData <- as.seqData(e_data = edata,
                          f_data = fdata,
                          edata_cname = 'ID_REF',
                          fdata_cname = 'Samples'
  )
  
  err_1 <- "lkpm is not a valid option for 'transformation'"
  err_2 <- "log is not a valid option for 'transformation'. "
  err_3 <- "unused argument" ## maybe change to be more specific?
  
  testthat::expect_error(plot(myseqData, transformation = "lkpm"), err_1)
  testthat::expect_error(plot(myseqData, transformation = "log"), err_2)
  testthat::expect_error(plot(pmartRdata::isobaric_object, transformation = "lcpm"), err_3)
  
})