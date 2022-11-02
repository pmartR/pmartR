context('test helper functions')

test_that('helper functions pull from correct attributes',{
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  myseqData <- as.seqData(e_data = edata,
                      f_data = fdata,
                      edata_cname = 'ID_REF',
                      fdata_cname = 'Samples'
  )
  
  expect_identical(get_check_names(myseqData), attr(myseqData, "check.names"))
  expect_identical(get_data_info(myseqData), attr(myseqData, "data_info"))
  expect_identical(get_data_norm(myseqData), attr(myseqData, "data_info")$norm_info$is_normalized)
  expect_identical(get_data_scale(myseqData), attr(myseqData, "data_info")$data_scale)
  expect_identical(get_edata_cname(myseqData), attr(myseqData, "cnames")$edata_cname)
  expect_identical(get_emeta_cname(myseqData), attr(myseqData, "cnames")$emeta_cname)
  expect_identical(get_fdata_cname(myseqData), attr(myseqData, "cnames")$fdata_cname)
  expect_identical(get_meta_info(myseqData), attr(myseqData, "meta_info"))
  expect_true(get_check_names(set_check_names(myseqData)))
  expect_false(get_check_names(set_check_names(myseqData, FALSE)))
  
  mes_1 <- "No filters have been applied"
  ## check after processing as well ##
  expect_null(expect_message(get_filters(myseqData), mes_1))
  expect_null(get_group_DF(myseqData))
  expect_null(get_group_table(myseqData))
  
  err_1 <- "object must be of class 'statRes' or 'trellData'"
  err_2 <- "dcObj object must be of class 'statRes' or 'trellData'"
  err_3 <- "omicsData must be of class 'isobaricpepData' and 'pepData'"
  err_4 <- "omicsData must be of class 'pepData', 'proData' or 'isobaricpepData'"
  err_5 <- "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'"
  err_6 <- "omicsData must be of class 'nmrData'"
  err_7 <- "No methods were selected for scoring, there is no 'best' set of parameters to return."
  
  ## Error helpers
  testthat::expect_error(get_comparisons(myseqData))
  testthat::expect_error(get_data_class(myseqData))
  testthat::expect_error(get_isobaric_info(myseqData))
  testthat::expect_error(get_isobaric_norm(myseqData))
  testthat::expect_error(get_nmr_info(myseqData))
  testthat::expect_error(get_nmr_norm(myseqData))
  testthat::expect_error(get_spans_params(myseqData))
  
})