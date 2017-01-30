
#### Check plot.pepData ####

context("plot as pepData")

test_that("plot.pepData fails on non-specified objects", {
  
  ## Create some counterfeits ##
  data(metab_object)
  data(lipid_object)
  e_data = e_dat = f_data = e_meta = c(1, 2, 3, 4, 5)
  
  # Object with correct names but wrong setup #
  almost_object = data.frame(e_data, f_data, e_meta)
  
  # Object with wrong name, dataframe setup #
  almost_object2 = data.frame(e_dat, f_data, e_meta)
  
  # Correct object, corrupted e_data #
  almost_object3 = pep_object
  almost_object3$e_data = NULL
  
  
  # Check the following: summary object, NULL pointer, non-pep dataframes, #
  #   and corrupted dataframes #
  expect_that(plot.pepData(s_object), throws_error())
  expect_that(plot.pepData(NULL), throws_error())
  expect_that(plot.pepData(almost_object2), throws_error())
  expect_that(plot.pepData(almost_object3), throws_error())
  
  ## This one's a concern, it's supposed to throw an error ##
  #   It does in the console, but not here? Odd.   #
  expect_that(plot.pepData(almost_object), throws_error())
  
  ## These guys are known errors, and aren't expected to be a problem ##
  # expect_that(plot.pepData(metab_object), throws_error())
  # expect_that(plot.pepData(lipid_object), throws_error())
})