library(pmartRdata)
library(pmartRqc)
library(testthat)


#### check summary.pepData() ####
context("pep_Object testing")

data(pep_object)

# Check pep_object's a pep_Data class #
expect_is(pep_object, "pepData")

# For ease of reading, some comparison tools #
n_col_e = ncol(pep_object$e_data)
n_row_e = nrow(pep_object$e_data)
n_row_f = nrow(pep_object$f_data)
sumMiss = sum(is.na(pep_object$e_data))
uniq_Sigs = unique(pep_object$e_meta$Ref_ID)

# The summary object to be tested #
s_object = summary(pep_object)

test_that("summary results are of correct type", {
  # Test type and size of s_object #
  expect_is(s_object, "list")
  expect_equal(length(s_object), 5)
})
test_that("summary results arecalculated correctly", {
  # Check all 5 objects in s_object #
  expect_equal(s_object$num_samps, n_row_f)
  expect_equal(s_object$num_edata, n_row_e)
  expect_equal(s_object$num_emeta, length(unique(pep_object$e_meta$Ref_ID)))

  expect_equal(s_object$num_miss_obs, sumMiss)
  expect_equal(s_object$prop_missing, round(sumMiss / (n_row_e * n_col_e), 4))
})


#### Check plot.pepData ####

context("plot.pepData errors")

test_that("plot.pepData fails on non-specified objects", {
  
  ## Create some counterfeits ##
  data(metab_object)
  data(lipid_object)
  e_data = e_data = f_data = e_meta = c(1, 2, 3, 4, 5)
  
  # Object with correct names but wrong setup #
  almost_object = data.frame(e_data, f_data, e_meta)
  
  # Object with wrong name, dataframe setup #
  almost_object2 = data.frame(e_dat, f_data, e_meta)

  # Correct object, corrupted e_data #
  almost_object3 = pep_object
  almost_object3$e_data = NULL


  # Check the following: summary object, NULL pointer, and non-pep_object dataframes #
  expect_that(plot.pepData(s_object), throws_error())
  expect_that(plot.pepData(NULL), throws_error())
  expect_that(plot.pepData(almost_object2), throws_error())
  expect_that(plot.pepData(almost_object3), throws_error())

  ## This one's a concern, it's supposed to throw an error ##
  expect_that(plot.pepData(almost_object), throws_error())
  
  ## These guys are known errors, and aren't expected to be a problem ##
  expect_that(plot.pepData(metab_object), throws_error())
  expect_that(plot.pepData(lipid_object), throws_error())
})

