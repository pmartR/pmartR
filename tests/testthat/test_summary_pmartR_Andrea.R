library(pmartRdata)
library(pmartRqc)
library(testthat)


#### check summary.pepData() ####
context("summary as pep_object")

data(pep_object)

# Check pep_object's a pep_Data class #
test_that("pepObject is the right type", {
  expect_is(pep_object, "pepData")})

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
test_that("summary results are calculated correctly", {
  # Check all 5 objects in s_object #
  expect_equal(s_object$num_samps, n_row_f)
  expect_equal(s_object$num_edata, n_row_e)
  expect_equal(s_object$num_emeta, length(unique(pep_object$e_meta$Ref_ID)))

  expect_equal(s_object$num_miss_obs, sumMiss)
  expect_equal(s_object$prop_missing, round(sumMiss / (n_row_e * n_col_e), 4))
})



