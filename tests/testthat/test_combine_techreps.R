context("aggregate_tech_reps")
library(testthat)
library(pmartR)
library(pmartRdata)

techrep_pepData <- pmartRdata::techrep_pep_object

pep_techrep_avg <- combine_techreps(techrep_pepData)
techrep_cname = attr(pep_techrep_avg, "cnames")$techrep_cname

test_that("bad input throws error", {
  # not a column
  expect_error(combine_techreps(techrep_pepData, bio_sample_names = "asdf"))
  # names do not correspond to biological sample assignment
  expect_error(combine_techreps(techrep_pepData, bio_sample_names = "DILUTION"))
  expect_error(combine_techreps(techrep_pepData, bio_sample_names = "FACTOR"))
  expect_error(combine_techreps(techrep_pepData, bio_sample_names = "RunID"))
})

test_that("columns were correctly aggregated", {
  # number of biological samples....
  n_groups = length(unique(techrep_pepData$f_data[,techrep_cname]))
  
  # ...should equal the number of columns in the collapsed e_data - 1 
  expect_equal(n_groups, length(colnames(pep_techrep_avg$e_data)) - 1)
  # all rows should still have at least 1 nonmissing value
  expect_true(all(rowSums(!is.na(pep_techrep_avg$e_data[,-which(colnames(pep_techrep_avg$e_data) == get_edata_cname(pep_techrep_avg))])) > 0))
  # averaging should result in an e_data with a number of columns strictly less than that of the original data
  expect_true(ncol(techrep_pepData$e_data) > ncol(pep_techrep_avg$e_data))
})

test_that("attributes correctly set", {
  expect_false(get_fdata_cname(techrep_pepData) == get_fdata_cname(pep_techrep_avg))
  expect_true(get_fdata_cname(pep_techrep_avg) == techrep_cname)
  expect_equal(length(setdiff(names(attributes(pep_techrep_avg)$tech_rep_info$tech_reps_by_sample), as.character(pep_techrep_avg$f_data[, get_fdata_cname(pep_techrep_avg)]))), 0)
  expect_true(setdiff(colnames(techrep_pepData$e_data), attributes(pep_techrep_avg)$tech_rep_info$tech_reps_by_sample %>% unlist()) == get_edata_cname(techrep_pepData))
})


