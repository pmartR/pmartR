# testing function for ppp_rip() in subset_funcs.R 

library(testthat)
library(pmartRdata)
data("pep_object")

#if(is.null(attr(pep_object,"group_DF"))){omicsData<- group_designation(omicsData = pep_object, main_effects = "Condition")}

omicsData<- group_designation(omicsData = pep_object, main_effects = "Condition")

result<- pmartRqc:::ppp_rip(e_data = omicsData$e_data, edata_id = attr(omicsData,"cnames")$edata_cname, fdata_id = attr(omicsData,"cnames")$fdata_cname, groupDF = attr(omicsData, "group_DF"))

#hardcoded results are from running ppp_rip() with "pep_object" 

hardcoded_result <- c("1110", "11078", "216191", "6709133", "6753571")
hardcoded_result_len <- 2812



e_data = omicsData$e_data
edata_id = attr(omicsData,"cnames")$edata_cname
fdata_id = attr(omicsData,"cnames")$fdata_cname
groupDF = attr(omicsData, "group_DF")

context("output tests for ppp_rip()")


test_that("some output tests",{
  #checking result for correct class
  expect_that(result, is_a("character"))
  
  #checking that result is smaller than e_data, implying that peptides were filtered
  expect_that(length(result) < nrow(omicsData$e_data), is_true())
 
  #checking that more than two peptides remained after applying filter
  expect_that(length(result) > 2 , is_true())

  #expect_that(length(result) == length(omicsData$e_data[,edata_id]), throws_error())
 
})

mat<- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

context("input tests for ppp_rip()")


test_that("some input tests",{
  
  # one NULL argument
 
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = NULL, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = NULL, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = NULL), throws_error())
  
  # incorrect arguments
  expect_that(pmartRqc:::ppp_rip(e_data = c(1,2,3),edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = 11, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = mat, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  
  #incorrect arguments for alpha and proportion
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = 50), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = "three"), throws_error())
  
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = 50), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = "one"), throws_error())
  
  
  
})

context("tests using hard coded results of ppp_rip()")

test_that("output of ppp_rip() matches hard coded results",{
  
  expect_that(result[1:5], equals(hardcoded_result))
  expect_that(length(result), equals(hardcoded_result_len))
})