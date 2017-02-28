# testing function for rip() in subset_funcs.R

library(testthat)
library(pmartRdata)
data("pep_object")

omicsData<- group_designation(omicsData = pro_object, main_effects = "Condition")

e_data = omicsData$e_data
edata_id = attr(omicsData,"cnames")$edata_cname
groupDF = attr(omicsData, "group_DF")
mat<- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

result<- pmartRqc:::rip(e_data = omicsData$e_data, edata_id = attr(omicsData,"cnames")$edata_cname, groupDF = attr(omicsData, "group_DF"))

#hardcoded results are from running rip() with "pep_object" 

hardcoded_rip_result <- c("6PGL_HUMAN", "AATM_HUMAN", "ACAD9_HUMAN", "ACON_HUMAN", "ACOT1_HUMAN")
hardcoded_result_len <- 223



context("output tests for rip()")

test_that("some output tests",{
  
  #checking result for correct class
  expect_that(result,is_a("character"))
  
  #checking that result is smaller than e_data, implying that peptides were filtered
  expect_that(length(result) < nrow(omicsData$e_data) ,is_true())
 
  #checking that more than two peptides remained after applying filter
  expect_that(length(result) > 2 ,is_true())
  
})



context("input tests for rip()")


test_that("some input tests",{
  
  # one NULL argument
  expect_that(pmartRqc:::rip(e_data = NULL, edata_id = edata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = NULL, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = NULL), throws_error())
  
  # incorrect arguments for e_data
  expect_that(pmartRqc:::rip(e_data = c(1,2,3), edata_id = edata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::rip(edata_id = 11, edata_id = edata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::rip(groupDF = mat, edata_id = edata_id, groupDF = groupDF), throws_error())

  #incorrect arguments for alpha and proportion
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, alpha = 50), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, alpha = "three"), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, alpha = c(1,2,3)), throws_error())
  
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, proportion = 1), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, proportion = "half"), throws_error())
  expect_that(pmartRqc:::rip(e_data = e_data, edata_id = edata_id, groupDF = group_DF, propotion = c(1,2,3)), throws_error())
  
})

context("tests using hard coded results of rip()")

test_that("output of rip() matches hard coded results",{
  
  expect_that(result[1:5], equals(hardcoded_rip_result))
  expect_that(length(result), equals(hardcoded_result_len))
  
})