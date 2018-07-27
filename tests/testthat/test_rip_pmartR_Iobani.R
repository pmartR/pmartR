context("output tests for rip()")
# testing function for rip() in subset_funcs.R
library(pmartR)
library(testthat)
library(pmartRdata)
data("pep_object")

omicsData<- group_designation(omicsData = pep_object, main_effects = "Condition")

e_data = omicsData$e_data
edata_id = attr(omicsData,"cnames")$edata_cname
fdata_id = attr(omicsData,"cnames")$fdata_cname
groupDF = attr(omicsData, "group_DF")
mat<- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

result<- pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF)

#hardcoded results are from running rip() with "pep_object" 

hardcoded_rip_result <- c("11078", "216191", "6709133", "6753571", "6781312")
hardcoded_result_len <- 1158

test_that("results are of appropriate class and length",{ 
  expect_that(result, is_a("character"))
  expect_that(length(result) >= 2, is_true())    
})


context("input tests for rip()")

test_that("one NULL argument throws an error",{     
  expect_that(pmartR:::rip(e_data = NULL, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartR:::rip(e_data = e_data, edata_id = NULL, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = NULL, groupDF = groupDF), throws_error())
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = NULL), throws_error())
  
})

test_that("incorrect input for e_data argument throws error",{    
  expect_that(pmartR:::rip(e_data = c(1,2,3), edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())   
  expect_that(pmartR:::rip(edata_id = 11, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())   
  expect_that(pmartR:::rip(groupDF = mat, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
})

test_that("incorrect input for edata_id argument throws error",{  
  expect_that(pmartR:::rip(e_data = e_data, edata_id = "mass_tag_id", fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartR:::rip(e_data = e_data, edata_id = "abcdefgh", fdata_id = fdata_id, groupDF = groupDF), throws_error())
})

test_that("incorrect input for fdata_id argument throws error",{  
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = "mass_tag_id", groupDF = groupDF), throws_error())
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = "abcdefgh", groupDF = groupDF), throws_error())
})

test_that("incorrect input for alpha argument throws error",{   
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = group_DF, alpha = 50), throws_error())  
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = group_DF, alpha = "three"), throws_error())   
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = group_DF, alpha = c(1,2,3)), throws_error())    
  expect_that(pmartR:::rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = group_DF, alpha = -2) , throws_error()) 
})


context("tests using hard coded results of rip()")

test_that("output of rip() matches hard coded results",{
  expect_that(result[1:5], equals(hardcoded_rip_result))
  expect_that(length(result), equals(hardcoded_result_len))
})