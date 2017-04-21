# testing function for ppp_rip() in subset_funcs.R 
library(pmartRqc)
library(testthat)
library(pmartRdata)
data("pep_object")

omicsData<- group_designation(omicsData = pep_object, main_effects = "Condition")

result<- pmartRqc:::ppp_rip(e_data = omicsData$e_data, edata_id = attr(omicsData,"cnames")$edata_cname, fdata_id = attr(omicsData,"cnames")$fdata_cname, groupDF = attr(omicsData, "group_DF"))

e_data = omicsData$e_data
f_data = omicsData$f_data
edata_id = attr(omicsData,"cnames")$edata_cname
fdata_id = attr(omicsData,"cnames")$fdata_cname
groupDF = attr(omicsData, "group_DF")
mat<- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

#hardcoded results are from running ppp_rip() with "pep_object" 

hardcoded_result <- c("1110", "11078", "216191", "6709133", "6753571")
hardcoded_result_len <- 2812

context("output tests for ppp_rip()")


test_that("results are of appropriate class and length",{  
  expect_that(result, is_a("character"))
  expect_that(length(result) >= 2 , is_true()) 
})

test_that("sample names in e_data match those in f_data",{
  # Check whether the column names of omicsData$e_data and omicsData$f_data$Sample_Name are the same   
  expect_that(all(colnames(e_data)[-which(colnames(e_data) == edata_id)] == as.character(f_data[[fdata_id]])) == TRUE, is_true())
})

context("input tests for ppp_rip()")


test_that("one NULL argument throws an error",{ 
  expect_that(pmartRqc:::ppp_rip(e_data = NULL, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())   
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = NULL, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = NULL, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = NULL), throws_error())
})

test_that("invalid input for e_data argument throws error",{      
  expect_that(pmartRqc:::ppp_rip(e_data = c(1,2,3),edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = 11, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
  expect_that(pmartRqc:::ppp_rip(e_data = mat, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF), throws_error())
})

test_that("invalid input for edata_id argument throws error",{     
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = "mass_tag_id", fdata_id = fdata_id, groupDF = groupDF), throws_error())   
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = "abcdefg", fdata_id = fdata_id, groupDF = groupDF), throws_error())
})

test_that("invalid input for fdata_id argument throws error",{     
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = "sampleid", groupDF = groupDF), throws_error())   
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = "abcdefg", groupDF = groupDF), throws_error())
})

test_that("invalid input for alpha argument throws error",{
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = 50), throws_error())
 # expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = -50), throws_error())   
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = "three"), throws_error()) 
 # expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, alpha = c(.1,.2)), throws_error()) 
})

test_that("invalid input for proportion argument throws error",{    
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = 50), throws_error()) 
  #expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = -50), throws_error())   
  expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = "one"), throws_error())  
  #expect_that(pmartRqc:::ppp_rip(e_data = e_data, edata_id = edata_id, fdata_id = fdata_id, groupDF = groupDF, proportion = c(.1,.2)), throws_error())                  
})
  
  
context("tests using hard coded results of ppp_rip()")

test_that("output of ppp_rip() matches hard coded results",{
  expect_that(result[1:5], equals(hardcoded_result))
  expect_that(length(result), equals(hardcoded_result_len))
})