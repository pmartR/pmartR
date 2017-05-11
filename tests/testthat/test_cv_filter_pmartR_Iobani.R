# testing function for cv_filter() 
library(pmartRqc)
library(testthat)
library(pmartRdata)
data("pep_object")

omicsData <- group_designation(omicsData = pep_object, main_effects = "Condition")
result <- pmartRqc:::cv_filter(omicsData)

samp_id <- attr(omicsData, "cnames")$fdata_cname
edata_id <- attr(omicsData, "cnames")$edata_cname
mat <- matrix(c(1,2,3,4,5,6),nrow = 2, ncol = 3)
vec <- c(1,2,3)

#hardcoded results are from running cv_filter with "pep_object" 
hardcoded_cv_result <- c(NaN, 36.01185, 32.76790, 31.93044, 47.89102)
hardcoded_result_dim <- c(17407, 2)
hardcoded_max  <- 91.95324
hardcoded_tot_nas <- 2173


context("output tests for cv_filter()")


test_that("result is of appropriate class and length",{ 
  expect_that(result, is_a(c("cvFilt","data.frame")))      
  expect_that(length(result), equals(2)) 
})

test_that("components of result match attributes of omicsData",{     
  #checking that groupDF in result is the same as groupDF in omicsData
  expect_that(attr(result,"group_DF"), equals(attr(omicsData,"group_DF")))
 
  #checking that number of rows in result matches the rows of e_data
  expect_that(nrow(result), equals(nrow(omicsData$e_data)))
 
  #checking that edata_cname  in result match those in e_data
  expect_that(as.character(sort(result[[edata_id]])), is_equivalent_to(as.character(sort(omicsData$e_data[,which(names(omicsData$e_data)%in%edata_id)]))))
 
  #checking for correct class
  expect_that(result$CV_pooled, is_a("numeric"))
  
  #checking tha sample_names of result match those in e_data
  expect_that(attr(result,"sample_names"), equals(names(omicsData$e_data)[-which(names(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]))
  
})

context("input tests for cv_filter()")

test_that("invalid input for omicsData argument throws error",{     
  expect_that(pmartRqc:::cv_filter(mat), throws_error())  
  expect_that(pmartRqc:::cv_filter(vec), throws_error())
})

context("tests using hard coded results of cv_filter()")

test_that("output of cv_filter() matches hard coded results",{
  
  expect_that(round(result$CV_pooled[1:5], 4), equals(round(hardcoded_cv_result, 4)))
  expect_that(dim(result), equals(hardcoded_result_dim))
  expect_that(round(attr(result,"max_x_val"), 4), equals(round(hardcoded_max, 4)))
  expect_that(attr(result,"tot_nas"), equals(hardcoded_tot_nas))
  
})
