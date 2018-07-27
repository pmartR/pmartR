context("output tests for nonmissing_per_group()")
# Testing function for nonmissing_per_group 
library(pmartR)
library(pmartRdata)
library(dplyr)
data("pep_object")
omicsData <- group_designation(omicsData = pep_object, main_effects = "Condition")

e_data <- omicsData$e_data
groupDF <- attr(omicsData, "group_DF")
samp_id <- attr(omicsData, "cnames")$fdata_cname
cname_id <- attr(omicsData, "cnames")$edata_cname
mat <- matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3)

#here we apply nonmissing_per_group to omicsData and store in result
result<- pmartR:::nonmissing_per_group(omicsData)
result_df <- result$nonmiss_totals[1:5,] 

#hardcoded results are from running nonmissing_per_group with "pep_object" 
#first we will hardcode a data.frame with nonmiss_totals
Mass_Tag_ID <- c("1024","1047","1055","1104","1110")
Infection <- c(0,9,8,9,8)
Mock <- c(1,3,3,3,3)

hardcoded_nonmissdf <- data.frame(Mass_Tag_ID, Infection, Mock, stringsAsFactors = FALSE)

test_that("results are of appropriate class and length",{ 
  expect_that(result, is_a("list"))   
  expect_that(length(result), equals(2))   
  expect_that(result$group_sizes, is_a("data.frame"))   
  expect_that(result$nonmiss_totals, is_a("data.frame"))
})
  
test_that("components of result match attributes of omicsData",{
  #here we double check that result$group_sizes is correct using dplyr::summarise   
  expect_that(result$group_sizes, equals(data.frame(dplyr::summarise(dplyr::group_by(groupDF, Group), n_group = n()))))     
  
  #here we make sure the number of rows of nonmiss_totals in result is the same as the number of rows in e_data of omicsData   
  expect_that(nrow(result$nonmiss_totals), equals(nrow(omicsData$e_data)))      
  
  #here we check that the peptide column in result is the same as in e_data of omicsData   
  expect_that(result$nonmiss_totals[,cname_id], equals(as.character(omicsData$e_data[,which(names(omicsData$e_data)==cname_id)])))      
  
  #double checking that the names(result$nonmiss_totals) is correct by comparing to levels of (groupDF$Group)   
  expect_that(names(result$nonmiss_totals)[match(levels(groupDF$Group),names(result$nonmiss_totals))], equals(levels(groupDF$Group)))    
})

context("input tests for nonmissing_per_group()")


test_that("invalid input for omicsData argument throws error",{     
  expect_that(pmartR::nonmissing_per_group(omicsData = mat),throws_error())
}) 
 
test_that("leaving all arguments null throws error",{     
  expect_that(pmartR::nonmissing_per_group(omicsData = NULL, e_data = NULL, groupDF = NULL, cname_id = NULL, samp_id = NULL),throws_error())
})  

test_that("inputing a single argument and leaving the rest null throws error",{ 
  expect_that(pmartR::nonmissing_per_group(e_data = e_data),throws_error())   
  expect_that(pmartR::nonmissing_per_group(groupDF = groupDF),throws_error())   
  expect_that(pmartR::nonmissing_per_group(cname_id = cname_id),throws_error())   
  expect_that(pmartR::nonmissing_per_group(samp_id = samp_id),throws_error())      
})

test_that("e_data and groupDF are specified but cname_id and samp_id are not, throws error ",{
  expect_that(pmartR::nonmissing_per_group(e_data = e_data, groupDF = groupDF, cname_id = NULL, samp_id = NULL),throws_error())   
})


context("tests using hard coded results of nonmissing_per_group()")

test_that("output of nonmissing_per_group() matches hard coded results",{      
  expect_that(result_df, equals(hardcoded_nonmissdf))   
})
