#testthat function for summary.proData()

library(pmartRdata)
data("pro_object")
omicsData <- pro_object
result <- summary(omicsData)


context("output tests for summary.proData()")

test_that("summary results are of correct type", {
  expect_that(result, is_a("list"))
  expect_that(length(result), equals(5))
})

test_that("result matches attributes of omicsData",{ 
  expect_that(result$num_samps, equals(attr(omicsData, "data_info")$num_samps))
  expect_that(result$num_edata, equals(attr(omicsData, "data_info")$num_edata))
  expect_that(result$num_emeta, equals(attr(omicsData, "data_info")$num_emeta))
  expect_that(result$num_miss_obs, equals(attr(omicsData, "data_info")$num_miss_obs))
  expect_that(result$prop_missing, equals(round(attr(omicsData, "data_info")$num_miss_obs/(nrow(omicsData$e_data)*ncol(omicsData$e_data)), 3)))
})