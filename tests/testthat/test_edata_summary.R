context("test edata_summary")

library(testthat)
library(pmartR)
library(pmartRdata)

mypepData <- pmartRdata::pep_object

mypepData$f_data$fakegroup1 <- 1:12
mypepData$f_data$fakegroup2 <- c(1:9, 1:3)
mypepData$f_data$fakegroup3 <- c(1:8, 1, 9:11)

bad_columns = list(
  bad_col = c("___COLNAME___"),
  bad_col_1 = c("Condition", "___COLNAME__"),
  bad_col_all = c("__COLNAME1__", "___COLNAME2__"),
  bad_col_NA = c("Condition", NA)
)

test_that("by sample throws errors", {
  expect_error(edata_summary(mypepData, by = "sample", groupvar = "Condition"))
  lapply(bad_columns, function(col){
    expect_error(edata_summary(mypepData, by = "sample", groupvar = col))
  })
})

test_that("by molecule throws errors", {
  expect_error(edata_summary(mypepData, by = "molecule", groupvar = "fakegroup1"))
  expect_error(edata_summary(mypepData, by = "molecule", groupvar = c("Condition", "fakegroup2")))
  
  lapply(bad_columns, function(col){
    expect_error(edata_summary(mypepData, by = "molecule", groupvar = col))
  })
  
})

test_that("test correctness of some output", {
  suppressWarnings({
    summary1 <- edata_summary(mypepData, by = "molecule", groupvar = "fakegroup3")
    summary2 <- edata_summary(mypepData, by = "molecule", groupvar = "fakegroup2")
    summary3 <- edata_summary(mypepData, by = "molecule", groupvar = c("Condition", "fakegroup3"))
    summary4 <- edata_summary(mypepData, by = "sample")
  })
  
  expect_true(nrow(summary1$n_per_grp) == 1)
  expect_true(nrow(summary2$n_per_grp) == 3)
  expect_true(nrow(summary3$n_per_grp) == 1)
  expect_true(all(lapply(summary4, nrow) == 12))
  
  lapply(list(summary1, summary2, summary3, summary4), function(summary){
    expect_true(!any(unlist(lapply(summary, is.null))))
    expect_true(all(unlist(lapply(summary, inherits, "data.frame"))))
  })
  
  
  
})
