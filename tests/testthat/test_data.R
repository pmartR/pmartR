context("Test that data package loads")
library(testthat)
library(pmartRdata)
library(pmartR)

data(pep_object)

expect_is(pep_object, "pepData")
