library(pmartRdata)
library(pmartRqc)

context("Test that data package loads")

data(pep_object)

expect_is(pep_object, "pepData")