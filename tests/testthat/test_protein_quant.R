#protein_quant testthat 
library(testthat)
library(pmartRqc)
library(pmartRdata)

data("pep_object")
data("lipid_object")
mypepData = edata_transform(pep_object, "log2")

edata_cname = attr(mypepData, "cnames")$edata_cname
#for some reason attr(pep_object, "cnames")$emeta_cname is NULL
attr(mypepData, "cnames")$emeta_cname<- "Protein"

zrollup_result = protein_quant(mypepData, method = 'zrollup', combine_fn = 'median', single_pep = TRUE, single_observation = TRUE)

vector = c(1,2,3)

#hardcoding results for zrollup



context("output tests for protein_quant(method = 'zrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(zrollup_result), equals("proData"))
  expect_that(length(zrollup_result), equals(3))
  expect_that(length(attr(zrollup_result, "filters")), equals(2))
})

context("input tests for protein_quant(method = 'zrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', single_pep = FALSE, single_observation = FALSE), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', isoformRes = vector, single_pep = TRUE, single_observation = TRUE, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', single_pep = TRUE, single_observation = TRUE, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'zrollup', single_pep = TRUE, single_observation = TRUE), throws_error())
  })


##################### tests for protein_quant(method = 'qrollup') #########################
#qrollup_result = protein_quant(mypepData, method = 'qrollup', qrollup_thresh = 2, combine_fn = 'mean')
qrollup_thresh = 2

#hardcoding results for qrollup 
#
#
#

context("output tests for protein_quant(method = 'qrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(qrollup_result), equals("proData"))
  expect_that(length(qrollup_result), equals(3))
})

context("input tests for protein_quant(method = 'qrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = "abcdef", combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = 2, isoformRes = vector, combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = 2, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'qrollup', qrollup_thresh = 2, combine_fn = "mean"), throws_error())
})


##################### tests for protein_quant(method = 'qrollup') #########################
#rrollup_result = protein_quant(mypepData, method = 'rrollup', combine_fn = 'median')

#hardcoding results for qrollup of peptides
#
#
#

context("output tests for protein_quant(method = 'rrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(rrollup_result), equals("proData"))
  expect_that(length(rrollup_result), equals(3))
})

context("input tests for protein_quant(method = 'rrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'rrollup', isoformRes = vector, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'rrollup', combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'rrollup', combine_fn = "median"), throws_error())
})

