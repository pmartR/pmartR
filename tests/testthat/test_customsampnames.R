context("Test custom_sampnames")
library(pmartR)
library(pmartRdata)
library(testthat)

fdata_cname = attr(pep_object, "cnames")$fdata_cname

# largest value of sample names
max_firstn <- min(nchar(as.character(pep_object$f_data[,fdata_cname])))

# longest sample name value
max_sampname <- max(nchar(as.character(pep_object$f_data[,fdata_cname])))

test_that("test custom_sampnames() for firstn input", {
  # object and plot made from custom_sampnames with known good parameters
  firstn <- sample(1:max_firstn, 1)
  good_obj <- custom_sampnames(pep_object, firstn = firstn)
  p <- plot(good_obj, use_VizSampNames = TRUE)
  
  # throw error on value larger than max_firstn and for character input
  expect_warning(custom_sampnames(pep_object, firstn = max_firstn + 1))
  expect_error(custom_sampnames(pep_object, firstn = "0111"))
  
  # all values should be trimmed to length firstn
  expect_true(all(nchar(as.character(good_obj$f_data[,"VizSampNames"])) == firstn))
  
  # trimmed values should be equal to a substring from index 1 to firstn of the sample ID column
  expect_setequal(substr(good_obj$f_data[,fdata_cname], 1, firstn), good_obj$f_data$VizSampNames)
  
  # plot has correct labels
  expect_setequal(good_obj$f_data$VizSampNames, p$scales$scales[[1]]$labels)
})

test_that("test custom_sampnames() for from-to input",{ 
  ## errors should be thrown for:
  
  # range completely excludes some or all sample names
  expect_error(custom_sampnames(pep_object, from = max_firstn + 1, to = max_firstn + sample(1:(max_sampname-max_firstn), 1)))
  expect_error(custom_sampnames(pep_object, from = max_sampname + 1, to = max_sampname + 2))
  
  # non-numeric inputs
  expect_error(custom_sampnames(pep_object, from = "2", to = "10"))
  
  # to is less than from
  from_1 <- sample(1:max_firstn, 1)
  expect_error(custom_sampnames(pep_object, from = from_1, to = from_1 - 1))
  
  ## test correct subsetting and plots have subsetted sample names ##
  from <- sample(1:max_firstn, 1) 
  to <- sample(from:max_sampname, 1) 
  
  # good parameters and plot objects
  obj_1 <- custom_sampnames(pep_object, from = from, to = to)
  obj_2 <- custom_sampnames(pep_object, from = from, to = 1000)
  p1 <- plot(obj_1, use_VizSampNames = TRUE)
  p2 <- plot(obj_2, use_VizSampNames = TRUE)
  
  # correct subsetting
  expect_setequal(
    unlist(lapply(strsplit(as.character(pep_object$f_data$SampleID), split = ""), 
                  function(x){x <- x[from:to] 
                  paste(x[!is.na(x)], collapse = "")
                  })),
    obj_1$f_data[,"VizSampNames"]
  )
  
  expect_setequal(
    unlist(lapply(strsplit(as.character(pep_object$f_data$SampleID), split = ""), 
                  function(x){x <- x[from:1000] 
                  paste(x[!is.na(x)], collapse = "")
                  })),
    obj_2$f_data[,"VizSampNames"]
  )
  
  # plot has correct labels
  expect_setequal(obj_1$f_data$VizSampNames, p1$scales$scales[[1]]$labels)
  expect_setequal(obj_2$f_data$VizSampNames, p2$scales$scales[[1]]$labels)
})

test_that("test custom_sampnames() for delim+components",{
  
  select_components <- c(1, sample(2:max_firstn, 2))
  
  # object made from custom_sampnames with known good parameters
  delim_c <- custom_sampnames(pep_object, delim = "c", components = 1)
  delim_empty <- custom_sampnames(pep_object, delim = "", components = select_components)
  
  # errors should be thrown for:
  # none of the delimiters specify an index in the split sample
  expect_error(custom_sampnames(pep_object, delim = "c", components = c(3,4)))
  
  # correct subsetting
  expect_equal(
    unlist(lapply(strsplit(as.character(pep_object$f_data$SampleID), split = "c"), function(x){x[1]})),
    delim_c$f_data[,"VizSampNames"]
  )
  expect_equal(
    unlist(lapply(strsplit(as.character(pep_object$f_data$SampleID), split = ""), function(x){paste(x[select_components], collapse = "")})),
    delim_empty$f_data[,"VizSampNames"]
  )
  
  # plot has correct labels
  plot_c <- plot(delim_c, use_VizSampNames = TRUE)
  plot_empty <- plot(delim_empty, use_VizSampNames = TRUE)
  
  expect_setequal(delim_c$f_data$VizSampNames, plot_c$scales$scales[[1]]$labels)
  expect_setequal(delim_empty$f_data$VizSampNames, plot_empty$scales$scales[[1]]$labels)
})


