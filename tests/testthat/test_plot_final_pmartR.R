context("input tests for plot.proData()")
# testing function for plot.proData()

library(testthat) 
library(pmartRdata) 
data("pro_object") 
omicsData<- pro_object

#Here we create a proData object that has e_data = NULL
almost_object <- omicsData
almost_object$e_data <- NULL

input<- c(1, 2, 3)
mat<- matrix(1:6, nrow = 2, ncol = 3)

test_that("invalid input for order_by argument throws error",{
  expect_that(plot(omicsData, order_by = input), throws_error())   
  expect_that(plot(omicsData, order_by = mat), throws_error())   
  expect_that(plot(omicsData, order_by = 11), throws_error())  
  expect_that(plot(omicsData, order_by = 1.223), throws_error())
  expect_that(plot(omicsData, order_by = -2), throws_error()) 
  expect_that(plot(omicsData, order_by = c("blue")), throws_error())   
  expect_that(plot(omicsData, order_by = c("Condition", "Status")), throws_error()) 
})

test_that("invalid input for color_by argument throws error",{   
  expect_that(plot(omicsData, color_by = input), throws_error())   
  expect_that(plot(omicsData, color_by = mat), throws_error())   
  expect_that(plot(omicsData, color_by = 11), throws_error()) 
  expect_that(plot(omicsData, color_by = 1.54), throws_error())
  expect_that(plot(omicsData, color_by = -4), throws_error())
  expect_that(plot(omicsData, color_by = c("blue")), throws_error())   
  expect_that(plot(omicsData, color_by = c("Condition", "Status")), throws_error())
}) 

test_that("invalid input for facet_by argument throws error",{          
  expect_that(plot(omicsData, facet_by = input), throws_error())   
  expect_that(plot(omicsData, facet_by = mat), throws_error())   
  expect_that(plot(omicsData, facet_by = 11), throws_error()) 
  expect_that(plot(omicsData, facet_by = 1.75), throws_error()) 
  expect_that(plot(omicsData, facet_by = -1), throws_error())   
  expect_that(plot(omicsData, facet_by = c("blue")), throws_error())   
  expect_that(plot(omicsData, facet_by = c("blue","green")), throws_error())  
})

test_that("invalid input for facet_cols argument throws error",{
  expect_that(plot(omicsData, facet_by = NULL,facet_cols = 1), throws_error())   
  expect_that(plot(omicsData, facet_by = names(omicsData$f_data[2]),facet_cols = 0), throws_error())      
})

test_that("invalid input for omicsData argument throws error",{  
  expect_that(pmartR::plot.proData(input), throws_error())   
  expect_that(pmartR::plot.proData(mat), throws_error()) 
})

test_that("invalid input for legend_position argument throws error",{  
  expect_that(plot(omicsData, legend_position = FALSE), throws_error())
})

test_that("invalid input for bw_theme argument throws error",{    
  expect_that(plot(omicsData, bw_theme = NULL), throws_error())
})

test_that("invalid input for title_size argument throws error",{    
  expect_that(plot(omicsData, title_size = "five"), throws_error())
})

test_that("invalid input for x_lab_size argument throws error",{    
  expect_that(plot(omicsData, x_lab_size = "five"), throws_error())
})

test_that("invalid input for y_lab_size argument throws error",{    
  expect_that(plot(omicsData, y_lab_size = "five"), throws_error())
})

test_that("invalid input for ylimit argument throws error",{    
  expect_that(plot(omicsData, ylimit = "one"), throws_error())
  expect_that(plot(omicsData, ylimit = 15), throws_error())
  expect_that(plot(omicsData, ylimit = c(1,2,3)), throws_error())
  expect_that(plot(omicsData, ylimit = c("one", "two")), throws_error())
})

