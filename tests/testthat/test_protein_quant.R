#protein_quant testthat 
library(testthat)
library(pmartRqc)
library(pmartRdata)

data("pep_object")
data("lipid_object")
mypepData = edata_transform(pep_object, "log2")

zrollup_median = protein_quant(myobject, method = 'zrollup', combine_fn = 'median', single_observation = T)
zrollup_mean = protein_quant(myobject, method = 'zrollup', combine_fn = 'mean', single_observation = T)

vector = c(1,2,3)

#hardcoding results for zrollup
zroll_median_arpc4human = c(0.2169473, 0.0000000, -1.14799757, -0.76259778, -0.15596116, -0.08445067, -0.45365500, -1.4851329, -1.5295777, 1.24065259, 0.63670702, -0.6028498)
zroll_median_rt05human = c(0.3215534, 0.4297530, 0.70710678, -1.06047369, -0.13417357, 0.63022510, -0.40965528, -2.3200769, NA, -0.32155337, NA, NA)

zroll_mean_arpc4human = c(0.37681456, 0.001846431, -1.14799757, -0.7625978, -0.1559612, -0.08445067, -0.45365500, -1.4851329, -1.5295777, 1.2406526, 0.5739629, -0.2640281)
zroll_mean_rt05human = c(0.32155337, 0.429752990, 0.70710678, -1.0604737, -0.1341736, 0.63022510, -0.40965528, -2.3200769, NA, -0.3215534, NA, NA)
                         
context("output tests for protein_quant(method = 'zrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(zrollup_median), equals("proData"))
  expect_that(length(zrollup_median), equals(3))
  expect_that(length(attr(zrollup_median, "filters")), equals(1))
})

context("input tests for protein_quant(method = 'zrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = myobject, method = 'zrollup', single_observation = FALSE), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'zrollup', isoformRes = vector, single_observation = TRUE, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'zrollup', single_observation = TRUE, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'zrollup', single_pep = TRUE, single_observation = TRUE), throws_error())
  })

context("tests using hard code results of protein_quant(method = 'zrollup')")

test_that("output of protein_quant(combine_fn = 'median') matches hard code results",{
  expect_that(round(as.numeric(zrollup_median$e_data[which(zrollup_median$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(zroll_median_arpc4human, 4)))
  expect_that(round(as.numeric(zrollup_median$e_data[which(zrollup_median$e_data$Protein == "RT05_HUMAN"), -1]), 4) , equals(round(zroll_median_rt05human, 4)))
})

test_that("output of protein_quant(combine_fn = 'mean') matches hard code results",{
  expect_that(round(as.numeric(zrollup_mean$e_data[which(zrollup_mean$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(zroll_mean_arpc4human, 4)))
  expect_that(round(as.numeric(zrollup_mean$e_data[which(zrollup_mean$e_data$Protein == "RT05_HUMAN"), -1]), 4) , equals(round(zroll_mean_rt05human, 4)))
})


##################### tests for protein_quant(method = 'qrollup') #########################
qrollup_median = protein_quant(myobject, method = 'qrollup', qrollup_thresh = .4, combine_fn = 'median')
qrollup_mean = protein_quant(myobject, method = 'qrollup', qrollup_thresh = .4, combine_fn = 'mean')

#hardcoding results for qrollup 
qroll_median_arpc4human = c(0.2474660, 0.1517315, -0.6131726, 0.33727198, -0.6803733, 0.7763576, 0.5373058, -0.9645733, -0.1593309, 1.1989259, -0.3182641, 1.1663487)   
qroll_median_sucb2human = c(-0.4816177, -0.7681180, -1.0095532, -1.41528146, -0.1800738, -0.7017922, -0.8099622, -1.4231692, -1.1513830, -1.5620134, -2.0310164, -1.5592478)    

qroll_mean_arpc4human = c(0.24746604, 0.1517315, -0.6131726,  0.337271983, -0.6803733, 0.7763576, 0.53730575, -0.9645733, -0.1593309, 1.1989259, -0.3898157, 1.1663487)
qroll_mean_sucb2human = c(-0.59685702, -0.8711415, -1.0966185, -1.414361405, -0.4728151, -0.4873553, -0.83881699, -1.4231692, -1.8523838, -1.5793305, -1.5602705, -1.3551594)    


context("output tests for protein_quant(method = 'qrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(qrollup_median), equals("proData"))
  expect_that(length(qrollup_median), equals(3))
})

context("input tests for protein_quant(method = 'qrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', qrollup_thresh = "abcdef", combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', qrollup_thresh = .4, isoformRes = vector, combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', qrollup_thresh = .4, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', qrollup_thresh = 2, combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', qrollup_thresh = 'abcdef', combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'qrollup', combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'qrollup', qrollup_thresh = .4, combine_fn = "mean"), throws_error())
})

context("tests using hard code results of protein_quant(method = 'qrollup')")

test_that("output of protein_quant(combine_fn = 'median') matches hard code results",{
  expect_that(round(as.numeric(qrollup_median$e_data[which(qrollup_median$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(qroll_median_arpc4human, 4)))
  expect_that(round(as.numeric(qrollup_median$e_data[which(qrollup_median$e_data$Protein == "SUCB2_HUMAN"), -1]), 4) , equals(round(qroll_median_sucb2human, 4)))
 })

test_that("output of protein_quant(combine_fn = 'mean') matches hard code results",{
  expect_that(round(as.numeric(qrollup_mean$e_data[which(qrollup_mean$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(qroll_mean_arpc4human, 4)))
  expect_that(round(as.numeric(qrollup_mean$e_data[which(qrollup_mean$e_data$Protein == "SUCB2_HUMAN"), -1]), 4) , equals(round(qroll_mean_sucb2human, 4)))
})

##################### tests for protein_quant(method = 'rollup') #########################
rollup_median = protein_quant(myobject, method = 'rollup', combine_fn = "median")
rollup_mean = protein_quant(myobject, method = 'rollup', combine_fn = "mean")

#hardcoding results for rollup 
roll_median_arpc4human = c(-0.4765741, -0.5822546, -0.6131726, 0.3372720, -0.680373272, 0.7763576, 0.5373058, -2.605045, -0.1593309, -0.9349052, -1.1835781, -4.47183371)      
roll_median_sucb2human = c(-0.8871225, -0.9424352, -1.2426063, -1.6586355, -0.799533249, -0.8742998, -0.9773790, -1.423169, -2.9598858, -1.5620134, -3.2954900, -1.64914103)   

roll_mean_arpc4human = c(-1.0584259, -1.277392, -0.6131726, 0.3372720, -0.680373272, 0.7763576, 0.5373058, -2.605045, -0.1593309, -0.9349052, -1.4795607, -3.15941063)     
roll_mean_sucb2human = c(-1.0101463, -1.166656, -1.6978494, -1.8028265, -1.596926353, -1.3936373, -0.9915361, -1.423169, -2.7693071, -1.5793305, -3.1176134, -2.91821324)   

context("output tests for protein_quant(method = 'rollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(rollup_median), equals("proData"))
  expect_that(length(rollup_median), equals(3))
})                          

context("input tests for protein_quant(method = 'rollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = myobject, method = 'rollup', isoformRes = vector, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'rollup', combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = myobject, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'rollup', combine_fn = "median"), throws_error())
})                           
                           
context("tests using hard code results of protein_quant(method = 'rollup')")

test_that("output of protein_quant(combine_fn = 'median') matches hard code results",{
  expect_that(round(as.numeric(rollup_median$e_data[which(rollup_median$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(roll_median_arpc4human, 4)))
  expect_that(round(as.numeric(rollup_median$e_data[which(rollup_median$e_data$Protein == "SUCB2_HUMAN"), -1]), 4) , equals(round(roll_median_sucb2human, 4)))
})

test_that("output of protein_quant(combine_fn = 'mean') matches hard code results",{
  expect_that(round(as.numeric(rollup_mean$e_data[which(rollup_mean$e_data$Protein == "ARPC4_HUMAN"), -1]), 4) , equals(round(roll_mean_arpc4human, 4)))
  expect_that(round(as.numeric(rollup_mean$e_data[which(rollup_mean$e_data$Protein == "SUCB2_HUMAN"), -1]), 4) , equals(round(roll_mean_sucb2human, 4)))
})                           
                           
                    
##################### tests for protein_quant(method = 'rrollup') #########################
rrollup_result = protein_quant(myobject, method = 'rrollup', combine_fn = 'median')

context("output tests for protein_quant(method = 'rrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(rrollup_result), equals("proData"))
  expect_that(length(rrollup_result), equals(3))
})

context("input tests for protein_quant(method = 'rrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = myobject, method = 'rrollup', isoformRes = vector, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = myobject, method = 'rrollup', combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'rrollup', combine_fn = "median"), throws_error())
})
