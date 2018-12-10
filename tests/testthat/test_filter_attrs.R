context("attributes proteomics filter")
# testing function for cv_filter() 
library(pmartR)
library(testthat)
library(pmartRdata)
library(dplyr)

sink(file = "logs_test_filter_attrs")

myfilt <- proteomics_filter(pep_object)
prefilt_n_edata <- attributes(pep_object)$data_info$num_edata
prefilt_n_emeta <- attributes(pep_object)$data_info$num_emeta

group_pep <- group_designation(pep_object, main_effects = "Condition")
groupfilt <- proteomics_filter(group_pep)
group_pep_filtered <- applyFilt(groupfilt, group_pep, min_num_peps = 2)
# summary(groupfilt, min_num_peps = 2)

test_that("apply proteomics filter and check summary", {
  for(i in 2:4){
      degen_peps = as.logical(i%%2)
      pepsummary <- summary(myfilt, min_num_peps = i)
      pep_filtered <- applyFilt(myfilt, pep_object, min_num_peps = i, degen_peps = degen_peps)
      
      expect_equal(prefilt_n_edata - pepsummary$num_pep_filtered, attributes(pep_filtered)$data_info$num_edata)
      expect_true((prefilt_n_emeta - pepsummary$num_pro_filtered) >= attributes(pep_filtered)$data_info$num_emeta)
      
      na_sums <- pep_filtered$e_data %>% dplyr::select(-one_of(get_edata_cname(pep_filtered))) %>% {!is.na(.)} %>% rowSums()
      
      expect_false(any(na_sums == 0))
    }    
  })

sink()

omicsData <- pro_object
omicsData$f_data$fakegroup <- c(rep("A",3), rep("B",4), rep("C",4))

# Weird groups
# omicsData$f_data$fakegroup <- 1:12
# omicsData$f_data$fakegroup2 <- c(1:8, 1, 9:11)

omicsData_fakegroup <- group_designation(omicsData = omicsData, main_effects = c("fakegroup"))
omicsData_conditiongroup <- group_designation(omicsData = omicsData, main_effects = c("Condition"))

imdanovafilt_fakegroup <- imdanova_filter(omicsData = omicsData_fakegroup)
imdanovafilt_conditiongroup <- imdanova_filter(omicsData = omicsData_conditiongroup)

nonmiss_params <- list(list(2, NULL), list(NULL, 3), list(2, 3))

lapply(nonmiss_params, function(params){
    omicsData_fakegroup <- applyFilt(filter_object = imdanovafilt_fakegroup, omicsData = omicsData_fakegroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
    omicsData_conditiongroup <- applyFilt(filter_object = imdanovafilt_conditiongroup, omicsData = omicsData_conditiongroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
    
    sink(file = "logs_statres_filter_summary")
    fakegroup_summary <- summary(imdanovafilt_fakegroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
    conditiongroup_summary <- summary(imdanovafilt_conditiongroup, min_nonmiss_anova = params[[1]], min_nonmiss_gtest = params[[2]])
    sink()
    
    #test correct number filtered
    expect_equal(attr(omicsData_fakegroup, "data_info")$num_edata, fakegroup_summary$num_not_filtered)
    expect_equal(attr(omicsData_conditiongroup, "data_info")$num_edata, conditiongroup_summary$num_not_filtered)
    
    expect_equal(length(attr(omicsData_fakegroup, "filters")$imdanovaFilt$filtered), fakegroup_summary$num_filtered)
    expect_equal(length(attr(omicsData_conditiongroup, "filters")$imdanovaFilt$filtered), conditiongroup_summary$num_filtered)
    
})
  

