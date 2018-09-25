context("attributes proteomics filter")
# testing function for cv_filter() 
library(pmartR)
library(testthat)
library(pmartRdata)

sink(file = "logs_test_filter_attrs")

data("pep_object")
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