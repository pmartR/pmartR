context("Test normalization")

library(testthat)
library(pmartR)
library(pmartRdata)

pep_object <- pmartRdata::pep_object
pep_object <- edata_transform(pep_object, "log2") %>%
  group_designation(main_effects = "Condition")

subset_funs <- c("all", "los", "ppp", "rip", "ppp_rip", "complete")
norm_funs <- c("median", "mean", "zscore", "mad")

test_that("attribute integrity of normalize_global output", {
  for(norm in norm_funs){
    for(subset in subset_funs){
      # test attribute setting
      tempres <- normalize_global(pep_object, subset_fn = subset, norm_fn = norm, apply_norm = TRUE)
      tempnormRes <- normalize_global(pep_object, subset_fn = subset, norm_fn = norm, apply_norm = FALSE)
      
      expect_true(attributes(tempres)$data_info$norm_info$subset_fn == subset)
      expect_true(attributes(tempres)$data_info$norm_info$norm_fn == norm)
      expect_true(attributes(tempres)$data_info$norm_info$is_normalized == TRUE)
      expect_false(is.null(attributes(tempres)$data_info$norm_info$params$norm_location))
      
      # test test_normRes
      norm_test <- test_normRes(tempres)
      normRes_test <- test_normRes(tempnormRes)
      
      # should always be a scale p-value
      expect_true(!is.null(norm_test$p_location) & !is.null(normRes_test$p_location))
      
      # if there was a scale parameter, should be a scale p-value
      if(!is.null(attributes(tempres)$data_info$norm_info$params$norm_scale)){
        expect_true(!is.null(norm_test$p_scale))
      }
      if(!is.null(attributes(tempres)$data_info$norm_info$params$norm_scale)){
        expect_true(!is.null(normRes_test$p_scale))
      }
      
      # test backtransformation
      compare_edata <- tempres$e_data
      compare_edata[-which(colnames(compare_edata) == get_edata_cname(tempres))] <- sapply(1:(ncol(compare_edata)-1), function(i){
        if(!is.null(attributes(tempres)$data_info$norm_info$params$norm_scale)){
          (compare_edata[,i+1] * attributes(tempres)$data_info$norm_info$params$norm_scale[i]) + attributes(tempres)$data_info$norm_info$params$norm_location[i]
        }
        else compare_edata[,i+1] + attributes(tempres)$data_info$norm_info$params$norm_location[i]
      })
      
      expect_equal(pep_object$e_data, compare_edata)
    }
  }
})

