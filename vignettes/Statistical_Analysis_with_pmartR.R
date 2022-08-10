## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message = FALSE, fig.width = 10, fig.height = 8, eval = FALSE)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup,include=FALSE,warning=FALSE----------------------------------------
#  library(pmartR)
#  library(pmartRdata)
#  
#  rm(list = ls())
#  
#  data(pep_object)
#  
#  pep_object <- edata_replace(pep_object, 0, NA)
#  pep_object <- edata_transform(pep_object, "log2")
#  
#  # create some fake sample information to use in examples #
#  pep_object$f_data$Condition2 <- c("A","A","A","B","B","B","C","C","C","D","D","D")
#  pep_object$f_data$Condition3 <- c(rep("Plus",6), rep("Minus",6))

## ---- warning=FALSE-----------------------------------------------------------
#  pep_object <- group_designation(omicsData = pep_object, main_effects = c("Condition"))
#  myfilt <- imdanova_filter(omicsData = pep_object)
#  pep_object <- applyFilt(filter_object = myfilt, omicsData = pep_object, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)
#  all_pairwise_results <- anova_test(omicsData = pep_object)

## ---- warning=FALSE-----------------------------------------------------------
#  all_pairwise_results_adjusted <- anova_test(omicsData = pep_object, pval_adjust = "Tukey")

## ---- warning=FALSE-----------------------------------------------------------
#  pep_object$f_data$Condition2 <- c(rep("A",3), rep("B", 3), rep("C", 3), rep("D", 3))
#  pep_object <- group_designation(pep_object, main_effects = "Condition2")
#  one_vs_all <- data.frame(Control=rep("A",3),Test=c("B","C","D"))
#  one_vs_all_results <- anova_test(omicsData = pep_object, comparisons = one_vs_all)

## ---- warning=FALSE-----------------------------------------------------------
#  pep_object <- group_designation(omicsData = pep_object, main_effects = c("Condition"))
#  all_pairwise_results <- anova_test(omicsData = pep_object)

## ----statres------------------------------------------------------------------
#  stat_results <- imd_anova(omicsData = pep_object, test_method = "combined")
#  summary(stat_results)
#  print(stat_results)
#  plot(stat_results)
#  plot(stat_results, plot_type = "volcano")
#  

