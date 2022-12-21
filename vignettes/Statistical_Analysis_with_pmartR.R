## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message = FALSE, fig.width = 6, fig.height = 5, eval = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, warning=FALSE-----------------------------------------------------
library(pmartR)
library(pmartRdata)

pep_object <- edata_transform(omicsData = pep_object, 
                              data_scale = "log2")

## ---- warning=FALSE-----------------------------------------------------------
pep_object <- group_designation(omicsData = pep_object, 
                                main_effects = c("Phenotype"))
myfilt <- imdanova_filter(omicsData = pep_object)
pep_object <- applyFilt(filter_object = myfilt, 
                        omicsData = pep_object, 
                        min_nonmiss_anova = 2, 
                        min_nonmiss_gtest = 3)
all_pairwise_results <- imd_anova(omicsData = pep_object, 
                                  test_method = "anova")

## ---- warning=FALSE-----------------------------------------------------------
all_pairwise_results_adjusted <- imd_anova(omicsData = pep_object, 
                                           test_method = "anova", 
                                           pval_adjust_a = "Tukey")

## ---- warning=FALSE-----------------------------------------------------------
one_vs_all <- data.frame(Control=rep("Phenotype1",2), 
                         Test=c("Phenotype2","Phenotype3"))
one_vs_all_results <- imd_anova(omicsData = pep_object, 
                                test_method = "anova", 
                                comparisons = one_vs_all)

## ---- warning=FALSE-----------------------------------------------------------
pep_object <- group_designation(omicsData = pep_object, 
                                main_effects = c("SecondPhenotype"))
all_pairwise_results <- imd_anova(omicsData = pep_object, 
                                  test_method = "anova")

## -----------------------------------------------------------------------------
# add fake pair info to f_data using SecondPhenotype
pep_object$f_data$Pairing <- pep_object$f_data$SecondPhenotype
pep_object$f_data$Pair_ID <- c(rep(1:4, 2), rep(5:8, 2), rep(9:12,2))
mypaired <- group_designation(omicsData = pep_object,
                              pair_group = "Pairing",
                              pair_id = "Pair_ID",
                              pair_denom = "B")

myfilt <- imdanova_filter(omicsData = mypaired)
mypaired <- applyFilt(filter_object = myfilt, 
                      omicsData = mypaired, 
                      min_nonmiss_anova = 2,
                      min_nonmiss_gtest = 3)

paired_stats <- imd_anova(omicsData = mypaired,
                          test_method = "anova")

## -----------------------------------------------------------------------------
mypep_covariates <- group_designation(pep_object, main_effects = "Phenotype", covariates = "Characteristic")
imd_covariates <- imd_anova(omicsData = mypep_covariates,
                            test_method = "gtest")

## -----------------------------------------------------------------------------
# add fake pair info to f_data using SecondPhenotype
mypaired <- group_designation(omicsData = pep_object,
                              main_effects = "Phenotype",
                              pair_group = "Pairing",
                              pair_id = "Pair_ID",
                              pair_denom = "B")

myfilt <- imdanova_filter(omicsData = mypaired)
mypaired <- applyFilt(filter_object = myfilt, 
                      omicsData = mypaired, 
                      min_nonmiss_anova = 2,
                      min_nonmiss_gtest = 3)

paired_stats <- imd_anova(omicsData = mypaired,
                          test_method = "gtest")

## ----statres------------------------------------------------------------------
stat_results <- imd_anova(omicsData = pep_object, 
                          test_method = "combined")
summary(stat_results)
print(stat_results)
plot(stat_results)
plot(stat_results, plot_type = "volcano")

