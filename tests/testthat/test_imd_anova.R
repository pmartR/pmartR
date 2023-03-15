context('IMD-ANOVA tests')

test_that('all tests conform to the decrees of the God of Stats',{

  # Load IMD-ANOVA filter objects ----------------------------------------------

  load(system.file('testdata',
                   'standards_filter.RData',
                   package = 'pmartR'))

  # Load IMD-ANOVA standards ---------------------------------------------------

  load(system.file('testdata',
                   'standards_imd_anova.RData',
                   package = 'pmartR'))

  # Compare IMD-ANOVAly --------------------------------------------------------

  # Default comparisons ---------------

  afruit_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova")
  gfruit_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest")
  cfruit_1_0_2 <- imd_anova(cfilta_1_0_2, test_method = "combined")

  afruit_1_1_3 <- imd_anova(afilta_1_1_3,
                            test_method = "anova")
  gfruit_1_1_3 <- imd_anova(gfilta_1_1_3,
                            test_method = "gtest",
                            parallel = FALSE)
  cfruit_1_1_3 <- imd_anova(cfilta_1_1_3,
                            test_method = "combined")

  afruit_1_2_3 <- imd_anova(afilta_1_2_3,
                            test_method = "anova")
  gfruit_1_2_3 <- imd_anova(gfilta_1_2_3,
                            test_method = "gtest",
                            parallel = FALSE)
  cfruit_1_2_3 <- imd_anova(cfilta_1_2_3,
                            test_method = "combined")

  afruit_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova")
  gfruit_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest")
  cfruit_2_0_3 <- imd_anova(cfilta_2_0_3, test_method = "combined")

  afruit_2_1_4 <- imd_anova(afilta_2_1_4,
                            test_method = "anova")
  gfruit_2_1_4 <- imd_anova(gfilta_2_1_4,
                            test_method = "gtest",
                            parallel = FALSE)
  cfruit_2_1_4 <- imd_anova(cfilta_2_1_4,
                            test_method = "combined")

  afruit_2_2_4 <- imd_anova(afilta_2_2_4,
                            test_method = "anova")
  gfruit_2_2_4 <- imd_anova(gfilta_2_2_4,
                            test_method = "gtest",
                            parallel = FALSE)
  cfruit_2_2_4 <- imd_anova(cfilta_2_2_4,
                            test_method = "combined")

  # Custom comparisons ---------------

  afruit_cus_2_0_3 <- imd_anova(omicsData = afilta_2_0_3,
                                comparisons = data.frame(
                                  Control = c("Infection_low",
                                              "Mock_none",
                                              "Mock_none"),
                                  Test = c("Infection_high",
                                           "Infection_high",
                                           "Infection_low")
                                ),
                                test_method = "anova")

  gfruit_cus_2_0_3 <- imd_anova(omicsData = gfilta_2_0_3,
                                comparisons = data.frame(
                                  Control = c("Infection_low",
                                              "Mock_none",
                                              "Mock_none"),
                                  Test = c("Infection_high",
                                           "Infection_high",
                                           "Infection_low")
                                ),
                                test_method = "gtest")

  # ANOVA: Adjusted p-values ---------------

  afruit_bon_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  afruit_dun_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")

  afruit_bon_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(3)
  afruit_dun_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")

  afruit_bon_1_2_3 <- imd_anova(afilta_1_2_3, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_1_2_3 <- imd_anova(afilta_1_2_3, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_1_2_3 <- imd_anova(afilta_1_2_3, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(5)
  afruit_dun_1_2_3 <- imd_anova(afilta_1_2_3, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")

  afruit_bon_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(4)
  afruit_dun_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")
  # Still gives an error even though all case-vs-control comparisons are being
  # made. Hmmmmmmmm. I took out the if statement that checks if a Dunnett
  # correction can be used because it throws an error when the conditions to be
  # used are met.

  afruit_bon_2_1_4 <- imd_anova(afilta_2_1_4, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_2_1_4 <- imd_anova(afilta_2_1_4, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_2_1_4 <- imd_anova(afilta_2_1_4, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(8)
  afruit_dun_2_1_4 <- imd_anova(afilta_2_1_4, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")

  afruit_bon_2_2_4 <- imd_anova(afilta_2_2_4, test_method = "anova",
                                pval_adjust_a_multcomp = "bonferroni")
  afruit_holm_2_2_4 <- imd_anova(afilta_2_2_4, test_method = "anova",
                                 pval_adjust_a_multcomp = "holm")
  afruit_tuk_2_2_4 <- imd_anova(afilta_2_2_4, test_method = "anova",
                                pval_adjust_a_multcomp = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(11)
  afruit_dun_2_2_4 <- imd_anova(afilta_2_2_4, test_method = "anova",
                                pval_adjust_a_multcomp = "dunnett")

  # G-Test: Adjusted p-values ---------------

  gfruit_bon_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest",
                                pval_adjust_g_multcomp = "bonferroni")
  gfruit_holm_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest",
                                 pval_adjust_g_multcomp = "holm")
  gfruit_bon_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                pval_adjust_g_multcomp = "bonferroni")
  gfruit_holm_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                 pval_adjust_g_multcomp = "holm")
  gfruit_bon_2_1_4 <- imd_anova(gfilta_2_1_4, test_method = "gtest",
                                pval_adjust_g_multcomp = "bonferroni")
  gfruit_holm_2_1_4 <- imd_anova(gfilta_2_1_4, test_method = "gtest",
                                 pval_adjust_g_multcomp = "holm")

  # Holy IMD-ANOVA unit tests, Statman! ----------------------------------------

  # ANOVA: Unadjusted p-values ---------------

  expect_equal(afruit_1_0_2, astan_1_0_2)
  expect_equal(afruit_1_1_3, astan_1_1_3)
  expect_equal(afruit_1_2_3, astan_1_2_3)
  expect_equal(afruit_2_0_3, astan_2_0_3)
  expect_equal(afruit_2_1_4, astan_2_1_4)
  expect_equal(afruit_2_2_4, astan_2_2_4)

  # Make certain the custom comparisons create same output as the default when
  # all possible comparisons are used.
  expect_equal(afruit_2_0_3, afruit_cus_2_0_3)

  # Assure adjusted p-values are the same as unadjusted p-values when there
  # is only one test between two groups.
  expect_equal(afruit_1_0_2$P_value_G_Infection_vs_Mock,
               afruit_bon_1_0_2$P_value_G_Infection_vs_Mock)
  expect_equal(afruit_1_0_2$P_value_G_Infection_vs_Mock,
               afruit_holm_1_0_2$P_value_G_Infection_vs_Mock)
  expect_equal(afruit_1_0_2$P_value_G_Infection_vs_Mock,
               afruit_tuk_1_0_2$P_value_G_Infection_vs_Mock)
  expect_equal(afruit_1_0_2$P_value_G_Infection_vs_Mock,
               afruit_dun_1_0_2$P_value_G_Infection_vs_Mock)

  # ANOVA: Adjusted p-values ---------------

  # NOTE: Currently commented out because the current version of imd_anova
  # produces p-values greater than 1. Holy inconceivable p-values, Statman!
  expect_equal(
    data.matrix(afruit_bon_1_1_3[, 11:13]),
    pmin(data.matrix(astan_1_1_3[, 11:13] * 3), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  expect_equal(
    unclass(afruit_holm_1_1_3[, 11:13]),
    unclass(
      data.frame(t(apply(astan_1_1_3[, 11:13],
                         1,
                         p.adjust,
                         method = "holm")))
    )
  )
  expect_equal(
    data.frame(afruit_tuk_1_1_3[, 11:13]),
    tukey_pval_1_1_3
  )
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_1_1_3[, 11:13]) -
            data.matrix(dunnett_1_1_3)),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # NOTE: Currently commented out because the current version of imd_anova
  # produces p-values greater than 1. Holy inconceivable p-values, Statman!
  expect_equal(
    data.matrix(afruit_bon_1_2_3[, 11:13]),
    pmin(data.matrix(astan_1_2_3[, 11:13] * 3), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  expect_equal(
    unclass(afruit_holm_1_2_3[, 11:13]),
    unclass(
      data.frame(t(apply(astan_1_2_3[, 11:13],
                         1,
                         p.adjust,
                         method = "holm")))
    )
  )
  expect_equal(
    data.frame(afruit_tuk_1_2_3[, 11:13]),
    tukey_pval_1_2_3
  )
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_1_2_3[, 11:13]) -
            data.matrix(dunnett_1_2_3)),
      na.rm = TRUE
    ) < 0.0001132809
  )

  expect_equal(
    data.matrix(afruit_bon_2_0_3[, 11:13]),
    pmin(data.matrix(astan_2_0_3[, 11:13] * 3), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  expect_equal(
    unclass(afruit_holm_2_0_3[, 11:13]),
    unclass(
      data.frame(t(apply(astan_2_0_3[, 11:13],
                         1,
                         p.adjust,
                         method = "holm")))
    )
  )
  expect_equal(
    data.frame(afruit_tuk_2_0_3[, 11:13]),
    tukey_pval_2_0_3
  )
  # Because of the random process in the mvtnorm::pmvt function we test the
  # Dunnett adjusted p-values using a threshold. The threshold is 5 standard
  # deviations above the mean difference between runs of the mvtnorm:pmvt
  # function. To create the threshold the Dunnett adjusted p-values were
  # calculated 10,001 times (the first run is used as a reference) using the
  # following code: imd_anova(afilta_2_0_3, test_method = "anova", 
  # pval_adjust_a_multcomp = "dunnett"). We calculated the mean of the absolute
  # difference between the p-values from the initial run and each subsequent run.
  # The mean and standard deviation of the mean absolute difference were then
  # calculated for the 10,000 runs.
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_2_0_3[, 11:13]) -
            data.matrix(dunnett_2_0_3)),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # NOTE: Currently commented out because the current version of imd_anova
  # produces p-values greater than 1. Holy inconceivable p-values, Statman!
  expect_equal(
    data.matrix(afruit_bon_2_1_4[, 16:21]),
    pmin(data.matrix(astan_2_1_4[, 16:21] * 6), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  expect_equal(
    unclass(afruit_holm_2_1_4[, 16:21]),
    unclass(
      data.frame(t(apply(astan_2_1_4[, 16:21],
                         1,
                         p.adjust,
                         method = "holm")))
    )
  )
  expect_equal(
    data.frame(afruit_tuk_2_1_4[, 16:21]),
    tukey_pval_2_1_4
  )
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_2_1_4[, 16:21]) -
            data.matrix(dunnett_2_1_4)),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # NOTE: Currently commented out because the current version of imd_anova
  # produces p-values greater than 1. Holy inconceivable p-values, Statman!
  expect_equal(
    data.matrix(afruit_bon_2_2_4[, 16:21]),
    pmin(data.matrix(astan_2_2_4[, 16:21] * 6), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  expect_equal(
    unclass(afruit_holm_2_2_4[, 16:21]),
    unclass(
      data.frame(t(apply(astan_2_2_4[, 16:21],
                         1,
                         p.adjust,
                         method = "holm")))
    )
  )
  expect_equal(
    data.frame(afruit_tuk_2_2_4[, 16:21]),
    tukey_pval_2_2_4
  )
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_2_2_4[, 16:21]) -
            data.matrix(dunnett_2_2_4)),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # ANOVA: FDR adjust ---------------------------
  afruit_fdr_bon_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova", pval_adjust_a_fdr = "bonferroni")
  expect_equal(afruit_fdr_bon_1_0_2$P_value_A_Infection_vs_Mock,
               c(apply(data.frame(afruit_1_0_2$P_value_A_Infection_vs_Mock), 2, p.adjust, method = "bonferroni")))
  
  afruit_fdr_bon_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova", pval_adjust_a_fdr = "bonferroni")
  expect_equal(afruit_fdr_bon_1_1_3$P_value_A_zombie_vs_human,
               c(apply(data.frame(afruit_1_1_3$P_value_A_zombie_vs_human), 2, p.adjust, method = "bonferroni")))
  
  afruit_fdr_bon_1_2_3 <- imd_anova(afilta_1_2_3, test_method = "anova", pval_adjust_a_fdr = "bonferroni")
  expect_equal(afruit_fdr_bon_1_2_3$P_value_A_zombie_vs_human,
               c(apply(data.frame(afruit_1_2_3$P_value_A_zombie_vs_human), 2, p.adjust, method = "bonferroni")))

  afruit_fdr_bon_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova", pval_adjust_a_fdr = "bonferroni")
  expect_equal(afruit_fdr_bon_2_0_3$P_value_A_Infection_high_vs_Mock_none,
               c(apply(data.frame(afruit_2_0_3$P_value_A_Infection_high_vs_Mock_none), 2, p.adjust, method = "bonferroni")))

  afruit_fdr_bon_2_1_4 <- imd_anova(afilta_2_1_4, test_method = "anova", pval_adjust_a_fdr = "bonferroni")
  expect_equal(afruit_fdr_bon_2_1_4$P_value_A_Infection_high_vs_Infection_low,
               c(apply(data.frame(afruit_2_1_4$P_value_A_Infection_high_vs_Infection_low), 2, p.adjust, method = "bonferroni")))

  # G-Test: Unadjusted p-values ---------------

  expect_equal(gfruit_1_0_2, gstan_1_0_2)
  expect_equal(gfruit_1_1_3, gstan_1_1_3)
  expect_equal(gfruit_1_2_3, gstan_1_2_3)
  expect_equal(gfruit_2_0_3, gstan_2_0_3)
  expect_equal(gfruit_2_1_4, gstan_2_1_4)
  expect_equal(gfruit_2_2_4, gstan_2_2_4)

  # Ensure the custom comparisons create same output as the default when all
  # possible comparisons are used.
  expect_equal(gfruit_2_0_3, gfruit_cus_2_0_3)

  # Make sure adjusted p-values are the same as unadjusted p-values when there
  # is only one test between two groups.
  expect_equal(gfruit_1_0_2$P_value_G_Infection_vs_Mock,
               gfruit_bon_1_0_2$P_value_G_Infection_vs_Mock)
  expect_equal(gfruit_1_0_2$P_value_G_Infection_vs_Mock,
               gfruit_holm_1_0_2$P_value_G_Infection_vs_Mock)

  # G-Test: Adjusted p-values ---------------

  expect_equal(
    data.frame(gfruit_bon_2_0_3[, 11:13]),
    data.frame(pmin(data.matrix(gstan_2_0_3[, 11:13] * 3), 1))
  )
  expect_equal(
    data.frame(gfruit_holm_2_0_3[, 11:13]),
    data.frame(t(apply(gstan_2_0_3[, 11:13],
                       1,
                       p.adjust,
                       method = "holm")))
  )

  expect_equal(
    data.frame(gfruit_bon_2_1_4[, 16:21]),
    data.frame(pmin(data.matrix(gstan_2_1_4[, 16:21] * 6), 1))
  )
  expect_equal(
    data.frame(gfruit_holm_2_1_4[, 16:21]),
    data.frame(t(apply(gstan_2_1_4[, 16:21],
                       1,
                       p.adjust,
                       method = "holm")))
  )
  
  # G-Test: FDR adjust --------------------------
  gfruit_fdr_bon_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest", pval_adjust_g_fdr = "bonferroni")
  expect_equal(gfruit_fdr_bon_1_0_2$P_value_G_Infection_vs_Mock,
               c(apply(data.frame(gfruit_1_0_2$P_value_G_Infection_vs_Mock), 2, p.adjust, method = "bonferroni")))
  
  gfruit_fdr_bon_1_1_3 <- imd_anova(gfilta_1_1_3, test_method = "gtest", pval_adjust_g_fdr = "bonferroni")
  expect_equal(gfruit_fdr_bon_1_1_3$P_value_G_zombie_vs_human,
               c(apply(data.frame(gfruit_1_1_3$P_value_G_zombie_vs_human), 2, p.adjust, method = "bonferroni")))
  
  gfruit_fdr_bon_1_2_3 <- imd_anova(gfilta_1_2_3, test_method = "gtest", pval_adjust_g_fdr = "bonferroni")
  expect_equal(gfruit_fdr_bon_1_2_3$P_value_G_zombie_vs_human,
               c(apply(data.frame(gfruit_1_2_3$P_value_G_zombie_vs_human), 2, p.adjust, method = "bonferroni")))
  
  gfruit_fdr_bon_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest", pval_adjust_g_fdr = "bonferroni")
  expect_equal(gfruit_fdr_bon_2_0_3$P_value_G_Infection_high_vs_Mock_none,
               c(apply(data.frame(gfruit_2_0_3$P_value_G_Infection_high_vs_Mock_none), 2, p.adjust, method = "bonferroni")))
  
  gfruit_fdr_bon_2_1_4 <- imd_anova(gfilta_2_1_4, test_method = "gtest", pval_adjust_g_fdr = "bonferroni")
  expect_equal(gfruit_fdr_bon_2_1_4$P_value_G_Infection_high_vs_Infection_low,
               c(apply(data.frame(gfruit_2_1_4$P_value_G_Infection_high_vs_Infection_low), 2, p.adjust, method = "bonferroni")))
  
  # Combined: Unadjusted p-values ---------------

  expect_equal(cfruit_1_0_2, cstan_1_0_2)
  expect_equal(cfruit_1_1_3, cstan_1_1_3)
  expect_equal(cfruit_1_2_3, cstan_1_2_3)
  expect_equal(cfruit_2_0_3, cstan_2_0_3)
  expect_equal(cfruit_2_1_4, cstan_2_1_4)
  expect_equal(cfruit_2_2_4, cstan_2_2_4)
  
  # Test argument checks -------------------------------------------------------

  # Maybe add later (depending on the position of the planets and stars).

})
