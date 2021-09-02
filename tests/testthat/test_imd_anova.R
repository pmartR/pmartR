context('IMD-ANOVA tests')

test_that('all tests conform to the decrees of the God of Stats',{
  
  # Load data and prepare omicsData objects ------------------------------------
  
  # Paired data ---------------
  
  load(system.file('testdata',
                   'little_pairdata.RData',
                   package = 'pmartR'))
  
  pairdata <- as.pepData(e_data = edata,
                         f_data = fdata,
                         e_meta = emeta,
                         edata_cname = "Mass_Tag_ID",
                         fdata_cname = "Name",
                         emeta_cname = "Protein")
  
  # Throw some natural logs on the data because the following functions require
  # the data to be logged.
  pairdata <- edata_transform(pairdata, "log")
  
  # Not so paired data ---------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")
  
  # Create additional f_data objects with different main effects and covariates.
  fdata2 <- fdata
  set.seed(22250)
  fdata2$Dose <- c(sample(c("high", "low"), 9, replace = TRUE),
                   "none", "none", "none")
  fdata2$Gender <- sample(0:1, 12, replace = TRUE)
  
  # Create a vector that has three different levels for Condition.
  fdataZ <- fdata
  set.seed(720)
  fdataZ$Condition <- c(sample(c("zombie", "mutant"), 9, replace = TRUE),
                        "human", "human", "human")
  fdataZ$Gender <- sample(0:1, 12, replace = TRUE)
  
  # Create pepData objects, log the data, add group_DF attributes --------------
  
  # NOTE: There are three numbers that follow each object's name. The first
  # number represents the number of main effects. The second number represents
  # the number of covariates. The third represents the number of groups (across
  # all main effects).
  
  # One main effect, zero covariates, two groups
  pdata_1_0_2 <- pdata
  pdata_1_0_2 <- edata_transform(pdata_1_0_2, "log")
  pdata_1_0_2 <- group_designation(omicsData = pdata_1_0_2,
                                   main_effects = "Condition")
  groupDF_1_0_2 <- attr(pdata_1_0_2, "group_DF")
  
  # Two main effects, zero covariates, three groups
  pdata_2_0_3 <- as.pepData(e_data = edata,
                            f_data = fdata2,
                            edata_cname = "Mass_Tag_ID",
                            fdata_cname = "SampleID")
  pdata_2_0_3 <- edata_transform(pdata_2_0_3, "log")
  pdata_2_0_3 <- group_designation(omicsData = pdata_2_0_3,
                                   main_effects = c("Condition", "Dose"))
  groupDF_2_0_3 <- attr(pdata_2_0_3, "group_DF")
  
  # One main effect, one covariate, three groups
  pdata_1_1_3 <- as.pepData(e_data = edata,
                            f_data = fdataZ,
                            e_meta = emeta,
                            edata_cname = "Mass_Tag_ID",
                            fdata_cname = "SampleID",
                            emeta_cname = "Protein")
  pdata_1_1_3 <- edata_transform(pdata_1_1_3, "log")
  pdata_1_1_3 <- group_designation(omicsData = pdata_1_1_3,
                                   main_effects = "Condition",
                                   covariates = "Gender")
  groupDF_1_1_3 <- attr(pdata_1_1_3, "group_DF")
  
  # Filter IMD-ANOVAly ---------------------------------------------------------
  
  # Create filtas ---------------
  
  filta_1_0_2 <- imdanova_filter(pdata_1_0_2)
  afilta_1_0_2 <- applyFilt(filta_1_0_2, pdata_1_0_2,
                            min_nonmiss_anova = 2,
                            remove_singleton_groups = FALSE)
  gfilta_1_0_2 <- applyFilt(filta_1_0_2, pdata_1_0_2,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)
  cfilta_1_0_2 <- applyFilt(filta_1_0_2, pdata_1_0_2,
                            min_nonmiss_anova = 2,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)
  
  filta_2_0_3 <- imdanova_filter(pdata_2_0_3)
  afilta_2_0_3 <- applyFilt(filta_2_0_3, pdata_2_0_3,
                            min_nonmiss_anova = 2,
                            remove_singleton_groups = FALSE)
  gfilta_2_0_3 <- applyFilt(filta_2_0_3, pdata_2_0_3,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)
  cfilta_2_0_3 <- applyFilt(filta_2_0_3, pdata_2_0_3,
                            min_nonmiss_anova = 2,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)
  
  # ANOVA and G-test functions -------------------------------------------------
  
  # A function for calculating the g-test statistic for two groups.
  g <- function (obs1, obs2, abs1, abs2, n1, n2, n_total) {
    
    samurai_obs <- obs1 * log(obs1 / ((obs1 + obs2) * n1 / n_total))
    samurai_obs[is.na(samurai_obs)] <- 0
    samurai_abs <- abs1 * log(abs1 / ((abs1 + abs2) * n1 / n_total))
    samurai_abs[is.na(samurai_abs)] <- 0
    
    samurai <- samurai_obs + samurai_abs
    
    katana_obs <- obs2 * log(obs2 / ((obs1 + obs2) * n2 / n_total))
    katana_obs[is.na(katana_obs)] <- 0
    katana_abs <- abs2 * log(abs2 / ((abs1 + abs2) * n2 / n_total))
    katana_abs[is.na(katana_abs)] <- 0
    
    katana <- katana_obs + katana_abs
    
    fierce <- unname(samurai + katana)
    
    return (2 * fierce)
    
  }
  
  gflag <- function (obs1, obs2, abs1, abs2, pvals, cutoff) {
    
    maiden <- obs1 / (obs1 + abs1)
    
    dragon <- obs2 / (obs2 + abs2)
    
    knight <- sign(maiden - dragon)
    
    knight[which(pvals >= cutoff)] <- 0
    
    return (2 * unname(knight))
    
  }
  
  # Construct count matrices ---------------------------------------------------
  
  # Objects that start with obs represent counts of non-missing (observed)
  # values and objects that start with abs represent counts of missing (absent)
  # values.
  
  obs_inf_1_0_2 <- rowSums(!is.na(gfilta_1_0_2$e_data[, 2:10]))
  obs_mock_1_0_2 <- rowSums(!is.na(gfilta_1_0_2$e_data[, 11:13]))
  abs_inf_1_0_2 <- rowSums(is.na(gfilta_1_0_2$e_data[, 2:10]))
  abs_mock_1_0_2 <- rowSums(is.na(gfilta_1_0_2$e_data[, 11:13]))
  
  obs_inf_high_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, c(2:4, 7:8)]))
  obs_inf_low_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, c(5:6, 9:10)]))
  obs_mock_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, 11:13]))
  abs_inf_high_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, c(2:4, 7:8)]))
  abs_inf_low_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, c(5:6, 9:10)]))
  abs_mock_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, 11:13]))
  
  # Assemble standards ---------------------------------------------------------
  
  pval_g_1_0_2 <- data.frame(
    Infection_vs_Mock = pchisq(
      q = g(obs1 = obs_inf_1_0_2,
            obs2 = obs_mock_1_0_2,
            abs1 = abs_inf_1_0_2,
            abs2 = abs_mock_1_0_2,
            n1 = 9,
            n2 = 3,
            n_total = 9 + 3),
      df = 1,
      lower.tail = FALSE
    )
  )
  
  flag_g_1_0_2 <- data.frame(
    Infection_vs_Mock = gflag(
      obs1 = obs_inf_1_0_2,
      obs2 = obs_mock_1_0_2,
      abs1 = abs_inf_1_0_2,
      abs2 = abs_mock_1_0_2,
      pvals = pval_g_1_0_2[, 1],
      cutoff = 0.05
    )
  )
  
  mean_1_0_2 <- data.frame(
    Mean_Infection = rowMeans(gfilta_1_0_2$e_data[, 2:10],
                              na.rm = TRUE),
    Mean_Mock = rowMeans(gfilta_1_0_2$e_data[, 11:13],
                         na.rm = TRUE)
  )
  
  gstan_1_0_2 <- list(
    Full_results = data.frame(
      Mass_Tag_ID = gfilta_1_0_2$e_data$Mass_Tag_ID,
      Count_Infection = unname(obs_inf_1_0_2),
      Count_Mock = unname(obs_mock_1_0_2),
      P_value_G_Infection_vs_Mock = pval_g_1_0_2[, 1],
      Flag_Infection_vs_Mock = flag_g_1_0_2[, 1],
      mean_1_0_2,
      Fold_change_Infection_vs_Mock = (
        mean_1_0_2[, 1] - mean_1_0_2[, 2]
      ),
      row.names = NULL
    ),
    
    Flags = data.frame(
      Mass_Tag_ID = gfilta_1_0_2$e_data$Mass_Tag_ID,
      flag_g_1_0_2,
      row.names = NULL
    ),
    
    P_values = data.frame(
      Mass_Tag_ID = gfilta_1_0_2$e_data$Mass_Tag_ID,
      pval_g_1_0_2,
      row.names = NULL
    )
  )
  
  class(gstan_1_0_2) <- "statRes"
  
  attr(gstan_1_0_2, "group_DF") <- groupDF_1_0_2
  attr(gstan_1_0_2, "comparisons") <- c("Infection_vs_Mock")
  attr(gstan_1_0_2, "number_significant") <- data.frame(
    Comparison = c("Infection_vs_Mock"),
    Up_total = c(length(which(flag_g_1_0_2[, 1] == 2))),
    Down_total = c(length(which(flag_g_1_0_2[, 1] == -2))),
    Up_anova = c(0),
    Down_anova = c(0),
    Up_gtest = c(length(which(flag_g_1_0_2[, 1] == 2))),
    Down_gtest = c(length(which(flag_g_1_0_2[, 1] == -2))),
    row.names = c("Infection_vs_Mock")
  )
  attr(gstan_1_0_2, "statistical_test") <- "gtest"
  attr(gstan_1_0_2, "adjustment_method") <- "none"
  attr(gstan_1_0_2, "pval_thresh") <- 0.05
  attr(gstan_1_0_2, "data_info") <- list(
    data_scale_orig = "abundance",
    data_scale = "log",
    norm_info = list(is_normalized = FALSE),
    num_edata = length(unique(gfilta_1_0_2$e_data$Mass_Tag_ID)),
    num_miss_obs = sum(is.na(gfilta_1_0_2$e_data[, -1])),
    prop_missing = (sum(is.na(gfilta_1_0_2$e_data[, -1])) /
                      prod(dim(gfilta_1_0_2$e_data[, -1]))),
    num_samps = dim(gfilta_1_0_2$f_data)[1],
    data_types = NULL
  )
  attr(gstan_1_0_2, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = "Protein",
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(gstan_1_0_2, "data_class") <- "pepData"
  
  pval_g_2_0_3 <- data.frame(
    Infection_high_vs_Infection_low = pchisq(
      q = g(obs1 = obs_inf_high_2_0_3,
            obs2 = obs_inf_low_2_0_3,
            abs1 = abs_inf_high_2_0_3,
            abs2 = abs_inf_low_2_0_3,
            n1 = 5,
            n2 = 4,
            n_total = 5 + 4),
      df = 1,
      lower.tail = FALSE
    ),
    Infection_high_vs_Mock_none = pchisq(
      q = g(obs1 = obs_inf_high_2_0_3,
            obs2 = obs_mock_2_0_3,
            abs1 = abs_inf_high_2_0_3,
            abs2 = abs_mock_2_0_3,
            n1 = 5,
            n2 = 3,
            n_total = 5 + 3),
      df = 1,
      lower.tail = FALSE
    ),
    Infection_low_vs_Mock_none = pchisq(
      q = g(obs1 = obs_inf_low_2_0_3,
            obs2 = obs_mock_2_0_3,
            abs1 = abs_inf_low_2_0_3,
            abs2 = abs_mock_2_0_3,
            n1 = 4,
            n2 = 3,
            n_total = 4 + 3),
      df = 1,
      lower.tail = FALSE
    )
  )
  
  flag_g_2_0_3 <- data.frame(
    Infection_high_vs_Infection_low = gflag(
      obs1 = obs_inf_high_2_0_3,
      obs2 = obs_inf_low_2_0_3,
      abs1 = abs_inf_high_2_0_3,
      abs2 = abs_inf_low_2_0_3,
      pvals = pval_g_2_0_3[, 1],
      cutoff = 0.05
    ),
    Infection_high_vs_Mock_none = gflag(
      obs1 = obs_inf_high_2_0_3,
      obs2 = obs_mock_2_0_3,
      abs1 = abs_inf_high_2_0_3,
      abs2 = abs_mock_2_0_3,
      pvals = pval_g_2_0_3[, 2],
      cutoff = 0.05
    ),
    Infection_low_vs_Mock_none = gflag(
      obs1 = obs_inf_low_2_0_3,
      obs2 = obs_mock_2_0_3,
      abs1 = abs_inf_low_2_0_3,
      abs2 = abs_mock_2_0_3,
      pvals = pval_g_2_0_3[, 3],
      cutoff = 0.05
    )
  )
  
  mean_2_0_3 <- data.frame(
    Mean_Infection_high = rowMeans(gfilta_2_0_3$e_data[, c(2:4, 7:8)],
                                   na.rm = TRUE),
    Mean_Infection_low = rowMeans(gfilta_2_0_3$e_data[, c(5:6, 9:10)],
                                  na.rm = TRUE),
    Mean_Mock_none = rowMeans(gfilta_2_0_3$e_data[, 11:13],
                              na.rm = TRUE)
  )
  
  gstan_2_0_3 <- list(
    Full_results = data.frame(
      Mass_Tag_ID = gfilta_2_0_3$e_data$Mass_Tag_ID,
      Count_Infection_high = unname(obs_inf_high_2_0_3),
      Count_Infection_low = unname(obs_inf_low_2_0_3),
      Count_Mock_none = unname(obs_mock_2_0_3),
      P_value_G_Infection_high_vs_Infection_low = pval_g_2_0_3[, 1],
      P_value_G_Infection_high_vs_Mock_none = pval_g_2_0_3[, 2],
      P_value_G_Infection_low_vs_Mock_none = pval_g_2_0_3[, 3],
      Flag_Infection_high_vs_Infection_low = flag_g_2_0_3[, 1],
      Flag_Infection_high_vs_Mock_none = flag_g_2_0_3[, 2],
      Flag_Infection_low_vs_Mock_none = flag_g_2_0_3[, 3],
      mean_2_0_3,
      Fold_change_Infection_high_vs_Infection_low = (
        mean_2_0_3[, 1] - mean_2_0_3[, 2]
      ),
      Fold_change_Infection_high_vs_Mock_none = (
        mean_2_0_3[, 1] - mean_2_0_3[, 3]
      ),
      Fold_change_Infection_low_vs_Mock_none = (
        mean_2_0_3[, 2] - mean_2_0_3[, 3]
      ),
      row.names = NULL
    ),
    
    Flags = data.frame(
      Mass_Tag_ID = gfilta_2_0_3$e_data$Mass_Tag_ID,
      flag_g_2_0_3,
      row.names = NULL
    ),
    
    P_values = data.frame(
      Mass_Tag_ID = gfilta_2_0_3$e_data$Mass_Tag_ID,
      pval_g_2_0_3,
      row.names = NULL
    )
  )
  
  class(gstan_2_0_3) <- "statRes"
  
  attr(gstan_2_0_3, "group_DF") <- groupDF_2_0_3
  attr(gstan_2_0_3, "comparisons") <- c("Infection_high_vs_Infection_low",
                                        "Infection_high_vs_Mock_none",
                                        "Infection_low_vs_Mock_none")
  attr(gstan_2_0_3, "number_significant") <- data.frame(
    Comparison = c("Infection_high_vs_Infection_low",
                   "Infection_high_vs_Mock_none",
                   "Infection_low_vs_Mock_none"),
    Up_total = c(length(which(flag_g_2_0_3[, 1] == 2)),
                 length(which(flag_g_2_0_3[, 2] == 2)),
                 length(which(flag_g_2_0_3[, 3] == 2))),
    Down_total = c(length(which(flag_g_2_0_3[, 1] == -2)),
                   length(which(flag_g_2_0_3[, 2] == -2)),
                   length(which(flag_g_2_0_3[, 3] == -2))),
    Up_anova = c(0, 0, 0),
    Down_anova = c(0, 0, 0),
    Up_gtest = c(length(which(flag_g_2_0_3[, 1] == 2)),
                 length(which(flag_g_2_0_3[, 2] == 2)),
                 length(which(flag_g_2_0_3[, 3] == 2))),
    Down_gtest = c(length(which(flag_g_2_0_3[, 1] == -2)),
                   length(which(flag_g_2_0_3[, 2] == -2)),
                   length(which(flag_g_2_0_3[, 3] == -2))),
    row.names = c("Infection_high_vs_Infection_low",
                  "Infection_high_vs_Mock_none",
                  "Infection_low_vs_Mock_none")
  )
  attr(gstan_2_0_3, "statistical_test") <- "gtest"
  attr(gstan_2_0_3, "adjustment_method") <- "none"
  attr(gstan_2_0_3, "pval_thresh") <- 0.05
  attr(gstan_2_0_3, "data_info") <- list(
    data_scale_orig = "abundance",
    data_scale = "log",
    norm_info = list(is_normalized = FALSE),
    num_edata = length(unique(gfilta_2_0_3$e_data$Mass_Tag_ID)),
    num_miss_obs = sum(is.na(gfilta_2_0_3$e_data[, -1])),
    prop_missing = (sum(is.na(gfilta_2_0_3$e_data[, -1])) /
                      prod(dim(gfilta_2_0_3$e_data[, -1]))),
    num_samps = dim(gfilta_2_0_3$f_data)[1],
    data_types = NULL
  )
  attr(gstan_2_0_3, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = NULL,
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(gstan_2_0_3, "data_class") <- "pepData"
  
  # Compare IMD-ANOVAly --------------------------------------------------------
  
  # Default comparisons ---------------
  
  afruit_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova")
  gfruit_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest")
  cfruit_1_0_2 <- imd_anova(cfilta_1_0_2, test_method = "combined")
  
  afruit_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova")
  gfruit_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest")
  cfruit_2_0_3 <- imd_anova(cfilta_2_0_3, test_method = "combined")
  
  # Custom comparisons ---------------
  
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
  
  # Adjusted p-values ---------------
  
  gfruit_bon_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest",
                                pval_adjust = "bonferroni")
  
  gfruit_bon_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                pval_adjust = "bonferroni")
  gfruit_holm_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                 pval_adjust = "holm")
  
  # Holy IMD-ANOVA unit tests, Statman! ----------------------------------------
  
  expect_equal(gfruit_1_0_2, gstan_1_0_2)
  expect_equal(gfruit_2_0_3, gstan_2_0_3)
  
  # Make sure the 
  expect_equal(gfruit_2_0_3, gfruit_cus_2_0_3)
  
  # Make sure adjusted p-values are the same as unadjusted p-values when there
  # is only one test between two groups.
  expect_equal(gfruit_1_0_2$Full_results$P_value_G_Infection_vs_Mock,
               gfruit_bon_1_0_2$Full_results$P_value_G_Infection_vs_Mock)
  expect_equal(attr(gfruit_bon_1_0_2, "adjustment_method"),
               "bonferroni")
  
  # Have a looksie at the p-values when more than one test is applied per
  # biomolecule. Check both the Full_results and P_values data frames.
  expect_equal(
    gfruit_bon_2_0_3$Full_results[, 5:7],
    data.frame(pmin(data.matrix(gstan_2_0_3$Full_results[, 5:7] * 2), 1))
  )
  expect_equal(
    gfruit_bon_2_0_3$P_values[, 2:4],
    data.frame(pmin(data.matrix(gstan_2_0_3$P_values[, 2:4] * 2), 1))
  )
  expect_equal(
    gfruit_holm_2_0_3$Full_results[, 5:7],
    data.frame(t(apply(gstan_2_0_3$Full_results[, 5:7],
                       1,
                       p.adjust,
                       method = "holm")))
  )
  expect_equal(
    gfruit_holm_2_0_3$P_values[, 2:4],
    data.frame(t(apply(gstan_2_0_3$P_values[, 2:4],
                       1,
                       p.adjust,
                       method = "holm")))
  )
  
  # Test argument checks -------------------------------------------------------
  
  # Maybe add later (depending on the position of the planets and stars).
  
})
