# Load functions for calculating IMD-ANOVA standards ---------------------------

source (system.file("testdata",
                    "imd_anova_standard_fns.R",
                    package = "pmartR"))

# Load data and prepare omicsData objects --------------------------------------

load(system.file('testdata',
                 'little_pairdata.RData',
                 package = 'pmartR'))

# Create additional f_data objects with different main effects and covariates.
fdata2 <- fdata
set.seed(40)
# Because this is paired data we only need to simulate 15 values for Gender and
# repeat them according to the three sample value (Mock, FM, and AM).
fdata2$Gender <- c(rep(sample(c("F", "M"), 5, replace = TRUE), 2),
                   rep(sample(c("F", "M"), 5, replace = TRUE), 2),
                   rep(sample(c("F", "M"), 5, replace = TRUE), 2))

# Create pepData objects, log the data, add group_DF attributes ----------------

# Paired data, zero main effects, zero covariates, one group
pairdata_0_0_3 <- as.pepData(e_data = edata,
                             f_data = fdata,
                             e_meta = emeta,
                             edata_cname = 'Mass_Tag_ID',
                             fdata_cname = 'Name',
                             emeta_cname = 'Protein')
pairdata_0_0_3 <- edata_transform(pairdata_0_0_3,
                                  data_scale = "log")
pairdata_0_0_3 <- group_designation(pairdata_0_0_3,
                                    pair_id = "PairID",
                                    pair_group = "Time",
                                    pair_denom = "18")

# Paired data, one main effect, zero covariates, three groups
pairdata_1_0_3 <- as.pepData(e_data = edata,
                             f_data = fdata,
                             e_meta = emeta,
                             edata_cname = 'Mass_Tag_ID',
                             fdata_cname = 'Name',
                             emeta_cname = 'Protein')
pairdata_1_0_3 <- edata_transform(pairdata_1_0_3,
                                  data_scale = "log")
pairdata_1_0_3 <- group_designation(pairdata_1_0_3,
                                    main_effects = "Virus",
                                    pair_id = "PairID",
                                    pair_group = "Time",
                                    pair_denom = "18")

# Paired data, one main effect, one covariate, three groups
pairdata_1_1_3 <- as.pepData(e_data = edata,
                             f_data = fdata2,
                             e_meta = emeta,
                             edata_cname = 'Mass_Tag_ID',
                             fdata_cname = 'Name',
                             emeta_cname = 'Protein')
pairdata_1_1_3 <- edata_transform(pairdata_1_1_3,
                                  data_scale = "log")
pairdata_1_1_3 <- group_designation(pairdata_1_1_3,
                                    main_effects = "Virus",
                                    covariates = "Gender",
                                    pair_id = "PairID",
                                    pair_group = "Time",
                                    pair_denom = "18")

# Filter IMD-ANOVAly -----------------------------------------------------------

filta_0_0_3 <- imdanova_filter(pairdata_0_0_3)
afilta_0_0_3 <- applyFilt(filta_0_0_3, pairdata_0_0_3,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_0_0_3 <- applyFilt(filta_0_0_3, pairdata_0_0_3,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_0_0_3 <- applyFilt(filta_0_0_3, pairdata_0_0_3,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

filta_1_0_3 <- imdanova_filter(pairdata_1_0_3)
afilta_1_0_3 <- applyFilt(filta_1_0_3, pairdata_1_0_3,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_1_0_3 <- applyFilt(filta_1_0_3, pairdata_1_0_3,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_1_0_3 <- applyFilt(filta_1_0_3, pairdata_1_0_3,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

filta_1_1_3 <- imdanova_filter(pairdata_1_1_3)
afilta_1_1_3 <- applyFilt(filta_1_1_3, pairdata_1_1_3,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_1_1_3 <- applyFilt(filta_1_1_3, pairdata_1_1_3,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_1_1_3 <- applyFilt(filta_1_1_3, pairdata_1_1_3,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

# Save filter objects  ---------------------------------------------------------

# save(
#
#   afilta_0_0_3,
#   gfilta_0_0_3,
#   cfilta_0_0_3,
#
#   afilta_1_0_3,
#   gfilta_1_0_3,
#   cfilta_1_0_3,
#
#   afilta_1_1_3,
#   gfilta_1_1_3,
#   cfilta_1_1_3,
#
#   file = file.path("/Users/mart077/Documents/r_packages/pmartR",
#                    "inst/testdata/standards_filter_paired.RData")
# )

# Prepare paired objects -------------------------------------------------------

# Use 18hr as the control because this is how the original functions to take the
# difference were written (the second half of a pair was always subtracted from
# the first half). The order of a pair was determined by the order of the
# samples in fdata. For example, if the 0hr sample occurred in row one and the
# 18hr sample in row two the difference would be take by 0hr - 18hr. If the
# order was mixed up (18hr occurred before 0hr in fdata) then for that pair 0hr
# would be subtracted from 18hr.
diff_a_1_0_3 <- data.frame(
  Mass_Tag_ID = afilta_1_0_3$e_data$Mass_Tag_ID,
  Pair_1 = afilta_1_0_3$e_data$Mock_0hr_1 - afilta_1_0_3$e_data$Mock_18hr_1,
  Pair_2 = afilta_1_0_3$e_data$Mock_0hr_2 - afilta_1_0_3$e_data$Mock_18hr_2,
  Pair_3 = afilta_1_0_3$e_data$Mock_0hr_3 - afilta_1_0_3$e_data$Mock_18hr_3,
  Pair_4 = afilta_1_0_3$e_data$Mock_0hr_4 - afilta_1_0_3$e_data$Mock_18hr_4,
  Pair_5 = afilta_1_0_3$e_data$Mock_0hr_5 - afilta_1_0_3$e_data$Mock_18hr_5,
  Pair_6 = afilta_1_0_3$e_data$FM_0hr_1 - afilta_1_0_3$e_data$FM_18hr_1,
  Pair_7 = afilta_1_0_3$e_data$FM_0hr_2 - afilta_1_0_3$e_data$FM_18hr_2,
  Pair_8 = afilta_1_0_3$e_data$FM_0hr_3 - afilta_1_0_3$e_data$FM_18hr_3,
  Pair_9 = afilta_1_0_3$e_data$FM_0hr_4 - afilta_1_0_3$e_data$FM_18hr_4,
  Pair_10 = afilta_1_0_3$e_data$FM_0hr_5 - afilta_1_0_3$e_data$FM_18hr_5,
  Pair_11 = afilta_1_0_3$e_data$AM_0hr_1 - afilta_1_0_3$e_data$AM_18hr_1,
  Pair_12 = afilta_1_0_3$e_data$AM_0hr_2 - afilta_1_0_3$e_data$AM_18hr_2,
  Pair_13 = afilta_1_0_3$e_data$AM_0hr_3 - afilta_1_0_3$e_data$AM_18hr_3,
  Pair_14 = afilta_1_0_3$e_data$AM_0hr_4 - afilta_1_0_3$e_data$AM_18hr_4,
  Pair_15 = afilta_1_0_3$e_data$AM_0hr_5 - afilta_1_0_3$e_data$AM_18hr_5,
  row.names = NULL
)

diff_a_0_0_3 <- data.frame(
  Mass_Tag_ID = afilta_0_0_3$e_data$Mass_Tag_ID,
  Pair_1 = afilta_0_0_3$e_data$Mock_0hr_1 - afilta_0_0_3$e_data$Mock_18hr_1,
  Pair_2 = afilta_0_0_3$e_data$Mock_0hr_2 - afilta_0_0_3$e_data$Mock_18hr_2,
  Pair_3 = afilta_0_0_3$e_data$Mock_0hr_3 - afilta_0_0_3$e_data$Mock_18hr_3,
  Pair_4 = afilta_0_0_3$e_data$Mock_0hr_4 - afilta_0_0_3$e_data$Mock_18hr_4,
  Pair_5 = afilta_0_0_3$e_data$Mock_0hr_5 - afilta_0_0_3$e_data$Mock_18hr_5,
  Pair_6 = afilta_0_0_3$e_data$FM_0hr_1 - afilta_0_0_3$e_data$FM_18hr_1,
  Pair_7 = afilta_0_0_3$e_data$FM_0hr_2 - afilta_0_0_3$e_data$FM_18hr_2,
  Pair_8 = afilta_0_0_3$e_data$FM_0hr_3 - afilta_0_0_3$e_data$FM_18hr_3,
  Pair_9 = afilta_0_0_3$e_data$FM_0hr_4 - afilta_0_0_3$e_data$FM_18hr_4,
  Pair_10 = afilta_0_0_3$e_data$FM_0hr_5 - afilta_0_0_3$e_data$FM_18hr_5,
  Pair_11 = afilta_0_0_3$e_data$AM_0hr_1 - afilta_0_0_3$e_data$AM_18hr_1,
  Pair_12 = afilta_0_0_3$e_data$AM_0hr_2 - afilta_0_0_3$e_data$AM_18hr_2,
  Pair_13 = afilta_0_0_3$e_data$AM_0hr_3 - afilta_0_0_3$e_data$AM_18hr_3,
  Pair_14 = afilta_0_0_3$e_data$AM_0hr_4 - afilta_0_0_3$e_data$AM_18hr_4,
  Pair_15 = afilta_0_0_3$e_data$AM_0hr_5 - afilta_0_0_3$e_data$AM_18hr_5,
  row.names = NULL
)

diff_g <- data.frame(
  Mass_Tag_ID = gfilta_1_0_3$e_data$Mass_Tag_ID,
  Pair_1 = gfilta_1_0_3$e_data$Mock_0hr_1 - gfilta_1_0_3$e_data$Mock_18hr_1,
  Pair_2 = gfilta_1_0_3$e_data$Mock_0hr_2 - gfilta_1_0_3$e_data$Mock_18hr_2,
  Pair_3 = gfilta_1_0_3$e_data$Mock_0hr_3 - gfilta_1_0_3$e_data$Mock_18hr_3,
  Pair_4 = gfilta_1_0_3$e_data$Mock_0hr_4 - gfilta_1_0_3$e_data$Mock_18hr_4,
  Pair_5 = gfilta_1_0_3$e_data$Mock_0hr_5 - gfilta_1_0_3$e_data$Mock_18hr_5,
  Pair_6 = gfilta_1_0_3$e_data$FM_0hr_1 - gfilta_1_0_3$e_data$FM_18hr_1,
  Pair_7 = gfilta_1_0_3$e_data$FM_0hr_2 - gfilta_1_0_3$e_data$FM_18hr_2,
  Pair_8 = gfilta_1_0_3$e_data$FM_0hr_3 - gfilta_1_0_3$e_data$FM_18hr_3,
  Pair_9 = gfilta_1_0_3$e_data$FM_0hr_4 - gfilta_1_0_3$e_data$FM_18hr_4,
  Pair_10 = gfilta_1_0_3$e_data$FM_0hr_5 - gfilta_1_0_3$e_data$FM_18hr_5,
  Pair_11 = gfilta_1_0_3$e_data$AM_0hr_1 - gfilta_1_0_3$e_data$AM_18hr_1,
  Pair_12 = gfilta_1_0_3$e_data$AM_0hr_2 - gfilta_1_0_3$e_data$AM_18hr_2,
  Pair_13 = gfilta_1_0_3$e_data$AM_0hr_3 - gfilta_1_0_3$e_data$AM_18hr_3,
  Pair_14 = gfilta_1_0_3$e_data$AM_0hr_4 - gfilta_1_0_3$e_data$AM_18hr_4,
  Pair_15 = gfilta_1_0_3$e_data$AM_0hr_5 - gfilta_1_0_3$e_data$AM_18hr_5,
  row.names = NULL
)

groupie <- data.frame(
  Name = paste("Pair", 1:15, sep = "_"),
  Group = c(rep("Mock", 5), rep("FM", 5), rep("AM", 5))
)

Xmatrix_1_1_3 = structure(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1,
                            0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                            0, 1, 0, 0, 1, 1, 0),
                          .Dim = c(15L, 5L))

# Assemble ANOVA standards -----------------------------------------------------

# main effects: 0; covariates: 0; groups: 3 ---------------

pval_a_0_0_3 <- data.frame(
  P_value_A_paired_diff = diff_a_0_0_3[, -1] %>%
    apply( 1, t.test) %>%
    lapply(`[[`, "p.value") %>%
    unlist()
)

flag_a_0_0_3 <- data.frame(
  Flag_A_paired_diff = aflag_diff(
    diff = rowMeans(diff_a_0_0_3[, 2:16],
                    na.rm = TRUE),
    pvals = pval_a_0_0_3[, 1],
    cutoff = 0.05
  )
)

astan_0_0_3 <- data.frame(
  Mass_Tag_ID = diff_a_0_0_3$Mass_Tag_ID,
  Count_paired_diff = unname(rowSums(!is.na(afilta_0_0_3$e_data[, 2:31]))),
  Mean_paired_diff = rowMeans(diff_a_0_0_3[, 2:16],
                       na.rm = TRUE),
  Fold_change_paired_diff = rowMeans(diff_a_0_0_3[, 2:16],
                              na.rm = TRUE), # Usually mean(grp1) - mean(grp2)
  pval_a_0_0_3,
  flag_a_0_0_3,
  row.names = NULL
)

class(astan_0_0_3) <- c("statRes", "data.frame")

attr(astan_0_0_3, "group_DF") <- attr(afilta_0_0_3, "group_DF")
attr(astan_0_0_3, "comparisons") <- c("paired_diff")
attr(astan_0_0_3, "number_significant") <- data.frame(
  Comparison = c("paired_diff"),
  Up_total = c(length(which(flag_a_0_0_3[, 1] == 1))),
  Down_total = c(length(which(flag_a_0_0_3[, 1] == -1))),
  Up_anova = c(length(which(flag_a_0_0_3[, 1] == 1))),
  Down_anova = c(length(which(flag_a_0_0_3[, 1] == -1))),
  Up_gtest = c(0),
  Down_gtest = c(0),
  row.names = NULL
)
attr(astan_0_0_3, "statistical_test") <- "anova"
attr(astan_0_0_3, "adjustment_method_a_multcomp") <- "none"
attr(astan_0_0_3, "adjustment_method_g_multcomp") <- "none"
attr(astan_0_0_3, "adjustment_method_a_fdr") <- "none"
attr(astan_0_0_3, "adjustment_method_g_fdr") <- "none"
attr(astan_0_0_3, "pval_thresh") <- 0.05
attr(astan_0_0_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(diff_a_0_0_3$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_0_0_3$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_0_0_3$e_data[, -1])) /
                    prod(dim(afilta_0_0_3$e_data[, -1]))),
  num_samps = dim(afilta_0_0_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(astan_0_0_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_0_0_3$e_data$Mass_Tag_ID,
  paired_diff = flag_a_0_0_3$Flag_A_paired_diff
)
attr(astan_0_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(astan_0_0_3, "data_class") <- "pepData"

# main effects: 1; covariates: 0; groups: 3 ---------------

mean_a_1_0_3 <- data.frame(
  Mean_Mock = rowMeans(diff_a_1_0_3[, 2:6],
                       na.rm = TRUE),
  Mean_FM = rowMeans(diff_a_1_0_3[, 7:11],
                     na.rm = TRUE),
  Mean_AM = rowMeans(diff_a_1_0_3[, 12:16],
                     na.rm = TRUE)
)

group_counts_1_0_3 <- data.frame(
  nona_Mock = rowSums(!is.na(diff_a_1_0_3[, 2:6])),
  nona_FM = rowSums(!is.na(diff_a_1_0_3[, 7:11])),
  nona_AM = rowSums(!is.na(diff_a_1_0_3[, 12:16]))
)

nona_grps_1_0_3 <- rowSums(group_counts_1_0_3 != 0)

sigma_1_0_3 <- dplyr::rowwise(diff_a_1_0_3[, -1]) %>%
  dplyr::mutate(
    stdev = tryCatch (
      summary(lm(
        dplyr::c_across(Pair_1:Pair_15) ~ groupie$Group
      ))$sigma,
      error = function (e) {NA}
    )
  ) %>%
  dplyr::pull(stdev)

nona_counts_1_0_3 <- rowSums(!is.na(diff_a_1_0_3[, -1]))

diffs_1_0_3 <- mean_a_1_0_3 %>%
  dplyr::mutate(
    diff_Mock_FM = Mean_Mock - Mean_FM,
    diff_Mock_AM = Mean_Mock - Mean_AM,
    diff_FM_AM = Mean_FM - Mean_AM
  ) %>%
  dplyr::select(diff_Mock_FM, diff_Mock_AM, diff_FM_AM)

test_stat_1_0_3 <- diffs_1_0_3 %>%
  dplyr::mutate(
    stat_Mock_FM = (diff_Mock_FM /
                      (sigma_1_0_3 * sqrt(1/group_counts_1_0_3[, 1] +
                                            1/group_counts_1_0_3[, 2]))),
    stat_Mock_AM = (diff_Mock_AM /
                      (sigma_1_0_3 * sqrt(1/group_counts_1_0_3[, 1] +
                                            1/group_counts_1_0_3[, 3]))),
    stat_FM_AM = (diff_FM_AM /
                    (sigma_1_0_3 * sqrt(1/group_counts_1_0_3[, 2] +
                                          1/group_counts_1_0_3[, 3])))
  ) %>%
  dplyr::select(stat_Mock_FM, stat_Mock_AM, stat_FM_AM)

pval_a_1_0_3 <- test_stat_1_0_3 %>%
  dplyr::mutate(
    P_value_A_Mock_vs_FM = pt(
      q = abs(stat_Mock_FM),
      df = nona_counts_1_0_3 - nona_grps_1_0_3,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_vs_AM = pt(
      q = abs(stat_Mock_AM),
      df = nona_counts_1_0_3 - nona_grps_1_0_3,
      lower.tail = FALSE
    ) * 2,
    P_value_A_FM_vs_AM = pt(
      q = abs(stat_FM_AM),
      df = nona_counts_1_0_3 - nona_grps_1_0_3,
      lower.tail = FALSE
    ) * 2
  ) %>%
  dplyr::select(P_value_A_Mock_vs_FM,
                P_value_A_Mock_vs_AM,
                P_value_A_FM_vs_AM)

flag_a_1_0_3 <- data.frame(
  Mock_vs_FM = aflag(
    grp1 = mean_a_1_0_3$Mean_Mock,
    grp2 = mean_a_1_0_3$Mean_FM,
    pvals = pval_a_1_0_3[, 1],
    cutoff = 0.05
  ),
  Mock_vs_AM = aflag(
    grp1 = mean_a_1_0_3$Mean_Mock,
    grp2 = mean_a_1_0_3$Mean_AM,
    pvals = pval_a_1_0_3[, 2],
    cutoff = 0.05
  ),
  FM_vs_AM = aflag(
    grp1 = mean_a_1_0_3$Mean_FM,
    grp2 = mean_a_1_0_3$Mean_AM,
    pvals = pval_a_1_0_3[, 3],
    cutoff = 0.05
  )
)

astan_1_0_3 <- data.frame(
  Mass_Tag_ID = diff_a_1_0_3$Mass_Tag_ID,
  Count_Mock = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 2:11]))),
  Count_FM = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 12:21]))),
  Count_AM = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 22:31]))),
  mean_a_1_0_3,
  Fold_change_Mock_vs_FM = (
    mean_a_1_0_3[, 1] - mean_a_1_0_3[, 2]
  ),
  Fold_change_Mock_vs_AM = (
    mean_a_1_0_3[, 1] - mean_a_1_0_3[, 3]
  ),
  Fold_change_FM_vs_AM = (
    mean_a_1_0_3[, 2] - mean_a_1_0_3[, 3]
  ),
  pval_a_1_0_3,
  Flag_A_Mock_vs_FM = flag_a_1_0_3[, 1],
  Flag_A_Mock_vs_AM = flag_a_1_0_3[, 2],
  Flag_A_FM_vs_AM = flag_a_1_0_3[, 3],
  row.names = NULL
)

class(astan_1_0_3) <- c("statRes", "data.frame")

attr(astan_1_0_3, "group_DF") <- attr(afilta_1_0_3, "group_DF")
attr(astan_1_0_3, "comparisons") <- c("Mock_vs_FM",
                                      "Mock_vs_AM",
                                      "FM_vs_AM")
attr(astan_1_0_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM",
                 "Mock_vs_AM",
                 "FM_vs_AM"),
  Up_total = c(length(which(flag_a_1_0_3[, 1] == 1)),
               length(which(flag_a_1_0_3[, 2] == 1)),
               length(which(flag_a_1_0_3[, 3] == 1))),
  Down_total = c(length(which(flag_a_1_0_3[, 1] == -1)),
                 length(which(flag_a_1_0_3[, 2] == -1)),
                 length(which(flag_a_1_0_3[, 3] == -1))),
  Up_anova = c(length(which(flag_a_1_0_3[, 1] == 1)),
               length(which(flag_a_1_0_3[, 2] == 1)),
               length(which(flag_a_1_0_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_0_3[, 1] == -1)),
                 length(which(flag_a_1_0_3[, 2] == -1)),
                 length(which(flag_a_1_0_3[, 3] == -1))),
  Up_gtest = c(0, 0, 0),
  Down_gtest = c(0, 0, 0),
  row.names = NULL
)
attr(astan_1_0_3, "statistical_test") <- "anova"
attr(astan_1_0_3, "adjustment_method_a_multcomp") <- "none"
attr(astan_1_0_3, "adjustment_method_g_multcomp") <- "none"
attr(astan_1_0_3, "adjustment_method_a_fdr") <- "none"
attr(astan_1_0_3, "adjustment_method_g_fdr") <- "none"
attr(astan_1_0_3, "pval_thresh") <- 0.05
attr(astan_1_0_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(diff_a_1_0_3$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_1_0_3$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_1_0_3$e_data[, -1])) /
                    prod(dim(afilta_1_0_3$e_data[, -1]))),
  num_samps = dim(afilta_1_0_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(astan_1_0_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_1_0_3$e_data$Mass_Tag_ID,
  flag_a_1_0_3
)
attr(astan_1_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(astan_1_0_3, "data_class") <- "pepData"

# main effects: 1; covariates: 1; groups: 3 ---------------

Xmatrix = Xmatrix_1_1_3

gp <- factor(c(rep(1,5), rep(2,5), rep(3,5)),labels=1:3,levels=c(1,2,3))
Betas = compute_betas(data_mat = data.matrix(diff_a_1_0_3[, -1]), Xmatrix = Xmatrix, gp = gp)

covariate_effects = Xmatrix[,4:5] %*% t(Betas[,4:5])

# Adjust the means to remove effect of covariates.
adj_data_1_1_3 <- data.matrix(diff_a_1_0_3[, -1]) - t(covariate_effects)
adj_data_1_1_3 <- as.data.frame(adj_data_1_1_3)

mean_a_1_1_3 <- data.frame(
  Mean_Mock = rowMeans(adj_data_1_1_3[, 1:5],
                       na.rm = TRUE),
  Mean_FM = rowMeans(adj_data_1_1_3[, 6:10],
                     na.rm = TRUE),
  Mean_AM = rowMeans(adj_data_1_1_3[, 11:15],
                     na.rm = TRUE)
)

group_counts_1_1_3 <- data.frame(
  nona_Mock = rowSums(!is.na(adj_data_1_1_3[, 1:5])),
  nona_FM = rowSums(!is.na(adj_data_1_1_3[, 6:10])),
  nona_AM = rowSums(!is.na(adj_data_1_1_3[, 11:15]))
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n_grp = sum(dplyr::c_across(nona_Mock:nona_AM) == 0)) %>%
  dplyr::ungroup()

nona_grps_1_1_3 <- rowSums(group_counts_1_1_3[, 1:3] != 0)

mymiss <- diff_a_1_0_3[, -1] %>% apply(1, function(x) !is.na(x))
ranks <- sapply(1:ncol(mymiss), function(i) {
  col = mymiss[,i]
  rank = Matrix::rankMatrix(Xmatrix[col,])
})

XTXs <- lapply(1:ncol(mymiss), function(i) {
  col = mymiss[,i]
  rank = MASS::ginv(t(Xmatrix[col,]) %*% Xmatrix[col,])
})

SEs <- lapply(1:length(XTXs), function(i) {
  nona_idx = mymiss[,i]
  y = unlist(diff_a_1_0_3[i,-1][nona_idx])
  H = Xmatrix[nona_idx,] %*% XTXs[[i]] %*% t(Xmatrix[nona_idx,])
  resids = y - (H %*% y)
  SSE = sum(resids^2)
  denom = sum(nona_idx) - ranks[i] + group_counts_1_1_3$n_grp[1]
  return(sqrt(SSE/denom))
})


nona_counts_1_1_3 <- rowSums(!is.na(adj_data_1_1_3))

diffs_1_1_3 <- mean_a_1_1_3 %>%
  dplyr::mutate(
    diff_m_f = Mean_Mock - Mean_FM,
    diff_m_a = Mean_Mock - Mean_AM,
    diff_f_a = Mean_FM - Mean_AM
  ) %>%
  dplyr::select(diff_m_f, diff_m_a, diff_f_a)

cmat = rbind(c(1, -1, 0, 0, 0), c(1, 0, -1, 0, 0), c(0, 1, -1, 0, 0))
diff_denoms <- lapply(1:length(XTXs), function(i) {
  sqrt(diag(cmat %*% XTXs[[i]] %*% t(cmat))) * SEs[[i]]
})

diff_denoms <- do.call(rbind, diff_denoms) %>% `colnames<-`(c("C1", "C2", "C3"))

test_stat_1_1_3 <- diffs_1_1_3 %>%
  cbind(group_counts_1_1_3) %>%
  cbind(diff_denoms) %>%
  dplyr::mutate(
    stat_m_f = (diff_m_f / C1),
    stat_m_a = (diff_m_a / C2),
    stat_f_a = (diff_f_a / C3)
  ) %>%
  dplyr::select(stat_m_f, stat_m_a, stat_f_a) %>%
  dplyr::ungroup() %>%
  `row.names<-`(NULL)

pval_a_1_1_3 <- test_stat_1_1_3 %>%
  dplyr::mutate(
    lg = group_counts_1_1_3$n_grp,
    ranks = ranks,
    P_value_A_Mock_vs_FM = pt(
      q = abs(stat_m_f),
      df = nona_counts_1_1_3 - ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_vs_AM = pt(
      q = abs(stat_m_a),
      df = nona_counts_1_1_3 - ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_FM_vs_AM = pt(
      q = abs(stat_f_a),
      df = nona_counts_1_1_3 - ranks,
      lower.tail = FALSE
    ) * 2
  ) %>%
  dplyr::select(P_value_A_Mock_vs_FM,
                P_value_A_Mock_vs_AM,
                P_value_A_FM_vs_AM) %>%
  `row.names<-`(NULL)

flag_a_1_1_3 <- data.frame(
  Mock_vs_FM = aflag(
    grp1 = mean_a_1_1_3$Mean_Mock,
    grp2 = mean_a_1_1_3$Mean_FM,
    pvals = pval_a_1_1_3[, 1],
    cutoff = 0.05
  ),
  Mock_vs_AM = aflag(
    grp1 = mean_a_1_1_3$Mean_Mock,
    grp2 = mean_a_1_1_3$Mean_AM,
    pvals = pval_a_1_1_3[, 2],
    cutoff = 0.05
  ),
  FM_vs_AM = aflag(
    grp1 = mean_a_1_1_3$Mean_FM,
    grp2 = mean_a_1_1_3$Mean_AM,
    pvals = pval_a_1_1_3[, 3],
    cutoff = 0.05
  )
)

astan_1_1_3 <- data.frame(
  Mass_Tag_ID = afilta_1_1_3$e_data$Mass_Tag_ID,
  Count_Mock = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 2:11]))),
  Count_FM = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 12:21]))),
  Count_AM = unname(rowSums(!is.na(afilta_1_0_3$e_data[, 22:31]))),
  mean_a_1_1_3,
  Fold_change_Mock_vs_FM = (
    mean_a_1_1_3[, 1] - mean_a_1_1_3[, 2]
  ),
  Fold_change_Mock_vs_AM = (
    mean_a_1_1_3[, 1] - mean_a_1_1_3[, 3]
  ),
  Fold_change_FM_vs_AM = (
    mean_a_1_1_3[, 2] - mean_a_1_1_3[, 3]
  ),
  pval_a_1_1_3,
  Flag_A_Mock_vs_FM = flag_a_1_1_3[, 1],
  Flag_A_Mock_vs_AM = flag_a_1_1_3[, 2],
  Flag_A_FM_vs_AM = flag_a_1_1_3[, 3],
  row.names = NULL
)

class(astan_1_1_3) <- c("statRes", "data.frame")

attr(astan_1_1_3, "group_DF") <- attr(afilta_1_1_3, "group_DF")
attr(astan_1_1_3, "comparisons") <- c("Mock_vs_FM",
                                      "Mock_vs_AM",
                                      "FM_vs_AM")
attr(astan_1_1_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM",
                 "Mock_vs_AM",
                 "FM_vs_AM"),
  Up_total = c(length(which(flag_a_1_1_3[, 1] == 1)),
               length(which(flag_a_1_1_3[, 2] == 1)),
               length(which(flag_a_1_1_3[, 3] == 1))),
  Down_total = c(length(which(flag_a_1_1_3[, 1] == -1)),
                 length(which(flag_a_1_1_3[, 2] == -1)),
                 length(which(flag_a_1_1_3[, 3] == -1))),
  Up_anova = c(length(which(flag_a_1_1_3[, 1] == 1)),
               length(which(flag_a_1_1_3[, 2] == 1)),
               length(which(flag_a_1_1_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_1_3[, 1] == -1)),
                 length(which(flag_a_1_1_3[, 2] == -1)),
                 length(which(flag_a_1_1_3[, 3] == -1))),
  Up_gtest = c(0, 0, 0),
  Down_gtest = c(0, 0, 0),
  row.names = NULL
)
attr(astan_1_1_3, "statistical_test") <- "anova"
attr(astan_1_1_3, "adjustment_method_a_multcomp") <- "none"
attr(astan_1_1_3, "adjustment_method_g_multcomp") <- "none"
attr(astan_1_1_3, "adjustment_method_a_fdr") <- "none"
attr(astan_1_1_3, "adjustment_method_g_fdr") <- "none"
attr(astan_1_1_3, "pval_thresh") <- 0.05
attr(astan_1_1_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(afilta_1_1_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_1_1_3$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_1_1_3$e_data[, -1])) /
                    prod(dim(afilta_1_1_3$e_data[, -1]))),
  num_samps = dim(afilta_1_1_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(astan_1_1_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_1_1_3$e_data$Mass_Tag_ID,
  flag_a_1_1_3
)
attr(astan_1_1_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(astan_1_1_3, "data_class") <- "pepData"

# Generate G-Test standards ----------------------------------------------------

# Construct count matrices: Objects that start with obs represent counts of
# non-missing (observed) values and objects that start with abs represent
# counts of missing (absent) values.

obs_mock <- rowSums(!is.na(gfilta_1_0_3$e_data[, 2:11]))
obs_fm <- rowSums(!is.na(gfilta_1_0_3$e_data[, 12:21]))
obs_am <- rowSums(!is.na(gfilta_1_0_3$e_data[, 22:31]))
abs_mock <- rowSums(is.na(gfilta_1_0_3$e_data[, 2:11]))
abs_fm <- rowSums(is.na(gfilta_1_0_3$e_data[, 12:21]))
abs_am <- rowSums(is.na(gfilta_1_0_3$e_data[, 22:31]))

# main effects: 1; covariates: 0; groups: 3 ---------------

pval_g_1_0_3 <- data.frame(
  P_value_G_Mock_vs_FM = rep(0, nrow(gfilta_1_0_3$e_data)),
  P_value_G_Mock_vs_AM = rep(0, nrow(gfilta_1_0_3$e_data)),
  P_value_G_FM_vs_AM = rep(0, nrow(gfilta_1_0_3$e_data))
)

for (e in 1:nrow(gfilta_1_0_3$e_data)) {

  pval_g_1_0_3[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_0_3$e_data[e, 2:21])) ~
        gfilta_1_0_3$f_data$PairID[1:20] +
        attr(gfilta_1_0_3, "group_DF")$Group[1:20],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  pval_g_1_0_3[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_0_3$e_data[e, c(2:11, 22:31)])) ~
        gfilta_1_0_3$f_data$PairID[c(1:10, 21:30)] +
        attr(gfilta_1_0_3, "group_DF")$Group[c(1:10, 21:30)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  pval_g_1_0_3[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_0_3$e_data[e, 12:31])) ~
        gfilta_1_0_3$f_data$PairID[11:30] +
        attr(gfilta_1_0_3, "group_DF")$Group[11:30],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]

}

flag_g_1_0_3 <- data.frame(
  Mock_vs_FM = gflag(
    obs1 = obs_mock,
    obs2 = obs_fm,
    abs1 = abs_mock,
    abs2 = abs_fm,
    pvals = pval_g_1_0_3[, 1],
    cutoff = 0.05
  ),
  Mock_vs_AM = gflag(
    obs1 = obs_mock,
    obs2 = obs_am,
    abs1 = abs_mock,
    abs2 = abs_am,
    pvals = pval_g_1_0_3[, 2],
    cutoff = 0.05
  ),
  FM_vs_AM = gflag(
    obs1 = obs_fm,
    obs2 = obs_am,
    abs1 = abs_fm,
    abs2 = abs_am,
    pvals = pval_g_1_0_3[, 3],
    cutoff = 0.05
  )
)

mean_g_1_0_3 <- data.frame(
  Mean_Mock = rowMeans(diff_g[, 2:6],
                       na.rm = TRUE),
  Mean_FM = rowMeans(diff_g[, 7:11],
                     na.rm = TRUE),
  Mean_AM = rowMeans(diff_g[, 12:16],
                     na.rm = TRUE)
)

gstan_1_0_3 <- data.frame(
  Mass_Tag_ID = gfilta_1_0_3$e_data$Mass_Tag_ID,
  Count_Mock = unname(obs_mock),
  Count_FM = unname(obs_fm),
  Count_AM = unname(obs_am),
  mean_g_1_0_3,
  Fold_change_Mock_vs_FM = (
    mean_g_1_0_3[, 1] - mean_g_1_0_3[, 2]
  ),
  Fold_change_Mock_vs_AM = (
    mean_g_1_0_3[, 1] - mean_g_1_0_3[, 3]
  ),
  Fold_change_FM_vs_AM = (
    mean_g_1_0_3[, 2] - mean_g_1_0_3[, 3]
  ),
  P_value_G_Mock_vs_FM = pval_g_1_0_3[, 1],
  P_value_G_Mock_vs_AM = pval_g_1_0_3[, 2],
  P_value_G_FM_vs_AM = pval_g_1_0_3[, 3],
  Flag_G_Mock_vs_FM = flag_g_1_0_3[, 1],
  Flag_G_Mock_vs_AM = flag_g_1_0_3[, 2],
  Flag_G_FM_vs_AM = flag_g_1_0_3[, 3],
  row.names = NULL
)

class(gstan_1_0_3) <- c("statRes", "data.frame")

attr(gstan_1_0_3, "group_DF") <- attr(gfilta_1_0_3, "group_DF")
attr(gstan_1_0_3, "comparisons") <- c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM")
attr(gstan_1_0_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM"),
  Up_total = c(length(which(flag_g_1_0_3[, 1] == 1)),
               length(which(flag_g_1_0_3[, 2] == 1)),
               length(which(flag_g_1_0_3[, 3] == 1))),
  Down_total = c(length(which(flag_g_1_0_3[, 1] == -1)),
                 length(which(flag_g_1_0_3[, 2] == -1)),
                 length(which(flag_g_1_0_3[, 3] == -1))),
  Up_anova = c(0, 0, 0),
  Down_anova = c(0, 0, 0),
  Up_gtest = c(length(which(flag_g_1_0_3[, 1] == 1)),
               length(which(flag_g_1_0_3[, 2] == 1)),
               length(which(flag_g_1_0_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_0_3[, 1] == -1)),
                 length(which(flag_g_1_0_3[, 2] == -1)),
                 length(which(flag_g_1_0_3[, 3] == -1))),
  row.names = NULL
)
attr(gstan_1_0_3, "statistical_test") <- "gtest"
attr(gstan_1_0_3, "adjustment_method_a_multcomp") <- "none"
attr(gstan_1_0_3, "adjustment_method_g_multcomp") <- "none"
attr(gstan_1_0_3, "adjustment_method_a_fdr") <- "none"
attr(gstan_1_0_3, "adjustment_method_g_fdr") <- "none"
attr(gstan_1_0_3, "pval_thresh") <- 0.05
attr(gstan_1_0_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_0_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_0_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_0_3$e_data[, -1])) /
                    prod(dim(gfilta_1_0_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_0_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(gstan_1_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(gstan_1_0_3, "data_class") <- "pepData"

# main effects: 1; covariates: 1; groups: 3 ---------------

pval_g_1_1_3 <- data.frame(
  P_value_G_Mock_vs_FM = rep(0, nrow(gfilta_1_1_3$e_data)),
  P_value_G_Mock_vs_AM = rep(0, nrow(gfilta_1_1_3$e_data)),
  P_value_G_FM_vs_AM = rep(0, nrow(gfilta_1_1_3$e_data))
)

for (e in 1:nrow(gfilta_1_1_3$e_data)) {

  pval_g_1_1_3[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, 2:21])) ~
        gfilta_1_1_3$f_data$PairID[1:20] +
        attr(attr(gfilta_1_1_3, "group_DF"), "covariates")$Gender[1:20] +
        attr(gfilta_1_1_3, "group_DF")$Group[1:20],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  pval_g_1_1_3[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, c(2:11, 22:31)])) ~
        gfilta_1_1_3$f_data$PairID[c(1:10, 21:30)] +
        attr(attr(gfilta_1_1_3, "group_DF"),
             "covariates")$Gender[c(1:10, 21:30)] +
        attr(gfilta_1_1_3, "group_DF")$Group[c(1:10, 21:30)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  pval_g_1_1_3[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, 12:31])) ~
        gfilta_1_1_3$f_data$PairID[11:30] +
        attr(attr(gfilta_1_1_3, "group_DF"), "covariates")$Gender[11:30] +
        attr(gfilta_1_1_3, "group_DF")$Group[11:30],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]

}

flag_g_1_1_3 <- data.frame(
  Mock_vs_FM = gflag(
    obs1 = obs_mock,
    obs2 = obs_fm,
    abs1 = abs_mock,
    abs2 = abs_fm,
    pvals = pval_g_1_1_3[, 1],
    cutoff = 0.05
  ),
  Mock_vs_AM = gflag(
    obs1 = obs_mock,
    obs2 = obs_am,
    abs1 = abs_mock,
    abs2 = abs_am,
    pvals = pval_g_1_1_3[, 2],
    cutoff = 0.05
  ),
  FM_vs_AM = gflag(
    obs1 = obs_fm,
    obs2 = obs_am,
    abs1 = abs_fm,
    abs2 = abs_am,
    pvals = pval_g_1_1_3[, 3],
    cutoff = 0.05
  )
)

Xmatrix = Xmatrix_1_1_3
gp <- factor(c(rep(1,5), rep(2,5), rep(3,5)),labels=1:3,levels=c(1,2,3))
Betas = compute_betas(data_mat = data.matrix(diff_g[, -1]), Xmatrix = Xmatrix, gp = gp)

covariate_effects = Xmatrix[,4:5] %*% t(Betas[,4:5])

# Adjust the means to remove effect of covariates.
adj_data_g_1_1_3 <- data.matrix(diff_g[, -1]) - t(covariate_effects)
adj_data_g_1_1_3 <- as.data.frame(adj_data_g_1_1_3)

mean_g_1_1_3 <- data.frame(
  Mean_Mock = rowMeans(adj_data_g_1_1_3[, 1:5],
                       na.rm = TRUE),
  Mean_FM = rowMeans(adj_data_g_1_1_3[, 6:10],
                     na.rm = TRUE),
  Mean_AM = rowMeans(adj_data_g_1_1_3[, 11:15],
                     na.rm = TRUE)
)

gstan_1_1_3 <- data.frame(
  Mass_Tag_ID = gfilta_1_1_3$e_data$Mass_Tag_ID,
  Count_Mock = unname(obs_mock),
  Count_FM = unname(obs_fm),
  Count_AM = unname(obs_am),
  mean_g_1_1_3,
  Fold_change_Mock_vs_FM = (
    mean_g_1_1_3[, 1] - mean_g_1_1_3[, 2]
  ),
  Fold_change_Mock_vs_AM = (
    mean_g_1_1_3[, 1] - mean_g_1_1_3[, 3]
  ),
  Fold_change_FM_vs_AM = (
    mean_g_1_1_3[, 2] - mean_g_1_1_3[, 3]
  ),
  P_value_G_Mock_vs_FM = pval_g_1_1_3[, 1],
  P_value_G_Mock_vs_AM = pval_g_1_1_3[, 2],
  P_value_G_FM_vs_AM = pval_g_1_1_3[, 3],
  Flag_G_Mock_vs_FM = flag_g_1_1_3[, 1],
  Flag_G_Mock_vs_AM = flag_g_1_1_3[, 2],
  Flag_G_FM_vs_AM = flag_g_1_1_3[, 3],
  row.names = NULL
)

class(gstan_1_1_3) <- c("statRes", "data.frame")

attr(gstan_1_1_3, "group_DF") <- attr(gfilta_1_1_3, "group_DF")
attr(gstan_1_1_3, "comparisons") <- c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM")
attr(gstan_1_1_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM"),
  Up_total = c(length(which(flag_g_1_1_3[, 1] == 1)),
               length(which(flag_g_1_1_3[, 2] == 1)),
               length(which(flag_g_1_1_3[, 3] == 1))),
  Down_total = c(length(which(flag_g_1_1_3[, 1] == -1)),
                 length(which(flag_g_1_1_3[, 2] == -1)),
                 length(which(flag_g_1_1_3[, 3] == -1))),
  Up_anova = c(0, 0, 0),
  Down_anova = c(0, 0, 0),
  Up_gtest = c(length(which(flag_g_1_1_3[, 1] == 1)),
               length(which(flag_g_1_1_3[, 2] == 1)),
               length(which(flag_g_1_1_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_1_3[, 1] == -1)),
                 length(which(flag_g_1_1_3[, 2] == -1)),
                 length(which(flag_g_1_1_3[, 3] == -1))),
  row.names = NULL
)
attr(gstan_1_1_3, "statistical_test") <- "gtest"
attr(gstan_1_1_3, "adjustment_method_a_multcomp") <- "none"
attr(gstan_1_1_3, "adjustment_method_g_multcomp") <- "none"
attr(gstan_1_1_3, "adjustment_method_a_fdr") <- "none"
attr(gstan_1_1_3, "adjustment_method_g_fdr") <- "none"
attr(gstan_1_1_3, "pval_thresh") <- 0.05
attr(gstan_1_1_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_1_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_1_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_1_3$e_data[, -1])) /
                    prod(dim(gfilta_1_1_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_1_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(gstan_1_1_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(gstan_1_1_3, "data_class") <- "pepData"

# Create combined standards ----------------------------------------------------

# main effects: 1; covariates: 0; groups: 3 ---------------

# Combine the G-test and ANOVA results and place columns in correct order.
cstan_1_0_3 <- dplyr::full_join(gstan_1_0_3[, c(1:4, 11:16)],
                                astan_1_0_3) %>%
  dplyr::relocate(
    dplyr::starts_with("Count_", vars = colnames(.data)),
    .after = "Mass_Tag_ID"
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Mean_", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Fold_change_", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("P_value_A", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("P_value_G", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Flag_A", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Flag_G", vars = colnames(.data)),
    .after = dplyr::last_col()
  )

# Replace all NaN values with NA.
cstan_1_0_3[is.nan(data.matrix(cstan_1_0_3))] <- NA

class(cstan_1_0_3) <- c("statRes", "data.frame")

attr(cstan_1_0_3, "group_DF") <- attr(cfilta_1_0_3, "group_DF")
attr(cstan_1_0_3, "comparisons") <- c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM")
attr(cstan_1_0_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM"),
  Up_total = c(sum(length(which(flag_g_1_0_3[, 1] == 1)),
                   length(which(flag_a_1_0_3[, 1] == 1))),
               sum(length(which(flag_g_1_0_3[, 2] == 1)),
                   length(which(flag_a_1_0_3[, 2] == 1))),
               sum(length(which(flag_g_1_0_3[, 3] == 1)),
                   length(which(flag_a_1_0_3[, 3] == 1)))),
  Down_total = c(sum(length(which(flag_g_1_0_3[, 1] == -1)),
                     length(which(flag_a_1_0_3[, 1] == -1))),
                 sum(length(which(flag_g_1_0_3[, 2] == -1)),
                     length(which(flag_a_1_0_3[, 2] == -1))),
                 sum(length(which(flag_g_1_0_3[, 3] == -1)),
                     length(which(flag_a_1_0_3[, 3] == -1)))),
  Up_anova = c(length(which(flag_a_1_0_3[, 1] == 1)),
               length(which(flag_a_1_0_3[, 2] == 1)),
               length(which(flag_a_1_0_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_0_3[, 1] == -1)),
                 length(which(flag_a_1_0_3[, 2] == -1)),
                 length(which(flag_a_1_0_3[, 3] == -1))),
  Up_gtest = c(length(which(flag_g_1_0_3[, 1] == 1)),
               length(which(flag_g_1_0_3[, 2] == 1)),
               length(which(flag_g_1_0_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_0_3[, 1] == -1)),
                 length(which(flag_g_1_0_3[, 2] == -1)),
                 length(which(flag_g_1_0_3[, 3] == -1))),
  row.names = NULL
)
attr(cstan_1_0_3, "statistical_test") <- "combined"
attr(cstan_1_0_3, "adjustment_method_a_multcomp") <- "none"
attr(cstan_1_0_3, "adjustment_method_g_multcomp") <- "none"
attr(cstan_1_0_3, "adjustment_method_a_fdr") <- "none"
attr(cstan_1_0_3, "adjustment_method_g_fdr") <- "none"
attr(cstan_1_0_3, "pval_thresh") <- 0.05
attr(cstan_1_0_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_0_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_0_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_0_3$e_data[, -1])) /
                    prod(dim(gfilta_1_0_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_0_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(cstan_1_0_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_1_0_3$e_data$Mass_Tag_ID,
  Mock_vs_FM = dplyr::case_when(
    is.na(cstan_1_0_3$Flag_A_Mock_vs_FM) ~
      cstan_1_0_3$Flag_G_Mock_vs_FM,
    (cstan_1_0_3$P_value_A_Mock_vs_FM > 0.05 &
       cstan_1_0_3$P_value_G_Mock_vs_FM < 0.05) ~
      cstan_1_0_3$Flag_G_Mock_vs_FM,
    !is.na(cstan_1_0_3$Flag_A_Mock_vs_FM) ~
      cstan_1_0_3$Flag_A_Mock_vs_FM
  ),
  Mock_vs_AM = dplyr::case_when(
    is.na(cstan_1_0_3$Flag_A_Mock_vs_AM) ~
      cstan_1_0_3$Flag_G_Mock_vs_AM,
    (cstan_1_0_3$P_value_A_Mock_vs_AM > 0.05 &
       cstan_1_0_3$P_value_G_Mock_vs_AM < 0.05) ~
      cstan_1_0_3$Flag_G_Mock_vs_AM,
    !is.na(cstan_1_0_3$Flag_A_Mock_vs_AM) ~
      cstan_1_0_3$Flag_A_Mock_vs_AM
  ),
  FM_vs_AM = dplyr::case_when(
    is.na(cstan_1_0_3$Flag_A_FM_vs_AM) ~
      cstan_1_0_3$Flag_G_FM_vs_AM,
    (cstan_1_0_3$P_value_A_FM_vs_AM > 0.05 &
       cstan_1_0_3$P_value_G_FM_vs_AM < 0.05) ~
      cstan_1_0_3$Flag_G_FM_vs_AM,
    !is.na(cstan_1_0_3$Flag_A_FM_vs_AM) ~
      cstan_1_0_3$Flag_A_FM_vs_AM
  )
)
attr(cstan_1_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(cstan_1_0_3, "data_class") <- "pepData"

# main effects: 1; covariates: 1; groups: 3 ---------------

# Combine the G-test and ANOVA results and place columns in correct order.
cstan_1_1_3 <- dplyr::full_join(gstan_1_1_3[, c(1:4, 11:16)],
                                astan_1_1_3) %>%
  dplyr::relocate(
    dplyr::starts_with("Count_", vars = colnames(.data)),
    .after = "Mass_Tag_ID"
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Mean_", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Fold_change_", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("P_value_A", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("P_value_G", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Flag_A", vars = colnames(.data)),
    .after = dplyr::last_col()
  ) %>%
  dplyr::relocate(
    dplyr::starts_with("Flag_G", vars = colnames(.data)),
    .after = dplyr::last_col()
  )

# Replace all NaN values with NA.
cstan_1_1_3[is.nan(data.matrix(cstan_1_1_3))] <- NA

class(cstan_1_1_3) <- c("statRes", "data.frame")

attr(cstan_1_1_3, "group_DF") <- attr(cfilta_1_1_3, "group_DF")
attr(cstan_1_1_3, "comparisons") <- c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM")
attr(cstan_1_1_3, "number_significant") <- data.frame(
  Comparison = c("Mock_vs_FM", "Mock_vs_AM", "FM_vs_AM"),
  Up_total = c(sum(length(which(flag_g_1_1_3[, 1] == 1)),
                   length(which(flag_a_1_1_3[, 1] == 1))),
               sum(length(which(flag_g_1_1_3[, 2] == 1)),
                   length(which(flag_a_1_1_3[, 2] == 1))),
               sum(length(which(flag_g_1_1_3[, 3] == 1)),
                   length(which(flag_a_1_1_3[, 3] == 1)))),
  Down_total = c(sum(length(which(flag_g_1_1_3[, 1] == -1)),
                     length(which(flag_a_1_1_3[, 1] == -1))),
                 sum(length(which(flag_g_1_1_3[, 2] == -1)),
                     length(which(flag_a_1_1_3[, 2] == -1))),
                 sum(length(which(flag_g_1_1_3[, 3] == -1)),
                     length(which(flag_a_1_1_3[, 3] == -1)))),
  Up_anova = c(length(which(flag_a_1_1_3[, 1] == 1)),
               length(which(flag_a_1_1_3[, 2] == 1)),
               length(which(flag_a_1_1_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_1_3[, 1] == -1)),
                 length(which(flag_a_1_1_3[, 2] == -1)),
                 length(which(flag_a_1_1_3[, 3] == -1))),
  Up_gtest = c(length(which(flag_g_1_1_3[, 1] == 1)),
               length(which(flag_g_1_1_3[, 2] == 1)),
               length(which(flag_g_1_1_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_1_3[, 1] == -1)),
                 length(which(flag_g_1_1_3[, 2] == -1)),
                 length(which(flag_g_1_1_3[, 3] == -1))),
  row.names = NULL
)
attr(cstan_1_1_3, "statistical_test") <- "combined"
attr(cstan_1_1_3, "adjustment_method_a_multcomp") <- "none"
attr(cstan_1_1_3, "adjustment_method_g_multcomp") <- "none"
attr(cstan_1_1_3, "adjustment_method_a_fdr") <- "none"
attr(cstan_1_1_3, "adjustment_method_g_fdr") <- "none"
attr(cstan_1_1_3, "pval_thresh") <- 0.05
attr(cstan_1_1_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_1_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_1_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_1_3$e_data[, -1])) /
                    prod(dim(gfilta_1_1_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_1_3$f_data)[1],
  data_types = NULL,
  batch_info = list(is_bc = FALSE)
)
attr(cstan_1_1_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_1_1_3$e_data$Mass_Tag_ID,
  Mock_vs_FM = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_Mock_vs_FM) ~
      cstan_1_1_3$Flag_G_Mock_vs_FM,
    (cstan_1_1_3$P_value_A_Mock_vs_FM > 0.05 &
       cstan_1_1_3$P_value_G_Mock_vs_FM < 0.05) ~
      cstan_1_1_3$Flag_G_Mock_vs_FM,
    !is.na(cstan_1_1_3$Flag_A_Mock_vs_FM) ~
      cstan_1_1_3$Flag_A_Mock_vs_FM
  ),
  Mock_vs_AM = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_Mock_vs_AM) ~
      cstan_1_1_3$Flag_G_Mock_vs_AM,
    (cstan_1_1_3$P_value_A_Mock_vs_AM > 0.05 &
       cstan_1_1_3$P_value_G_Mock_vs_AM < 0.05) ~
      cstan_1_1_3$Flag_G_Mock_vs_AM,
    !is.na(cstan_1_1_3$Flag_A_Mock_vs_AM) ~
      cstan_1_1_3$Flag_A_Mock_vs_AM
  ),
  FM_vs_AM = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_FM_vs_AM) ~
      cstan_1_1_3$Flag_G_FM_vs_AM,
    (cstan_1_1_3$P_value_A_FM_vs_AM > 0.05 &
       cstan_1_1_3$P_value_G_FM_vs_AM < 0.05) ~
      cstan_1_1_3$Flag_G_FM_vs_AM,
    !is.na(cstan_1_1_3$Flag_A_FM_vs_AM) ~
      cstan_1_1_3$Flag_A_FM_vs_AM
  )
)
attr(cstan_1_1_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "Name",
  techrep_cname = NULL
)
attr(cstan_1_1_3, "data_class") <- "pepData"

# Save standards for paired IMD-ANOVA tests ------------------------------------

save(

  astan_0_0_3,

  astan_1_0_3,
  gstan_1_0_3,
  cstan_1_0_3,

  astan_1_1_3,
  gstan_1_1_3,
  cstan_1_1_3,

  file = file.path("inst/testdata/standards_imd_anova_paired.RData")

)
