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

pairdata_1_0_3 <- as.pepData(e_data = edata,
                             f_data = fdata,
                             e_meta = emeta,
                             edata_cname = 'Mass_Tag_ID',
                             fdata_cname = 'Name',
                             emeta_cname = 'Protein')
pairdata_1_0_3 <- edata_transform(pairdata_1_0_3,
                                  data_scale = "log")
pairdata_1_0_3 <- group_designation(pairdata_1_0_3,
                                    main_effects = "Virus")

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
                                    covariates = "Gender")

# Filter IMD-ANOVAly -----------------------------------------------------------

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
#   afilta_1_0_3,
#   gfilta_1_0_3,
#   cfilta_1_0_3,
#
#   file = "/Users/mart077/pmartR/inst/testdata/standards_filter_paired.RData"
# )

# Prepare paired objects -------------------------------------------------------

differential <- data.frame(
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

groupie <- data.frame(
  Name = paste("Pair", 1:15, sep = "_"),
  Group = c(rep("Mock", 5), rep("FM", 5), rep("AM", 5))
)

# Assemble ANOVA standards -----------------------------------------------------

# main effects: 1; covariates: 0; groups: 3 ---------------

mean_a_1_0_3 <- data.frame(
  Mean_Mock = rowMeans(differential[, 2:6],
                            na.rm = TRUE),
  Mean_FM = rowMeans(differential[, 7:11],
                       na.rm = TRUE),
  Mean_AM = rowMeans(differential[, 12:16],
                     na.rm = TRUE)
)

group_counts_1_0_3 <- data.frame(
  nona_Mock = rowSums(!is.na(differential[, 2:6])),
  nona_FM = rowSums(!is.na(differential[, 7:11])),
  nona_AM = rowSums(!is.na(differential[, 12:16]))
)

nona_grps_1_0_3 <- rowSums(group_counts_1_0_3 != 0)

sigma_1_0_3 <- dplyr::rowwise(differential[, -1]) %>%
  dplyr::mutate(
    stdev = tryCatch (
      summary(lm(
        dplyr::c_across(Pair_1:Pair_15) ~ groupie$Group
      ))$sigma,
      error = function (e) {NA}
    )
  ) %>%
  dplyr::pull(stdev)

nona_counts_1_0_3 <- rowSums(!is.na(differential[, -1]))

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
  Mass_Tag_ID = differential$Mass_Tag_ID,
  Count_Mock = unname(
    rowSums(!is.na(differential[, 2:6]))
  ),
  Count_FM = unname(
    rowSums(!is.na(differential[, 7:11]))
  ),
  Count_AM = unname(rowSums(!is.na(differential[, 12:16]))),
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
attr(astan_1_0_3, "adjustment_method_a") <- "none"
attr(astan_1_0_3, "adjustment_method_g") <- "none"
attr(astan_1_0_3, "pval_thresh") <- 0.05
attr(astan_1_0_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(differential$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(differential[, -1])),
  prop_missing = (sum(is.na(differential[, -1])) /
                    prod(dim(differential[, -1]))),
  num_samps = dim(afilta_1_0_3$f_data)[1],
  data_types = NULL
)
attr(astan_1_0_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_1_0_3$e_data$Mass_Tag_ID,
  flag_a_1_0_3
)
attr(astan_1_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(astan_1_0_3, "data_class") <- "pepData"

# main effects: 1; covariates: 1; groups: 3 ---------------

# Adjust the means to remove effect of covariates.
adj_data_1_1_3 <- project_to_null(
  data_mat = data.matrix(differential[, -1]),
  Xmatrix = structure(c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1,
                        0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0,
                        0, 1, 0, 0, 1, 1, 0),
                      .Dim = c(15L, 5L)),
  ngroups = 3
)

cov_df_1_2_3 <- structure(
  c(1, 2, 2, 2, 2, 2, 1, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2,
    2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 0, 2, 2, 0, 2,
    2, 2, 1, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 0, 2, 1, 0,
    2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1,
    2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 0),
  .Dim = c(150L, 1L)
)
colnames(cov_df_1_2_3) <- "df"

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

sigma_1_1_3 <- adj_data_1_1_3 %>%
  dplyr::mutate(mMock = mean_a_1_1_3$Mean_Mock,
                mFM = mean_a_1_1_3$Mean_FM,
                mAM = mean_a_1_1_3$Mean_AM,
                lg = group_counts_1_1_3$n_grp) %>%
  cbind(cov_df_1_2_3) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    sse = sum(
      c((dplyr::c_across(Pair_1:Pair_5) - mMock)^2,
        (dplyr::c_across(Pair_6:Pair_10) - mFM)^2,
        (dplyr::c_across(Pair_11:Pair_15) - mAM)^2),
      na.rm = TRUE
    ),
    vari = sse / (sum(!is.na(dplyr::c_across(Pair_1:Pair_15))) -
                    # first number: number groups
                    # lg: number groups with all missing data
                    # cov_df_1_2_3: degrees of freedom lost due to covariates
                    (3 - lg) - df)
  ) %>%
  dplyr::pull(vari)

nona_counts_1_1_3 <- rowSums(!is.na(adj_data_1_1_3))

diffs_1_1_3 <- mean_a_1_1_3 %>%
  dplyr::mutate(
    diff_m_f = Mean_Mock - Mean_FM,
    diff_m_a = Mean_Mock - Mean_AM,
    diff_f_a = Mean_FM - Mean_AM
  ) %>%
  dplyr::select(diff_m_f, diff_m_a, diff_f_a)

test_stat_1_1_3 <- diffs_1_1_3 %>%
  cbind(group_counts_1_1_3) %>%
  dplyr::mutate(
    stat_m_f = (diff_m_f /
                  sqrt((1/nona_Mock +
                          1/nona_FM) * sigma_1_1_3)),
    stat_m_a = (diff_m_a /
                  sqrt((1/nona_Mock +
                          1/nona_AM) * sigma_1_1_3)),
    stat_f_a = (diff_f_a /
                  sqrt((1/nona_FM +
                          1/nona_AM) * sigma_1_1_3))
  ) %>%
  dplyr::select(stat_m_f, stat_m_a, stat_f_a) %>%
  dplyr::ungroup() %>%
  `row.names<-`(NULL)

pval_a_1_1_3 <- test_stat_1_1_3 %>%
  dplyr::mutate(
    lg = group_counts_1_1_3$n_grp,
    P_value_A_Mock_vs_FM = pt(
      q = abs(stat_m_f),
      df = nona_counts_1_1_3 - (3 - lg),
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_vs_AM = pt(
      q = abs(stat_m_a),
      df = nona_counts_1_1_3 - (3 - lg),
      lower.tail = FALSE
    ) * 2,
    P_value_A_FM_vs_AM = pt(
      q = abs(stat_f_a),
      df = nona_counts_1_1_3 - (3 - lg),
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
  Count_Mock = group_counts_1_1_3$nona_Mock,
  Count_FM = group_counts_1_1_3$nona_FM,
  Count_AM = group_counts_1_1_3$nona_AM,
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
attr(astan_1_1_3, "adjustment_method_a") <- "none"
attr(astan_1_1_3, "adjustment_method_g") <- "none"
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
  data_types = NULL
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
