context('IMD-ANOVA tests')

test_that('all tests conform to the decrees of the God of Stats',{

  # Load data and prepare omicsData objects ------------------------------------

  # Paired data ---------------

  # load(system.file('testdata',
  #                  'little_pairdata.RData',
  #                  package = 'pmartR'))
  #
  # pairdata <- as.pepData(e_data = edata,
  #                        f_data = fdata,
  #                        e_meta = emeta,
  #                        edata_cname = "Mass_Tag_ID",
  #                        fdata_cname = "Name",
  #                        emeta_cname = "Protein")
  #
  # # Throw some natural logs on the data because the following functions require
  # # the data to be logged.
  # pairdata <- edata_transform(pairdata, "log")

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
  fdataZ$Gender <- sample(c("F", "M"), 12, replace = TRUE)

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

  filta_1_1_3 <- imdanova_filter(pdata_1_1_3)
  afilta_1_1_3 <- applyFilt(filta_1_1_3, pdata_1_1_3,
                            min_nonmiss_anova = 2,
                            remove_singleton_groups = FALSE)
  gfilta_1_1_3 <- applyFilt(filta_1_1_3, pdata_1_1_3,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)
  cfilta_1_1_3 <- applyFilt(filta_1_1_3, pdata_1_1_3,
                            min_nonmiss_anova = 2,
                            min_nonmiss_gtest = 2,
                            remove_singleton_groups = FALSE)

  # ANOVA and G-test functions -------------------------------------------------

  # Function to construct the null space projection matrix.
  proj_mat <- function(X, ngroups){
    #If the X matrix has atleast two rows, find projection matrix
    #into the null space corresponding to X
    if(!is.null(nrow(X))){
      Imat <- diag(1,nrow(X))
      Px <- MASS::ginv(t(X)%*%X)%*%t(X)
      # Zero out the first ngroups rows of Px because we only want to remove the
      # effect of the covariates (which appear in the last n_covariates rows of
      # Px).
      Px[1:ngroups, ] <- 0
      Px <- X%*%Px
      return(Imat-Px)
    }

    #If only one datapoint is left then just return 1
    return(1)
  }


  # Project each row of the data.matrix into X's null space.
  project_to_null <- function(data_mat, Xmatrix, ngroups){

    data_no_x <- data_mat

    for(i in 1:nrow(data_mat)){

      to_rem <- which(is.na(data_mat[i,]))
      if(length(to_rem)>0){
        roi <-  data_mat[i,-to_rem]
        IPxi <- proj_mat(Xmatrix[-to_rem,], ngroups)
        data_no_x[i,-to_rem] <- IPxi%*%matrix(roi,ncol=1)
      }else{
        roi <- data_mat[i,]
        IPxi <- proj_mat(Xmatrix, ngroups)
        data_no_x[i,] <- IPxi%*%matrix(roi,ncol=1)
      }

    }

    return(data.frame(data_no_x))

  }

  a <- function (data, groups) {

    # nas <- !is.na(data)

    nova <- vector(length = nrow(data))

    for (e in 1:nrow(data)) {

      # nova[[e]] <- tryCatch (
      #   summary(
      #     aov(as.numeric(data[e, nas[e, ]]) ~ groups[nas[e, ]])
      #   )[[1]]$`Pr(>F)`[[1]],
      #   error = function (e) {NaN}
      # )

      nova[[e]] <- tryCatch (
        anova(
          lm(as.numeric(data[e, ]) ~ groups)
        )$`Pr(>F)`[[1]],
        error = function (e) {NaN}
      )

    }

    return (nova)

  }

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

    return (unname(knight))

  }

  aflag <- function (grp1, grp2, pvals, cutoff) {

    ninja <- sign(grp1 - grp2)

    ninja[which(pvals >= cutoff)] <- 0

    ninja[is.na(ninja)] <- 0

    return (unname(ninja))

  }

  # Assemble ANOVA standards ---------------------------------------------------

  # main effects: 1; covariates: 0; groups: 2 ---------------

  mean_a_1_0_2 <- data.frame(
    Mean_Infection = rowMeans(afilta_1_0_2$e_data[, 2:10],
                              na.rm = TRUE),
    Mean_Mock = rowMeans(afilta_1_0_2$e_data[, 11:13],
                         na.rm = TRUE)
  )

  pval_a_1_0_2 <- data.frame(
    Infection_vs_Mock = a(data = afilta_1_0_2$e_data[, -1],
                          groups = attr(afilta_1_0_2, "group_DF")$Group)
  )

  flag_a_1_0_2 <- data.frame(
    Infection_vs_Mock = aflag(
      grp1 = mean_a_1_0_2$Mean_Infection,
      grp2 = mean_a_1_0_2$Mean_Mock,
      pvals = pval_a_1_0_2[, 1],
      cutoff = 0.05
    )
  )

  astan_1_0_2 <- data.frame(
    Mass_Tag_ID = afilta_1_0_2$e_data$Mass_Tag_ID,
    Count_Infection = unname(rowSums(!is.na(afilta_1_0_2$e_data[, 2:10]))),
    Count_Mock = unname(rowSums(!is.na(afilta_1_0_2$e_data[, 11:13]))),
    mean_a_1_0_2,
    Fold_change_Infection_vs_Mock = (
      mean_a_1_0_2[, 1] - mean_a_1_0_2[, 2]
    ),
    P_value_A_Infection_vs_Mock = pval_a_1_0_2$Infection_vs_Mock,
    Flag_A_Infection_vs_Mock = flag_a_1_0_2$Infection_vs_Mock,
    row.names = NULL
  )

  class(astan_1_0_2) <- c("statRes", "data.frame")

  attr(astan_1_0_2, "group_DF") <- groupDF_1_0_2
  attr(astan_1_0_2, "comparisons") <- c("Infection_vs_Mock")
  attr(astan_1_0_2, "number_significant") <- data.frame(
    Comparison = c("Infection_vs_Mock"),
    Up_total = c(length(which(flag_a_1_0_2[, 1] == 1))),
    Down_total = c(length(which(flag_a_1_0_2[, 1] == -1))),
    Up_anova = c(length(which(flag_a_1_0_2[, 1] == 1))),
    Down_anova = c(length(which(flag_a_1_0_2[, 1] == -1))),
    Up_gtest = c(0),
    Down_gtest = c(0),
    row.names = NULL
  )
  attr(astan_1_0_2, "statistical_test") <- "anova"
  attr(astan_1_0_2, "adjustment_method_a") <- "none"
  attr(astan_1_0_2, "adjustment_method_g") <- "none"
  attr(astan_1_0_2, "pval_thresh") <- 0.05
  attr(astan_1_0_2, "data_info") <- list(
    data_scale_orig = "abundance",
    data_scale = "log",
    norm_info = list(is_normalized = FALSE),
    num_edata = length(unique(afilta_1_0_2$e_data$Mass_Tag_ID)),
    num_miss_obs = sum(is.na(afilta_1_0_2$e_data[, -1])),
    prop_missing = (sum(is.na(afilta_1_0_2$e_data[, -1])) /
                      prod(dim(afilta_1_0_2$e_data[, -1]))),
    num_samps = dim(afilta_1_0_2$f_data)[1],
    data_types = NULL
  )
  attr(astan_1_0_2, "bpFlags") <- data.frame(
    Mass_Tag_ID = afilta_1_0_2$e_data$Mass_Tag_ID,
    flag_a_1_0_2
  )
  attr(astan_1_0_2, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = "Protein",
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(astan_1_0_2, "data_class") <- "pepData"

  # main effects: 2; covariates: 0; groups: 3 ---------------

  mean_a_2_0_3 <- data.frame(
    Mean_Infection_high = rowMeans(afilta_2_0_3$e_data[, c(2:4, 7:8)],
                                   na.rm = TRUE),
    Mean_Infection_low = rowMeans(afilta_2_0_3$e_data[, c(5:6, 9:10)],
                                  na.rm = TRUE),
    Mean_Mock_none = rowMeans(afilta_2_0_3$e_data[, 11:13],
                              na.rm = TRUE)
  )

  group_counts_2_0_3 <- data.frame(
    nona_iHigh = rowSums(!is.na(afilta_2_0_3$e_data[, c(2:4, 7:8)])),
    nona_iLow = rowSums(!is.na(afilta_2_0_3$e_data[, c(5:6, 9:10)])),
    nona_mNone = rowSums(!is.na(afilta_2_0_3$e_data[, 11:13]))
  )

  nona_grps_2_0_3 <- rowSums(group_counts_2_0_3 != 0)

  sigma_2_0_3 <- dplyr::rowwise(afilta_2_0_3$e_data[, -1]) %>%
    dplyr::mutate(
      stdev = summary(lm(
        dplyr::c_across(Infection1:Mock3) ~ groupDF_2_0_3$Group
      ))$sigma
    ) %>%
    dplyr::pull(stdev)

  nona_counts_2_0_3 <- rowSums(!is.na(afilta_2_0_3$e_data[, -1]))

  diffs_2_0_3 <- mean_a_2_0_3 %>%
    dplyr::mutate(
      diff_high_low = Mean_Infection_high - Mean_Infection_low,
      diff_high_none = Mean_Infection_high - Mean_Mock_none,
      diff_low_none = Mean_Infection_low - Mean_Mock_none
    ) %>%
    dplyr::select(diff_high_low, diff_high_none, diff_low_none)

  test_stat_2_0_3 <- diffs_2_0_3 %>%
    dplyr::mutate(
      stat_high_low = (diff_high_low /
                         (sigma_2_0_3 * sqrt(1/group_counts_2_0_3[, 1] +
                                               1/group_counts_2_0_3[, 2]))),
      stat_high_none = (diff_high_none /
                          (sigma_2_0_3 * sqrt(1/group_counts_2_0_3[, 1] +
                                                1/group_counts_2_0_3[, 3]))),
      stat_low_none = (diff_low_none /
                         (sigma_2_0_3 * sqrt(1/group_counts_2_0_3[, 2] +
                                               1/group_counts_2_0_3[, 3])))
    ) %>%
    dplyr::select(stat_high_low, stat_high_none, stat_low_none)

  pval_a_2_0_3 <- test_stat_2_0_3 %>%
    dplyr::mutate(
      P_value_A_Infection_high_vs_Infection_low = pt(
        q = abs(stat_high_low),
        df = nona_counts_2_0_3 - nona_grps_2_0_3,
        lower.tail = FALSE
      ) * 2,
      P_value_A_Infection_high_vs_Mock_none = pt(
        q = abs(stat_high_none),
        df = nona_counts_2_0_3 - nona_grps_2_0_3,
        lower.tail = FALSE
      ) * 2,
      P_value_A_Infection_low_vs_Mock_none = pt(
        q = abs(stat_low_none),
        df = nona_counts_2_0_3 - nona_grps_2_0_3,
        lower.tail = FALSE
      ) * 2
    ) %>%
    dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                  P_value_A_Infection_high_vs_Mock_none,
                  P_value_A_Infection_low_vs_Mock_none)

  flag_a_2_0_3 <- data.frame(
    Infection_high_vs_Infection_low = aflag(
      grp1 = mean_a_2_0_3$Mean_Infection_high,
      grp2 = mean_a_2_0_3$Mean_Infection_low,
      pvals = pval_a_2_0_3[, 1],
      cutoff = 0.05
    ),
    Infection_high_vs_Mock_none = aflag(
      grp1 = mean_a_2_0_3$Mean_Infection_high,
      grp2 = mean_a_2_0_3$Mean_Mock_none,
      pvals = pval_a_2_0_3[, 2],
      cutoff = 0.05
    ),
    Infection_low_vs_Mock_none = aflag(
      grp1 = mean_a_2_0_3$Mean_Infection_low,
      grp2 = mean_a_2_0_3$Mean_Mock_none,
      pvals = pval_a_2_0_3[, 3],
      cutoff = 0.05
    )
  )

  astan_2_0_3 <- data.frame(
    Mass_Tag_ID = afilta_2_0_3$e_data$Mass_Tag_ID,
    Count_Infection_high = unname(
      rowSums(!is.na(afilta_2_0_3$e_data[, c(2:4, 7:8)]))
    ),
    Count_Infection_low = unname(
      rowSums(!is.na(afilta_2_0_3$e_data[, c(5:6, 9:10)]))
    ),
    Count_Mock_none = unname(rowSums(!is.na(afilta_2_0_3$e_data[, 11:13]))),
    mean_a_2_0_3,
    Fold_change_Infection_high_vs_Infection_low = (
      mean_a_2_0_3[, 1] - mean_a_2_0_3[, 2]
    ),
    Fold_change_Infection_high_vs_Mock_none = (
      mean_a_2_0_3[, 1] - mean_a_2_0_3[, 3]
    ),
    Fold_change_Infection_low_vs_Mock_none = (
      mean_a_2_0_3[, 2] - mean_a_2_0_3[, 3]
    ),
    pval_a_2_0_3,
    Flag_A_Infection_high_vs_Infection_low = flag_a_2_0_3[, 1],
    Flag_A_Infection_high_vs_Mock_none = flag_a_2_0_3[, 2],
    Flag_A_Infection_low_vs_Mock_none = flag_a_2_0_3[, 3],
    row.names = NULL
  )

  class(astan_2_0_3) <- c("statRes", "data.frame")

  attr(astan_2_0_3, "group_DF") <- groupDF_2_0_3
  attr(astan_2_0_3, "comparisons") <- c("Infection_high_vs_Infection_low",
                                        "Infection_high_vs_Mock_none",
                                        "Infection_low_vs_Mock_none")
  attr(astan_2_0_3, "number_significant") <- data.frame(
    Comparison = c("Infection_high_vs_Infection_low",
                   "Infection_high_vs_Mock_none",
                   "Infection_low_vs_Mock_none"),
    Up_total = c(length(which(flag_a_2_0_3[, 1] == 1)),
                 length(which(flag_a_2_0_3[, 2] == 1)),
                 length(which(flag_a_2_0_3[, 3] == 1))),
    Down_total = c(length(which(flag_a_2_0_3[, 1] == -1)),
                   length(which(flag_a_2_0_3[, 2] == -1)),
                   length(which(flag_a_2_0_3[, 3] == -1))),
    Up_anova = c(length(which(flag_a_2_0_3[, 1] == 1)),
                 length(which(flag_a_2_0_3[, 2] == 1)),
                 length(which(flag_a_2_0_3[, 3] == 1))),
    Down_anova = c(length(which(flag_a_2_0_3[, 1] == -1)),
                   length(which(flag_a_2_0_3[, 2] == -1)),
                   length(which(flag_a_2_0_3[, 3] == -1))),
    Up_gtest = c(0, 0, 0),
    Down_gtest = c(0, 0, 0),
    row.names = NULL
  )
  attr(astan_2_0_3, "statistical_test") <- "anova"
  attr(astan_2_0_3, "adjustment_method_a") <- "none"
  attr(astan_2_0_3, "adjustment_method_g") <- "none"
  attr(astan_2_0_3, "pval_thresh") <- 0.05
  attr(astan_2_0_3, "data_info") <- list(
    data_scale_orig = "abundance",
    data_scale = "log",
    norm_info = list(is_normalized = FALSE),
    num_edata = length(unique(afilta_2_0_3$e_data$Mass_Tag_ID)),
    num_miss_obs = sum(is.na(afilta_2_0_3$e_data[, -1])),
    prop_missing = (sum(is.na(afilta_2_0_3$e_data[, -1])) /
                      prod(dim(afilta_2_0_3$e_data[, -1]))),
    num_samps = dim(afilta_2_0_3$f_data)[1],
    data_types = NULL
  )
  attr(astan_2_0_3, "bpFlags") <- data.frame(
    Mass_Tag_ID = afilta_2_0_3$e_data$Mass_Tag_ID,
    flag_a_2_0_3
  )
  attr(astan_2_0_3, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = NULL,
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(astan_2_0_3, "data_class") <- "pepData"

  tukey_stat_2_0_3 <- test_stat_2_0_3 * sqrt(2)
  suppressWarnings(
    tukey_pval_2_0_3 <- tukey_stat_2_0_3 %>%
      dplyr::mutate(
        P_value_A_Infection_high_vs_Infection_low = ptukey(
          abs(stat_high_low),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_2_0_3 - 3,
          lower.tail = FALSE
        ),
        P_value_A_Infection_high_vs_Mock_none = ptukey(
          abs(stat_high_none),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_2_0_3 - 3,
          lower.tail = FALSE
        ),
        P_value_A_Infection_low_vs_Mock_none = ptukey(
          abs(stat_low_none),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_2_0_3 - 3,
          lower.tail = FALSE
        )
      ) %>%
      dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                    P_value_A_Infection_high_vs_Mock_none,
                    P_value_A_Infection_low_vs_Mock_none) %>%
      `rownames<-`(NULL)
  )

  # Set a seed because the mvtnorm::pmvt function has a random process.
  set.seed(5)
  dunnett_2_0_3 <- pval_a_2_0_3 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      no_nas = sum(!is.na(dplyr::c_across(tidyselect::everything())))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dFree = nona_counts_2_0_3 - no_nas) %>%
    cbind(test_stat_2_0_3) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pval_high_low = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_high_low),
               no_nas),
          rep(abs(stat_high_low),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_high_none = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_high_none),
               no_nas),
          rep(abs(stat_high_none),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_low_none = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_low_none),
               no_nas),
          rep(abs(stat_low_none),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      )
    )

  # main effects: 1; covariates: 1; groups: 3 ---------------

  # Adjust the means to remove effect of covariates.
  adj_data_1_1_3 <- project_to_null(
    data_mat = data.matrix(afilta_1_1_3$e_data[, -1]),
    Xmatrix = structure(c(1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0,
                          0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                          1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
                          0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1),
                        .Dim = c(12L, 5L)),
    ngroups = 3
  )

  mean_a_1_1_3 <- data.frame(
    Mean_mutant = rowMeans(adj_data_1_1_3[, c(1, 3:4, 9)],
                           na.rm = TRUE),
    Mean_zombie = rowMeans(adj_data_1_1_3[, c(2, 5:8)],
                           na.rm = TRUE),
    Mean_human = rowMeans(adj_data_1_1_3[, 10:12],
                          na.rm = TRUE)
  )

  group_counts_1_1_3 <- data.frame(
    nona_mutant = rowSums(!is.na(adj_data_1_1_3[, c(1, 3:4, 9)])),
    nona_zombie = rowSums(!is.na(adj_data_1_1_3[, c(2, 5:8)])),
    nona_human = rowSums(!is.na(adj_data_1_1_3[, 10:12]))
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n_grp = sum(dplyr::c_across(nona_mutant:nona_human) == 0)) %>%
    dplyr::ungroup()

  vari_1_1_3 <- adj_data_1_1_3 %>%
    dplyr::mutate(mMutant = mean_a_1_1_3$Mean_mutant,
                  mZombie = mean_a_1_1_3$Mean_zombie,
                  mHuman = mean_a_1_1_3$Mean_human,
                  lg = group_counts_1_1_3$n_grp) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sse = sum(
        c((dplyr::c_across(c(Infection1, Infection3, Infection4, Infection9)) -
             mMutant)^2,
          (dplyr::c_across(c(Infection2, Infection5:Infection8)) - mZombie)^2,
          (dplyr::c_across(Mock1:Mock3) - mHuman)^2),
        na.rm = TRUE
      ),
      vari = sse / (sum(!is.na(dplyr::c_across(Infection1:Mock3))) -
                      (3 - lg) - 2)
    ) %>%
    dplyr::pull(vari)

  nona_counts_1_1_3 <- rowSums(!is.na(adj_data_1_1_3))

  diffs_1_1_3 <- mean_a_1_1_3 %>%
    dplyr::mutate(
      diff_m_z = Mean_mutant - Mean_zombie,
      diff_m_h = Mean_mutant - Mean_human,
      diff_z_h = Mean_zombie - Mean_human
    ) %>%
    dplyr::select(diff_m_z, diff_m_h, diff_z_h)

  test_stat_1_1_3 <- diffs_1_1_3 %>%
    cbind(group_counts_1_1_3) %>%
    dplyr::mutate(
      stat_m_z = (diff_m_z /
                    sqrt((1/nona_mutant +
                            1/nona_zombie) * vari_1_1_3)),
      stat_m_h = (diff_m_h /
                    sqrt((1/nona_mutant +
                            1/nona_human) * vari_1_1_3)),
      stat_z_h = (diff_z_h /
                    sqrt((1/nona_zombie +
                            1/nona_human) * vari_1_1_3))
    ) %>%
    dplyr::select(stat_m_z, stat_m_h, stat_z_h) %>%
    dplyr::ungroup() %>%
    `row.names<-`(NULL)

  pval_a_1_1_3 <- test_stat_1_1_3 %>%
    dplyr::mutate(
      lg = group_counts_1_1_3$n_grp,
      P_value_A_mutant_vs_zombie = pt(
        q = abs(stat_m_z),
        df = nona_counts_1_1_3 - (3 - lg),
        lower.tail = FALSE
      ) * 2,
      P_value_A_mutant_vs_human = pt(
        q = abs(stat_m_h),
        df = nona_counts_1_1_3 - (3 - lg),
        lower.tail = FALSE
      ) * 2,
      P_value_A_zombie_vs_human = pt(
        q = abs(stat_z_h),
        df = nona_counts_1_1_3 - (3 - lg),
        lower.tail = FALSE
      ) * 2
    ) %>%
    dplyr::select(P_value_A_mutant_vs_zombie,
                  P_value_A_mutant_vs_human,
                  P_value_A_zombie_vs_human) %>%
    `row.names<-`(NULL)

  flag_a_1_1_3 <- data.frame(
    mutant_vs_zombie = aflag(
      grp1 = mean_a_1_1_3$Mean_mutant,
      grp2 = mean_a_1_1_3$Mean_zombie,
      pvals = pval_a_1_1_3[, 1],
      cutoff = 0.05
    ),
    mutant_vs_human = aflag(
      grp1 = mean_a_1_1_3$Mean_mutant,
      grp2 = mean_a_1_1_3$Mean_human,
      pvals = pval_a_1_1_3[, 2],
      cutoff = 0.05
    ),
    zombie_vs_human = aflag(
      grp1 = mean_a_1_1_3$Mean_zombie,
      grp2 = mean_a_1_1_3$Mean_human,
      pvals = pval_a_1_1_3[, 3],
      cutoff = 0.05
    )
  )

  astan_1_1_3 <- data.frame(
    Mass_Tag_ID = afilta_1_1_3$e_data$Mass_Tag_ID,
    Count_mutant = group_counts_1_1_3$nona_mutant,
    Count_zombie = group_counts_1_1_3$nona_zombie,
    Count_human = group_counts_1_1_3$nona_human,
    mean_a_1_1_3,
    Fold_change_mutant_vs_zombie = (
      mean_a_1_1_3[, 1] - mean_a_1_1_3[, 2]
    ),
    Fold_change_mutant_vs_human = (
      mean_a_1_1_3[, 1] - mean_a_1_1_3[, 3]
    ),
    Fold_change_zombie_vs_human = (
      mean_a_1_1_3[, 2] - mean_a_1_1_3[, 3]
    ),
    pval_a_1_1_3,
    Flag_A_mutant_vs_zombie = flag_a_1_1_3[, 1],
    Flag_A_mutant_vs_human = flag_a_1_1_3[, 2],
    Flag_A_zombie_vs_human = flag_a_1_1_3[, 3],
    row.names = NULL
  )

  class(astan_1_1_3) <- c("statRes", "data.frame")

  attr(astan_1_1_3, "group_DF") <- groupDF_1_1_3
  attr(astan_1_1_3, "comparisons") <- c("mutant_vs_zombie",
                                        "mutant_vs_human",
                                        "zombie_vs_human")
  attr(astan_1_1_3, "number_significant") <- data.frame(
    Comparison = c("mutant_vs_zombie",
                   "mutant_vs_human",
                   "zombie_vs_human"),
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
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(astan_1_1_3, "data_class") <- "pepData"

  tukey_stat_1_1_3 <- test_stat_1_1_3 * sqrt(2)
  suppressWarnings(
    tukey_pval_1_1_3 <- tukey_stat_1_1_3 %>%
      dplyr::mutate(
        P_value_A_mutant_vs_zombie = ptukey(
          abs(stat_m_z),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_1_1_3 - 3,
          lower.tail = FALSE
        ),
        P_value_A_mutant_vs_human = ptukey(
          abs(stat_m_h),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_1_1_3 - 3,
          lower.tail = FALSE
        ),
        P_value_A_zombie_vs_human = ptukey(
          abs(stat_z_h),
          nranges = 1,
          nmeans = 3,
          df = nona_counts_1_1_3 - 3,
          lower.tail = FALSE
        )
      ) %>%
      dplyr::select(P_value_A_mutant_vs_zombie,
                    P_value_A_mutant_vs_human,
                    P_value_A_zombie_vs_human) %>%
      `rownames<-`(NULL)
  )

  # Set a seed because the mvtnorm::pmvt function has a random process.
  set.seed(4)
  dunnett_1_1_3 <- pval_a_1_1_3 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      no_nas = sum(!is.na(dplyr::c_across(tidyselect::everything())))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dFree = nona_counts_1_1_3 - no_nas) %>%
    cbind(test_stat_1_1_3) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pval_m_z = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_m_z),
               no_nas),
          rep(abs(stat_m_z),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_m_h = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_m_h),
               no_nas),
          rep(abs(stat_m_h),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_z_h = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_z_h),
               no_nas),
          rep(abs(stat_z_h),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      )
    )

  # Generate G-Test standards --------------------------------------------------

  # Construct count matrices: Objects that start with obs represent counts of
  # non-missing (observed) values and objects that start with abs represent
  # counts of missing (absent) values.

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

  # main effects: 1; covariates: 0; groups: 2 ---------------

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

  gstan_1_0_2 <- data.frame(
    Mass_Tag_ID = gfilta_1_0_2$e_data$Mass_Tag_ID,
    Count_Infection = unname(obs_inf_1_0_2),
    Count_Mock = unname(obs_mock_1_0_2),
    mean_1_0_2,
    Fold_change_Infection_vs_Mock = (
      mean_1_0_2[, 1] - mean_1_0_2[, 2]
    ),
    P_value_G_Infection_vs_Mock = pval_g_1_0_2[, 1],
    Flag_G_Infection_vs_Mock = flag_g_1_0_2[, 1],
    row.names = NULL
  )

  class(gstan_1_0_2) <- c("statRes", "data.frame")

  attr(gstan_1_0_2, "group_DF") <- groupDF_1_0_2
  attr(gstan_1_0_2, "comparisons") <- c("Infection_vs_Mock")
  attr(gstan_1_0_2, "number_significant") <- data.frame(
    Comparison = c("Infection_vs_Mock"),
    Up_total = c(length(which(flag_g_1_0_2[, 1] == 1))),
    Down_total = c(length(which(flag_g_1_0_2[, 1] == -1))),
    Up_anova = c(0),
    Down_anova = c(0),
    Up_gtest = c(length(which(flag_g_1_0_2[, 1] == 1))),
    Down_gtest = c(length(which(flag_g_1_0_2[, 1] == -1))),
    row.names = NULL
  )
  attr(gstan_1_0_2, "statistical_test") <- "gtest"
  attr(gstan_1_0_2, "adjustment_method_a") <- "none"
  attr(gstan_1_0_2, "adjustment_method_g") <- "none"
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

  # main effects: 2; covariates: 0; groups: 3 ---------------

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

  gstan_2_0_3 <- data.frame(
    Mass_Tag_ID = gfilta_2_0_3$e_data$Mass_Tag_ID,
    Count_Infection_high = unname(obs_inf_high_2_0_3),
    Count_Infection_low = unname(obs_inf_low_2_0_3),
    Count_Mock_none = unname(obs_mock_2_0_3),
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
    P_value_G_Infection_high_vs_Infection_low = pval_g_2_0_3[, 1],
    P_value_G_Infection_high_vs_Mock_none = pval_g_2_0_3[, 2],
    P_value_G_Infection_low_vs_Mock_none = pval_g_2_0_3[, 3],
    Flag_G_Infection_high_vs_Infection_low = flag_g_2_0_3[, 1],
    Flag_G_Infection_high_vs_Mock_none = flag_g_2_0_3[, 2],
    Flag_G_Infection_low_vs_Mock_none = flag_g_2_0_3[, 3],
    row.names = NULL
  )

  class(gstan_2_0_3) <- c("statRes", "data.frame")

  attr(gstan_2_0_3, "group_DF") <- groupDF_2_0_3
  attr(gstan_2_0_3, "comparisons") <- c("Infection_high_vs_Infection_low",
                                        "Infection_high_vs_Mock_none",
                                        "Infection_low_vs_Mock_none")
  attr(gstan_2_0_3, "number_significant") <- data.frame(
    Comparison = c("Infection_high_vs_Infection_low",
                   "Infection_high_vs_Mock_none",
                   "Infection_low_vs_Mock_none"),
    Up_total = c(length(which(flag_g_2_0_3[, 1] == 1)),
                 length(which(flag_g_2_0_3[, 2] == 1)),
                 length(which(flag_g_2_0_3[, 3] == 1))),
    Down_total = c(length(which(flag_g_2_0_3[, 1] == -1)),
                   length(which(flag_g_2_0_3[, 2] == -1)),
                   length(which(flag_g_2_0_3[, 3] == -1))),
    Up_anova = c(0, 0, 0),
    Down_anova = c(0, 0, 0),
    Up_gtest = c(length(which(flag_g_2_0_3[, 1] == 1)),
                 length(which(flag_g_2_0_3[, 2] == 1)),
                 length(which(flag_g_2_0_3[, 3] == 1))),
    Down_gtest = c(length(which(flag_g_2_0_3[, 1] == -1)),
                   length(which(flag_g_2_0_3[, 2] == -1)),
                   length(which(flag_g_2_0_3[, 3] == -1))),
    row.names = NULL
  )
  attr(gstan_2_0_3, "statistical_test") <- "gtest"
  attr(gstan_2_0_3, "adjustment_method_a") <- "none"
  attr(gstan_2_0_3, "adjustment_method_g") <- "none"
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

  # Create combined standards --------------------------------------------------

  # main effects: 1; covariates: 0; groups: 2 ---------------

  # Combine the G-test and ANOVA results and place columns in correct order.
  cstan_1_0_2 <- dplyr::full_join(gstan_1_0_2[, c(1:3, 7:8)],
                                  astan_1_0_2) %>%
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
  cstan_1_0_2[is.nan(data.matrix(cstan_1_0_2))] <- NA

  class(cstan_1_0_2) <- c("statRes", "data.frame")

  attr(cstan_1_0_2, "group_DF") <- groupDF_1_0_2
  attr(cstan_1_0_2, "comparisons") <- c("Infection_vs_Mock")
  attr(cstan_1_0_2, "number_significant") <- data.frame(
    Comparison = c("Infection_vs_Mock"),
    Up_total = sum(length(which(flag_g_1_0_2[, 1] == 1)),
                   length(which(flag_a_1_0_2[, 1] == 1))),
    Down_total = sum(length(which(flag_g_1_0_2[, 1] == -1)),
                     length(which(flag_a_1_0_2[, 1] == -1))),
    Up_anova = c(length(which(flag_a_1_0_2[, 1] == 1))),
    Down_anova = c(length(which(flag_a_1_0_2[, 1] == -1))),
    Up_gtest = c(length(which(flag_g_1_0_2[, 1] == 1))),
    Down_gtest = c(length(which(flag_g_1_0_2[, 1] == -1))),
    row.names = NULL
  )
  attr(cstan_1_0_2, "statistical_test") <- "combined"
  attr(cstan_1_0_2, "adjustment_method_a") <- "none"
  attr(cstan_1_0_2, "adjustment_method_g") <- "none"
  attr(cstan_1_0_2, "pval_thresh") <- 0.05
  attr(cstan_1_0_2, "data_info") <- list(
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
  attr(cstan_1_0_2, "bpFlags") <- data.frame(
    Mass_Tag_ID = gfilta_1_0_2$e_data$Mass_Tag_ID,
    Infection_vs_Mock = dplyr::case_when(
      is.na(cstan_1_0_2$Flag_A_Infection_vs_Mock) ~
        cstan_1_0_2$Flag_G_Infection_vs_Mock,
      (cstan_1_0_2$P_value_A_Infection_vs_Mock > 0.05 &
         cstan_1_0_2$P_value_G_Infection_vs_Mock < 0.05) ~
        cstan_1_0_2$Flag_G_Infection_vs_Mock,
      !is.na(cstan_1_0_2$Flag_A_Infection_vs_Mock) ~
        cstan_1_0_2$Flag_A_Infection_vs_Mock
    )
  )
  attr(cstan_1_0_2, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = "Protein",
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(cstan_1_0_2, "data_class") <- "pepData"

  # main effects: 2; covariates: 0; groups: 3 ---------------

  # Combine the G-test and ANOVA results and place columns in correct order.
  cstan_2_0_3 <- dplyr::full_join(gstan_2_0_3[, c(1:4, 11:16)],
                                  astan_2_0_3) %>%
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
  cstan_2_0_3[is.nan(data.matrix(cstan_2_0_3))] <- NA

  class(cstan_2_0_3) <- c("statRes", "data.frame")

  attr(cstan_2_0_3, "group_DF") <- groupDF_2_0_3
  attr(cstan_2_0_3, "comparisons") <- c(
    "Infection_high_vs_Infection_low",
    "Infection_high_vs_Mock_none",
    "Infection_low_vs_Mock_none"
  )
  attr(cstan_2_0_3, "number_significant") <- data.frame(
    Comparison = c("Infection_high_vs_Infection_low",
                   "Infection_high_vs_Mock_none",
                   "Infection_low_vs_Mock_none"),
    Up_total = c(sum(length(which(flag_g_2_0_3[, 1] == 1)),
                     length(which(flag_a_2_0_3[, 1] == 1))),
                 sum(length(which(flag_g_2_0_3[, 2] == 1)),
                     length(which(flag_a_2_0_3[, 2] == 1))),
                 sum(length(which(flag_g_2_0_3[, 3] == 1)),
                     length(which(flag_a_2_0_3[, 3] == 1)))),
    Down_total = c(sum(length(which(flag_g_2_0_3[, 1] == -1)),
                       length(which(flag_a_2_0_3[, 1] == -1))),
                   sum(length(which(flag_g_2_0_3[, 2] == -1)),
                       length(which(flag_a_2_0_3[, 2] == -1))),
                   sum(length(which(flag_g_2_0_3[, 3] == -1)),
                       length(which(flag_a_2_0_3[, 3] == -1)))),
    Up_anova = c(length(which(flag_a_2_0_3[, 1] == 1)),
                 length(which(flag_a_2_0_3[, 2] == 1)),
                 length(which(flag_a_2_0_3[, 3] == 1))),
    Down_anova = c(length(which(flag_a_2_0_3[, 1] == -1)),
                   length(which(flag_a_2_0_3[, 2] == -1)),
                   length(which(flag_a_2_0_3[, 3] == -1))),
    Up_gtest = c(length(which(flag_g_2_0_3[, 1] == 1)),
                 length(which(flag_g_2_0_3[, 2] == 1)),
                 length(which(flag_g_2_0_3[, 3] == 1))),
    Down_gtest = c(length(which(flag_g_2_0_3[, 1] == -1)),
                   length(which(flag_g_2_0_3[, 2] == -1)),
                   length(which(flag_g_2_0_3[, 3] == -1))),
    row.names = NULL
  )
  attr(cstan_2_0_3, "statistical_test") <- "combined"
  attr(cstan_2_0_3, "adjustment_method_a") <- "none"
  attr(cstan_2_0_3, "adjustment_method_g") <- "none"
  attr(cstan_2_0_3, "pval_thresh") <- 0.05
  attr(cstan_2_0_3, "data_info") <- list(
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
  attr(cstan_2_0_3, "bpFlags") <- data.frame(
    Mass_Tag_ID = gfilta_2_0_3$e_data$Mass_Tag_ID,
    Infection_high_vs_Infection_low = dplyr::case_when(
      is.na(cstan_2_0_3$Flag_A_Infection_high_vs_Infection_low) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Infection_low,
      is.na(cstan_2_0_3$P_value_A_Infection_high_vs_Infection_low) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Infection_low,
      (cstan_2_0_3$P_value_A_Infection_high_vs_Infection_low > 0.05 &
         cstan_2_0_3$P_value_G_Infection_high_vs_Infection_low < 0.05) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Infection_low,
      !is.na(cstan_2_0_3$Flag_A_Infection_high_vs_Infection_low) ~
        cstan_2_0_3$Flag_A_Infection_high_vs_Infection_low
    ),
    Infection_high_vs_Mock_none = dplyr::case_when(
      is.na(cstan_2_0_3$Flag_A_Infection_high_vs_Mock_none) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Mock_none,
      is.na(cstan_2_0_3$P_value_A_Infection_high_vs_Mock_none) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Mock_none,
      (cstan_2_0_3$P_value_A_Infection_high_vs_Mock_none > 0.05 &
         cstan_2_0_3$P_value_G_Infection_high_vs_Mock_none < 0.05) ~
        cstan_2_0_3$Flag_G_Infection_high_vs_Mock_none,
      !is.na(cstan_2_0_3$Flag_A_Infection_high_vs_Mock_none) ~
        cstan_2_0_3$Flag_A_Infection_high_vs_Mock_none
    ),
    Infection_low_vs_Mock_none = dplyr::case_when(
      is.na(cstan_2_0_3$Flag_A_Infection_low_vs_Mock_none) ~
        cstan_2_0_3$Flag_G_Infection_low_vs_Mock_none,
      is.na(cstan_2_0_3$P_value_A_Infection_low_vs_Mock_none) ~
        cstan_2_0_3$Flag_G_Infection_low_vs_Mock_none,
      (cstan_2_0_3$P_value_A_Infection_low_vs_Mock_none > 0.05 &
         cstan_2_0_3$P_value_G_Infection_low_vs_Mock_none < 0.05) ~
        cstan_2_0_3$Flag_G_Infection_low_vs_Mock_none,
      !is.na(cstan_2_0_3$Flag_A_Infection_low_vs_Mock_none) ~
        cstan_2_0_3$Flag_A_Infection_low_vs_Mock_none
    )
  )
  attr(cstan_2_0_3, "cnames") <- list(
    edata_cname = "Mass_Tag_ID",
    emeta_cname = NULL,
    fdata_cname = "SampleID",
    techrep_cname = NULL
  )
  attr(cstan_2_0_3, "data_class") <- "pepData"

  # Compare IMD-ANOVAly --------------------------------------------------------

  # Default comparisons ---------------

  afruit_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova")
  gfruit_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest")
  cfruit_1_0_2 <- imd_anova(cfilta_1_0_2, test_method = "combined")

  afruit_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova")
  gfruit_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest")
  cfruit_2_0_3 <- imd_anova(cfilta_2_0_3, test_method = "combined")

  afruit_1_1_3 <- imd_anova(afilta_1_1_3,
                            test_method = "anova",
                            covariates = fdataZ[, c(1, 3)])
  gfruit_1_1_3 <- imd_anova(gfilta_1_1_3, test_method = "gtest")
  cfruit_1_1_3 <- imd_anova(cfilta_1_1_3, test_method = "combined")

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

  # Adjusted p-values ---------------

  afruit_bon_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a = "bonferroni")
  afruit_holm_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                 pval_adjust_a = "holm")
  afruit_tuk_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a = "tukey")
  afruit_dun_1_0_2 <- imd_anova(afilta_1_0_2, test_method = "anova",
                                pval_adjust_a = "dunnett")

  afruit_bon_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a = "bonferroni")
  afruit_holm_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                 pval_adjust_a = "holm")
  afruit_tuk_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a = "tukey")
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(4)
  afruit_dun_2_0_3 <- imd_anova(afilta_2_0_3, test_method = "anova",
                                pval_adjust_a = "dunnett")
  # Still gives an error even though all case-vs-control comparisons are being
  # made. Hmmmmmmmm. I took out the if statement that checks if a Dunnett
  # correction can be used because it throws an error when the conditions to be
  # used are met.

  afruit_bon_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a = "bonferroni",
                                covariates = fdataZ[, c(1, 3)])
  afruit_holm_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                 pval_adjust_a = "holm",
                                 covariates = fdataZ[, c(1, 3)])
  afruit_tuk_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a = "tukey",
                                covariates = fdataZ[, c(1, 3)])
  # Set a seed because the mvtnorm::pmvt function--which is called when doing a
  # Dunnett p-value correction--has a random process.
  set.seed(3)
  afruit_dun_1_1_3 <- imd_anova(afilta_1_1_3, test_method = "anova",
                                pval_adjust_a = "dunnett",
                                covariates = fdataZ[, c(1, 3)])

  gfruit_bon_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest",
                                pval_adjust_g = "bonferroni")
  gfruit_holm_1_0_2 <- imd_anova(gfilta_1_0_2, test_method = "gtest",
                                 pval_adjust_g = "holm")
  gfruit_bon_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                pval_adjust_g = "bonferroni")
  gfruit_holm_2_0_3 <- imd_anova(gfilta_2_0_3, test_method = "gtest",
                                 pval_adjust_g = "holm")

  # Holy IMD-ANOVA unit tests, Statman! ----------------------------------------

  # ANOVA: Unadjusted p-values ---------------

  expect_equal(afruit_1_0_2, astan_1_0_2)
  expect_equal(afruit_2_0_3, astan_2_0_3)
  expect_equal(afruit_1_1_3, astan_1_1_3)

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

  expect_equal(
    data.matrix(afruit_bon_2_0_3[, 11:13]),
    pmin(data.matrix(astan_2_0_3[, 11:13] * 3), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  # expect_equal(
  #   afruit_holm_2_0_3[, 11:13],
  #   data.frame(t(apply(astan_2_0_3[, 11:13],
  #                      1,
  #                      p.adjust,
  #                      method = "holm")))
  # )
  expect_equal(
    data.frame(afruit_tuk_2_0_3[, 11:13]),
    tukey_pval_2_0_3
  )
  # Because of the random process in the mvtnorm::pmvt function we test the
  # Dunnett adjusted p-values using a threshold. The threshold is 5 standard
  # deviations above the mean difference between runs of the mvtnorm:pmvt
  # function. To create the threshold the Dunnett adjusted p-values were
  # calculated 10,001 times (the first run is used as a reference) using the
  # following code: imd_anova(afilta_2_0_3, test_method = "anova", pval_adjust =
  # "dunnett"). We calculated the mean of the absolute difference between the
  # p-values from the initial run and each subsequent run. The mean and standard
  # deviation of the mean absolute difference were then calculated for the
  # 10,000 runs.
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_2_0_3[, 11:13]) -
            data.matrix(dunnett_2_0_3[, 9:11])),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # NOTE: Currently commented out because the current version of imd_anova
  # produces p-values greater than 1. Holy inconceivable p-values, Statman!
  expect_equal(
    data.matrix(afruit_bon_1_1_3[, 11:13]),
    pmin(data.matrix(astan_1_1_3[, 11:13] * 3), 1)
  )
  # Calculate the Holm adjusted p-values.
  # Current implementation adjusts p-values even if only one test is performed.
  # expect_equal(
  #   afruit_holm_1_1_3[, 11:13],
  #   data.frame(t(apply(astan_1_1_3[, 11:13],
  #                      1,
  #                      p.adjust,
  #                      method = "holm")))
  # )
  expect_equal(
    data.frame(afruit_tuk_1_1_3[, 11:13]),
    tukey_pval_1_1_3
  )
  expect_true(
    mean(
      abs(data.matrix(afruit_dun_1_1_3[, 11:13]) -
            data.matrix(dunnett_1_1_3[, 9:11])),
      na.rm = TRUE
    ) < 0.0001132809
  )

  # G-Test: Unadjusted p-values ---------------

  expect_equal(gfruit_1_0_2, gstan_1_0_2)
  expect_equal(gfruit_2_0_3, gstan_2_0_3)

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

  # Combined output ---------------

  expect_equal(cfruit_1_0_2, cstan_1_0_2)
  expect_equal(cfruit_2_0_3, cstan_2_0_3)

  # Test argument checks -------------------------------------------------------

  # Maybe add later (depending on the position of the planets and stars).

})
