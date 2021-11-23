# Load data and prepare omicsData objects --------------------------------------

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
fdata2$Level <- c(sample(c("high", "low"), 9, replace = TRUE),
                  "none", "none", "none")

# Create a vector that has three different levels for Condition.
fdataZ <- fdata
set.seed(720)
fdataZ$Condition <- c(sample(c("zombie", "mutant"), 9, replace = TRUE),
                      "human", "human", "human")
fdataZ$Gender <- sample(c("F", "M"), 12, replace = TRUE)
fdataZ$Age <- round(runif(12, min = 19, max = 89), 2)

# Copy edata so the names of the samples can be changed.
edata2_2 <- edata

# Change some of the Infection samples to Mock samples.
names(edata2_2) <- c("Mass_Tag_ID",
                     paste0("Infection", 1:6),
                     paste0("Mock", 1:6))

# Create additional f_data objects with different main effects and covariates.
fdata2_2 <- fdata

# Update the sample names in f_data.
fdata2_2$SampleID <- c(paste0("Infection", 1:6),
                       paste0("Mock", 1:6))

# Update the first main effect to account for changing some infection samples to
# mock samples.
fdata2_2$Condition <- c(rep("Infection", 6),
                        rep("Mock", 6))

fdata2_2$Level <- c("high", "low", "high", "low", "high", "low", "high",
                    "high", "low", "low", "low", "high")
fdata2_2$Gender <- fdataZ$Gender
fdata2_2$Age <- fdataZ$Age

# Create pepData objects, log the data, add group_DF attributes ----------------

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

# One main effect, two covariates, three groups
pdata_1_2_3 <- as.pepData(e_data = edata,
                          f_data = fdataZ,
                          e_meta = emeta,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID",
                          emeta_cname = "Protein")
pdata_1_2_3 <- edata_transform(pdata_1_2_3, "log")
pdata_1_2_3 <- group_designation(omicsData = pdata_1_2_3,
                                 main_effects = "Condition",
                                 covariates = c("Gender", "Age"))
groupDF_1_2_3 <- attr(pdata_1_2_3, "group_DF")

# Two main effects, zero covariates, three groups
pdata_2_0_3 <- as.pepData(e_data = edata,
                          f_data = fdata2,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID")
pdata_2_0_3 <- edata_transform(pdata_2_0_3, "log")
pdata_2_0_3 <- group_designation(omicsData = pdata_2_0_3,
                                 main_effects = c("Condition", "Level"))
groupDF_2_0_3 <- attr(pdata_2_0_3, "group_DF")

# Two main effects, one covariate, four groups
pdata_2_1_4 <- as.pepData(e_data = edata2_2,
                          f_data = fdata2_2,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID")
pdata_2_1_4 <- edata_transform(pdata_2_1_4, "log")
pdata_2_1_4 <- group_designation(omicsData = pdata_2_1_4,
                                 main_effects = c("Condition", "Level"),
                                 covariates = "Gender")
groupDF_2_1_4 <- attr(pdata_2_1_4, "group_DF")

# Two main effects, two covariates, four groups
pdata_2_2_4 <- as.pepData(e_data = edata2_2,
                          f_data = fdata2_2,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID")
pdata_2_2_4 <- edata_transform(pdata_2_2_4, "log")
pdata_2_2_4 <- group_designation(omicsData = pdata_2_2_4,
                                 main_effects = c("Condition", "Level"),
                                 covariates = c("Gender", "Age"))
groupDF_2_2_4 <- attr(pdata_2_2_4, "group_DF")

# Filter IMD-ANOVAly -----------------------------------------------------------

# One main effect, zero covariates, two groups
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

# One main effect, one covariate, three groups
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

# One main effect, two covariates, three groups
filta_1_2_3 <- imdanova_filter(pdata_1_2_3)
afilta_1_2_3 <- applyFilt(filta_1_2_3, pdata_1_2_3,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_1_2_3 <- applyFilt(filta_1_2_3, pdata_1_2_3,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_1_2_3 <- applyFilt(filta_1_2_3, pdata_1_2_3,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

# Two main effects, zero covariates, three groups
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

# Two main effects, one covariate, three groups
filta_2_1_4 <- imdanova_filter(pdata_2_1_4)
afilta_2_1_4 <- applyFilt(filta_2_1_4, pdata_2_1_4,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_2_1_4 <- applyFilt(filta_2_1_4, pdata_2_1_4,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_2_1_4 <- applyFilt(filta_2_1_4, pdata_2_1_4,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

# Two main effects, two covariates, three groups
filta_2_2_4 <- imdanova_filter(pdata_2_2_4)
afilta_2_2_4 <- applyFilt(filta_2_2_4, pdata_2_2_4,
                          min_nonmiss_anova = 2,
                          remove_singleton_groups = FALSE)
gfilta_2_2_4 <- applyFilt(filta_2_2_4, pdata_2_2_4,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)
cfilta_2_2_4 <- applyFilt(filta_2_2_4, pdata_2_2_4,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

# Save filter objects  ---------------------------------------------------------

# save(
#
#   afilta_1_0_2,
#   gfilta_1_0_2,
#   cfilta_1_0_2,
#
#   afilta_1_1_3,
#   gfilta_1_1_3,
#   cfilta_1_1_3,
#
#   afilta_1_2_3,
#   gfilta_1_2_3,
#   cfilta_1_2_3,
#
#   afilta_2_0_3,
#   gfilta_2_0_3,
#   cfilta_2_0_3,
#
#   afilta_2_1_4,
#   gfilta_2_1_4,
#   cfilta_2_1_4,
#
#   afilta_2_2_4,
#   gfilta_2_2_4,
#   cfilta_2_2_4,
#
#   file = "/Users/mart077/pmartR/inst/testdata/standards_filter.RData"
# )

# ANOVA and G-test functions ---------------------------------------------------

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

# This function and the following function, two_factor_anova_r, are used for two
# main effects (with or without covariates) when there is a significant
# interaction between the main effects for some of the biomolecules (rows of
# e_data). These functions are not used for examples when there is no
# significant interaction between the main effects, for any of the biomolecules,
# because the mean and variance computations are the same for every row.
run_two_factor <- function (data, gpData, red_df) {

  #Create design matrix for reduced model, i.e., no interaction effect
  colnames(gpData)[-c(1,2)] <- c("Factor1","Factor2")
  gpData <- cbind(gpData,y=1:nrow(gpData))

  #Create design matrix for the full model, i.e., all first order and
  #interaction effects
  Xred <- unname(model.matrix(lm(y ~ Factor1 + Factor2 - 1, data = gpData)))
  Xfull <- unname(model.matrix(lm(y ~ Factor1*Factor2 - 1, data = gpData)))

  res <- two_factor_anova_r(
    y = data,
    X_full = Xfull,
    X_red = Xred,
    red_df = red_df,
    group_ids = as.numeric(factor(gpData$Group,
                                  levels = unique(gpData$Group)))
  )

  #Get the unique group levels to translate ANOVA parameters into group means
  red_gpData <- dplyr::distinct(gpData, Group, .keep_all=TRUE)
  red_gpData <- dplyr::arrange(red_gpData, y)

  res$par_estimates <- res$par_estimates[, 1:length(red_gpData$Group)]
  colnames(res$par_estimates) <- red_gpData$Group

  return (res)

}

two_factor_anova_r <- function (y,
                                X_full,
                                X_red,
                                red_df,
                                group_ids) {

  n <- nrow(y)
  p_red <- ncol(X_red)
  p_full <- ncol(X_full)
  pval <- vector(mode = "numeric", length = n)
  Fstat <- vector(mode = "numeric", length = n)
  sig_est <- vector(mode = "numeric", length = n)
  par_ests <- vector(mode = "numeric", length = p_full)

  parmat <- matrix(nrow = n,
                   ncol = p_full)
  group_sizes <- matrix(nrow = n,
                        ncol = p_full)

  #Loop over rows in y
  for (i in 1:n) {

    yrowi <- y[i, ]
    to_remove <- which(!is.finite(yrowi))

    yrowi_nona <- yrowi
    X_red_nona <- X_red
    X_full_nona <- X_full
    num_to_remove <- length(to_remove)
    group_ids_nona <- group_ids # -1 #Make group_ids zero indexed
    # Indices in R start at 1. There is no need to subtract 1 from group_ids to
    # correctly subset rows/columns by group.

    if (num_to_remove > 0) {

      # Subset matrices by the to_remove object.
      yrowi_nona <- yrowi_nona[-to_remove]
      X_red_nona <- X_red_nona[-to_remove, ]
      X_full_nona <- X_full_nona[-to_remove, ]
      group_ids_nona <- group_ids_nona[-to_remove]

    }

    # Remove completely empty columns.
    csums <- colSums(X_full_nona)
    zero_cols <- which(csums == 0)
    non_zero_cols <- X_full_nona[, -zero_cols]

    #Subtract off df spent elsewhere, e.g., on covariates
    df_red <- (nrow(X_red_nona) -
                 Matrix::rankMatrix(X_red_nona)[[1]] -
                 red_df[i, ])
    df_full <- (nrow(X_full_nona) -
                  Matrix::rankMatrix(X_full_nona)[[1]] -
                  red_df[i, ])

    PxRed <- X_red_nona %*% MASS::ginv(
      t(X_red_nona) %*% X_red_nona
    ) %*% t(X_red_nona)

    # Create an identity matrix with the same dimension as PxRed
    diag_mat <- diag(nrow(PxRed))

    sigma2_red <- yrowi_nona %*% (diag_mat - PxRed) %*% yrowi_nona / df_red


    if ((df_red - df_full) <= 0) {

      #If interaction can't be estimated, automatically select smaller model
      pval[[i]] <- 1
      Fstat[[i]] <- 0

    } else {

      PxFull <- (X_full_nona %*% MASS::ginv(t(X_full_nona) %*% X_full_nona) %*%
                   t(X_full_nona))

      # Create an identity matrix with the same dimension as PxFull
      diag_mat <- diag(nrow(PxFull))

      sigma2_full <- (yrowi_nona %*% (diag_mat - PxFull) %*% yrowi_nona /
                        df_full)

      Fstat[[i]] <- (sigma2_red  * df_red - sigma2_full * df_full) / sigma2_full
      pval[[i]] <- pf(q = Fstat[[i]],
                      df1 = df_red - df_full,
                      df2 = df_full,
                      lower.tail = FALSE,
                      log.p = FALSE)

      if (!is.finite(pval[[i]])) {

        pval[[i]] <- 1

      }

    }

    if (pval[[i]] < 0.05) {

      #Reject null hypothesis that reduced model is good enough, use full model
      XFinal <- X_full_nona
      sig_est[[i]] <- sigma2_full

    } else {

      XFinal <- X_full_nona

      # Zero out the interaction terms if they're insignificant.
      XFinal[, (p_red + 1):p_full] <- 0
      sig_est[[i]] <- sigma2_red

    }

    #"Parameter estimates" are the group means: Xbeta_hat=X(XpX)^{-1}XpY
    par_ests_temp <- (XFinal %*% MASS::ginv(t(XFinal) %*% XFinal) %*%
                        t(XFinal) %*% yrowi_nona)

    #Find groups that had at least one non-missing value
    group_ids_nona_unq <- unique(group_ids_nona)

    # Fill in the par_ests vector, put NaN if group was missing or the average
    # effect in groups with no missing data
    if (length(group_ids_nona_unq) < p_full) {

      for (j in 1:p_full) {

        missing_gp <- which(group_ids_nona_unq == j)

        if (length(missing_gp) > 0) {

          par_ests[[j]] <- mean(par_ests_temp[which(group_ids_nona == j)])

        } else {

          par_ests[[j]] <- NA

        }

      }

    } else {

      for (k in 1:p_full) {

        par_ests[[k]] <- mean(par_ests_temp[which(group_ids_nona == k)])

      }

    }

    parmat[i, ] <- par_ests

    gsizes <- rep(0, p_full)

    # Compute group sizes after accounting for NaNs
    for (j in 1:p_full) {

      size_j <- which(group_ids_nona == j)

      if (length(size_j) > 0) {

        gsizes[[j]] <- length(size_j)

      }

    }

    group_sizes[i, ] <- gsizes # For now don't return the interaction groups

  }

  return (list(par_estimates = parmat,
               group_sizes = group_sizes,
               Sigma2 = sig_est,
               Fstats = Fstat,
               pvalue = pval))

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

  ninja <- rep(0, length(pvals))

  star <- which(pvals < cutoff)

  ninja[star] <- sign(grp1 - grp2)[star]

  return (ninja)

}

# Assemble ANOVA standards -----------------------------------------------------

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

cov_df_1_2_3 <- 2

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

nona_grps_1_1_3 <- rowSums(group_counts_1_1_3[, 1:3] != 0)

sigma_1_1_3 <- adj_data_1_1_3 %>%
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
                    # first number: number groups
                    # lg: number groups with all missing data
                    # cov_df_1_2_3: degrees of freedom lost due to covariates
                    (3 - lg) - cov_df_1_2_3)
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
                          1/nona_zombie) * sigma_1_1_3)),
    stat_m_h = (diff_m_h /
                  sqrt((1/nona_mutant +
                          1/nona_human) * sigma_1_1_3)),
    stat_z_h = (diff_z_h /
                  sqrt((1/nona_zombie +
                          1/nona_human) * sigma_1_1_3))
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

n_tests_1_1_3 <- rowSums(!is.na(tukey_stat_1_1_3))

suppressWarnings(
  tukey_pval_1_1_3 <- tukey_stat_1_1_3 %>%
    dplyr::mutate(
      P_value_A_mutant_vs_zombie = ptukey(
        abs(stat_m_z),
        nranges = 1,
        nmeans = n_tests_1_1_3,
        df = nona_counts_1_1_3 - n_tests_1_1_3,
        lower.tail = FALSE
      ),
      P_value_A_mutant_vs_human = ptukey(
        abs(stat_m_h),
        nranges = 1,
        nmeans = n_tests_1_1_3,
        df = nona_counts_1_1_3 - n_tests_1_1_3,
        lower.tail = FALSE
      ),
      P_value_A_zombie_vs_human = ptukey(
        abs(stat_z_h),
        nranges = 1,
        nmeans = n_tests_1_1_3,
        df = nona_counts_1_1_3 - n_tests_1_1_3,
        lower.tail = FALSE
      )
    ) %>%
    dplyr::select(P_value_A_mutant_vs_zombie,
                  P_value_A_mutant_vs_human,
                  P_value_A_zombie_vs_human) %>%
    `rownames<-`(NULL) %>%
    data.matrix()
)

# Find rows where just one test occurred.
just_one_1_1_3 <- which(rowSums(!is.na(tukey_stat_1_1_3)) == 1)

# Nab index where each test occurred in the subsetted p-value data frame.
not_na_1_1_3 <- which(!is.na(pval_a_1_1_3[just_one_1_1_3, ]))

# Replace NAs (from only one test being performed) with the original p-value.
tukey_pval_1_1_3[just_one_1_1_3, ][not_na_1_1_3] <- data.matrix(
  pval_a_1_1_3[just_one_1_1_3, ]
)[not_na_1_1_3]

# Convert back to a data frame.
tukey_pval_1_1_3 <- data.frame(tukey_pval_1_1_3)

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
  ) %>%
  dplyr::select(pval_m_z, pval_m_h, pval_z_h)

# main effects: 1; covariates: 2; groups: 3 ---------------

# Adjust the means to remove effect of covariates.
adj_data_1_2_3 <- project_to_null(
  data_mat = data.matrix(afilta_1_2_3$e_data[, -1]),
  Xmatrix = structure(
    c(1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1,
      1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0,
      1, 0, 1, 0, 0, 1, 0, 0, 1, 59.68, 37.16, 69.4, 73.12,
      69.27, 47.68, 53.79, 19.45, 27.06, 30.07, 22.04, 31.57),
    .Dim = c(12L, 6L)
  ),
  ngroups = 3
)

cov_df_1_2_3 <- 3

mean_a_1_2_3 <- data.frame(
  Mean_mutant = rowMeans(adj_data_1_2_3[, c(1, 3:4, 9)],
                         na.rm = TRUE),
  Mean_zombie = rowMeans(adj_data_1_2_3[, c(2, 5:8)],
                         na.rm = TRUE),
  Mean_human = rowMeans(adj_data_1_2_3[, 10:12],
                        na.rm = TRUE)
)

group_counts_1_2_3 <- data.frame(
  nona_mutant = rowSums(!is.na(adj_data_1_2_3[, c(1, 3:4, 9)])),
  nona_zombie = rowSums(!is.na(adj_data_1_2_3[, c(2, 5:8)])),
  nona_human = rowSums(!is.na(adj_data_1_2_3[, 10:12]))
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n_grp = sum(dplyr::c_across(nona_mutant:nona_human) == 0)) %>%
  dplyr::ungroup()

nona_grps_1_2_3 <- rowSums(group_counts_1_2_3[, 1:3] != 0)

sigma_1_2_3 <- adj_data_1_2_3 %>%
  dplyr::mutate(mMutant = mean_a_1_2_3$Mean_mutant,
                mZombie = mean_a_1_2_3$Mean_zombie,
                mHuman = mean_a_1_2_3$Mean_human,
                lg = group_counts_1_2_3$n_grp) %>%
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
                    # first number: number groups
                    # lg: number groups with all missing data
                    # cov_df_1_2_3: degrees of freedom lost due to covariates
                    (3 - lg) - cov_df_1_2_3)
  ) %>%
  dplyr::pull(vari)

nona_counts_1_2_3 <- rowSums(!is.na(adj_data_1_2_3))

diffs_1_2_3 <- mean_a_1_2_3 %>%
  dplyr::mutate(
    diff_m_z = Mean_mutant - Mean_zombie,
    diff_m_h = Mean_mutant - Mean_human,
    diff_z_h = Mean_zombie - Mean_human
  ) %>%
  dplyr::select(diff_m_z, diff_m_h, diff_z_h)

suppressWarnings(
  test_stat_1_2_3 <- diffs_1_2_3 %>%
    cbind(group_counts_1_2_3) %>%
    dplyr::mutate(
      stat_m_z = (diff_m_z /
                    sqrt((1/nona_mutant +
                            1/nona_zombie) * sigma_1_2_3)),
      stat_m_h = (diff_m_h /
                    sqrt((1/nona_mutant +
                            1/nona_human) * sigma_1_2_3)),
      stat_z_h = (diff_z_h /
                    sqrt((1/nona_zombie +
                            1/nona_human) * sigma_1_2_3))
    ) %>%
    dplyr::select(stat_m_z, stat_m_h, stat_z_h) %>%
    dplyr::ungroup() %>%
    `row.names<-`(NULL)
)

pval_a_1_2_3 <- test_stat_1_2_3 %>%
  dplyr::mutate(
    lg = group_counts_1_2_3$n_grp,
    P_value_A_mutant_vs_zombie = pt(
      q = abs(stat_m_z),
      df = nona_counts_1_2_3 - (3 - lg),
      lower.tail = FALSE
    ) * 2,
    P_value_A_mutant_vs_human = pt(
      q = abs(stat_m_h),
      df = nona_counts_1_2_3 - (3 - lg),
      lower.tail = FALSE
    ) * 2,
    P_value_A_zombie_vs_human = pt(
      q = abs(stat_z_h),
      df = nona_counts_1_2_3 - (3 - lg),
      lower.tail = FALSE
    ) * 2
  ) %>%
  dplyr::select(P_value_A_mutant_vs_zombie,
                P_value_A_mutant_vs_human,
                P_value_A_zombie_vs_human) %>%
  `row.names<-`(NULL)

flag_a_1_2_3 <- data.frame(
  mutant_vs_zombie = aflag(
    grp1 = mean_a_1_2_3$Mean_mutant,
    grp2 = mean_a_1_2_3$Mean_zombie,
    pvals = pval_a_1_2_3[, 1],
    cutoff = 0.05
  ),
  mutant_vs_human = aflag(
    grp1 = mean_a_1_2_3$Mean_mutant,
    grp2 = mean_a_1_2_3$Mean_human,
    pvals = pval_a_1_2_3[, 2],
    cutoff = 0.05
  ),
  zombie_vs_human = aflag(
    grp1 = mean_a_1_2_3$Mean_zombie,
    grp2 = mean_a_1_2_3$Mean_human,
    pvals = pval_a_1_2_3[, 3],
    cutoff = 0.05
  )
)

astan_1_2_3 <- data.frame(
  Mass_Tag_ID = afilta_1_2_3$e_data$Mass_Tag_ID,
  Count_mutant = group_counts_1_2_3$nona_mutant,
  Count_zombie = group_counts_1_2_3$nona_zombie,
  Count_human = group_counts_1_2_3$nona_human,
  mean_a_1_2_3,
  Fold_change_mutant_vs_zombie = (
    mean_a_1_2_3[, 1] - mean_a_1_2_3[, 2]
  ),
  Fold_change_mutant_vs_human = (
    mean_a_1_2_3[, 1] - mean_a_1_2_3[, 3]
  ),
  Fold_change_zombie_vs_human = (
    mean_a_1_2_3[, 2] - mean_a_1_2_3[, 3]
  ),
  pval_a_1_2_3,
  Flag_A_mutant_vs_zombie = flag_a_1_2_3[, 1],
  Flag_A_mutant_vs_human = flag_a_1_2_3[, 2],
  Flag_A_zombie_vs_human = flag_a_1_2_3[, 3],
  row.names = NULL
)

class(astan_1_2_3) <- c("statRes", "data.frame")

attr(astan_1_2_3, "group_DF") <- groupDF_1_2_3
attr(astan_1_2_3, "comparisons") <- c("mutant_vs_zombie",
                                      "mutant_vs_human",
                                      "zombie_vs_human")
attr(astan_1_2_3, "number_significant") <- data.frame(
  Comparison = c("mutant_vs_zombie",
                 "mutant_vs_human",
                 "zombie_vs_human"),
  Up_total = c(length(which(flag_a_1_2_3[, 1] == 1)),
               length(which(flag_a_1_2_3[, 2] == 1)),
               length(which(flag_a_1_2_3[, 3] == 1))),
  Down_total = c(length(which(flag_a_1_2_3[, 1] == -1)),
                 length(which(flag_a_1_2_3[, 2] == -1)),
                 length(which(flag_a_1_2_3[, 3] == -1))),
  Up_anova = c(length(which(flag_a_1_2_3[, 1] == 1)),
               length(which(flag_a_1_2_3[, 2] == 1)),
               length(which(flag_a_1_2_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_2_3[, 1] == -1)),
                 length(which(flag_a_1_2_3[, 2] == -1)),
                 length(which(flag_a_1_2_3[, 3] == -1))),
  Up_gtest = c(0, 0, 0),
  Down_gtest = c(0, 0, 0),
  row.names = NULL
)
attr(astan_1_2_3, "statistical_test") <- "anova"
attr(astan_1_2_3, "adjustment_method_a") <- "none"
attr(astan_1_2_3, "adjustment_method_g") <- "none"
attr(astan_1_2_3, "pval_thresh") <- 0.05
attr(astan_1_2_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(afilta_1_2_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_1_2_3$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_1_2_3$e_data[, -1])) /
                    prod(dim(afilta_1_2_3$e_data[, -1]))),
  num_samps = dim(afilta_1_2_3$f_data)[1],
  data_types = NULL
)
attr(astan_1_2_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_1_2_3$e_data$Mass_Tag_ID,
  flag_a_1_2_3
)
attr(astan_1_2_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(astan_1_2_3, "data_class") <- "pepData"

tukey_stat_1_2_3 <- test_stat_1_2_3 * sqrt(2)

n_tests_1_2_3 <- rowSums(!is.na(tukey_stat_1_2_3))

suppressWarnings(
  tukey_pval_1_2_3 <- tukey_stat_1_2_3 %>%
    dplyr::mutate(
      P_value_A_mutant_vs_zombie = ptukey(
        abs(stat_m_z),
        nranges = 1,
        nmeans = n_tests_1_2_3,
        df = nona_counts_1_2_3 - n_tests_1_2_3,
        lower.tail = FALSE
      ),
      P_value_A_mutant_vs_human = ptukey(
        abs(stat_m_h),
        nranges = 1,
        nmeans = n_tests_1_2_3,
        df = nona_counts_1_2_3 - n_tests_1_2_3,
        lower.tail = FALSE
      ),
      P_value_A_zombie_vs_human = ptukey(
        abs(stat_z_h),
        nranges = 1,
        nmeans = n_tests_1_2_3,
        df = nona_counts_1_2_3 - n_tests_1_2_3,
        lower.tail = FALSE
      )
    ) %>%
    dplyr::select(P_value_A_mutant_vs_zombie,
                  P_value_A_mutant_vs_human,
                  P_value_A_zombie_vs_human) %>%
    `rownames<-`(NULL) %>%
    data.matrix()
)

# Find rows where just one test occurred.
just_one_1_2_3 <- which(rowSums(!is.na(tukey_stat_1_2_3)) == 1)

# Nab index where each test occurred in the subsetted p-value data frame.
not_na_1_2_3 <- which(!is.na(pval_a_1_2_3[just_one_1_2_3, ]))

# Replace NAs (from only one test being performed) with the original p-value.
tukey_pval_1_2_3[just_one_1_2_3, ][not_na_1_2_3] <- data.matrix(
  pval_a_1_2_3[just_one_1_2_3, ]
)[not_na_1_2_3]

# Convert back to a data frame.
tukey_pval_1_2_3 <- data.frame(tukey_pval_1_2_3)

# Set a seed because the mvtnorm::pmvt function has a random process.
set.seed(4)
suppressWarnings(
  dunnett_1_2_3 <- pval_a_1_2_3 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      no_nas = sum(!is.na(dplyr::c_across(tidyselect::everything())))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dFree = nona_counts_1_2_3 - no_nas) %>%
    cbind(test_stat_1_2_3) %>%
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
    ) %>%
    dplyr::select(pval_m_z, pval_m_h, pval_z_h)
)

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

n_tests_2_0_3 <- rowSums(!is.na(tukey_stat_2_0_3))

suppressWarnings(
  tukey_pval_2_0_3 <- tukey_stat_2_0_3 %>%
    dplyr::mutate(
      P_value_A_Infection_high_vs_Infection_low = ptukey(
        abs(stat_high_low),
        nranges = 1,
        nmeans = n_tests_2_0_3,
        df = nona_counts_2_0_3 - n_tests_2_0_3,
        lower.tail = FALSE
      ),
      P_value_A_Infection_high_vs_Mock_none = ptukey(
        abs(stat_high_none),
        nranges = 1,
        nmeans = n_tests_2_0_3,
        df = nona_counts_2_0_3 - n_tests_2_0_3,
        lower.tail = FALSE
      ),
      P_value_A_Infection_low_vs_Mock_none = ptukey(
        abs(stat_low_none),
        nranges = 1,
        nmeans = n_tests_2_0_3,
        df = nona_counts_2_0_3 - n_tests_2_0_3,
        lower.tail = FALSE
      )
    ) %>%
    dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                  P_value_A_Infection_high_vs_Mock_none,
                  P_value_A_Infection_low_vs_Mock_none) %>%
    `rownames<-`(NULL) %>%
    data.matrix()
)

# Find rows where just one test occurred.
just_one_2_0_3 <- which(rowSums(!is.na(tukey_stat_2_0_3)) == 1)

# Nab index where each test occurred in the subsetted p-value data frame.
not_na_2_0_3 <- which(!is.na(pval_a_2_0_3[just_one_2_0_3, ]))

# Replace NAs (from only one test being performed) with the original p-value.
tukey_pval_2_0_3[just_one_2_0_3, ][not_na_2_0_3] <- data.matrix(
  pval_a_2_0_3[just_one_2_0_3, ]
)[not_na_2_0_3]

# Convert back to a data frame.
tukey_pval_2_0_3 <- data.frame(tukey_pval_2_0_3)

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
  ) %>%
  dplyr::select(pval_high_low, pval_high_none, pval_low_none)

# main effects: 2; covariates: 1; groups: 4 ---------------

# Adjust the means to remove effect of covariates.
adj_data_2_1_4 <- project_to_null(
  data_mat = data.matrix(afilta_2_1_4$e_data[, -1]),
  Xmatrix = structure(c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
                        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
                        0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0,
                        1, 0, 0, 1),
                      .Dim = c(12L, 6L)),
  ngroups = 4
)

cov_df_2_1_4 <- 2

cobra <- run_two_factor(data = data.matrix(adj_data_2_1_4),
                        gpData = groupDF_2_1_4,
                        red_df = matrix(cov_df_2_1_4,
                                        nrow = nrow(adj_data_2_1_4),
                                        ncol = 1))

group_counts_2_1_4 <- data.frame(
  nona_Infection_high = rowSums(!is.na(adj_data_2_1_4[, c(1, 3, 5)])),
  nona_Infection_low = rowSums(!is.na(adj_data_2_1_4[, c(2, 4, 6)])),
  nona_Mock_high = rowSums(!is.na(adj_data_2_1_4[, c(7, 8, 12)])),
  nona_Mock_low = rowSums(!is.na(adj_data_2_1_4[, c(9, 10, 11)]))
) %>%
  dplyr::ungroup()

nona_grps_2_1_4 <- unname(rowSums(group_counts_2_1_4 != 0))

nona_counts_2_1_4 <- unname(rowSums(!is.na(adj_data_2_1_4)))

mean_a_2_1_4 <- data.frame(
  Mean_Infection_high = cobra$par_estimates[, 1],
  Mean_Infection_low = cobra$par_estimates[, 2],
  Mean_Mock_high = cobra$par_estimates[, 3],
  Mean_Mock_low = cobra$par_estimates[, 4]
)

sigma_2_1_4 <- cobra$Sigma2

diffs_2_1_4 <- mean_a_2_1_4 %>%
  dplyr::mutate(
    diff_ih_il = Mean_Infection_high - Mean_Infection_low,
    diff_ih_mh = Mean_Infection_high - Mean_Mock_high,
    diff_ih_ml = Mean_Infection_high - Mean_Mock_low,
    diff_il_mh = Mean_Infection_low - Mean_Mock_high,
    diff_il_ml = Mean_Infection_low - Mean_Mock_low,
    diff_mh_ml = Mean_Mock_high - Mean_Mock_low
  ) %>%
  dplyr::select(diff_ih_il, diff_ih_mh, diff_ih_ml,
                diff_il_mh, diff_il_ml, diff_mh_ml)

test_stat_2_1_4 <- diffs_2_1_4 %>%
  cbind(group_counts_2_1_4) %>%
  dplyr::mutate(
    stat_ih_il = (diff_ih_il /
                    sqrt((1/nona_Infection_high +
                            1/nona_Infection_low) * sigma_2_1_4)),
    stat_ih_mh = (diff_ih_mh /
                    sqrt((1/nona_Infection_high +
                            1/nona_Mock_high) * sigma_2_1_4)),
    stat_ih_ml = (diff_ih_ml /
                    sqrt((1/nona_Infection_high +
                            1/nona_Mock_low) * sigma_2_1_4)),
    stat_il_mh = (diff_il_mh /
                    sqrt((1/nona_Infection_low +
                            1/nona_Mock_high) * sigma_2_1_4)),
    stat_il_ml = (diff_il_ml /
                    sqrt((1/nona_Infection_low +
                            1/nona_Mock_low) * sigma_2_1_4)),
    stat_mh_ml = (diff_mh_ml /
                    sqrt((1/nona_Mock_high +
                            1/nona_Mock_low) * sigma_2_1_4))
  ) %>%
  dplyr::select(stat_ih_il, stat_ih_mh, stat_ih_ml,
                stat_il_mh, stat_il_ml, stat_mh_ml) %>%
  dplyr::ungroup() %>%
  `row.names<-`(NULL)

pval_a_2_1_4 <- test_stat_2_1_4 %>%
  dplyr::mutate(
    lg = group_counts_2_1_4$n_grp,
    P_value_A_Infection_high_vs_Infection_low = pt(
      q = abs(stat_ih_il),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_high = pt(
      q = abs(stat_ih_mh),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_low = pt(
      q = abs(stat_ih_ml),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_high = pt(
      q = abs(stat_il_mh),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_low = pt(
      q = abs(stat_il_ml),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_high_vs_Mock_low = pt(
      q = abs(stat_mh_ml),
      df = nona_counts_2_1_4 - nona_grps_2_1_4,
      lower.tail = FALSE
    ) * 2
  ) %>%
  dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                P_value_A_Infection_high_vs_Mock_high,
                P_value_A_Infection_high_vs_Mock_low,
                P_value_A_Infection_low_vs_Mock_high,
                P_value_A_Infection_low_vs_Mock_low,
                P_value_A_Mock_high_vs_Mock_low) %>%
  `row.names<-`(NULL)

flag_a_2_1_4 <- data.frame(
  Infection_high_vs_Infection_low = aflag(
    grp1 = mean_a_2_1_4$Mean_Infection_high,
    grp2 = mean_a_2_1_4$Mean_Infection_low,
    pvals = pval_a_2_1_4[, 1],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_high = aflag(
    grp1 = mean_a_2_1_4$Mean_Infection_high,
    grp2 = mean_a_2_1_4$Mean_Mock_high,
    pvals = pval_a_2_1_4[, 2],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4$Mean_Infection_high,
    grp2 = mean_a_2_1_4$Mean_Mock_low,
    pvals = pval_a_2_1_4[, 3],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_high = aflag(
    grp1 = mean_a_2_1_4$Mean_Infection_low,
    grp2 = mean_a_2_1_4$Mean_Mock_high,
    pvals = pval_a_2_1_4[, 4],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4$Mean_Infection_low,
    grp2 = mean_a_2_1_4$Mean_Mock_low,
    pvals = pval_a_2_1_4[, 5],
    cutoff = 0.05
  ),
  Mock_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4$Mean_Mock_high,
    grp2 = mean_a_2_1_4$Mean_Mock_low,
    pvals = pval_a_2_1_4[, 6],
    cutoff = 0.05
  )
)

astan_2_1_4 <- data.frame(
  Mass_Tag_ID = afilta_2_1_4$e_data$Mass_Tag_ID,
  Count_Infection_high = group_counts_2_1_4$nona_Infection_high,
  Count_Infection_low = group_counts_2_1_4$nona_Infection_low,
  Count_Mock_high = group_counts_2_1_4$nona_Mock_high,
  Count_Mock_low = group_counts_2_1_4$nona_Mock_low,
  mean_a_2_1_4,
  Fold_change_Infection_high_vs_Infection_low = (
    mean_a_2_1_4[, 1] - mean_a_2_1_4[, 2]
  ),
  Fold_change_Infection_high_vs_Mock_high = (
    mean_a_2_1_4[, 1] - mean_a_2_1_4[, 3]
  ),
  Fold_change_Infection_high_vs_Mock_low = (
    mean_a_2_1_4[, 1] - mean_a_2_1_4[, 4]
  ),
  Fold_change_Infection_low_vs_Mock_high = (
    mean_a_2_1_4[, 2] - mean_a_2_1_4[, 3]
  ),
  Fold_change_Infection_low_vs_Mock_low = (
    mean_a_2_1_4[, 2] - mean_a_2_1_4[, 4]
  ),
  Fold_change_Mock_high_vs_Mock_low = (
    mean_a_2_1_4[, 3] - mean_a_2_1_4[, 4]
  ),
  pval_a_2_1_4,
  Flag_A_Infection_high_vs_Infection_low = flag_a_2_1_4[, 1],
  Flag_A_Infection_high_vs_Mock_high = flag_a_2_1_4[, 2],
  Flag_A_Infection_high_vs_Mock_low = flag_a_2_1_4[, 3],
  Flag_A_Infection_low_vs_Mock_high = flag_a_2_1_4[, 4],
  Flag_A_Infection_low_vs_Mock_low = flag_a_2_1_4[, 5],
  Flag_A_Mock_high_vs_Mock_low = flag_a_2_1_4[, 6],
  row.names = NULL
)

class(astan_2_1_4) <- c("statRes", "data.frame")

attr(astan_2_1_4, "group_DF") <- groupDF_2_1_4
attr(astan_2_1_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(astan_2_1_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(length(which(flag_a_2_1_4[, 1] == 1)),
               length(which(flag_a_2_1_4[, 2] == 1)),
               length(which(flag_a_2_1_4[, 3] == 1)),
               length(which(flag_a_2_1_4[, 4] == 1)),
               length(which(flag_a_2_1_4[, 5] == 1)),
               length(which(flag_a_2_1_4[, 6] == 1))),
  Down_total = c(length(which(flag_a_2_1_4[, 1] == -1)),
                 length(which(flag_a_2_1_4[, 2] == -1)),
                 length(which(flag_a_2_1_4[, 3] == -1)),
                 length(which(flag_a_2_1_4[, 4] == -1)),
                 length(which(flag_a_2_1_4[, 5] == -1)),
                 length(which(flag_a_2_1_4[, 6] == -1))),
  Up_anova = c(length(which(flag_a_2_1_4[, 1] == 1)),
               length(which(flag_a_2_1_4[, 2] == 1)),
               length(which(flag_a_2_1_4[, 3] == 1)),
               length(which(flag_a_2_1_4[, 4] == 1)),
               length(which(flag_a_2_1_4[, 5] == 1)),
               length(which(flag_a_2_1_4[, 6] == 1))),
  Down_anova = c(length(which(flag_a_2_1_4[, 1] == -1)),
                 length(which(flag_a_2_1_4[, 2] == -1)),
                 length(which(flag_a_2_1_4[, 3] == -1)),
                 length(which(flag_a_2_1_4[, 4] == -1)),
                 length(which(flag_a_2_1_4[, 5] == -1)),
                 length(which(flag_a_2_1_4[, 6] == -1))),
  Up_gtest = c(0, 0, 0, 0, 0, 0),
  Down_gtest = c(0, 0, 0, 0, 0, 0),
  row.names = NULL
)
attr(astan_2_1_4, "statistical_test") <- "anova"
attr(astan_2_1_4, "adjustment_method_a") <- "none"
attr(astan_2_1_4, "adjustment_method_g") <- "none"
attr(astan_2_1_4, "pval_thresh") <- 0.05
attr(astan_2_1_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(afilta_2_1_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_2_1_4$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_2_1_4$e_data[, -1])) /
                    prod(dim(afilta_2_1_4$e_data[, -1]))),
  num_samps = dim(afilta_2_1_4$f_data)[1],
  data_types = NULL
)
attr(astan_2_1_4, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_2_1_4$e_data$Mass_Tag_ID,
  flag_a_2_1_4
)
attr(astan_2_1_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(astan_2_1_4, "data_class") <- "pepData"

tukey_stat_2_1_4 <- test_stat_2_1_4 * sqrt(2)

n_tests_2_1_4 <- rowSums(!is.na(tukey_stat_2_1_4))

suppressWarnings(
  tukey_pval_2_1_4 <- tukey_stat_2_1_4 %>%
    dplyr::mutate(
      P_value_A_Infection_high_vs_Infection_low = ptukey(
        abs(stat_ih_il),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_high_vs_Mock_high = ptukey(
        abs(stat_ih_mh),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_high_vs_Mock_low = ptukey(
        abs(stat_ih_ml),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_low_vs_Mock_high = ptukey(
        abs(stat_il_mh),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_low_vs_Mock_low = ptukey(
        abs(stat_il_ml),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      ),
      P_value_A_Mock_high_vs_Mock_low = ptukey(
        abs(stat_mh_ml),
        nranges = 1,
        nmeans = n_tests_2_1_4,
        df = nona_counts_2_1_4 - n_tests_2_1_4,
        lower.tail = FALSE
      )
    ) %>%
    dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                  P_value_A_Infection_high_vs_Mock_high,
                  P_value_A_Infection_high_vs_Mock_low,
                  P_value_A_Infection_low_vs_Mock_high,
                  P_value_A_Infection_low_vs_Mock_low,
                  P_value_A_Mock_high_vs_Mock_low,) %>%
    `rownames<-`(NULL) %>%
    data.matrix()
)

# Find rows where just one test occurred.
just_one_2_1_4 <- which(rowSums(!is.na(tukey_stat_2_1_4)) == 1)

# Nab index where each test occurred in the subsetted p-value data frame.
not_na_2_1_4 <- which(!is.na(pval_a_2_1_4[just_one_2_1_4, ]))

# Replace NAs (from only one test being performed) with the original p-value.
tukey_pval_2_1_4[just_one_2_1_4, ][not_na_2_1_4] <- data.matrix(
  pval_a_2_1_4[just_one_2_1_4, ]
)[not_na_2_1_4]

# Convert back to a data frame.
tukey_pval_2_1_4 <- data.frame(tukey_pval_2_1_4)

# Set a seed because the mvtnorm::pmvt function has a random process.
set.seed(4)
dunnett_2_1_4 <- pval_a_2_1_4 %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    no_nas = sum(!is.na(dplyr::c_across(tidyselect::everything())))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dFree = nona_counts_2_1_4 - no_nas) %>%
  cbind(test_stat_2_1_4) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    pval_ih_il = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_ih_il),
             no_nas),
        rep(abs(stat_ih_il),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    ),
    pval_ih_mh = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_ih_mh),
             no_nas),
        rep(abs(stat_ih_mh),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    ),
    pval_ih_ml = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_ih_ml),
             no_nas),
        rep(abs(stat_ih_ml),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    ),
    pval_il_mh = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_il_mh),
             no_nas),
        rep(abs(stat_il_mh),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    ),
    pval_il_ml = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_il_ml),
             no_nas),
        rep(abs(stat_il_ml),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    ),
    pval_mh_ml = tryCatch (
      1 - mvtnorm::pmvt(
        -rep(abs(stat_mh_ml),
             no_nas),
        rep(abs(stat_mh_ml),
            no_nas),
        df = dFree
      ),
      error = function (e) {NA}
    )
  ) %>%
  dplyr::select(pval_ih_il, pval_ih_mh, pval_ih_ml,
                pval_il_mh, pval_il_ml, pval_mh_ml)

# main effects: 2; covariates: 2; groups: 4 ---------------

# Adjust the means to remove effect of covariates.
adj_data_2_2_4 <- project_to_null(
  data_mat = data.matrix(afilta_2_2_4$e_data[, -1]),
  Xmatrix = structure(c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
                        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
                        0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0,
                        1, 0, 0, 1, 59.68, 37.16, 69.4, 73.12, 69.27, 47.68,
                        53.79, 19.45, 27.06, 30.07, 22.04, 31.57),
                      .Dim = c(12L, 7L)),
  ngroups = 4
)

cov_df_2_2_4 <- 3

kungfu <- run_two_factor(data = data.matrix(adj_data_2_2_4),
                         gpData = groupDF_2_2_4,
                         red_df = matrix(cov_df_2_2_4,
                                         nrow = nrow(adj_data_2_2_4),
                                         ncol = 1))

mean_a_2_2_4 <- data.frame(
  Mean_Infection_high = kungfu$par_estimates[, 1],
  Mean_Infection_low = kungfu$par_estimates[, 2],
  Mean_Mock_high = kungfu$par_estimates[, 3],
  Mean_Mock_low = kungfu$par_estimates[, 4]
)

sigma_2_2_4 <- kungfu$Sigma2

group_counts_2_2_4 <- data.frame(
  nona_Infection_high = rowSums(!is.na(adj_data_2_2_4[, c(1, 3, 5)])),
  nona_Infection_low = rowSums(!is.na(adj_data_2_2_4[, c(2, 4, 6)])),
  nona_Mock_high = rowSums(!is.na(adj_data_2_2_4[, c(7, 8, 12)])),
  nona_Mock_low = rowSums(!is.na(adj_data_2_2_4[, c(9, 10, 11)]))
) %>%
  dplyr::ungroup()

nona_grps_2_2_4 <- rowSums(group_counts_2_2_4 != 0)

nona_counts_2_2_4 <- rowSums(!is.na(adj_data_2_2_4))

diffs_2_2_4 <- mean_a_2_2_4 %>%
  dplyr::mutate(
    diff_ih_il = Mean_Infection_high - Mean_Infection_low,
    diff_ih_mh = Mean_Infection_high - Mean_Mock_high,
    diff_ih_ml = Mean_Infection_high - Mean_Mock_low,
    diff_il_mh = Mean_Infection_low - Mean_Mock_high,
    diff_il_ml = Mean_Infection_low - Mean_Mock_low,
    diff_mh_ml = Mean_Mock_high - Mean_Mock_low
  ) %>%
  dplyr::select(diff_ih_il, diff_ih_mh, diff_ih_ml,
                diff_il_mh, diff_il_ml, diff_mh_ml)

test_stat_2_2_4 <- diffs_2_2_4 %>%
  cbind(group_counts_2_2_4) %>%
  dplyr::mutate(
    stat_ih_il = (diff_ih_il /
                    sqrt((1/nona_Infection_high +
                            1/nona_Infection_low) * sigma_2_2_4)),
    stat_ih_mh = (diff_ih_mh /
                    sqrt((1/nona_Infection_high +
                            1/nona_Mock_high) * sigma_2_2_4)),
    stat_ih_ml = (diff_ih_ml /
                    sqrt((1/nona_Infection_high +
                            1/nona_Mock_low) * sigma_2_2_4)),
    stat_il_mh = (diff_il_mh /
                    sqrt((1/nona_Infection_low +
                            1/nona_Mock_high) * sigma_2_2_4)),
    stat_il_ml = (diff_il_ml /
                    sqrt((1/nona_Infection_low +
                            1/nona_Mock_low) * sigma_2_2_4)),
    stat_mh_ml = (diff_mh_ml /
                    sqrt((1/nona_Mock_high +
                            1/nona_Mock_low) * sigma_2_2_4))
  ) %>%
  dplyr::select(stat_ih_il, stat_ih_mh, stat_ih_ml,
                stat_il_mh, stat_il_ml, stat_mh_ml) %>%
  dplyr::ungroup() %>%
  `row.names<-`(NULL)

pval_a_2_2_4 <- test_stat_2_2_4 %>%
  dplyr::mutate(
    lg = group_counts_2_2_4$n_grp,
    P_value_A_Infection_high_vs_Infection_low = pt(
      q = abs(stat_ih_il),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_high = pt(
      q = abs(stat_ih_mh),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_low = pt(
      q = abs(stat_ih_ml),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_high = pt(
      q = abs(stat_il_mh),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_low = pt(
      q = abs(stat_il_ml),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_high_vs_Mock_low = pt(
      q = abs(stat_mh_ml),
      df = nona_counts_2_2_4 - nona_grps_2_2_4,
      lower.tail = FALSE
    ) * 2
  ) %>%
  dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                P_value_A_Infection_high_vs_Mock_high,
                P_value_A_Infection_high_vs_Mock_low,
                P_value_A_Infection_low_vs_Mock_high,
                P_value_A_Infection_low_vs_Mock_low,
                P_value_A_Mock_high_vs_Mock_low) %>%
  `row.names<-`(NULL)

flag_a_2_2_4 <- data.frame(
  Infection_high_vs_Infection_low = aflag(
    grp1 = mean_a_2_2_4$Mean_Infection_high,
    grp2 = mean_a_2_2_4$Mean_Infection_low,
    pvals = pval_a_2_2_4[, 1],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_high = aflag(
    grp1 = mean_a_2_2_4$Mean_Infection_high,
    grp2 = mean_a_2_2_4$Mean_Mock_high,
    pvals = pval_a_2_2_4[, 2],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4$Mean_Infection_high,
    grp2 = mean_a_2_2_4$Mean_Mock_low,
    pvals = pval_a_2_2_4[, 3],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_high = aflag(
    grp1 = mean_a_2_2_4$Mean_Infection_low,
    grp2 = mean_a_2_2_4$Mean_Mock_high,
    pvals = pval_a_2_2_4[, 4],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4$Mean_Infection_low,
    grp2 = mean_a_2_2_4$Mean_Mock_low,
    pvals = pval_a_2_2_4[, 5],
    cutoff = 0.05
  ),
  Mock_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4$Mean_Mock_high,
    grp2 = mean_a_2_2_4$Mean_Mock_low,
    pvals = pval_a_2_2_4[, 6],
    cutoff = 0.05
  )
)

astan_2_2_4 <- data.frame(
  Mass_Tag_ID = afilta_2_2_4$e_data$Mass_Tag_ID,
  Count_Infection_high = group_counts_2_2_4$nona_Infection_high,
  Count_Infection_low = group_counts_2_2_4$nona_Infection_low,
  Count_Mock_high = group_counts_2_2_4$nona_Mock_high,
  Count_Mock_low = group_counts_2_2_4$nona_Mock_low,
  mean_a_2_2_4,
  Fold_change_Infection_high_vs_Infection_low = (
    mean_a_2_2_4[, 1] - mean_a_2_2_4[, 2]
  ),
  Fold_change_Infection_high_vs_Mock_high = (
    mean_a_2_2_4[, 1] - mean_a_2_2_4[, 3]
  ),
  Fold_change_Infection_high_vs_Mock_low = (
    mean_a_2_2_4[, 1] - mean_a_2_2_4[, 4]
  ),
  Fold_change_Infection_low_vs_Mock_high = (
    mean_a_2_2_4[, 2] - mean_a_2_2_4[, 3]
  ),
  Fold_change_Infection_low_vs_Mock_low = (
    mean_a_2_2_4[, 2] - mean_a_2_2_4[, 4]
  ),
  Fold_change_Mock_high_vs_Mock_low = (
    mean_a_2_2_4[, 3] - mean_a_2_2_4[, 4]
  ),
  pval_a_2_2_4,
  Flag_A_Infection_high_vs_Infection_low = flag_a_2_2_4[, 1],
  Flag_A_Infection_high_vs_Mock_high = flag_a_2_2_4[, 2],
  Flag_A_Infection_high_vs_Mock_low = flag_a_2_2_4[, 3],
  Flag_A_Infection_low_vs_Mock_high = flag_a_2_2_4[, 4],
  Flag_A_Infection_low_vs_Mock_low = flag_a_2_2_4[, 5],
  Flag_A_Mock_high_vs_Mock_low = flag_a_2_2_4[, 6],
  row.names = NULL
)

class(astan_2_2_4) <- c("statRes", "data.frame")

attr(astan_2_2_4, "group_DF") <- groupDF_2_2_4
attr(astan_2_2_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(astan_2_2_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(length(which(flag_a_2_2_4[, 1] == 1)),
               length(which(flag_a_2_2_4[, 2] == 1)),
               length(which(flag_a_2_2_4[, 3] == 1)),
               length(which(flag_a_2_2_4[, 4] == 1)),
               length(which(flag_a_2_2_4[, 5] == 1)),
               length(which(flag_a_2_2_4[, 6] == 1))),
  Down_total = c(length(which(flag_a_2_2_4[, 1] == -1)),
                 length(which(flag_a_2_2_4[, 2] == -1)),
                 length(which(flag_a_2_2_4[, 3] == -1)),
                 length(which(flag_a_2_2_4[, 4] == -1)),
                 length(which(flag_a_2_2_4[, 5] == -1)),
                 length(which(flag_a_2_2_4[, 6] == -1))),
  Up_anova = c(length(which(flag_a_2_2_4[, 1] == 1)),
               length(which(flag_a_2_2_4[, 2] == 1)),
               length(which(flag_a_2_2_4[, 3] == 1)),
               length(which(flag_a_2_2_4[, 4] == 1)),
               length(which(flag_a_2_2_4[, 5] == 1)),
               length(which(flag_a_2_2_4[, 6] == 1))),
  Down_anova = c(length(which(flag_a_2_2_4[, 1] == -1)),
                 length(which(flag_a_2_2_4[, 2] == -1)),
                 length(which(flag_a_2_2_4[, 3] == -1)),
                 length(which(flag_a_2_2_4[, 4] == -1)),
                 length(which(flag_a_2_2_4[, 5] == -1)),
                 length(which(flag_a_2_2_4[, 6] == -1))),
  Up_gtest = c(0, 0, 0, 0, 0, 0),
  Down_gtest = c(0, 0, 0, 0, 0, 0),
  row.names = NULL
)
attr(astan_2_2_4, "statistical_test") <- "anova"
attr(astan_2_2_4, "adjustment_method_a") <- "none"
attr(astan_2_2_4, "adjustment_method_g") <- "none"
attr(astan_2_2_4, "pval_thresh") <- 0.05
attr(astan_2_2_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(afilta_2_2_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(afilta_2_2_4$e_data[, -1])),
  prop_missing = (sum(is.na(afilta_2_2_4$e_data[, -1])) /
                    prod(dim(afilta_2_2_4$e_data[, -1]))),
  num_samps = dim(afilta_2_2_4$f_data)[1],
  data_types = NULL
)
attr(astan_2_2_4, "bpFlags") <- data.frame(
  Mass_Tag_ID = afilta_2_2_4$e_data$Mass_Tag_ID,
  flag_a_2_2_4
)
attr(astan_2_2_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(astan_2_2_4, "data_class") <- "pepData"

tukey_stat_2_2_4 <- test_stat_2_2_4 * sqrt(2)

n_tests_2_2_4 <- rowSums(!is.na(tukey_stat_2_2_4))

suppressWarnings(
  tukey_pval_2_2_4 <- tukey_stat_2_2_4 %>%
    dplyr::mutate(
      P_value_A_Infection_high_vs_Infection_low = ptukey(
        abs(stat_ih_il),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_high_vs_Mock_high = ptukey(
        abs(stat_ih_mh),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_high_vs_Mock_low = ptukey(
        abs(stat_ih_ml),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_low_vs_Mock_high = ptukey(
        abs(stat_il_mh),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      ),
      P_value_A_Infection_low_vs_Mock_low = ptukey(
        abs(stat_il_ml),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      ),
      P_value_A_Mock_high_vs_Mock_low = ptukey(
        abs(stat_mh_ml),
        nranges = 1,
        nmeans = n_tests_2_2_4,
        df = nona_counts_2_2_4 - n_tests_2_2_4,
        lower.tail = FALSE
      )
    ) %>%
    dplyr::select(P_value_A_Infection_high_vs_Infection_low,
                  P_value_A_Infection_high_vs_Mock_high,
                  P_value_A_Infection_high_vs_Mock_low,
                  P_value_A_Infection_low_vs_Mock_high,
                  P_value_A_Infection_low_vs_Mock_low,
                  P_value_A_Mock_high_vs_Mock_low,) %>%
    `rownames<-`(NULL) %>%
    data.matrix()
)

# Find rows where just one test occurred.
just_one_2_2_4 <- which(rowSums(!is.na(tukey_stat_2_2_4)) == 1)

# Nab index where each test occurred in the subsetted p-value data frame.
not_na_2_2_4 <- which(!is.na(pval_a_2_2_4[just_one_2_2_4, ]))

# Replace NAs (from only one test being performed) with the original p-value.
tukey_pval_2_2_4[just_one_2_2_4, ][not_na_2_2_4] <- data.matrix(
  pval_a_2_2_4[just_one_2_2_4, ]
)[not_na_2_2_4]

# Convert back to a data frame.
tukey_pval_2_2_4 <- data.frame(tukey_pval_2_2_4)

# Set a seed because the mvtnorm::pmvt function has a random process.
set.seed(4)
suppressWarnings(
  dunnett_2_2_4 <- pval_a_2_2_4 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      no_nas = sum(!is.na(dplyr::c_across(tidyselect::everything())))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dFree = nona_counts_2_2_4 - no_nas) %>%
    cbind(test_stat_2_2_4) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pval_ih_il = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_ih_il),
               no_nas),
          rep(abs(stat_ih_il),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_ih_mh = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_ih_mh),
               no_nas),
          rep(abs(stat_ih_mh),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_ih_ml = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_ih_ml),
               no_nas),
          rep(abs(stat_ih_ml),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_il_mh = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_il_mh),
               no_nas),
          rep(abs(stat_il_mh),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_il_ml = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_il_ml),
               no_nas),
          rep(abs(stat_il_ml),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      ),
      pval_mh_ml = tryCatch (
        1 - mvtnorm::pmvt(
          -rep(abs(stat_mh_ml),
               no_nas),
          rep(abs(stat_mh_ml),
              no_nas),
          df = dFree
        ),
        error = function (e) {NA}
      )
    ) %>%
    dplyr::select(pval_ih_il, pval_ih_mh, pval_ih_ml,
                  pval_il_mh, pval_il_ml, pval_mh_ml)
)

# Generate G-Test standards ----------------------------------------------------

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

# Create combined standards ----------------------------------------------------

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

# Save standards for IMD-ANOVA tests -------------------------------------------

# save(
#   astan_1_0_2,
#   gstan_1_0_2,
#   cstan_1_0_2,
#
#   astan_1_1_3,
#   astan_1_2_3,
#
#   astan_2_0_3,
#   gstan_2_0_3,
#   cstan_2_0_3,
#
#   astan_2_1_4,
#   astan_2_2_4,
#
#   tukey_pval_1_1_3,
#   tukey_pval_1_2_3,
#   tukey_pval_2_0_3,
#   tukey_pval_2_1_4,
#   tukey_pval_2_2_4,
#
#   dunnett_1_1_3,
#   dunnett_1_2_3,
#   dunnett_2_0_3,
#   dunnett_2_1_4,
#   dunnett_2_2_4,
#
#   file = "/Users/mart077/pmartR/inst/testdata/standards_imd_anova.RData"
# )
