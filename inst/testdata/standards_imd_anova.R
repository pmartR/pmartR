# Load functions for calculating IMD-ANOVA standards ---------------------------

source(system.file("testdata",
                    "imd_anova_standard_fns.R",
                    package = "pmartR"))

# Load data and prepare omicsData objects --------------------------------------

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

# Two main effects, one covariate, four groups
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

# Two main effects, two covariates, four groups
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

# Assemble ANOVA standards -----------------------------------------------------

Xmatrix_1_1_3 <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 
1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0), dim = c(12L, 4L), dimnames = list(
    c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", 
    "12"), c("(Intercept)", "Groupzombie", "Grouphuman", "GenderM"
    )), assign = c(0L, 1L, 1L, 2L), contrasts = list(Group = "contr.treatment", 
    Gender = "contr.treatment"))

Xmatrix_1_2_3 <-  structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 
1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 59.68, 37.16, 69.4, 73.12, 69.27, 
47.68, 53.79, 19.45, 27.06, 30.07, 22.04, 31.57), dim = c(12L, 
5L))

Xmatrix_2_1_4  <-  structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0), dim = c(12L, 4L))

# with interaction effect
Xmatrix_2_1_4_full  <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 
0), dim = c(12L, 5L))

Xmatrix_2_2_4 <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 59.68, 37.16, 69.4, 73.12, 69.27, 
47.68, 53.79, 19.45, 27.06, 30.07, 22.04, 31.57), dim = c(12L, 
5L))

# with interaction effect
Xmatrix_2_2_4_full <- structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 
1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 59.68, 37.16, 69.4, 73.12, 69.27, 
47.68, 53.79, 19.45, 27.06, 30.07, 22.04, 31.57, 0, 0, 0, 0, 
0, 0, 0, 0, 1, 1, 1, 0), dim = c(12L, 6L))

# Selected Full/Reduced Model
which_X_2_1_4_a <- structure(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0))

which_X_2_2_4_a <- structure(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0))

which_X_2_1_4_g <- structure(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

which_X_2_2_4_g <- structure(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))


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
attr(astan_1_0_2, "adjustment_method_a_multcomp") <- "none"
attr(astan_1_0_2, "adjustment_method_g_multcomp") <- "none"
attr(astan_1_0_2, "adjustment_method_a_fdr") <- "none"
attr(astan_1_0_2, "adjustment_method_g_fdr") <- "none"
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
attr(astan_1_0_2, "which_X") <- rep(0, nrow(astan_1_0_2))
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
data_1_1_3 = afilta_1_1_3$e_data[, -1]

Betas = compute_betas(data_mat = data.matrix(data_1_1_3), Xmatrix = data.matrix(Xmatrix_1_1_3))
group_sampnames <- groupDF_1_1_3[,get_fdata_cname(afilta_1_1_3)]
groupData <- groupDF_1_1_3[group_sampnames %in% colnames(afilta_1_1_3$e_data),]
groupData <- groupData %>% 
  dplyr::left_join(afilta_1_1_3$f_data)
groupData[,"Condition"] <- lapply(groupData["Condition"], function(x) factor(x, levels=unique(x)))

pred_grid <- get_pred_grid(groupData, "Condition", covariate_names = "Gender")

mean_a_1_1_3 <- get_lsmeans(data = data_1_1_3, xmatrix = Xmatrix_1_1_3, pred_grid = pred_grid, Betas = Betas)

beta_to_mu = pred_grid
beta_to_mu[,4] = 0
beta_to_mu = unique(beta_to_mu)
cmat = rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
cmat = cmat %*% beta_to_mu

test_values <- get_test_values(afilta_1_1_3$e_data[, -1], Xmatrix_1_1_3, cmat)

group_counts_1_1_3 <- data.frame(
  nona_mutant = rowSums(!is.na(data_1_1_3[, c(1, 3:4, 9)])),
  nona_zombie = rowSums(!is.na(data_1_1_3[, c(2, 5:8)])),
  nona_human = rowSums(!is.na(data_1_1_3[, 10:12]))
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n_grp = sum(dplyr::c_across(nona_mutant:nona_human) == 0)) %>%
  dplyr::ungroup()

nona_grps_1_1_3 <- rowSums(group_counts_1_1_3[, 1:3] != 0)
nona_counts_1_1_3 <- rowSums(!is.na(data_1_1_3))

# set mean_a_1_1_3 to NA for groups with zero observations
mean_a_1_1_3[group_counts_1_1_3[,-4] == 0] <- NA

diffs_1_1_3 <- mean_a_1_1_3 %>%
  dplyr::mutate(
    diff_m_z = Mean_mutant - Mean_zombie,
    diff_m_h = Mean_mutant - Mean_human,
    diff_z_h = Mean_zombie - Mean_human
  ) %>%
  dplyr::select(diff_m_z, diff_m_h, diff_z_h)

diff_denoms <- test_values$diff_denoms
colnames(diff_denoms) <- c("C1", "C2", "C3")

test_stat_1_1_3 <- diffs_1_1_3 %>%
  cbind(group_counts_1_1_3) %>%
  cbind(diff_denoms) %>%
  dplyr::mutate(
    stat_m_z = (diff_m_z / C1),
    stat_m_h = (diff_m_h / C2),
    stat_z_h = (diff_z_h / C3)
  ) %>%
  dplyr::select(stat_m_z, stat_m_h, stat_z_h) %>%
  dplyr::ungroup() %>%
  `row.names<-`(NULL)

pval_a_1_1_3 <- test_stat_1_1_3 %>%
  dplyr::mutate(
    P_value_A_mutant_vs_zombie = pt(
      q = abs(stat_m_z),
      df = nona_counts_1_1_3 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_mutant_vs_human = pt(
      q = abs(stat_m_h),
      df = nona_counts_1_1_3 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_zombie_vs_human = pt(
      q = abs(stat_z_h),
      df = nona_counts_1_1_3 - test_values$ranks,
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
  data_types = NULL
)
attr(astan_1_1_3, "which_X") <- rep(0, nrow(astan_1_1_3))
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
    no_nas = sum(!is.na(dplyr::c_across(dplyr::everything())))
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

data_1_2_3 = afilta_1_2_3$e_data[, -1]

Betas = compute_betas(data_mat = data.matrix(data_1_2_3), Xmatrix = data.matrix(Xmatrix_1_2_3))
group_sampnames <- groupDF_1_2_3[,get_fdata_cname(afilta_1_2_3)]
groupData <- groupDF_1_2_3[group_sampnames %in% colnames(afilta_1_2_3$e_data),]
groupData <- groupData %>% 
  dplyr::left_join(afilta_1_2_3$f_data)
groupData[,"Condition"] <- lapply(groupData["Condition"], function(x) factor(x, levels=unique(x)))
pred_grid <- get_pred_grid(groupData, "Condition", c("Gender", "Age"))

mean_a_1_2_3 <- get_lsmeans(data = data_1_2_3, xmatrix = Xmatrix_1_2_3, pred_grid = pred_grid, Betas = Betas, continuous_covar_inds = 5)

beta_to_mu = pred_grid
beta_to_mu[,4] = 0
beta_to_mu = unique(beta_to_mu)
cmat = rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
cmat = cmat %*% beta_to_mu

group_counts_1_2_3 <- data.frame(
  nona_mutant = rowSums(!is.na(data_1_2_3[, c(1, 3:4, 9)])),
  nona_zombie = rowSums(!is.na(data_1_2_3[, c(2, 5:8)])),
  nona_human = rowSums(!is.na(data_1_2_3[, 10:12]))
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(n_grp = sum(dplyr::c_across(nona_mutant:nona_human) == 0)) %>%
  dplyr::ungroup()

nona_grps_1_2_3 <- rowSums(group_counts_1_2_3[, 1:3] != 0)
nona_counts_1_2_3 <- rowSums(!is.na(data_1_2_3))

mean_a_1_2_3[group_counts_1_2_3[,-4] == 0] <- NA

diffs_1_2_3 <- mean_a_1_2_3 %>%
  dplyr::mutate(
    diff_m_z = Mean_mutant - Mean_zombie,
    diff_m_h = Mean_mutant - Mean_human,
    diff_z_h = Mean_zombie - Mean_human
  ) %>%
  dplyr::select(diff_m_z, diff_m_h, diff_z_h)

test_values <- get_test_values(afilta_1_2_3$e_data[, -1], Xmatrix_1_2_3, cmat)
diff_denoms <- test_values$diff_denoms
colnames(diff_denoms) <- c("C1", "C2", "C3")

suppressWarnings(
  test_stat_1_2_3 <- diffs_1_2_3 %>%
    cbind(diff_denoms) %>%
    cbind(group_counts_1_2_3) %>%
    dplyr::mutate(
      stat_m_z = (diff_m_z / C1),
      stat_m_h = (diff_m_h / C2),
      stat_z_h = (diff_z_h / C3)
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
      df = nona_counts_1_2_3 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_mutant_vs_human = pt(
      q = abs(stat_m_h),
      df = nona_counts_1_2_3 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_zombie_vs_human = pt(
      q = abs(stat_z_h),
      df = nona_counts_1_2_3 - test_values$ranks,
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
attr(astan_1_2_3, "adjustment_method_a_multcomp") <- "none"
attr(astan_1_2_3, "adjustment_method_g_multcomp") <- "none"
attr(astan_1_2_3, "adjustment_method_a_fdr") <- "none"
attr(astan_1_2_3, "adjustment_method_g_fdr") <- "none"
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
attr(astan_1_2_3, "which_X") <- rep(0, nrow(astan_1_2_3))
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
      no_nas = sum(!is.na(dplyr::c_across(dplyr::everything())))
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

# no covariates? should just be the observed means.
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
attr(astan_2_0_3, "adjustment_method_a_multcomp") <- "none"
attr(astan_2_0_3, "adjustment_method_g_multcomp") <- "none"
attr(astan_2_0_3, "adjustment_method_a_fdr") <- "none"
attr(astan_2_0_3, "adjustment_method_g_fdr") <- "none"
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
attr(astan_2_0_3, "which_X") <- rep(0, nrow(astan_2_0_3))
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
    no_nas = sum(!is.na(dplyr::c_across(dplyr::everything())))
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

data_2_1_4 = data.matrix(afilta_2_1_4$e_data[, -1])

covariate_names = colnames(attr(attr(afilta_2_1_4, "group_DF"), "covariates"))[-1]
main_effect_names = attr(attr(afilta_2_1_4, "group_DF"), "main_effects")
group_sampnames <- groupDF_2_1_4[,get_fdata_cname(afilta_2_1_4)]
groupData <- groupDF_2_1_4[group_sampnames %in% colnames(afilta_2_1_4$e_data),]
groupData[,main_effect_names] <- lapply(groupData[main_effect_names], function(x) factor(x, levels=unique(x)))
groupData <- groupData %>% 
  dplyr::left_join(afilta_2_1_4$f_data)

# no continuous covariates
pred_grid_red_2_1_4 = get_pred_grid(groupData, c("Condition", "Level"), "Gender")
pred_grid_full_2_1_4 = get_pred_grid(groupData, c("Condition", "Level"), "Gender", as.formula("~Condition*Level+Gender"))

cobra <- run_twofactor_cpp(
  data = data.matrix(data_2_1_4),
  gpData = groupData[,c("Group", main_effect_names, covariate_names)], 
  Xfull=Xmatrix_2_1_4_full, Xred = Xmatrix_2_1_4,
  pred_grid_full = pred_grid_full_2_1_4, pred_grid_red = pred_grid_red_2_1_4,
  continuous_covar_inds = numeric(0)
)

group_counts_2_1_4 <- data.frame(
  nona_Infection_high = rowSums(!is.na(data_2_1_4[, c(1, 3, 5)])),
  nona_Infection_low = rowSums(!is.na(data_2_1_4[, c(2, 4, 6)])),
  nona_Mock_high = rowSums(!is.na(data_2_1_4[, c(7, 8, 12)])),
  nona_Mock_low = rowSums(!is.na(data_2_1_4[, c(9, 10, 11)]))
) %>%
  dplyr::ungroup()

nona_grps_2_1_4 <- unname(rowSums(group_counts_2_1_4 != 0))
nona_counts_2_1_4 <- unname(rowSums(!is.na(data_2_1_4)))

mean_a_2_1_4 <- data.frame(cobra$lsmeans) 
colnames(mean_a_2_1_4) <- paste0("Mean_", unique(groupData$Group))

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


cmat = rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1),
             c(0, 1, -1, 0), c(0, 1, 0, -1), c(0, 0, 1, -1))

beta_to_mu <- pred_grid_red_2_1_4
beta_to_mu[,4] <- 0
beta_to_mu <- unique(beta_to_mu)
cmat_red <- cmat %*% beta_to_mu

beta_to_mu_full <- pred_grid_full_2_1_4
beta_to_mu_full[,4] <- 0
beta_to_mu_full <- unique(beta_to_mu_full)
cmat_full <- cmat %*% beta_to_mu_full

test_values <- get_test_values_twofactor(data_2_1_4, Xmatrix_2_1_4, Xmatrix_2_1_4_full, cmat_red, cmat_full, cobra$which_X)

diff_denoms <- test_values$diff_denoms
colnames(diff_denoms) <- c("C1", "C2", "C3", "C4", "C5", "C6")

test_stat_2_1_4 <- diffs_2_1_4 %>%
  cbind(group_counts_2_1_4) %>%
  cbind(diff_denoms) %>%
  dplyr::mutate(
    stat_ih_il = (diff_ih_il / C1),
    stat_ih_mh = (diff_ih_mh / C2),
    stat_ih_ml = (diff_ih_ml / C3),
    stat_il_mh = (diff_il_mh / C4),
    stat_il_ml = (diff_il_ml / C5),
    stat_mh_ml = (diff_mh_ml / C6)
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
      df = nona_counts_2_1_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_high = pt(
      q = abs(stat_ih_mh),
      df = nona_counts_2_1_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_low = pt(
      q = abs(stat_ih_ml),
      df = nona_counts_2_1_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_high = pt(
      q = abs(stat_il_mh),
      df = nona_counts_2_1_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_low = pt(
      q = abs(stat_il_ml),
      df = nona_counts_2_1_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_high_vs_Mock_low = pt(
      q = abs(stat_mh_ml),
      df = nona_counts_2_1_4 - test_values$ranks,
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
    grp1 = mean_a_2_1_4[,1],
    grp2 = mean_a_2_1_4[,2],
    pvals = pval_a_2_1_4[, 1],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_high = aflag(
    grp1 = mean_a_2_1_4[,1],
    grp2 = mean_a_2_1_4[,3],
    pvals = pval_a_2_1_4[, 2],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4[,1],
    grp2 = mean_a_2_1_4[,4],
    pvals = pval_a_2_1_4[, 3],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_high = aflag(
    grp1 = mean_a_2_1_4[,2],
    grp2 = mean_a_2_1_4[,3],
    pvals = pval_a_2_1_4[, 4],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4[,2],
    grp2 = mean_a_2_1_4[,4],
    pvals = pval_a_2_1_4[, 5],
    cutoff = 0.05
  ),
  Mock_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_1_4[,3],
    grp2 = mean_a_2_1_4[,4],
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
  Fold_change_Infection_high_vs_Infection_low = diffs_2_1_4$diff_ih_il,
  Fold_change_Infection_high_vs_Mock_high = diffs_2_1_4$diff_ih_mh,
  Fold_change_Infection_high_vs_Mock_low = diffs_2_1_4$diff_ih_ml,
  Fold_change_Infection_low_vs_Mock_high = diffs_2_1_4$diff_il_mh,
  Fold_change_Infection_low_vs_Mock_low = diffs_2_1_4$diff_il_ml,
  Fold_change_Mock_high_vs_Mock_low = diffs_2_1_4$diff_mh_ml,
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
attr(astan_2_1_4, "adjustment_method_a_multcomp") <- "none"
attr(astan_2_1_4, "adjustment_method_g_multcomp") <- "none"
attr(astan_2_1_4, "adjustment_method_a_fdr") <- "none"
attr(astan_2_1_4, "adjustment_method_g_fdr") <- "none"
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
attr(astan_2_1_4, "which_X") <- which_X_2_1_4_a
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
    no_nas = sum(!is.na(dplyr::c_across(dplyr::everything())))
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
data_2_2_4 = data.matrix(afilta_2_2_4$e_data[, -1])

covariate_names = colnames(attr(attr(afilta_2_2_4, "group_DF"), "covariates"))[-1]
main_effect_names = attr(attr(afilta_2_2_4, "group_DF"), "main_effects")
group_sampnames <- groupDF_2_2_4[,get_fdata_cname(afilta_2_2_4)]
groupData <- groupDF_2_2_4[group_sampnames %in% colnames(afilta_2_2_4$e_data),]
groupData <- groupData %>% 
  dplyr::left_join(afilta_2_2_4$f_data)
groupData[,main_effect_names] <- lapply(groupData[main_effect_names], function(x) factor(x, levels=unique(x)))

pred_grid_red_2_2_4 = get_pred_grid(groupData, main_effect_names, covariate_names)
pred_grid_full_2_2_4 = get_pred_grid(groupData, main_effect_names, covariate_names, fspec = as.formula("~Condition*Level+Gender+Age"))

cobra <- run_twofactor_cpp(
  data = data.matrix(data_2_2_4),
  gpData = groupData[,c("Group", main_effect_names, covariate_names)], 
  Xfull=Xmatrix_2_2_4_full, Xred = Xmatrix_2_2_4,
  pred_grid_full = pred_grid_full_2_2_4, pred_grid_red = pred_grid_red_2_2_4,
  continuous_covar_inds = 5
)

group_counts_2_2_4 <- data.frame(
  nona_Infection_high = rowSums(!is.na(data_2_2_4[, c(1, 3, 5)])),
  nona_Infection_low = rowSums(!is.na(data_2_2_4[, c(2, 4, 6)])),
  nona_Mock_high = rowSums(!is.na(data_2_2_4[, c(7, 8, 12)])),
  nona_Mock_low = rowSums(!is.na(data_2_2_4[, c(9, 10, 11)]))
) %>%
  dplyr::ungroup()

nona_grps_2_2_4 <- unname(rowSums(group_counts_2_2_4 != 0))
nona_counts_2_2_4 <- unname(rowSums(!is.na(data_2_2_4)))

mean_a_2_2_4 <- data.frame(cobra$lsmeans)
colnames(mean_a_2_2_4) <- paste0("Mean_", unique(groupData$Group))

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

cmat = rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1),
             c(0, 1, -1, 0), c(0, 1, 0, -1), c(0, 0, 1, -1))

beta_to_mu = pred_grid_red_2_2_4
beta_to_mu[,4:5] <- 0
beta_to_mu <- unique(beta_to_mu)
cmat_red <- cmat %*% beta_to_mu

beta_to_mu_full = pred_grid_full_2_2_4
beta_to_mu_full[,4:5] <- 0
beta_to_mu_full <- unique(beta_to_mu_full)
cmat_full <- cmat %*% beta_to_mu_full

test_values <- get_test_values_twofactor(data_2_2_4, Xmatrix_2_2_4, Xmatrix_2_2_4_full, cmat_red, cmat_full, cobra$which_X)

diff_denoms <- test_values$diff_denoms
colnames(diff_denoms) <- c("C1", "C2", "C3", "C4", "C5", "C6")

test_stat_2_2_4 <- diffs_2_2_4 %>%
  cbind(group_counts_2_2_4) %>%
  cbind(diff_denoms) %>%
  dplyr::mutate(
    stat_ih_il = (diff_ih_il / C1),
    stat_ih_mh = (diff_ih_mh / C2),
    stat_ih_ml = (diff_ih_ml / C3),
    stat_il_mh = (diff_il_mh / C4),
    stat_il_ml = (diff_il_ml / C5),
    stat_mh_ml = (diff_mh_ml / C6)
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
      df = nona_counts_2_2_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_high = pt(
      q = abs(stat_ih_mh),
      df = nona_counts_2_2_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_high_vs_Mock_low = pt(
      q = abs(stat_ih_ml),
      df = nona_counts_2_2_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_high = pt(
      q = abs(stat_il_mh),
      df = nona_counts_2_2_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Infection_low_vs_Mock_low = pt(
      q = abs(stat_il_ml),
      df = nona_counts_2_2_4 - test_values$ranks,
      lower.tail = FALSE
    ) * 2,
    P_value_A_Mock_high_vs_Mock_low = pt(
      q = abs(stat_mh_ml),
      df = nona_counts_2_2_4 - test_values$ranks,
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
    grp1 = mean_a_2_2_4[,1],
    grp2 = mean_a_2_2_4[,2],
    pvals = pval_a_2_2_4[, 1],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_high = aflag(
    grp1 = mean_a_2_2_4[,1],
    grp2 = mean_a_2_2_4[,3],
    pvals = pval_a_2_2_4[, 2],
    cutoff = 0.05
  ),
  Infection_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4[,1],
    grp2 = mean_a_2_2_4[,4],
    pvals = pval_a_2_2_4[, 3],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_high = aflag(
    grp1 = mean_a_2_2_4[,2],
    grp2 = mean_a_2_2_4[,3],
    pvals = pval_a_2_2_4[, 4],
    cutoff = 0.05
  ),
  Infection_low_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4[,2],
    grp2 = mean_a_2_2_4[,4],
    pvals = pval_a_2_2_4[, 5],
    cutoff = 0.05
  ),
  Mock_high_vs_Mock_low = aflag(
    grp1 = mean_a_2_2_4[,3],
    grp2 = mean_a_2_2_4[,4],
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
  Fold_change_Infection_high_vs_Infection_low = diffs_2_2_4$diff_ih_il,
  Fold_change_Infection_high_vs_Mock_high = diffs_2_2_4$diff_ih_mh,
  Fold_change_Infection_high_vs_Mock_low = diffs_2_2_4$diff_ih_ml,
  Fold_change_Infection_low_vs_Mock_high = diffs_2_2_4$diff_il_mh,
  Fold_change_Infection_low_vs_Mock_low = diffs_2_2_4$diff_il_ml,
  Fold_change_Mock_high_vs_Mock_low = diffs_2_2_4$diff_mh_ml,
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
attr(astan_2_2_4, "adjustment_method_a_multcomp") <- "none"
attr(astan_2_2_4, "adjustment_method_g_multcomp") <- "none"
attr(astan_2_2_4, "adjustment_method_a_fdr") <- "none"
attr(astan_2_2_4, "adjustment_method_g_fdr") <- "none"
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
attr(astan_2_2_4, "which_X") <- which_X_2_2_4_a
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
      no_nas = sum(!is.na(dplyr::c_across(dplyr::everything())))
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

obs_mut_1_1_3 <- rowSums(!is.na(gfilta_1_1_3$e_data[, c(2, 4, 5, 10)]))
obs_zom_1_1_3 <- rowSums(!is.na(gfilta_1_1_3$e_data[, c(3, 6:9)]))
obs_hum_1_1_3 <- rowSums(!is.na(gfilta_1_1_3$e_data[, c(11:13)]))
abs_mut_1_1_3 <- rowSums(is.na(gfilta_1_1_3$e_data[, c(2, 4, 5, 10)]))
abs_zom_1_1_3 <- rowSums(is.na(gfilta_1_1_3$e_data[, c(3, 6:9)]))
abs_hum_1_1_3 <- rowSums(is.na(gfilta_1_1_3$e_data[, c(11:13)]))

obs_mut_1_2_3 <- rowSums(!is.na(gfilta_1_2_3$e_data[, c(2, 4, 5, 10)]))
obs_zom_1_2_3 <- rowSums(!is.na(gfilta_1_2_3$e_data[, c(3, 6:9)]))
obs_hum_1_2_3 <- rowSums(!is.na(gfilta_1_2_3$e_data[, c(11:13)]))
abs_mut_1_2_3 <- rowSums(is.na(gfilta_1_2_3$e_data[, c(2, 4, 5, 10)]))
abs_zom_1_2_3 <- rowSums(is.na(gfilta_1_2_3$e_data[, c(3, 6:9)]))
abs_hum_1_2_3 <- rowSums(is.na(gfilta_1_2_3$e_data[, c(11:13)]))

obs_inf_high_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, c(2:4, 7:8)]))
obs_inf_low_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, c(5:6, 9:10)]))
obs_mock_2_0_3 <- rowSums(!is.na(gfilta_2_0_3$e_data[, 11:13]))
abs_inf_high_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, c(2:4, 7:8)]))
abs_inf_low_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, c(5:6, 9:10)]))
abs_mock_2_0_3 <- rowSums(is.na(gfilta_2_0_3$e_data[, 11:13]))

obs_ih_2_1_4 <- rowSums(!is.na(gfilta_2_1_4$e_data[, c(2, 4, 6)]))
obs_il_2_1_4 <- rowSums(!is.na(gfilta_2_1_4$e_data[, c(3, 5, 7)]))
obs_mh_2_1_4 <- rowSums(!is.na(gfilta_2_1_4$e_data[, c(8, 9, 13)]))
obs_ml_2_1_4 <- rowSums(!is.na(gfilta_2_1_4$e_data[, 10:12]))
abs_ih_2_1_4 <- rowSums(is.na(gfilta_2_1_4$e_data[, c(2, 4, 6)]))
abs_il_2_1_4 <- rowSums(is.na(gfilta_2_1_4$e_data[, c(3, 5, 7)]))
abs_mh_2_1_4 <- rowSums(is.na(gfilta_2_1_4$e_data[, c(8, 9, 13)]))
abs_ml_2_1_4 <- rowSums(is.na(gfilta_2_1_4$e_data[, 10:12]))

obs_ih_2_2_4 <- rowSums(!is.na(gfilta_2_2_4$e_data[, c(2, 4, 6)]))
obs_il_2_2_4 <- rowSums(!is.na(gfilta_2_2_4$e_data[, c(3, 5, 7)]))
obs_mh_2_2_4 <- rowSums(!is.na(gfilta_2_2_4$e_data[, c(8, 9, 13)]))
obs_ml_2_2_4 <- rowSums(!is.na(gfilta_2_2_4$e_data[, 10:12]))
abs_ih_2_2_4 <- rowSums(is.na(gfilta_2_2_4$e_data[, c(2, 4, 6)]))
abs_il_2_2_4 <- rowSums(is.na(gfilta_2_2_4$e_data[, c(3, 5, 7)]))
abs_mh_2_2_4 <- rowSums(is.na(gfilta_2_2_4$e_data[, c(8, 9, 13)]))
abs_ml_2_2_4 <- rowSums(is.na(gfilta_2_2_4$e_data[, 10:12]))

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
attr(gstan_1_0_2, "adjustment_method_a_multcomp") <- "none"
attr(gstan_1_0_2, "adjustment_method_g_multcomp") <- "none"
attr(gstan_1_0_2, "adjustment_method_a_fdr") <- "none"
attr(gstan_1_0_2, "adjustment_method_g_fdr") <- "none"
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
attr(gstan_1_0_2, "which_X") <- rep(0, nrow(gstan_1_0_2))
attr(gstan_1_0_2, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_1_0_2, "data_class") <- "pepData"

# main effects: 1; covariates: 1; groups: 3 ---------------

pval_g_1_1_3 <- data.frame(
  P_value_G_mutant_vs_zombie = rep(0, nrow(gfilta_1_1_3$e_data)),
  P_value_G_mutant_vs_human = rep(0, nrow(gfilta_1_1_3$e_data)),
  P_value_G_zombie_vs_human = rep(0, nrow(gfilta_1_1_3$e_data))
)

for (e in 1:nrow(gfilta_1_1_3$e_data)) {

  pval_g_1_1_3[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, 2:10])) ~
        attr(groupDF_1_1_3, "covariates")$Gender[1:9] +
        groupDF_1_1_3$Group[1:9],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  pval_g_1_1_3[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, c(2, 4, 5, 10:13)])) ~
        attr(groupDF_1_1_3, "covariates")$Gender[c(1, 3, 4, 9:12)] +
        groupDF_1_1_3$Group[c(1, 3, 4, 9:12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  pval_g_1_1_3[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_1_3$e_data[e, c(3, 6:9, 11:13)])) ~
        attr(groupDF_1_1_3, "covariates")$Gender[c(2, 5:8, 10:12)] +
        groupDF_1_1_3$Group[c(2, 5:8, 10:12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]

}

flag_g_1_1_3 <- data.frame(
  Flag_G_mutant_vs_zombie = gflag(
    obs1 = obs_mut_1_1_3,
    obs2 = obs_zom_1_1_3,
    abs1 = abs_mut_1_1_3,
    abs2 = abs_zom_1_1_3,
    pvals = pval_g_1_1_3[, 1],
    cutoff = 0.05
  ),
  Flag_G_mutant_vs_human = gflag(
    obs1 = obs_mut_1_1_3,
    obs2 = obs_hum_1_1_3,
    abs1 = abs_mut_1_1_3,
    abs2 = abs_hum_1_1_3,
    pvals = pval_g_1_1_3[, 2],
    cutoff = 0.05
  ),
  Flag_G_zombie_vs_human = gflag(
    obs1 = obs_zom_1_1_3,
    obs2 = obs_hum_1_1_3,
    abs1 = abs_zom_1_1_3,
    abs2 = abs_hum_1_1_3,
    pvals = pval_g_1_1_3[, 3],
    cutoff = 0.05
  )
)

data_1_1_3 <- gfilta_1_1_3$e_data[, -1]

Betas = compute_betas(data_mat = data.matrix(data_1_1_3), Xmatrix = data.matrix(Xmatrix_1_1_3))
group_sampnames <- gfilta_1_1_3$f_data$SampleID
groupData <- groupDF_1_1_3[group_sampnames %in% colnames(gfilta_1_1_3$e_data),]
groupData <- groupData %>% 
  dplyr::left_join(gfilta_1_1_3$f_data)
groupData[,"Condition"] <- lapply(groupData["Condition"], function(x) factor(x, levels=unique(x)))

pred_grid <- get_pred_grid(groupData, "Condition", "Gender")
mean_1_1_3 <- get_lsmeans(data = data_1_1_3, xmatrix = Xmatrix_1_1_3, pred_grid = pred_grid, Betas = Betas)

counts_1_1_3 <- data.frame(
  "Count_mutant" = unname(obs_mut_1_1_3),
  "Count_zombie" = unname(obs_zom_1_1_3),
  "Count_human" = unname(obs_hum_1_1_3)
)

mean_1_1_3[counts_1_1_3 == 0] <- NA

cmat = rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
diffs_1_1_3 <- fold_change_diff(data.matrix(mean_1_1_3), cmat)

gstan_1_1_3 <- data.frame(
  Mass_Tag_ID = gfilta_1_1_3$e_data$Mass_Tag_ID,
  counts_1_1_3,
  mean_1_1_3,
  Fold_change_mutant_vs_zombie = diffs_1_1_3[,1],
  Fold_change_mutant_vs_human = diffs_1_1_3[,2],
  Fold_change_zombie_vs_human = diffs_1_1_3[,3],
  pval_g_1_1_3,
  flag_g_1_1_3,
  row.names = NULL
)

class(gstan_1_1_3) <- c("statRes", "data.frame")

attr(gstan_1_1_3, "group_DF") <- groupDF_1_1_3
attr(gstan_1_1_3, "comparisons") <- c("mutant_vs_zombie",
                                      "mutant_vs_human",
                                      "zombie_vs_human")
attr(gstan_1_1_3, "number_significant") <- data.frame(
  Comparison = c("mutant_vs_zombie",
                 "mutant_vs_human",
                 "zombie_vs_human"),
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
  data_types = NULL
)
attr(gstan_1_1_3, "which_X") <- rep(0, nrow(gstan_1_1_3))
attr(gstan_1_1_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_1_1_3, "data_class") <- "pepData"

# main effects: 1; covariates: 2; groups: 3 ---------------

pval_g_1_2_3 <- data.frame(
  P_value_G_mutant_vs_zombie = rep(0, nrow(gfilta_1_2_3$e_data)),
  P_value_G_mutant_vs_human = rep(0, nrow(gfilta_1_2_3$e_data)),
  P_value_G_zombie_vs_human = rep(0, nrow(gfilta_1_2_3$e_data))
)

for (e in 1:nrow(gfilta_1_2_3$e_data)) {

  pval_g_1_2_3[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_2_3$e_data[e, 2:10])) ~
        attr(groupDF_1_2_3, "covariates")$Gender[1:9] +
        attr(groupDF_1_2_3, "covariates")$Age[1:9] +
        groupDF_1_2_3$Group[1:9],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  pval_g_1_2_3[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_2_3$e_data[e, c(2, 4, 5, 10:13)])) ~
        attr(groupDF_1_2_3, "covariates")$Gender[c(1, 3, 4, 9:12)] +
        attr(groupDF_1_2_3, "covariates")$Age[c(1, 3, 4, 9:12)] +
        groupDF_1_2_3$Group[c(1, 3, 4, 9:12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  pval_g_1_2_3[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_1_2_3$e_data[e, c(3, 6:9, 11:13)])) ~
        attr(groupDF_1_2_3, "covariates")$Gender[c(2, 5:8, 10:12)] +
        attr(groupDF_1_2_3, "covariates")$Age[c(2, 5:8, 10:12)] +
        groupDF_1_2_3$Group[c(2, 5:8, 10:12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]

}

flag_g_1_2_3 <- data.frame(
  Flag_G_mutant_vs_zombie = gflag(
    obs1 = obs_mut_1_2_3,
    obs2 = obs_zom_1_2_3,
    abs1 = abs_mut_1_2_3,
    abs2 = abs_zom_1_2_3,
    pvals = pval_g_1_2_3[, 1],
    cutoff = 0.05
  ),
  Flag_G_mutant_vs_human = gflag(
    obs1 = obs_mut_1_2_3,
    obs2 = obs_hum_1_2_3,
    abs1 = abs_mut_1_2_3,
    abs2 = abs_hum_1_2_3,
    pvals = pval_g_1_2_3[, 2],
    cutoff = 0.05
  ),
  Flag_G_zombie_vs_human = gflag(
    obs1 = obs_zom_1_2_3,
    obs2 = obs_hum_1_2_3,
    abs1 = abs_zom_1_2_3,
    abs2 = abs_hum_1_2_3,
    pvals = pval_g_1_2_3[, 3],
    cutoff = 0.05
  )
)

data_1_2_3 = data.matrix(gfilta_1_2_3$e_data[, -1])
Betas = compute_betas(data_mat = data_1_2_3, Xmatrix = data.matrix(Xmatrix_1_2_3))
group_sampnames <- gfilta_1_2_3$f_data$SampleID
groupData <- groupDF_1_2_3[group_sampnames %in% colnames(gfilta_1_2_3$e_data),]
groupData <- groupData %>% 
  dplyr::left_join(gfilta_1_2_3$f_data)
groupData[,"Condition"] <- lapply(groupData["Condition"], function(x) factor(x, levels=unique(x)))

pred_grid <- get_pred_grid(groupData, "Condition", c("Gender", "Age"))
mean_1_2_3 <- get_lsmeans(data = data_1_2_3, xmatrix = Xmatrix_1_2_3, pred_grid = pred_grid, Betas = Betas, continuous_covar_inds = 5)

counts_1_2_3 <- data.frame(
  "Count_mutant" = unname(obs_mut_1_2_3),
  "Count_zombie" = unname(obs_zom_1_2_3),
  "Count_human" = unname(obs_hum_1_2_3)
)

mean_1_2_3[counts_1_2_3 == 0] <- NA

cmat = rbind(c(1, -1, 0), c(1, 0, -1), c(0, 1, -1))
diffs_1_2_3 <- fold_change_diff(data.matrix(mean_1_2_3), cmat)

gstan_1_2_3 <- data.frame(
  Mass_Tag_ID = gfilta_1_2_3$e_data$Mass_Tag_ID,
  counts_1_2_3,
  mean_1_2_3,
  Fold_change_mutant_vs_zombie = diffs_1_2_3[,1],
  Fold_change_mutant_vs_human = diffs_1_2_3[,2],
  Fold_change_zombie_vs_human = diffs_1_2_3[,3],
  pval_g_1_2_3,
  flag_g_1_2_3,
  row.names = NULL
)

class(gstan_1_2_3) <- c("statRes", "data.frame")

attr(gstan_1_2_3, "group_DF") <- groupDF_1_2_3
attr(gstan_1_2_3, "comparisons") <- c("mutant_vs_zombie",
                                      "mutant_vs_human",
                                      "zombie_vs_human")
attr(gstan_1_2_3, "number_significant") <- data.frame(
  Comparison = c("mutant_vs_zombie",
                 "mutant_vs_human",
                 "zombie_vs_human"),
  Up_total = c(length(which(flag_g_1_2_3[, 1] == 1)),
               length(which(flag_g_1_2_3[, 2] == 1)),
               length(which(flag_g_1_2_3[, 3] == 1))),
  Down_total = c(length(which(flag_g_1_2_3[, 1] == -1)),
                 length(which(flag_g_1_2_3[, 2] == -1)),
                 length(which(flag_g_1_2_3[, 3] == -1))),
  Up_anova = c(0, 0, 0),
  Down_anova = c(0, 0, 0),
  Up_gtest = c(length(which(flag_g_1_2_3[, 1] == 1)),
               length(which(flag_g_1_2_3[, 2] == 1)),
               length(which(flag_g_1_2_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_2_3[, 1] == -1)),
                 length(which(flag_g_1_2_3[, 2] == -1)),
                 length(which(flag_g_1_2_3[, 3] == -1))),
  row.names = NULL
)
attr(gstan_1_2_3, "statistical_test") <- "gtest"
attr(gstan_1_2_3, "adjustment_method_a_multcomp") <- "none"
attr(gstan_1_2_3, "adjustment_method_g_multcomp") <- "none"
attr(gstan_1_2_3, "adjustment_method_a_fdr") <- "none"
attr(gstan_1_2_3, "adjustment_method_g_fdr") <- "none"
attr(gstan_1_2_3, "pval_thresh") <- 0.05
attr(gstan_1_2_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_2_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_2_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_2_3$e_data[, -1])) /
                    prod(dim(gfilta_1_2_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_2_3$f_data)[1],
  data_types = NULL
)
attr(gstan_1_2_3, "which_X") <- rep(0, nrow(gstan_1_2_3))
attr(gstan_1_2_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_1_2_3, "data_class") <- "pepData"

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
attr(gstan_2_0_3, "adjustment_method_a_multcomp") <- "none"
attr(gstan_2_0_3, "adjustment_method_g_multcomp") <- "none"
attr(gstan_2_0_3, "adjustment_method_a_fdr") <- "none"
attr(gstan_2_0_3, "adjustment_method_g_fdr") <- "none"
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
attr(gstan_2_0_3, "which_X") <- rep(0, nrow(gstan_2_0_3))
attr(gstan_2_0_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_2_0_3, "data_class") <- "pepData"

# main effects: 2; covariates: 1; groups: 4 ---------------

pval_g_2_1_4 <- data.frame(
  P_value_G_Infection_high_vs_Infection_low = rep(0, nrow(gfilta_2_1_4$e_data)),
  P_value_G_Infection_high_vs_Mock_high = rep(0, nrow(gfilta_2_1_4$e_data)),
  P_value_G_Infection_high_vs_Mock_low = rep(0, nrow(gfilta_2_1_4$e_data)),
  P_value_G_Infection_low_vs_Mock_high = rep(0, nrow(gfilta_2_1_4$e_data)),
  P_value_G_Infection_low_vs_Mock_low = rep(0, nrow(gfilta_2_1_4$e_data)),
  P_value_G_Mock_high_vs_Mock_low = rep(0, nrow(gfilta_2_1_4$e_data))
)

for (e in 1:nrow(gfilta_2_1_4$e_data)) {

  # Compute the p-value for the Infection_high_vs_Infection_low test.
  pval_g_2_1_4[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, 2:7])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[1:6] +
        groupDF_2_1_4$Group[1:6],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  # Compute the p-value for the Infection_high_vs_Mock_high test.
  pval_g_2_1_4[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, c(2, 4, 6, 8, 9, 13)])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[c(1, 3, 5, 7, 8, 12)] +
        groupDF_2_1_4$Group[c(1, 3, 5, 7, 8, 12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  # Compute the p-value for the Infection_high_vs_Mock_low test.
  pval_g_2_1_4[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, c(2, 4, 6, 10:12)])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[c(1, 3, 5, 9:11)] +
        groupDF_2_1_4$Group[c(1, 3, 5, 9:11)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  # Compute the p-value for the Infection_low_vs_Mock_high test.
  pval_g_2_1_4[e, 4] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, c(3, 5, 7, 8, 9, 13)])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[c(2, 4, 6, 7, 8, 12)] +
        groupDF_2_1_4$Group[c(2, 4, 6, 7, 8, 12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  # Compute the p-value for the Infection_low_vs_Mock_low test.
  pval_g_2_1_4[e, 5] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, c(3, 5, 7, 10:12)])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[c(2, 4, 6, 9:11)] +
        groupDF_2_1_4$Group[c(2, 4, 6, 9:11)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]
  # Compute the p-value for the Mock_high_vs_Mock_low test.
  pval_g_2_1_4[e, 6] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_1_4$e_data[e, 8:13])) ~
        attr(groupDF_2_1_4, "covariates")$Gender[7:12] +
        groupDF_2_1_4$Group[7:12],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[3]]

}

flag_g_2_1_4 <- data.frame(
  Flag_G_Infection_high_vs_Infection_low = gflag(
    obs1 = obs_ih_2_1_4,
    obs2 = obs_il_2_1_4,
    abs1 = abs_ih_2_1_4,
    abs2 = abs_il_2_1_4,
    pvals = pval_g_2_1_4[, 1],
    cutoff = 0.05
  ),
  Flag_G_Infection_high_vs_Mock_high = gflag(
    obs1 = obs_ih_2_1_4,
    obs2 = obs_mh_2_1_4,
    abs1 = abs_ih_2_1_4,
    abs2 = abs_mh_2_1_4,
    pvals = pval_g_2_1_4[, 2],
    cutoff = 0.05
  ),
  Flag_G_Infection_high_vs_Mock_low = gflag(
    obs1 = obs_ih_2_1_4,
    obs2 = obs_ml_2_1_4,
    abs1 = abs_ih_2_1_4,
    abs2 = abs_ml_2_1_4,
    pvals = pval_g_2_1_4[, 3],
    cutoff = 0.05
  ),
  Flag_G_Infection_low_vs_Mock_high = gflag(
    obs1 = obs_il_2_1_4,
    obs2 = obs_mh_2_1_4,
    abs1 = abs_il_2_1_4,
    abs2 = abs_mh_2_1_4,
    pvals = pval_g_2_1_4[, 4],
    cutoff = 0.05
  ),
  Flag_G_Infection_low_vs_Mock_low = gflag(
    obs1 = obs_il_2_1_4,
    obs2 = obs_ml_2_1_4,
    abs1 = abs_il_2_1_4,
    abs2 = abs_ml_2_1_4,
    pvals = pval_g_2_1_4[, 5],
    cutoff = 0.05
  ),
  Flag_G_Mock_high_vs_Mock_low = gflag(
    obs1 = obs_mh_2_1_4,
    obs2 = obs_ml_2_1_4,
    abs1 = abs_mh_2_1_4,
    abs2 = abs_ml_2_1_4,
    pvals = pval_g_2_1_4[, 6],
    cutoff = 0.05
  )
)

gdf = dplyr::left_join(groupDF_2_1_4, attr(groupDF_2_1_4, "covariates")) %>%
  dplyr::select(-SampleID)

dragon <- run_twofactor_cpp(
  data = data.matrix(gfilta_2_1_4$e_data[, -1]),
  gpData = gdf,
  Xfull = Xmatrix_2_1_4_full,
  Xred = Xmatrix_2_1_4,
  pred_grid_full = pred_grid_full_2_1_4,
  pred_grid_red = pred_grid_red_2_1_4,
  continuous_covar_inds = numeric(0)
)

mean_2_1_4 <- data.frame(
  Mean_Infection_high = dragon$adj_group_means[, 1],
  Mean_Infection_low = dragon$adj_group_means[, 2],
  Mean_Mock_high = dragon$adj_group_means[, 3],
  Mean_Mock_low = dragon$adj_group_means[, 4]
)

cmat = rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1),
             c(0, 1, -1, 0), c(0, 1, 0, -1), c(0, 0, 1, -1))

diffs_2_1_4 <- fold_change_diff(data.matrix(mean_2_1_4), cmat)
diffs_2_1_4 <- as.data.frame(diffs_2_1_4) %>%
  `colnames<-`(c("Fold_change_Infection_high_vs_Infection_low", "Fold_change_Infection_high_vs_Mock_high", "Fold_change_Infection_high_vs_Mock_low", "Fold_change_Infection_low_vs_Mock_high", "Fold_change_Infection_low_vs_Mock_low", "Fold_change_Mock_high_vs_Mock_low"))

gstan_2_1_4 <- data.frame(
  Mass_Tag_ID = gfilta_2_1_4$e_data$Mass_Tag_ID,
  Count_Infection_high = unname(obs_ih_2_1_4),
  Count_Infection_low = unname(obs_il_2_1_4),
  Count_Mock_high = unname(obs_mh_2_1_4),
  Count_Mock_low = unname(obs_ml_2_1_4),
  mean_2_1_4,
  diffs_2_1_4,
  pval_g_2_1_4,
  flag_g_2_1_4,
  row.names = NULL
)

class(gstan_2_1_4) <- c("statRes", "data.frame")

attr(gstan_2_1_4, "group_DF") <- groupDF_2_1_4
attr(gstan_2_1_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(gstan_2_1_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(length(which(flag_g_2_1_4[, 1] == 1)),
               length(which(flag_g_2_1_4[, 2] == 1)),
               length(which(flag_g_2_1_4[, 3] == 1)),
               length(which(flag_g_2_1_4[, 4] == 1)),
               length(which(flag_g_2_1_4[, 5] == 1)),
               length(which(flag_g_2_1_4[, 6] == 1))),
  Down_total = c(length(which(flag_g_2_1_4[, 1] == -1)),
                 length(which(flag_g_2_1_4[, 2] == -1)),
                 length(which(flag_g_2_1_4[, 3] == -1)),
                 length(which(flag_g_2_1_4[, 4] == -1)),
                 length(which(flag_g_2_1_4[, 5] == -1)),
                 length(which(flag_g_2_1_4[, 6] == -1))),
  Up_anova = c(0, 0, 0, 0, 0, 0),
  Down_anova = c(0, 0, 0, 0, 0, 0),
  Up_gtest = c(length(which(flag_g_2_1_4[, 1] == 1)),
               length(which(flag_g_2_1_4[, 2] == 1)),
               length(which(flag_g_2_1_4[, 3] == 1)),
               length(which(flag_g_2_1_4[, 4] == 1)),
               length(which(flag_g_2_1_4[, 5] == 1)),
               length(which(flag_g_2_1_4[, 6] == 1))),
  Down_gtest = c(length(which(flag_g_2_1_4[, 1] == -1)),
                 length(which(flag_g_2_1_4[, 2] == -1)),
                 length(which(flag_g_2_1_4[, 3] == -1)),
                 length(which(flag_g_2_1_4[, 4] == -1)),
                 length(which(flag_g_2_1_4[, 5] == -1)),
                 length(which(flag_g_2_1_4[, 6] == -1))),
  row.names = NULL
)
attr(gstan_2_1_4, "statistical_test") <- "gtest"
attr(gstan_2_1_4, "adjustment_method_a_multcomp") <- "none"
attr(gstan_2_1_4, "adjustment_method_g_multcomp") <- "none"
attr(gstan_2_1_4, "adjustment_method_a_fdr") <- "none"
attr(gstan_2_1_4, "adjustment_method_g_fdr") <- "none"
attr(gstan_2_1_4, "pval_thresh") <- 0.05
attr(gstan_2_1_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_2_1_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_2_1_4$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_2_1_4$e_data[, -1])) /
                    prod(dim(gfilta_2_1_4$e_data[, -1]))),
  num_samps = dim(gfilta_2_1_4$f_data)[1],
  data_types = NULL
)
attr(gstan_2_1_4, "which_X") <- which_X_2_1_4_g
attr(gstan_2_1_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_2_1_4, "data_class") <- "pepData"

# main effects: 2; covariates: 2; groups: 4 ---------------

pval_g_2_2_4 <- data.frame(
  P_value_G_Infection_high_vs_Infection_low = rep(0, nrow(gfilta_2_2_4$e_data)),
  P_value_G_Infection_high_vs_Mock_high = rep(0, nrow(gfilta_2_2_4$e_data)),
  P_value_G_Infection_high_vs_Mock_low = rep(0, nrow(gfilta_2_2_4$e_data)),
  P_value_G_Infection_low_vs_Mock_high = rep(0, nrow(gfilta_2_2_4$e_data)),
  P_value_G_Infection_low_vs_Mock_low = rep(0, nrow(gfilta_2_2_4$e_data)),
  P_value_G_Mock_high_vs_Mock_low = rep(0, nrow(gfilta_2_2_4$e_data))
)

for (e in 1:nrow(gfilta_2_2_4$e_data)) {

  # Compute the p-value for the Infection_high_vs_Infection_low test.
  pval_g_2_2_4[e, 1] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, 2:7])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[1:6] +
        attr(groupDF_2_2_4, "covariates")$Age[1:6] +
        groupDF_2_2_4$Group[1:6],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  # Compute the p-value for the Infection_high_vs_Mock_high test.
  pval_g_2_2_4[e, 2] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, c(2, 4, 6, 8, 9, 13)])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[c(1, 3, 5, 7, 8, 12)] +
        attr(groupDF_2_2_4, "covariates")$Age[c(1, 3, 5, 7, 8, 12)] +
        groupDF_2_2_4$Group[c(1, 3, 5, 7, 8, 12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  # Compute the p-value for the Infection_high_vs_Mock_low test.
  pval_g_2_2_4[e, 3] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, c(2, 4, 6, 10:12)])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[c(1, 3, 5, 9:11)] +
        attr(groupDF_2_2_4, "covariates")$Age[c(1, 3, 5, 9:11)] +
        groupDF_2_2_4$Group[c(1, 3, 5, 9:11)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  # Compute the p-value for the Infection_low_vs_Mock_high test.
  pval_g_2_2_4[e, 4] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, c(3, 5, 7, 8, 9, 13)])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[c(2, 4, 6, 7, 8, 12)] +
        attr(groupDF_2_2_4, "covariates")$Age[c(2, 4, 6, 7, 8, 12)] +
        groupDF_2_2_4$Group[c(2, 4, 6, 7, 8, 12)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  # Compute the p-value for the Infection_low_vs_Mock_low test.
  pval_g_2_2_4[e, 5] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, c(3, 5, 7, 10:12)])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[c(2, 4, 6, 9:11)] +
        attr(groupDF_2_2_4, "covariates")$Age[c(2, 4, 6, 9:11)] +
        groupDF_2_2_4$Group[c(2, 4, 6, 9:11)],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]
  # Compute the p-value for the Mock_high_vs_Mock_low test.
  pval_g_2_2_4[e, 6] <- anova(
    glm(
      as.numeric(!is.na(gfilta_2_2_4$e_data[e, 8:13])) ~
        attr(groupDF_2_2_4, "covariates")$Gender[7:12] +
        attr(groupDF_2_2_4, "covariates")$Age[7:12] +
        groupDF_2_2_4$Group[7:12],
      family = binomial
    ),
    test = "Chisq"
  )$`Pr(>Chi)`[[4]]

}

flag_g_2_2_4 <- data.frame(
  Flag_G_Infection_high_vs_Infection_low = gflag(
    obs1 = obs_ih_2_2_4,
    obs2 = obs_il_2_2_4,
    abs1 = abs_ih_2_2_4,
    abs2 = abs_il_2_2_4,
    pvals = pval_g_2_2_4[, 1],
    cutoff = 0.05
  ),
  Flag_G_Infection_high_vs_Mock_high = gflag(
    obs1 = obs_ih_2_2_4,
    obs2 = obs_mh_2_2_4,
    abs1 = abs_ih_2_2_4,
    abs2 = abs_mh_2_2_4,
    pvals = pval_g_2_2_4[, 2],
    cutoff = 0.05
  ),
  Flag_G_Infection_high_vs_Mock_low = gflag(
    obs1 = obs_ih_2_2_4,
    obs2 = obs_ml_2_2_4,
    abs1 = abs_ih_2_2_4,
    abs2 = abs_ml_2_2_4,
    pvals = pval_g_2_2_4[, 3],
    cutoff = 0.05
  ),
  Flag_G_Infection_low_vs_Mock_high = gflag(
    obs1 = obs_il_2_2_4,
    obs2 = obs_mh_2_2_4,
    abs1 = abs_il_2_2_4,
    abs2 = abs_mh_2_2_4,
    pvals = pval_g_2_2_4[, 4],
    cutoff = 0.05
  ),
  Flag_G_Infection_low_vs_Mock_low = gflag(
    obs1 = obs_il_2_2_4,
    obs2 = obs_ml_2_2_4,
    abs1 = abs_il_2_2_4,
    abs2 = abs_ml_2_2_4,
    pvals = pval_g_2_2_4[, 5],
    cutoff = 0.05
  ),
  Flag_G_Mock_high_vs_Mock_low = gflag(
    obs1 = obs_mh_2_2_4,
    obs2 = obs_ml_2_2_4,
    abs1 = abs_mh_2_2_4,
    abs2 = abs_ml_2_2_4,
    pvals = pval_g_2_2_4[, 6],
    cutoff = 0.05
  )
)

gdf = dplyr::left_join(groupDF_2_2_4, attr(groupDF_2_2_4, "covariates")) %>%
  dplyr::select(-SampleID)

mustang <- run_twofactor_cpp(
  data = data.matrix(gfilta_2_2_4$e_data[, -1]),
  gpData = gdf,
  Xfull=Xmatrix_2_2_4_full, Xred = Xmatrix_2_2_4,
  pred_grid_full = pred_grid_full_2_2_4, pred_grid_red = pred_grid_red_2_2_4,
  continuous_covar_inds = 5
)

mean_2_2_4 <- data.frame(
  Mean_Infection_high = mustang$adj_group_means[, 1],
  Mean_Infection_low = mustang$adj_group_means[, 2],
  Mean_Mock_high = mustang$adj_group_means[, 3],
  Mean_Mock_low = mustang$adj_group_means[, 4]
)

cmat = rbind(c(1, -1, 0, 0), c(1, 0, -1, 0), c(1, 0, 0, -1),
             c(0, 1, -1, 0), c(0, 1, 0, -1), c(0, 0, 1, -1))

diffs_2_2_4 <- fold_change_diff(data.matrix(mean_2_2_4), cmat)
diffs_2_2_4 <- as.data.frame(diffs_2_2_4) %>%
  `colnames<-`(c("Fold_change_Infection_high_vs_Infection_low", "Fold_change_Infection_high_vs_Mock_high", "Fold_change_Infection_high_vs_Mock_low", "Fold_change_Infection_low_vs_Mock_high", "Fold_change_Infection_low_vs_Mock_low", "Fold_change_Mock_high_vs_Mock_low"))

gstan_2_2_4 <- data.frame(
  Mass_Tag_ID = gfilta_2_2_4$e_data$Mass_Tag_ID,
  Count_Infection_high = unname(obs_ih_2_2_4),
  Count_Infection_low = unname(obs_il_2_2_4),
  Count_Mock_high = unname(obs_mh_2_2_4),
  Count_Mock_low = unname(obs_ml_2_2_4),
  mean_2_2_4,
  diffs_2_2_4,
  pval_g_2_2_4,
  flag_g_2_2_4,
  row.names = NULL
)

class(gstan_2_2_4) <- c("statRes", "data.frame")

attr(gstan_2_2_4, "group_DF") <- groupDF_2_2_4
attr(gstan_2_2_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(gstan_2_2_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(length(which(flag_g_2_2_4[, 1] == 1)),
               length(which(flag_g_2_2_4[, 2] == 1)),
               length(which(flag_g_2_2_4[, 3] == 1)),
               length(which(flag_g_2_2_4[, 4] == 1)),
               length(which(flag_g_2_2_4[, 5] == 1)),
               length(which(flag_g_2_2_4[, 6] == 1))),
  Down_total = c(length(which(flag_g_2_2_4[, 1] == -1)),
                 length(which(flag_g_2_2_4[, 2] == -1)),
                 length(which(flag_g_2_2_4[, 3] == -1)),
                 length(which(flag_g_2_2_4[, 4] == -1)),
                 length(which(flag_g_2_2_4[, 5] == -1)),
                 length(which(flag_g_2_2_4[, 6] == -1))),
  Up_anova = c(0, 0, 0, 0, 0, 0),
  Down_anova = c(0, 0, 0, 0, 0, 0),
  Up_gtest = c(length(which(flag_g_2_2_4[, 1] == 1)),
               length(which(flag_g_2_2_4[, 2] == 1)),
               length(which(flag_g_2_2_4[, 3] == 1)),
               length(which(flag_g_2_2_4[, 4] == 1)),
               length(which(flag_g_2_2_4[, 5] == 1)),
               length(which(flag_g_2_2_4[, 6] == 1))),
  Down_gtest = c(length(which(flag_g_2_2_4[, 1] == -1)),
                 length(which(flag_g_2_2_4[, 2] == -1)),
                 length(which(flag_g_2_2_4[, 3] == -1)),
                 length(which(flag_g_2_2_4[, 4] == -1)),
                 length(which(flag_g_2_2_4[, 5] == -1)),
                 length(which(flag_g_2_2_4[, 6] == -1))),
  row.names = NULL
)
attr(gstan_2_2_4, "statistical_test") <- "gtest"
attr(gstan_2_2_4, "adjustment_method_a_multcomp") <- "none"
attr(gstan_2_2_4, "adjustment_method_g_multcomp") <- "none"
attr(gstan_2_2_4, "adjustment_method_a_fdr") <- "none"
attr(gstan_2_2_4, "adjustment_method_g_fdr") <- "none"
attr(gstan_2_2_4, "pval_thresh") <- 0.05
attr(gstan_2_2_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_2_2_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_2_2_4$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_2_2_4$e_data[, -1])) /
                    prod(dim(gfilta_2_2_4$e_data[, -1]))),
  num_samps = dim(gfilta_2_2_4$f_data)[1],
  data_types = NULL
)
attr(gstan_2_2_4, "which_X") <- which_X_2_2_4_g
attr(gstan_2_2_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(gstan_2_2_4, "data_class") <- "pepData"


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
attr(cstan_1_0_2, "adjustment_method_a_multcomp") <- "none"
attr(cstan_1_0_2, "adjustment_method_g_multcomp") <- "none"
attr(cstan_1_0_2, "adjustment_method_a_fdr") <- "none"
attr(cstan_1_0_2, "adjustment_method_g_fdr") <- "none"
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
attr(cstan_1_0_2, "which_X") <- 
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, NA, 0, 0, NA, NA, 
    0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, NA, 
    0, 0, NA, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, 0, 
    0, 0, NA, NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, NA, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 
    0, 0)
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

attr(cstan_1_1_3, "group_DF") <- groupDF_1_1_3
attr(cstan_1_1_3, "comparisons") <- c("mutant_vs_zombie",
                                      "mutant_vs_human",
                                      "zombie_vs_human")
attr(cstan_1_1_3, "number_significant") <- data.frame(
  Comparison = c("mutant_vs_zombie",
                 "mutant_vs_human",
                 "zombie_vs_human"),
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
  data_types = NULL
)
attr(cstan_1_1_3, "which_X") <- 
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, NA, 0, 0, NA, 0, 
    0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 
    0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, NA, 
    NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
attr(cstan_1_1_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_1_1_3$e_data$Mass_Tag_ID,
  mutant_vs_zombie = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_mutant_vs_zombie) ~
      cstan_1_1_3$Flag_G_mutant_vs_zombie,
    is.na(cstan_1_1_3$P_value_A_mutant_vs_zombie) ~
      cstan_1_1_3$Flag_G_mutant_vs_zombie,
    (cstan_1_1_3$P_value_A_mutant_vs_zombie > 0.05 &
       cstan_1_1_3$P_value_G_mutant_vs_zombie < 0.05) ~
      cstan_1_1_3$Flag_G_mutant_vs_zombie,
    !is.na(cstan_1_1_3$Flag_A_mutant_vs_zombie) ~
      cstan_1_1_3$Flag_A_mutant_vs_zombie
  ),
  mutant_vs_human = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_mutant_vs_human) ~
      cstan_1_1_3$Flag_G_mutant_vs_human,
    is.na(cstan_1_1_3$P_value_A_mutant_vs_human) ~
      cstan_1_1_3$Flag_G_mutant_vs_human,
    (cstan_1_1_3$P_value_A_mutant_vs_human > 0.05 &
       cstan_1_1_3$P_value_G_mutant_vs_human < 0.05) ~
      cstan_1_1_3$Flag_G_mutant_vs_human,
    !is.na(cstan_1_1_3$Flag_A_mutant_vs_human) ~
      cstan_1_1_3$Flag_A_mutant_vs_human
  ),
  zombie_vs_human = dplyr::case_when(
    is.na(cstan_1_1_3$Flag_A_zombie_vs_human) ~
      cstan_1_1_3$Flag_G_zombie_vs_human,
    is.na(cstan_1_1_3$P_value_A_zombie_vs_human) ~
      cstan_1_1_3$Flag_G_zombie_vs_human,
    (cstan_1_1_3$P_value_A_zombie_vs_human > 0.05 &
       cstan_1_1_3$P_value_G_zombie_vs_human < 0.05) ~
      cstan_1_1_3$Flag_G_zombie_vs_human,
    !is.na(cstan_1_1_3$Flag_A_zombie_vs_human) ~
      cstan_1_1_3$Flag_A_zombie_vs_human
  )
)
attr(cstan_1_1_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(cstan_1_1_3, "data_class") <- "pepData"

# main effects: 1; covariates: 2; groups: 3 ---------------

# Combine the G-test and ANOVA results and place columns in correct order.
cstan_1_2_3 <- dplyr::full_join(gstan_1_2_3[, c(1:4, 11:16)],
                                astan_1_2_3) %>%
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
cstan_1_2_3[is.nan(data.matrix(cstan_1_2_3))] <- NA

class(cstan_1_2_3) <- c("statRes", "data.frame")

attr(cstan_1_2_3, "group_DF") <- groupDF_1_2_3
attr(cstan_1_2_3, "comparisons") <- c("mutant_vs_zombie",
                                      "mutant_vs_human",
                                      "zombie_vs_human")
attr(cstan_1_2_3, "number_significant") <- data.frame(
  Comparison = c("mutant_vs_zombie",
                 "mutant_vs_human",
                 "zombie_vs_human"),
  Up_total = c(sum(length(which(flag_g_1_2_3[, 1] == 1)),
                   length(which(flag_a_1_2_3[, 1] == 1))),
               sum(length(which(flag_g_1_2_3[, 2] == 1)),
                   length(which(flag_a_1_2_3[, 2] == 1))),
               sum(length(which(flag_g_1_2_3[, 3] == 1)),
                   length(which(flag_a_1_2_3[, 3] == 1)))),
  Down_total = c(sum(length(which(flag_g_1_2_3[, 1] == -1)),
                     length(which(flag_a_1_2_3[, 1] == -1))),
                 sum(length(which(flag_g_1_2_3[, 2] == -1)),
                     length(which(flag_a_1_2_3[, 2] == -1))),
                 sum(length(which(flag_g_1_2_3[, 3] == -1)),
                     length(which(flag_a_1_2_3[, 3] == -1)))),
  Up_anova = c(length(which(flag_a_1_2_3[, 1] == 1)),
               length(which(flag_a_1_2_3[, 2] == 1)),
               length(which(flag_a_1_2_3[, 3] == 1))),
  Down_anova = c(length(which(flag_a_1_2_3[, 1] == -1)),
                 length(which(flag_a_1_2_3[, 2] == -1)),
                 length(which(flag_a_1_2_3[, 3] == -1))),
  Up_gtest = c(length(which(flag_g_1_2_3[, 1] == 1)),
               length(which(flag_g_1_2_3[, 2] == 1)),
               length(which(flag_g_1_2_3[, 3] == 1))),
  Down_gtest = c(length(which(flag_g_1_2_3[, 1] == -1)),
                 length(which(flag_g_1_2_3[, 2] == -1)),
                 length(which(flag_g_1_2_3[, 3] == -1))),
  row.names = NULL
)
attr(cstan_1_2_3, "statistical_test") <- "combined"
attr(cstan_1_2_3, "adjustment_method_a_multcomp") <- "none"
attr(cstan_1_2_3, "adjustment_method_g_multcomp") <- "none"
attr(cstan_1_2_3, "adjustment_method_a_fdr") <- "none"
attr(cstan_1_2_3, "adjustment_method_g_fdr") <- "none"
attr(cstan_1_2_3, "pval_thresh") <- 0.05
attr(cstan_1_2_3, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_1_2_3$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_1_2_3$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_1_2_3$e_data[, -1])) /
                    prod(dim(gfilta_1_2_3$e_data[, -1]))),
  num_samps = dim(gfilta_1_2_3$f_data)[1],
  data_types = NULL
)
attr(cstan_1_2_3, "which_X") <- 
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, NA, 0, 0, NA, 0, 
    0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 
    0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, NA, 
    NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
attr(cstan_1_2_3, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_1_2_3$e_data$Mass_Tag_ID,
  mutant_vs_zombie = dplyr::case_when(
    is.na(cstan_1_2_3$Flag_A_mutant_vs_zombie) ~
      cstan_1_2_3$Flag_G_mutant_vs_zombie,
    is.na(cstan_1_2_3$P_value_A_mutant_vs_zombie) ~
      cstan_1_2_3$Flag_G_mutant_vs_zombie,
    (cstan_1_2_3$P_value_A_mutant_vs_zombie > 0.05 &
       cstan_1_2_3$P_value_G_mutant_vs_zombie < 0.05) ~
      cstan_1_2_3$Flag_G_mutant_vs_zombie,
    !is.na(cstan_1_2_3$Flag_A_mutant_vs_zombie) ~
      cstan_1_2_3$Flag_A_mutant_vs_zombie
  ),
  mutant_vs_human = dplyr::case_when(
    is.na(cstan_1_2_3$Flag_A_mutant_vs_human) ~
      cstan_1_2_3$Flag_G_mutant_vs_human,
    is.na(cstan_1_2_3$P_value_A_mutant_vs_human) ~
      cstan_1_2_3$Flag_G_mutant_vs_human,
    (cstan_1_2_3$P_value_A_mutant_vs_human > 0.05 &
       cstan_1_2_3$P_value_G_mutant_vs_human < 0.05) ~
      cstan_1_2_3$Flag_G_mutant_vs_human,
    !is.na(cstan_1_2_3$Flag_A_mutant_vs_human) ~
      cstan_1_2_3$Flag_A_mutant_vs_human
  ),
  zombie_vs_human = dplyr::case_when(
    is.na(cstan_1_2_3$Flag_A_zombie_vs_human) ~
      cstan_1_2_3$Flag_G_zombie_vs_human,
    is.na(cstan_1_2_3$P_value_A_zombie_vs_human) ~
      cstan_1_2_3$Flag_G_zombie_vs_human,
    (cstan_1_2_3$P_value_A_zombie_vs_human > 0.05 &
       cstan_1_2_3$P_value_G_zombie_vs_human < 0.05) ~
      cstan_1_2_3$Flag_G_zombie_vs_human,
    !is.na(cstan_1_2_3$Flag_A_zombie_vs_human) ~
      cstan_1_2_3$Flag_A_zombie_vs_human
  )
)
attr(cstan_1_2_3, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = "Protein",
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(cstan_1_2_3, "data_class") <- "pepData"

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
attr(cstan_2_0_3, "adjustment_method_a_multcomp") <- "none"
attr(cstan_2_0_3, "adjustment_method_g_multcomp") <- "none"
attr(cstan_2_0_3, "adjustment_method_a_fdr") <- "none"
attr(cstan_2_0_3, "adjustment_method_g_fdr") <- "none"
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
attr(cstan_2_0_3, "which_X") <- 
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, NA, NA, 
    0, 0, 0, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 
    0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA, 0, 0, 0, 0, 
    0, 0, NA, NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
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

# main effects: 2; covariates: 1; groups: 4 ---------------

# Combine the G-test and ANOVA results and place columns in correct order.
cstan_2_1_4 <- dplyr::full_join(gstan_2_1_4[, c(1:5, 16:27)],
                                astan_2_1_4) %>%
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
cstan_2_1_4[is.nan(data.matrix(cstan_2_1_4))] <- NA

class(cstan_2_1_4) <- c("statRes", "data.frame")

attr(cstan_2_1_4, "group_DF") <- groupDF_2_1_4
attr(cstan_2_1_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(cstan_2_1_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(sum(length(which(flag_g_2_1_4[, 1] == 1)),
                   length(which(flag_a_2_1_4[, 1] == 1))),
               sum(length(which(flag_g_2_1_4[, 2] == 1)),
                   length(which(flag_a_2_1_4[, 2] == 1))),
               sum(length(which(flag_g_2_1_4[, 3] == 1)),
                   length(which(flag_a_2_1_4[, 3] == 1))),
               sum(length(which(flag_g_2_1_4[, 4] == 1)),
                   length(which(flag_a_2_1_4[, 4] == 1))),
               sum(length(which(flag_g_2_1_4[, 5] == 1)),
                   length(which(flag_a_2_1_4[, 5] == 1))),
               sum(length(which(flag_g_2_1_4[, 6] == 1)),
                   length(which(flag_a_2_1_4[, 6] == 1)))),
  Down_total = c(sum(length(which(flag_g_2_1_4[, 1] == -1)),
                     length(which(flag_a_2_1_4[, 1] == -1))),
                 sum(length(which(flag_g_2_1_4[, 2] == -1)),
                     length(which(flag_a_2_1_4[, 2] == -1))),
                 sum(length(which(flag_g_2_1_4[, 3] == -1)),
                     length(which(flag_a_2_1_4[, 3] == -1))),
                 sum(length(which(flag_g_2_1_4[, 4] == -1)),
                     length(which(flag_a_2_1_4[, 4] == -1))),
                 sum(length(which(flag_g_2_1_4[, 5] == -1)),
                     length(which(flag_a_2_1_4[, 5] == -1))),
                 sum(length(which(flag_g_2_1_4[, 6] == -1)),
                     length(which(flag_a_2_1_4[, 6] == -1)))),
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
  Up_gtest = c(length(which(flag_g_2_1_4[, 1] == 1)),
               length(which(flag_g_2_1_4[, 2] == 1)),
               length(which(flag_g_2_1_4[, 3] == 1)),
               length(which(flag_g_2_1_4[, 4] == 1)),
               length(which(flag_g_2_1_4[, 5] == 1)),
               length(which(flag_g_2_1_4[, 6] == 1))),
  Down_gtest = c(length(which(flag_g_2_1_4[, 1] == -1)),
                 length(which(flag_g_2_1_4[, 2] == -1)),
                 length(which(flag_g_2_1_4[, 3] == -1)),
                 length(which(flag_g_2_1_4[, 4] == -1)),
                 length(which(flag_g_2_1_4[, 5] == -1)),
                 length(which(flag_g_2_1_4[, 6] == -1))),
  row.names = NULL
)
attr(cstan_2_1_4, "statistical_test") <- "combined"
attr(cstan_2_1_4, "adjustment_method_a_multcomp") <- "none"
attr(cstan_2_1_4, "adjustment_method_g_multcomp") <- "none"
attr(cstan_2_1_4, "adjustment_method_a_fdr") <- "none"
attr(cstan_2_1_4, "adjustment_method_g_fdr") <- "none"
attr(cstan_2_1_4, "pval_thresh") <- 0.05
attr(cstan_2_1_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_2_1_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_2_1_4$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_2_1_4$e_data[, -1])) /
                    prod(dim(gfilta_2_1_4$e_data[, -1]))),
  num_samps = dim(gfilta_2_1_4$f_data)[1],
  data_types = NULL
)
attr(cstan_2_1_4, "which_X") <- 
  c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 
    1, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 1, 
    NA, 0, 0, 0, 1, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, 0, 0, NA, NA, 
    NA, 0, 0, 0, NA, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
attr(cstan_2_1_4, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_2_1_4$e_data$Mass_Tag_ID,
  Infection_high_vs_Infection_low = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Infection_low) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Infection_low,
    is.na(cstan_2_1_4$P_value_A_Infection_high_vs_Infection_low) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Infection_low,
    (cstan_2_1_4$P_value_A_Infection_high_vs_Infection_low > 0.05 &
       cstan_2_1_4$P_value_G_Infection_high_vs_Infection_low < 0.05) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Infection_low,
    !is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Infection_low) ~
      cstan_2_1_4$Flag_A_Infection_high_vs_Infection_low
  ),
  Infection_high_vs_Mock_high = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Mock_high) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_high,
    is.na(cstan_2_1_4$P_value_A_Infection_high_vs_Mock_high) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_high,
    (cstan_2_1_4$P_value_A_Infection_high_vs_Mock_high > 0.05 &
       cstan_2_1_4$P_value_G_Infection_high_vs_Mock_high < 0.05) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_high,
    !is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Mock_high) ~
      cstan_2_1_4$Flag_A_Infection_high_vs_Mock_high
  ),
  Infection_high_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_low,
    is.na(cstan_2_1_4$P_value_A_Infection_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_low,
    (cstan_2_1_4$P_value_A_Infection_high_vs_Mock_low > 0.05 &
       cstan_2_1_4$P_value_G_Infection_high_vs_Mock_low < 0.05) ~
      cstan_2_1_4$Flag_G_Infection_high_vs_Mock_low,
    !is.na(cstan_2_1_4$Flag_A_Infection_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_A_Infection_high_vs_Mock_low
  ),
  Infection_low_vs_Mock_high = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Infection_low_vs_Mock_high) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_high,
    is.na(cstan_2_1_4$P_value_A_Infection_low_vs_Mock_high) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_high,
    (cstan_2_1_4$P_value_A_Infection_low_vs_Mock_high > 0.05 &
       cstan_2_1_4$P_value_G_Infection_low_vs_Mock_high < 0.05) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_high,
    !is.na(cstan_2_1_4$Flag_A_Infection_low_vs_Mock_high) ~
      cstan_2_1_4$Flag_A_Infection_low_vs_Mock_high
  ),
  Infection_low_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Infection_low_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_low,
    is.na(cstan_2_1_4$P_value_A_Infection_low_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_low,
    (cstan_2_1_4$P_value_A_Infection_low_vs_Mock_low > 0.05 &
       cstan_2_1_4$P_value_G_Infection_low_vs_Mock_low < 0.05) ~
      cstan_2_1_4$Flag_G_Infection_low_vs_Mock_low,
    !is.na(cstan_2_1_4$Flag_A_Infection_low_vs_Mock_low) ~
      cstan_2_1_4$Flag_A_Infection_low_vs_Mock_low
  ),
  Mock_high_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_1_4$Flag_A_Mock_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Mock_high_vs_Mock_low,
    is.na(cstan_2_1_4$P_value_A_Mock_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_G_Mock_high_vs_Mock_low,
    (cstan_2_1_4$P_value_A_Mock_high_vs_Mock_low > 0.05 &
       cstan_2_1_4$P_value_G_Mock_high_vs_Mock_low < 0.05) ~
      cstan_2_1_4$Flag_G_Mock_high_vs_Mock_low,
    !is.na(cstan_2_1_4$Flag_A_Mock_high_vs_Mock_low) ~
      cstan_2_1_4$Flag_A_Mock_high_vs_Mock_low
  )
)
attr(cstan_2_1_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(cstan_2_1_4, "data_class") <- "pepData"

# main effects: 2; covariates: 2; groups: 4 ---------------

# Combine the G-test and ANOVA results and place columns in correct order.
cstan_2_2_4 <- dplyr::full_join(gstan_2_2_4[, c(1:5, 16:27)],
                                astan_2_2_4) %>%
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
cstan_2_2_4[is.nan(data.matrix(cstan_2_2_4))] <- NA

class(cstan_2_2_4) <- c("statRes", "data.frame")

attr(cstan_2_2_4, "group_DF") <- groupDF_2_2_4
attr(cstan_2_2_4, "comparisons") <- c("Infection_high_vs_Infection_low",
                                      "Infection_high_vs_Mock_high",
                                      "Infection_high_vs_Mock_low",
                                      "Infection_low_vs_Mock_high",
                                      "Infection_low_vs_Mock_low",
                                      "Mock_high_vs_Mock_low")
attr(cstan_2_2_4, "number_significant") <- data.frame(
  Comparison = c("Infection_high_vs_Infection_low",
                 "Infection_high_vs_Mock_high",
                 "Infection_high_vs_Mock_low",
                 "Infection_low_vs_Mock_high",
                 "Infection_low_vs_Mock_low",
                 "Mock_high_vs_Mock_low"),
  Up_total = c(sum(length(which(flag_g_2_2_4[, 1] == 1)),
                   length(which(flag_a_2_2_4[, 1] == 1))),
               sum(length(which(flag_g_2_2_4[, 2] == 1)),
                   length(which(flag_a_2_2_4[, 2] == 1))),
               sum(length(which(flag_g_2_2_4[, 3] == 1)),
                   length(which(flag_a_2_2_4[, 3] == 1))),
               sum(length(which(flag_g_2_2_4[, 4] == 1)),
                   length(which(flag_a_2_2_4[, 4] == 1))),
               sum(length(which(flag_g_2_2_4[, 5] == 1)),
                   length(which(flag_a_2_2_4[, 5] == 1))),
               sum(length(which(flag_g_2_2_4[, 6] == 1)),
                   length(which(flag_a_2_2_4[, 6] == 1)))),
  Down_total = c(sum(length(which(flag_g_2_2_4[, 1] == -1)),
                     length(which(flag_a_2_2_4[, 1] == -1))),
                 sum(length(which(flag_g_2_2_4[, 2] == -1)),
                     length(which(flag_a_2_2_4[, 2] == -1))),
                 sum(length(which(flag_g_2_2_4[, 3] == -1)),
                     length(which(flag_a_2_2_4[, 3] == -1))),
                 sum(length(which(flag_g_2_2_4[, 4] == -1)),
                     length(which(flag_a_2_2_4[, 4] == -1))),
                 sum(length(which(flag_g_2_2_4[, 5] == -1)),
                     length(which(flag_a_2_2_4[, 5] == -1))),
                 sum(length(which(flag_g_2_2_4[, 6] == -1)),
                     length(which(flag_a_2_2_4[, 6] == -1)))),
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
  Up_gtest = c(length(which(flag_g_2_2_4[, 1] == 1)),
               length(which(flag_g_2_2_4[, 2] == 1)),
               length(which(flag_g_2_2_4[, 3] == 1)),
               length(which(flag_g_2_2_4[, 4] == 1)),
               length(which(flag_g_2_2_4[, 5] == 1)),
               length(which(flag_g_2_2_4[, 6] == 1))),
  Down_gtest = c(length(which(flag_g_2_2_4[, 1] == -1)),
                 length(which(flag_g_2_2_4[, 2] == -1)),
                 length(which(flag_g_2_2_4[, 3] == -1)),
                 length(which(flag_g_2_2_4[, 4] == -1)),
                 length(which(flag_g_2_2_4[, 5] == -1)),
                 length(which(flag_g_2_2_4[, 6] == -1))),
  row.names = NULL
)
attr(cstan_2_2_4, "statistical_test") <- "combined"
attr(cstan_2_2_4, "adjustment_method_a_multcomp") <- "none"
attr(cstan_2_2_4, "adjustment_method_g_multcomp") <- "none"
attr(cstan_2_2_4, "adjustment_method_a_fdr") <- "none"
attr(cstan_2_2_4, "adjustment_method_g_fdr") <- "none"
attr(cstan_2_2_4, "pval_thresh") <- 0.05
attr(cstan_2_2_4, "data_info") <- list(
  data_scale_orig = "abundance",
  data_scale = "log",
  norm_info = list(is_normalized = FALSE),
  num_edata = length(unique(gfilta_2_2_4$e_data$Mass_Tag_ID)),
  num_miss_obs = sum(is.na(gfilta_2_2_4$e_data[, -1])),
  prop_missing = (sum(is.na(gfilta_2_2_4$e_data[, -1])) /
                    prod(dim(gfilta_2_2_4$e_data[, -1]))),
  num_samps = dim(gfilta_2_2_4$f_data)[1],
  data_types = NULL
)
attr(cstan_2_2_4, "which_X") <- 
  c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 
  1, 0, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, NA, 0, 0, NA, 0, 0, 1, 
  NA, 0, 0, 0, 1, 0, 0, 0, 0, 0, NA, NA, 0, 0, 0, 0, 0, NA, NA, 
  NA, 0, 0, 0, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
attr(cstan_2_2_4, "bpFlags") <- data.frame(
  Mass_Tag_ID = gfilta_2_2_4$e_data$Mass_Tag_ID,
  Infection_high_vs_Infection_low = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Infection_low) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Infection_low,
    is.na(cstan_2_2_4$P_value_A_Infection_high_vs_Infection_low) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Infection_low,
    (cstan_2_2_4$P_value_A_Infection_high_vs_Infection_low > 0.05 &
       cstan_2_2_4$P_value_G_Infection_high_vs_Infection_low < 0.05) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Infection_low,
    !is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Infection_low) ~
      cstan_2_2_4$Flag_A_Infection_high_vs_Infection_low
  ),
  Infection_high_vs_Mock_high = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Mock_high) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_high,
    is.na(cstan_2_2_4$P_value_A_Infection_high_vs_Mock_high) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_high,
    (cstan_2_2_4$P_value_A_Infection_high_vs_Mock_high > 0.05 &
       cstan_2_2_4$P_value_G_Infection_high_vs_Mock_high < 0.05) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_high,
    !is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Mock_high) ~
      cstan_2_2_4$Flag_A_Infection_high_vs_Mock_high
  ),
  Infection_high_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_low,
    is.na(cstan_2_2_4$P_value_A_Infection_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_low,
    (cstan_2_2_4$P_value_A_Infection_high_vs_Mock_low > 0.05 &
       cstan_2_2_4$P_value_G_Infection_high_vs_Mock_low < 0.05) ~
      cstan_2_2_4$Flag_G_Infection_high_vs_Mock_low,
    !is.na(cstan_2_2_4$Flag_A_Infection_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_A_Infection_high_vs_Mock_low
  ),
  Infection_low_vs_Mock_high = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Infection_low_vs_Mock_high) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_high,
    is.na(cstan_2_2_4$P_value_A_Infection_low_vs_Mock_high) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_high,
    (cstan_2_2_4$P_value_A_Infection_low_vs_Mock_high > 0.05 &
       cstan_2_2_4$P_value_G_Infection_low_vs_Mock_high < 0.05) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_high,
    !is.na(cstan_2_2_4$Flag_A_Infection_low_vs_Mock_high) ~
      cstan_2_2_4$Flag_A_Infection_low_vs_Mock_high
  ),
  Infection_low_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Infection_low_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_low,
    is.na(cstan_2_2_4$P_value_A_Infection_low_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_low,
    (cstan_2_2_4$P_value_A_Infection_low_vs_Mock_low > 0.05 &
       cstan_2_2_4$P_value_G_Infection_low_vs_Mock_low < 0.05) ~
      cstan_2_2_4$Flag_G_Infection_low_vs_Mock_low,
    !is.na(cstan_2_2_4$Flag_A_Infection_low_vs_Mock_low) ~
      cstan_2_2_4$Flag_A_Infection_low_vs_Mock_low
  ),
  Mock_high_vs_Mock_low = dplyr::case_when(
    is.na(cstan_2_2_4$Flag_A_Mock_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Mock_high_vs_Mock_low,
    is.na(cstan_2_2_4$P_value_A_Mock_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_G_Mock_high_vs_Mock_low,
    (cstan_2_2_4$P_value_A_Mock_high_vs_Mock_low > 0.05 &
       cstan_2_2_4$P_value_G_Mock_high_vs_Mock_low < 0.05) ~
      cstan_2_2_4$Flag_G_Mock_high_vs_Mock_low,
    !is.na(cstan_2_2_4$Flag_A_Mock_high_vs_Mock_low) ~
      cstan_2_2_4$Flag_A_Mock_high_vs_Mock_low
  )
)
attr(cstan_2_2_4, "cnames") <- list(
  edata_cname = "Mass_Tag_ID",
  emeta_cname = NULL,
  fdata_cname = "SampleID",
  techrep_cname = NULL
)
attr(cstan_2_2_4, "data_class") <- "pepData"

# Save standards for IMD-ANOVA tests -------------------------------------------
# save() the following objects in "inst/testdata/standards_imd_anova.RData"
# astan_1_0_2,
# gstan_1_0_2,
# cstan_1_0_2,

# astan_1_1_3,
# gstan_1_1_3,
# cstan_1_1_3,

# astan_1_2_3,
# gstan_1_2_3,
# cstan_1_2_3,

# astan_2_0_3,
# gstan_2_0_3,
# cstan_2_0_3,

# astan_2_1_4,
# gstan_2_1_4,
# cstan_2_1_4,

# astan_2_2_4,
# gstan_2_2_4,
# cstan_2_2_4,

# tukey_pval_1_1_3,
# tukey_pval_1_2_3,
# tukey_pval_2_0_3,
# tukey_pval_2_1_4,
# tukey_pval_2_2_4,

# dunnett_1_1_3,
# dunnett_1_2_3,
# dunnett_2_0_3,
# dunnett_2_1_4,
# dunnett_2_2_4

# DO NOT PUT CODE HERE TO SAVE THE ABOVE OBJECTS.  It may create problems with CRAN submission.