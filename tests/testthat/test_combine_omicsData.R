source(system.file('testdata', 'load_data.R', package = 'pmartR'), local = TRUE)

## Generate good lipid data ##
good_lipid_data <- edata_transform(ldata, "log2")

good_lipid_data <- normalize_global(good_lipid_data, "all", "median", apply_norm = TRUE)

fake_cov <- c(rep("A", 5), rep("B", 6))
fake_cov2 <- c(rep(LETTERS[1:4], 2), rep(LETTERS[5], 3))
good_lipid_data$f_data["cov"] <- fake_cov
good_lipid_data$f_data['cov_2'] <- fake_cov2

good_lipid_data <- group_designation(good_lipid_data, "Condition", covariates = "cov")
good_lipid_data <- applyFilt(molecule_filter(good_lipid_data), good_lipid_data, min_num = 2)

## Generate second good lipid data ##

good_lipid_data2 <- good_lipid_data
# Some fake edata ID's to make it unique
good_lipid_data2$e_data[, get_edata_cname(good_lipid_data2)] <- paste0(
  "good_lipid_data2_", good_lipid_data2$e_data[, get_edata_cname(good_lipid_data2)]
  )
good_lipid_data2$e_meta[, get_edata_cname(good_lipid_data2)] <- paste0(
  "good_lipid_data2_", good_lipid_data2$e_meta[, get_edata_cname(good_lipid_data2)]
  )

good_lipid_data2 <- applyFilt(cv_filter(good_lipid_data2), good_lipid_data2, cv_thresh = 60)

## Generate good metab data ##
good_metab_data <- edata_transform(mdata, "log2")
good_metab_data <- normalize_global(good_metab_data, "all", "median", apply_norm = TRUE)

fake_cov <- c(rep("A", 6), rep("B", 6))
fake_cov2 <- c(rep(LETTERS[1:4], 2), rep(LETTERS[5], 4))
good_metab_data$f_data["cov"] <- fake_cov
good_metab_data$f_data['cov_2'] <- fake_cov2

good_metab_data <- group_designation(good_metab_data, "Condition", covariates = "cov")

## Generate second good metab data ##

good_metab_data2 <- good_metab_data
# Some fake edata ID's to make it unique
good_metab_data2$e_data[, get_edata_cname(good_metab_data2)] <- paste0(
  "good_metab_data2_", good_metab_data2$e_data[, get_edata_cname(good_metab_data2)]
)
good_metab_data2$e_meta[, get_edata_cname(good_metab_data2)] <- paste0(
  "good_metab_data2_", good_metab_data2$e_meta[, get_edata_cname(good_metab_data2)]
)

## Generate data to sometimes throw errors due to differences to good_lipid_data ##

# object used to check if objects are have the same group_designation/normalization
only_trans_lipid_data <- edata_transform(ldata, "log2")

# objects to be used with extra covariates -- errors on retain_groups = T
lipid_diff_cov <- good_lipid_data
lipid_diff_cov <- group_designation(lipid_diff_cov, "Condition", covariates = c("cov", "cov_2"))

# object to be used with extra covariates
lipid_diff_cov2 <- lipid_diff_cov
lipid_diff_cov2$e_data[, get_edata_cname(lipid_diff_cov2)] <- paste0("lipid_diff_cov2_", lipid_diff_cov2$e_data[, get_edata_cname(lipid_diff_cov2)])
lipid_diff_cov2$e_meta[, get_edata_cname(lipid_diff_cov2)] <- paste0("lipid_diff_cov2_", lipid_diff_cov2$e_meta[, get_edata_cname(lipid_diff_cov2)])

# object to be used with extra main effects
lipid_extra_effects <- good_lipid_data
lipid_extra_effects <- group_designation(lipid_extra_effects, c("Condition", "cov"), covariates = "cov_2")

suppressWarnings({
  combn1 <- combine_omicsData(good_lipid_data, good_lipid_data2)
  combn2 <- combine_omicsData(good_lipid_data, good_lipid_data2, retain_groups = TRUE)
  combn3 <- combine_omicsData(good_lipid_data, good_lipid_data2, retain_groups = F, retain_filters = TRUE)
  combn4 <- combine_omicsData(good_lipid_data, good_lipid_data2, retain_groups = TRUE, retain_filters = TRUE)
  combn5 <- combine_omicsData(good_lipid_data, lipid_diff_cov)
  combn6 <- combine_omicsData(good_lipid_data, lipid_diff_cov, retain_filters = TRUE)
  combn7 <- combine_omicsData(good_lipid_data, lipid_extra_effects)
  combn8 <- combine_omicsData(good_lipid_data2, lipid_diff_cov, retain_filters = TRUE)
  combn9 <- combine_omicsData(good_metab_data, good_metab_data2, retain_groups = TRUE, retain_filters = TRUE)
})

test_that("bad class errors", {
  suppressWarnings({
    expect_error(combine_omicsData(5, good_lipid_data2, retain_groups = TRUE), 
                 regexp = "Objects must be of the same class")
    expect_error(combine_omicsData(5, 5, retain_groups = TRUE), 
                 regexp = "Currently only support lipidData or metabData")
  })
})

test_that("pipeline errors", {
  suppressWarnings({
    expect_error(combine_omicsData(
      edata_transform(good_lipid_data, "log10"),
      good_lipid_data2
    ), regexp = "Objects must be on the same scale")
    expect_error(combine_omicsData(good_lipid_data, only_trans_lipid_data),
      regexp = "Both objects must have the same normalization status"
    )
    expect_error(combine_omicsData(good_lipid_data, 
                                   normalize_global(only_trans_lipid_data, 
                                                    "all", "median", 
                                                    apply_norm = TRUE), 
                                   retain_groups = TRUE),
      regexp = "Both objects must be grouped."
    )
  })
})

test_that("sample errors", {
  suppressWarnings({
    expect_error(combine_omicsData(
      applyFilt(custom_filter(good_lipid_data, f_data_remove = "Mock2"), good_lipid_data),
      good_lipid_data2
    ), regexp = "Number of samples must be the same in both objects")
  })
})

test_that("bad group/covariate structures throw an error", {
  suppressWarnings({
    expect_error(combine_omicsData(lipid_diff_cov, good_lipid_data2, retain_groups = TRUE), 
                 regexp = "covariate structure")
    expect_error(combine_omicsData(good_lipid_data2, lipid_diff_cov2, retain_groups = TRUE), 
                 regexp = "covariate structure")
    expect_error(combine_omicsData(good_lipid_data, lipid_extra_effects, retain_groups = TRUE), 
                 regexp = "main effect")
    expect_error(combine_omicsData(lipid_diff_cov, lipid_extra_effects, retain_groups = TRUE), 
                 regexp = "main effect")
  })
})

test_that("warnings thrown for duplicate e_data/e_meta identifiers", {
  expect_warning(combine_omicsData(good_lipid_data, lipid_extra_effects), regexp = "Duplicate molecule identifiers")
  expect_warning(combine_omicsData(good_lipid_data, good_lipid_data2), regexp = "molecule identifiers")
})

test_that("attributes correctly stored", {
  # drop filters and grouping info
  expect_true(all(
    is.null(attr(combn1, "group_DF")),
    length(attr(combn1, "filters")) == 0
  ))

  expect_true(all(
    is.null(attr(combn7, "group_DF")),
    length(attr(combn7, "filters")) == 0
  ))

  expect_true(all(
    is.null(attr(combn5, "group_DF")),
    length(attr(combn5, "filters")) == 0
  ))

  # no filters, keep groups
  expect_true(all(
    !is.null(attr(combn2, "group_DF")),
    length(attr(combn2, "filters")) == 0
  ))

  expect_true(all(
    !is.null(attr(combn2, "group_DF")),
    length(attr(combn2, "filters")) == 0
  ))

  ### no groups, keep filters
  ftypes <- attr(combn3, "filters") %>%
    lapply(function(x) x$type)

  expect_true(all(
    is.null(attr(combn3, "group_DF")),
    length(ftypes) == 3,
    all(ftypes == c("moleculeFilt", "moleculeFilt", "cvFilt")) ### Is this intended Dan?
  ))

  ftypes <- attr(combn6, "filters") %>%
    lapply(function(x) x$type)

  expect_true(all(
    is.null(attr(combn6, "group_DF")),
    length(ftypes) == 2,
    all(ftypes == c("moleculeFilt", "moleculeFilt"))
  ))

  ftypes <- attr(combn8, "filters") %>%
    lapply(function(x) x$type)

  expect_true(all(
    is.null(attr(combn6, "group_DF")),
    length(ftypes) == 3,
    all(ftypes == c("moleculeFilt", "cvFilt", "moleculeFilt"))
  ))

  # keep both filters and groups
  ftypes <- attr(combn4, "filters") %>%
    lapply(function(x) x$type)

  expect_true(all(
    !is.null(attr(combn4, "group_DF")),
    length(ftypes) == 3,
    all(ftypes == c("moleculeFilt", "moleculeFilt", "cvFilt"))
  ))
  
  # keep both filters and groups
  ftypes <- attr(combn9, "filters") %>%
    lapply(function(x) x$type)
  
  expect_true(all(
    !is.null(attr(combn9, "group_DF")),
    length(ftypes) == 0)
  )
})
