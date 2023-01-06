source(system.file('testdata', 'load_data.R', package = 'pmartR'), local = TRUE)

obj1 <- edata_transform(ldata, "log2")
obj12 <- obj1
obj1 <- normalize_global(obj1, "all", "median", apply_norm = TRUE)

fake_cov <- c(rep("A", 5), rep("B", 6))
fake_cov2 <- c(rep(LETTERS[1:4],2), rep(LETTERS[5], 3))
obj1$f_data["cov"] <- fake_cov
obj1$f_data['cov_2'] <- fake_cov2

# object to be used with extra covariates
obj3 <- obj1

# object to be used with extra main effects
obj5 <- obj1

obj1 <- group_designation(obj1, "Condition", covariates = "cov")
obj2 <- obj1

obj3 <- group_designation(obj3, "Condition", covariates = c("cov", "cov_2"))
obj4 <- obj3

obj5 <- group_designation(obj5, c("Condition", "cov"), covariates = "cov_2")

# Some fake edata ID's to make it unique
obj2$e_data[,get_edata_cname(obj2)] <- paste0("obj2_", obj2$e_data[,get_edata_cname(obj2)])
obj2$e_meta[,get_edata_cname(obj2)] <- paste0("obj2_", obj2$e_meta[,get_edata_cname(obj2)])

obj4$e_data[,get_edata_cname(obj4)] <- paste0("obj4_", obj4$e_data[,get_edata_cname(obj4)])
obj4$e_meta[,get_edata_cname(obj4)] <- paste0("obj4_", obj4$e_meta[,get_edata_cname(obj4)])

obj1 <- applyFilt(molecule_filter(obj1),obj1, min_num = 2)
obj2 <- applyFilt(cv_filter(obj2),obj2, cv_thresh = 60)

suppressWarnings({
  combn1 <- combine_lipidData(obj1, obj2)
  combn2 <- combine_lipidData(obj1, obj2, retain_groups = TRUE)
  combn3 <- combine_lipidData(obj1, obj2, retain_groups = F, retain_filters = TRUE)
  combn4 <- combine_lipidData(obj1, obj2, retain_groups = TRUE, retain_filters = TRUE)
  combn5 <- combine_lipidData(obj1, obj3)
  combn6 <- combine_lipidData(obj1, obj3, retain_filters = TRUE)
  combn7 <- combine_lipidData(obj1, obj5)
  combn8 <- combine_lipidData(obj2, obj3, retain_filters = TRUE)
})

test_that("bad class errors", {
  suppressWarnings({
    expect_error(combine_lipidData(5, obj2, retain_groups = TRUE), regexp = "Objects must be of the same class")
    expect_error(combine_lipidData(5, 5, retain_groups = TRUE), regexp = "Currently only support lipidData")
  })
})

test_that("pipeline errors", {
  suppressWarnings({
    expect_error(combine_lipidData(edata_transform(obj1, "log10"), 
                                   obj2), regexp = "Objects must be on the same scale")
    expect_error(combine_lipidData(obj1, obj12), 
                 regexp = "Both objects must have the same normalization status")
    expect_error(combine_lipidData(obj1, normalize_global(obj12, "all", "median", apply_norm = TRUE), retain_groups = TRUE), 
                 regexp = "Both objects must be grouped.")
  })
})

test_that("sample errors", {
  suppressWarnings({
    expect_error(combine_lipidData(applyFilt(custom_filter(obj1, f_data_remove = "Mock2"), obj1), 
                                   obj2), regexp = "Number of samples must be the same in both objects")
  })
})

test_that("bad group/covariate structures throw an error", {
  suppressWarnings({
    expect_error(combine_lipidData(obj3, obj2, retain_groups = TRUE), regexp = "covariate structure")
    expect_error(combine_lipidData(obj2, obj4, retain_groups = TRUE), regexp = "covariate structure")
    expect_error(combine_lipidData(obj1, obj5, retain_groups = TRUE), regexp = "main effect")
    expect_error(combine_lipidData(obj3, obj5, retain_groups = TRUE), regexp = "main effect")
  })
})

test_that("warnings thrown for duplicate e_data/e_meta identifiers", {
  expect_warning(combine_lipidData(obj1, obj5), regexp = "Duplicate molecule identifiers")
  expect_warning(combine_lipidData(obj1, obj2), regexp = "molecule identifiers")
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
    length(ftypes) == 2,
    all(ftypes == c("moleculeFilt", "cvFilt"))
  ))
  
  ftypes <- attr(combn6, "filters") %>% 
    lapply(function(x) x$type)
  
  expect_true(all(
    is.null(attr(combn6, "group_DF")),
    length(ftypes) == 1,
    all(ftypes == c("moleculeFilt"))
  ))
  
  ftypes <- attr(combn8, "filters") %>% 
    lapply(function(x) x$type)
  
  expect_true(all(
    is.null(attr(combn6, "group_DF")),
    length(ftypes) == 1,
    all(ftypes == c("cvFilt"))
  ))
  
  # keep both filters and groups
  ftypes <- attr(combn4, "filters") %>% 
    lapply(function(x) x$type)
  
  expect_true(all(
    !is.null(attr(combn4, "group_DF")),
    length(ftypes) == 2,
    all(ftypes == c("moleculeFilt", "cvFilt"))
  ))
})

