source(system.file('testdata', 'load_data.R', package = 'pmartR'), local = T)

obj1 <- edata_transform(ldata, "log2")
obj1 <- normalize_global(obj1, "all", "median", apply_norm = T)

fake_cov <- c(rep("A", 5), rep("B", 6))
obj1$f_data["cov"] <- fake_cov

obj1 <- group_designation(obj1, "Condition", covariates = "cov")
obj2 <- obj1

# Some fake edata ID's to make it unique
obj2$e_data[,get_edata_cname(obj2)] <- paste0("obj2_", obj2$e_data[,get_edata_cname(obj2)])
obj2$e_meta[,get_edata_cname(obj2)] <- paste0("obj2_", obj2$e_meta[,get_edata_cname(obj2)])

obj1 <- applyFilt(molecule_filter(obj1),obj1, min_num = 2)
obj2 <- applyFilt(cv_filter(obj2),obj2, cv_thresh = 60)

suppressWarnings({
  combn1 <- combine_lipidData(obj1, obj2)
  combn2 <- combine_lipidData(obj1, obj2, retain_groups = T)
  combn3 <- combine_lipidData(obj1, obj2, retain_groups = F, retain_filters = T)
  combn4 <- combine_lipidData(obj1, obj2, retain_groups = T, retain_filters = T)
})

test_that("attributes correctly stored", {
  expect_true(all(
    is.null(attr(combn1, "group_DF")),
    length(attr(combn1, "filters")) == 0
  ))

  # drop filters and grouping info
  expect_true(all(
    !is.null(attr(combn2, "group_DF")),
    length(attr(combn2, "filters")) == 0
  ))

  # no filters, keep groups
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

  ftypes <- attr(combn4, "filters") %>%
    lapply(function(x) x$type)

  # keep both filters and groups
  expect_true(all(
    !is.null(attr(combn4, "group_DF")),
    length(ftypes) == 2,
    all(ftypes == c("moleculeFilt", "cvFilt"))
  ))
})

