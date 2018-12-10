context("test spans_procedure....will take a minute.....")
library(testthat)
library(pmartRdata)
library(pmartR)

valid_norm_fn <- c("median", "mean", "zscore", "mad")
valid_subset_fn = c("all", "los", "ppp", "rip", "ppp_rip")

rand_norm_fn <- c("zscore", "mean")
rand_subset_fn <- c("ppp_rip", "all")

# object must be log2 transformed and grouped
expect_error(spans_procedure(pep_object))

myobject <- edata_transform(pep_object, data_scale = "log2")
expect_error(spans_procedure(myobject))

myobject <- group_designation(myobject, main_effects = "Condition")

expect_error(spans_procedure(myobject, norm_fn = rand_norm_fn, subset_fn = "badname"))
expect_error(spans_procedure(myobject, norm_fn = "badname", subset_fn = rand_subset_fn))
expect_error(spans_procedure(myobject, norm_fn = rand_norm_fn, subset_fn = c("all", "ppp", "los"), params = list(ppp = c(1.1, 0.2), los = c(0.2, 0.4))))
expect_error(spans_procedure(myobject, norm_fn = rand_norm_fn, subset_fn = c("all", "ppp", "ppp_rip"), params = list(ppp = c(0.1, 0.2), los = c(0.2, 0.4))))
expect_error(spans_procedure(myobject, norm_fn = "mean", subset_fn = "ppp_rip", params = list(ppp_rip = list(c(0.1, 0.2)))))

# make a small spansres_object so testing isn't slowed too much
spansres_obj <- spans_procedure(myobject, norm_fn = rand_norm_fn, subset_fn = rand_subset_fn, verbose = FALSE)
# get the parameters of the highest scoring methods
spans_params <- get_spans_params(spansres_obj)

# randomly sample 1 of these sets of methods/parameters and see that it works when passed to normalize_global
i <- sample(1:length(spans_params), 1)
norm_fn <- spans_params[[i]]$norm_fn
subset_fn <- spans_params[[i]]$subset_fn
params <- spans_params[[i]]$params
norm_object <- normalize_global(myobject, norm_fn = norm_fn, subset_fn = subset_fn, params = params)

# get matrix values and group vector to pass to kw_rcpp
norm_params <- norm_object$parameters[[1]]
groupvec = as.character(attr(spansres_obj, "group"))

# output from kw_rcpp
p_location <- kw_rcpp(matrix(norm_params$location, nrow = 1), groupvec)
if(!is.null(norm_params$scale)){
  p_scale <- kw_rcpp(matrix(norm_params$scale, nrow = 1), groupvec)
}

# location and scale p-values from spansres object
ploc_comp <- attr(spansres_obj, "method_selection_pvals")[i, "location_p_value"]
pscale_comp <- attr(spansres_obj, "method_selection_pvals")[i, "scale_p_value"]

test_that("spansres output is consistent with input", {
  expect_true(norm_fn %in% rand_norm_fn)
  expect_true(subset_fn %in% rand_subset_fn)
  expect_equal(ploc_comp, p_location)
  if(!is.na(pscale_comp)) expect_equal(pscale_comp, p_scale)
})
