#' Calculate SPANS Score for a Number of Normalization Methods
#'
#' Ranks different combinations of subset and normalization methods based on
#' a score that captures how much bias a particular normalization procedure
#' introduces into the data. Higher score implies less bias.
#'
#' @param omicsData aobject of the class 'pepData' or 'proData' created by
#'   \code{\link{as.pepData}} or \code{\link{as.proData}} respectively. The data
#'   must be log transformed (using edata_transform()) and have a grouping
#'   structure, usually set by calling group_designation() on the object.
#' @param subset_fn character vector indicating which subset functions to test.
#'   See details for the current offerings.
#' @param norm_fn character vector indicating the normalization functions to
#'   test. See details for the current offerings.
#' @param params list of additional arguments passed to the chosen subset
#'   functions. See details for parameter specification and default values.
#' @param group character specifying a column name in f_data that gives the
#'   group assignment of the samples. Defaults to NULL, in which case the
#'   grouping structure given in \code{attr(omicsData, 'group_DF')} is used.
#' @param n_iter number of iterations used in calculating the background
#'   distribution in step 0 of SPANS. Defaults to 1000.
#' @param sig_thresh numeric value that specifies the maximum p-value for which
#'   a biomolecule can be considered highly significant based on a
#'   Kruskal-Wallis test. Defaults to 0.0001.
#' @param nonsig_thresh numeric value that specifies the minimum p-value for
#'   which a biomolecule can be considered non-significant based on a
#'   Kruskal-Wallis test. Defaults to 0.5.
#' @param min_sig integer value specifying the minimum number of highly
#'   significant biomolecules identified in step 0 of SPANS in order to proceed.
#'   sig_thresh will be adjusted to the minimum value that gives this many
#'   biomolecules.
#' @param min_nonsig integer value specifying the minimum number of
#'   non-significant biomolecules identified in step 0 of SPANS in order to
#'   proceed.  nonsig_thresh will be adjusted to the maximum value that gives
#'   this many biomolecules.
#' @param max_sig integer value specifying the maximum number of
#'   highly significant biomolecules identified in step 0 if SPANS in order to
#'   proceed.  Excesses of highly significant biomolecules will be randomly
#'   sampled down to these values.
#' @param max_nonsig integer value specifying the maximum number of
#'   non-significant biomolecules identified in step 0 if SPANS in order to
#'   proceed.  Excesses of non-significant biomolecules will be randomly sampled
#'   down to these values.
#'
#' @param ... Additional arguments \tabular{ll}{ \code{location_thresh,
#'   scale_thresh} The minimum p-value resulting from a Kruskal-Wallis test on
#'   the location and scale parameters resulting from a normalization method in
#'   order for that method to be considered a candidate for scoring.\cr
#'   \code{verbose} Logical specifying whether to print the completion of SPANS
#'   procedure steps to console. Defaults to TRUE.\cr \code{parallel} Logical
#'   specifying whether to use a parallel backend.  Depending on the size of
#'   your data, setting this to FALSE can cause the algorithm to be very slow.
#'   Defaults to TRUE. }
#' @details Below are details for specifying function and parameter options.
#' @section Subset Functions: Specifying a subset function indicates the subset
#'   of features (rows of \code{e_data}) that should be used for computing
#'   normalization factors. The following are valid options: "all", "los",
#'   "ppp", "rip", and "ppp_rip". \cr \tabular{ll}{ \tab "all" is the subset
#'   that includes all features (i.e. no subsetting is done). \cr \tab "los"
#'   identifies the subset of the features associated with the top \code{L},
#'   where \code{L} is a proportion between 0 and 1, order statistics.
#'   Specifically, the features with the top \code{L} proportion of highest
#'   absolute abundance are retained for each sample, and the union of these
#'   features is taken as the subset identified (Wang et al., 2006). \cr \tab
#'   "ppp" (orignally stands for percentage of peptides present) identifies the
#'   subset of features that are present/non-missing for a minimum
#'   \code{proportion} of samples (Karpievitch et al., 2009; Kultima et al.,
#'   2009). \cr \tab "complete" subset of features that have no missing data
#'   across all samples.  Equivalent to "ppp" with proportion = 1. \cr \tab
#'   "rip" identifies features with complete data that have a p-value greater
#'   than a defined threshold \code{alpha} (common values include 0.1 or 0.25)
#'   when subjected to a Kruskal-Wallis test based (non-parametric one-way
#'   ANOVA) on group membership (Webb-Robertson et al., 2011). \cr \tab
#'   "ppp_rip" is equivalent to "rip" however rather than requiring features
#'   with complete data, features with at least a \code{proportion} of
#'   non-missing values are subject to the Kruskal-Wallis test.\cr }
#' @section Normalization Functions: Specifying a normalization function
#'   indicates how normalization scale and location parameters should be
#'   calculated. The following are valid options: "median", "mean", "zscore",
#'   and "mad". Parameters for median centering are calculated if "median" is
#'   specified. The location estimates are the sample-wise medians of the subset
#'   data. There are no scale estimates for median centering. Parameters for
#'   mean centering are calculated if "mean" is specified. The location
#'   estimates are the sample-wise means of the subset data. There are no scale
#'   estimates for median centering. Parameters for z-score transformation are
#'   calculated if "zscore" is specified. The location estimates are the subset
#'   means for each sample. The scale estimates are the subset standard
#'   deviations for each sample. Parameters for median absolute deviation (MAD)
#'   transformation are calculated if "mad" is specified.
#'
#' @section Specifying Subset Parameters Using the \code{params} argument:
#'   Parameters for the chosen subset function should be specified in a list.
#'   The list elements should have names corresponding to the subset function
#'   inputs and contain a \emph{list} of numeric values.  The elements of
#'   ppp_rip will be length 2 numeric vectors, corresponding to the parameters
#'   for ppp and rip. See examples.
#'
#'   The following subset functions have parameters that can be specified:
#'   \tabular{ll}{ los \tab list of values between 0 and 1 indicating the top
#'   proportion of order statistics. Defaults to list(0.05,0.1,0.2,0.3) if
#'   unspecified. \cr \tab \cr ppp \tab list of values between 0 and 1
#'   specifying the proportion of samples that must have non-missing values for
#'   a feature to be retained. Defaults to list(0.1,0.25,0.50,0.75) if
#'   unspecified. \cr \tab \cr rip \tab list of values between 0 and 1
#'   specifying the p-value threshold for determining rank invariance. Defaults
#'   to list(0.1,0.15,0.2,0.25) if unspecified. \cr \tab \cr ppp_rip \tab list
#'   of length 2 numeric vectors corresponding to the RIP and PPP parameters
#'   above. Defaults list(c(0.1,0.1), c(0.25, 0.15), c(0.5, 0.2), c(0.75,0.25))
#'   if unspecified. \cr }
#'
#' @return An object of class 'SPANSRes', which is a dataframe containing
#'   columns for the subset method and normalization used, the parameters used
#'   in the subset method, and the corresponding SPANS score.  \cr
#'
#'   The column 'mols_used_in_norm' contains the number of molecules that were
#'   selected by the subset method and subsequently used to determine the
#'   location/scale parameters for normalization.  The column 'passed selection'
#'   is \code{TRUE} if the subset+normalization procedure was selected for
#'   scoring.\cr
#'
#'   The attribute 'method_selection_pvals' is a dataframe containing
#'   information on the p values used to determine if a method was selected for
#'   scoring (location_p_value, scale_p_value) as well as the probabilities
#'   (F_log_HSmPV, F_log_NSmPV) given by the empirical cdfs generated in the
#'   first step of SPANS.
#'
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' pep_object <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' pep_object <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
#'
#' ## default parameters
#' spans_res <- spans_procedure(omicsData = pep_object)
#'
#' ## specify only certain subset and normalization functions
#' spans_res <- spans_procedure(omicsData = pep_object, norm_fn = c("median", "zscore"), subset_fn = c("all", "los", "ppp"))
#'
#' ## specify parameters for supplied subset functions, notice ppp_rip takes a vector of two numeric arguments.
#' spans_res <- spans_procedure(omicsData = pep_object, subset_fn = c("all", "los", "ppp"), params = list(los = list(0.25, 0.5), ppp = list(0.15, 0.25)))
#' spans_res <- spans_procedure(omicsData = pep_object, subset_fn = c("all", "rip", "ppp_rip"), params = list(rip = list(0.3, 0.4), ppp_rip = list(c(0.15, 0.5), c(0.25, 0.5))))
#' }
#'
#' @author Daniel Claborne
#'
#' @references Webb-Robertson BJ, Matzke MM, Jacobs JM, Pounds JG, Waters KM. A
#'   statistical selection strategy for normalization procedures in LC-MS
#'   proteomics experiments through dataset-dependent ranking of normalization
#'   scaling factors. Proteomics. 2011;11(24):4736-41.
#'
#' @export
#'
spans_procedure <- function(omicsData,
                            norm_fn = c("median", "mean", "zscore", "mad"),
                            subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"),
                            params = NULL, group = NULL, n_iter = 1000,
                            sig_thresh = 0.0001, nonsig_thresh = 0.5,
                            min_nonsig = 20, min_sig = 20, max_nonsig = NULL,
                            max_sig = NULL, ...) {
  .spans_procedure(omicsData,
    norm_fn = norm_fn, subset_fn = subset_fn,
    params = params, group = group, n_iter = n_iter,
    sig_thresh = sig_thresh, nonsig_thresh = nonsig_thresh,
    min_nonsig = min_nonsig, min_sig = min_sig,
    max_sig = max_sig, max_nonsig = max_nonsig, ...
  )
}

.spans_procedure <- function(omicsData,
                             norm_fn = c("median", "mean", "zscore", "mad"),
                             subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"),
                             params = NULL, group = NULL, n_iter = 1000,
                             sig_thresh = 0.0001, nonsig_thresh = 0.5,
                             min_nonsig = 20, min_sig = 20, max_nonsig = NULL,
                             max_sig = NULL, location_thresh = 0.05,
                             scale_thresh = 0.05, verbose = TRUE,
                             parallel = TRUE) {
  edata_cname = get_edata_cname(omicsData)
  fdata_cname = get_fdata_cname(omicsData)
  nsamps = attributes(omicsData)$data_info$num_samps

  # error checks
  if (!inherits(omicsData, c("pepData", "proData"))) stop("omicsData must be of class 'pepData', or 'proData'")

  # params defaults if none specified
  if (is.null(params)) {
    params <- list(
      "los" = list(0.05, 0.1, 0.2, 0.3),
      "ppp" = list(0.1, 0.25, 0.50, 0.75),
      "rip" = list(0.1, 0.15, 0.2, 0.25),
      "ppp_rip" = list(c(0.1, 0.1), c(0.25, 0.15), c(0.5, 0.2), c(0.75, 0.25))
    )

    for (name in names(params)) {
      if (!(name %in% subset_fn)) {
        params[[name]] <- NULL
      }
    }
  }

  # simple function to check if an element of params contains all values between 0 and 1
  checkvals <- function(list) {
    sapply(list, function(x) {
      isTRUE(all(x >= 0) & all(x <= 1))
    })
  }

  if (!isTRUE(all(names(params) %in% c("los", "ppp", "rip", "ppp_rip")))) stop("params must be a named list with names of the normalization functions to be tested, one or more of: 'los', 'ppp', 'rip', 'ppp_rip'")

  # subset function parameter value checks
  if ("ppp_rip" %in% subset_fn) {
    if (is.null(params$ppp_rip)) stop("params must contain a list element named 'ppp_rip")
    if (!all(sapply(params$ppp_rip, length) == 2) | !all(checkvals(params$ppp_rip))) stop("each element of 'ppp_rip' in params must be a numeric vector of length 2 with entries between 0 and 1")
  }
  if ("los" %in% subset_fn) {
    if (is.null(params$los)) stop("params must contain a list element named 'los'")
    if (!all(checkvals(params$los))) stop("each element of 'los' in params must be a numeric value between 0 and 1")
  }
  if ("ppp" %in% subset_fn) {
    if (is.null(params$ppp)) stop("params must contain a list element named 'ppp'")
    if (!all(checkvals(params$ppp))) stop("each element of 'ppp' in params must be a numeric value between 0 and 1")
  }
  if ("rip" %in% subset_fn) {
    if (is.null(params$rip)) stop("params must contain a list element named 'rip'")
    if (!all(checkvals(params$rip))) stop("each element of 'rip' in params must be a numeric value between 0 and 1")
  }

  # assign group variable or internally call group_designation() if user specified a grouping column in f_data
  if (is.null(group)) {
    if (!inherits(attr(omicsData, "group_DF"), "data.frame")) stop("the omicsData object must have a grouping structure, usually set by calling group_designation() on the object")
    group_df = attr(omicsData, "group_DF")
    reorder = match(colnames(omicsData$e_data)[-which(colnames(omicsData$e_data) == edata_cname)], as.character(group_df[, fdata_cname]))
    group = group_df[reorder, ]$Group
  } else {
    if (!is.character(group) || !(group %in% colnames(omicsData$f_data[, -which(colnames(omicsData$f_data) == fdata_cname)]))) stop("group must be a string specifying a column in f_data by which to group by")
    omicsData <- group_designation(omicsData, main_effects = group)
    group_df = attr(omicsData, "group_DF")
    reorder = match(colnames(omicsData$e_data)[-which(colnames(omicsData$e_data) == edata_cname)], as.character(group_df[, fdata_cname]))
    group = group_df[reorder, ]$Group
  }

  # misc input checks
  if (!inherits(attr(omicsData, "group_DF"), "data.frame")) stop("the omicsData object must have a grouping structure, usually set by calling group_designation() on the object")
  if (!is.numeric(min_nonsig)) stop("min_nonsig must be numeric.")
  if (!is.numeric(min_sig)) stop("min_sig must be numeric.")
  if (any(c(min_nonsig, min_sig) < 1) | any(c(min_nonsig, min_sig) > nrow(omicsData$e_data))) stop("min_nonsig and min_sig must be an integer value no greater than the number of observed biomolecules")
  if (!is.null(max_sig)) {
    if (!is.numeric(max_sig) | max_sig <= min_sig) stop('max_sig must be an integer value greater than min_sig')
  }
  if (!is.null(max_nonsig)) {
    if (!is.numeric(max_nonsig) | max_nonsig <= min_nonsig) stop('max_nonsig must be an integer value greater than min_nonsig')
  }
  if (!is.numeric(sig_thresh)) stop("sig_thresh must be numeric.")
  if (!is.numeric(nonsig_thresh)) stop("nonsig_thresh must be numeric.")
  if (any(c(sig_thresh, nonsig_thresh) > 1) | any(c(sig_thresh, nonsig_thresh) < 0)) stop("sig_thresh and nonsig_thresh must be numeric values between 0 and 1")
  if (isTRUE(attributes(omicsData)$data_info$norm_info$norm_type == "global")) stop("a global normalization scheme has already been applied to this data, SPANS should be run on an unnormalized log-transformed pepData or proData object.")
  if (!isTRUE(grepl("log", attr(omicsData, "data_info")$data_scale))) stop("omicsData object must have been transformed to the log scale.  Either specify the attribute attr(omicsData, 'data_info')$data_scale or call edata_transform() on the omicsData object.")
  if (!all(subset_fn %in% c("all", "los", "ppp", "complete", "rip", "ppp_rip"))) stop("subset_fn must be a character vector containing more than one of the elements 'all', 'los', 'ppp', 'complete', 'rip', 'ppp_rip'")
  if (!all(norm_fn %in% c("median", "mean", "zscore", "mad"))) stop("norm_fn must be a character vector containing more than one of the elements 'median', 'mean', 'zscore', 'mad'")
  if (is.null(attr(omicsData, "group_DF"))) stop("omicsData object must have a grouping structure set by calling group_designation()")

  ### end main error checking ###

  # Calculate p-values on e_data by group --------------------------------------

  # get indices of significant and nonsignificant p-values
  kw_pvals <- kw_rcpp(omicsData$e_data %>% dplyr::select(-dplyr::all_of(edata_cname)) %>% as.matrix(), as.character(group))

  # initial storage of both vectors
  sig_inds <- (kw_pvals <= sig_thresh & !is.na(kw_pvals))
  nonsig_inds <- (kw_pvals >= nonsig_thresh & !is.na(kw_pvals))

  ## these two while loops change the p-value threshold until at least min_sig and min_nonsig indices are selected

  # low to high p-values from KW test
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>%
    unique() %>%
    sort()
  iter <- 1
  # while we dont have enough indices, add the indices corresponding to the 1:iter lowest p-values
  while (sum(sig_inds) < min_sig) {
    sig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }

  # same thing but for high p-values
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>%
    unique() %>%
    sort(decreasing = TRUE)
  iter <- 1
  while (sum(nonsig_inds) < min_nonsig) {
    nonsig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }

  # if user set a maximum for the number of significant and nonsignificant molecules and we are over the maximum
  # randomly sample max_sig indices from the indices of significant molecules ...
  if (!is.null(max_sig)) {
    if (sum(sig_inds) > max_sig) {
      make_false <- sample(which(sig_inds), sum(sig_inds) - max_sig)
      sig_inds[make_false] <- FALSE
    }
  }
  # ... and max_nonsig indices from the indices of non-significant molecules.
  if (!is.null(max_nonsig)) {
    if (sum(nonsig_inds) > max_nonsig) {
      make_false <- sample(which(nonsig_inds), sum(nonsig_inds) - max_nonsig)
      nonsig_inds[make_false] <- FALSE
    }
  }

  # get a vector of n_iter sample sizes for randomly selecting peptides to determine normalization factors
  scaling_factor <- sum(!is.na(omicsData$e_data %>% dplyr::select(-dplyr::all_of(edata_cname)))) / 100
  select_n <- ceiling(runif(n_iter, nsamps / scaling_factor, 100) * scaling_factor) - nsamps

  ### produce a list with all combinations of subset functions, normalization functions, and parameters ###
  all_calls <- list()

  # for each normalization/subset method combination...
  for (nf in norm_fn) {
    for (sf in subset_fn) {
      # ...if it is a normalization function that doesn't have parameters specified, just append a list of the normalization and subset function names...
      if (is.null(params[[sf]])) {
        all_calls[[length(all_calls) + 1]] <- list(norm_fn = nf, subset_fn = sf)
        # ...otherwise append the same list with a parameters variable as well
      } else {
        for (par in params[[sf]]) {
          if (sf == "rip") temp_par <- list("rip" = par)
          if (sf == "los") temp_par <- list("los" = par)
          if (sf == "ppp") temp_par <- list("ppp" = par)
          if (sf == "ppp_rip") temp_par <- list("ppp_rip" = list("ppp" = par[1], "rip" = par[2]))
          all_calls[[length(all_calls) + 1]] <- list(norm_fn = nf, subset_fn = sf, params = temp_par)
        }
      }
    }
  }

  if (length(all_calls) < 2) stop("Your input parameters did not result in more than 1 normalization method.")

  n_methods <- length(all_calls)

  # STEP 0: create random distribution ------------------------------------------

  # set up parallel backend
  if (parallel) {
    cores <- parallelly::availableCores()
    cl <- parallelly::makeClusterPSOCK(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  # get a median significant and non-significant p-value for n_iter iterations
  background_distribution <- foreach::foreach(
    i = 1:n_iter,
    .packages = "pmartR",
    .export = c(
      "spans_make_distribution",
      "kw_rcpp",
      "normalize_global_basic"
    )
  ) %dopar% {
    # function defined at bottom of this script
    spans_make_distribution(
      omicsData, group, norm_fn, sig_inds,
      nonsig_inds, select_n[i]
    )
  }

  # make empirical cdfs based on vectors of n_iter median p-values
  sig_cdf <- sapply(background_distribution, function(el) {
    el[[1]]
  }) %>% ecdf()
  nonsig_cdf <- sapply(background_distribution, function(el) {
    el[[2]]
  }) %>% ecdf()

  if (verbose) print("Finished creating background distribution, moving to method candidate selection")

  # STEP 1: select methods that will go on -------------------------------------

  # determine which methods (subset function + normalization function +
  # parameters combination) will be assessed in step 2. returned list contains
  # information on the method applied and a T/F value for whether it passed to
  # step 2.
  which_spans <- foreach::foreach(el = all_calls, .packages = "pmartR", .export = "kw_rcpp") %dopar% {
    norm_object <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params)
    params <- norm_object$parameters[[1]]

    p_location <- kw_rcpp(matrix(params$location, nrow = 1), group = as.character(group))

    if (!is.null(params$scale)) {
      p_scale <- kw_rcpp(matrix(params$scale, nrow = 1), group = as.character(group))
      if (any(c(p_location, p_scale) < c(location_thresh, scale_thresh))) res <- list(passfail = FALSE, step1_pvals = c(p_location, p_scale)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, p_scale))
    } else if (p_location < location_thresh) res <- list(passfail = FALSE, step1_pvals = c(p_location, NA)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, NA))

    res <- c(el, res, list(n_features_calc = norm_object$n_features_calc))

    return(res)
  }

  if (verbose) print("Finished method candidate selection, proceeding to score selected methods.")

  # STEP 2: score methods ------------------------------------------------------

  # Score each method that passed step 1 by normalizing the full data and
  # getting median Kruskal-Wallis p-values for significant and nonsignificant
  # peptides.
  scores <- foreach::foreach(el = which_spans, .packages = "pmartR", .export = c("kw_rcpp")) %dopar% {
    if (el$passfail) {
      norm_data <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params, apply_norm = TRUE)
      abundance_matrix <- norm_data$e_data %>%
        dplyr::select(-dplyr::all_of(edata_cname)) %>%
        as.matrix()

      sig_score <- -log10(median(kw_rcpp(abundance_matrix[sig_inds, ], group = as.character(group)), na.rm = TRUE))
      non_sig_score <- log10(median(kw_rcpp(abundance_matrix[nonsig_inds, ], group = as.character(group)), na.rm = TRUE))

      score <- (sig_cdf(sig_score) + nonsig_cdf(non_sig_score)) / 2

      return(list(score, sig_cdf(sig_score), nonsig_cdf(non_sig_score)))
    } else return(list(NA, NA, NA))
  }

  # The scores object is a lists of lists. The first element of each sub-list in
  # scores is the SPANS_score reported at the end of the function. Check if all
  # of these elements are NA. If they are the user will be smitten with an error
  # by the stats demigod Protector of R inputs.
  if (all(is.na(lapply(scores, `[[`, 1)))) {
    # stats demigod Protector of R inputs: "You shall not pass!"
    #
    # pmartR user: "Noooooooooooo! Why must it be like this Bad Data?"
    #
    # stats demigod Protector of R inputs: "You did not prepare your inputs
    # sufficiently. I am sorry but there is no help for you. BWAHAHAHAHAHA!"
    #
    # stats demigod Protector of R inputs smites down the pmartR user who then
    # dies dramatically.
    stop(paste("All SPANS scores are NA. This can be caused by specifying",
      "too many comparison groups with few samples.",
      sep = " "
    ))
  }

  if (verbose) print("Finished scoring selected methods")

  # Combine results into SPANS output ------------------------------------------

  # create dataframe with selected methods
  spansres_obj <- data.frame(
    "subset_method" = character(n_methods),
    "normalization_method" = character(n_methods),
    "SPANS_score" = numeric(n_methods),
    "parameters" = character(n_methods),
    "mols_used_in_norm" = numeric(n_methods),
    "passed_selection" = logical(n_methods),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  extra_info <- data.frame(
    "subset_method" = character(n_methods),
    "normalization_method" = character(n_methods),
    "parameters" = character(n_methods),
    "location_p_value" = numeric(n_methods),
    "scale_p_value" = numeric(n_methods),
    "F_log_HSmPV" = numeric(n_methods),
    "F_log_NSmPV" = numeric(n_methods),
    "SPANS_score" = numeric(n_methods),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # populate the dataframe from which_spans
  for (i in 1:n_methods) {
    ss <- which_spans[[i]]$subset_fn
    norm <- which_spans[[i]]$norm_fn
    score <- scores[[i]][[1]]
    num_mols <- which_spans[[i]]$n_features_calc
    params <- which_spans[[i]]$params %>%
      unlist() %>%
      as.character() %>%
      paste(collapse = ";")
    p_loc <- which_spans[[i]]$step1_pvals[1]
    p_scale <- which_spans[[i]]$step1_pvals[2]
    F_HSmPV <- scores[[i]][[2]]
    F_NSmPV <- scores[[i]][[3]]
    pass_fail <- which_spans[[i]]$passfail

    # store into row of df
    spansres_obj[i, ] <- list(ss, norm, score, params, num_mols, pass_fail)
    extra_info[i, ] <- list(ss, norm, params, p_loc, p_scale, F_HSmPV, F_NSmPV, score)
  }

  spansres_obj <- dplyr::arrange(spansres_obj, dplyr::desc(SPANS_score))
  extra_info <- dplyr::arrange(extra_info, dplyr::desc(SPANS_score)) %>% dplyr::select(-SPANS_score)

  attr(spansres_obj, "method_selection_pvals") <- extra_info
  attr(spansres_obj, "group_vector") = group
  attr(spansres_obj, "significant_thresh") = sig_thresh
  attr(spansres_obj, "nonsignificant_thresh") = nonsig_thresh
  attr(spansres_obj, "n_not_significant") = sum(nonsig_inds)
  attr(spansres_obj, "n_significant") = sum(sig_inds)
  attr(spansres_obj, "location_threshold") = location_thresh
  attr(spansres_obj, "scale_thresh") = scale_thresh
  class(spansres_obj) <- c("SPANSRes", "data.frame")

  return(spansres_obj)
}

#' Creates the list of median p-values used to make the background distribution
#' used to compute the SPANS score in step 2.
#'
#' @param omicsData an object of the class 'pepData' or 'proData' created by
#'   \code{\link{as.pepData}} or \code{\link{as.proData}}, respectively.
#' @param group_vector A character vector from the group_DF attribute specifying
#'   the order of the samples. This order is the same as the order of the
#'   samples (columns) in e_data.
#' @param norm_fn a character vector of normalization methods to choose from.
#'   Current options are 'mean', 'median', 'zscore', and 'mad'.
#' @param sig_inds significant peptide indices (row indices) based on a
#'   Kruskal-Wallis test on the un-normalized data
#' @param nonsig_inds non-significant peptide indices (row indices) based on a
#'   Kruskal-Wallis test on the un-normalized data
#' @param select_n number of peptide by sample indices in the data to randomly
#'   select to determine normalization parameters
#'
#' @return a list with 2 elements. The median of highly significant p-values,
#'   and the median of nonsignificant p-values. These are obtained from a SINGLE
#'   Kruskal-Wallis test on data normalized by scale/location factors determined
#'   from a randomly selected subset of peptides and normalization method
#'
spans_make_distribution <- function(omicsData, group_vector, norm_fn, sig_inds,
                                    nonsig_inds, select_n) {
  # Extract the name of the column containing the peptide IDs.
  edata_cname <- get_edata_cname(omicsData)

  # Determine the number of samples in the data. This will be used to correctly
  # subset the abundance_matrix when sampling non-NA values.
  nsamps <- attributes(omicsData)$data_info$num_samps

  # need a matrix to pass to kw_rcpp
  abundance_matrix <- omicsData$e_data %>%
    dplyr::select(-dplyr::all_of(edata_cname)) %>%
    as.matrix()

  # indices vector that will be used for subsetting
  inds <- NULL

  # Fish out the number of rows in abundance_matrix. This will be used to
  # subset this matrix by vector index instead of row/column indices.
  n_row <- nrow(abundance_matrix)

  # for each sample, randomly select an index to include for that sample, this
  # ensures each sample gets at least 1 observation remember matrices are stored
  # as vectors, we select elements using a single number
  for (j in 1:nsamps) {
    # The integer that is added to the end is to correctly subset the
    # abundance_matrix by its vector index (R numbers matrices down each row
    # then across each column).
    forced_ind <- (sample(which(!is.na(abundance_matrix[, j])), 3) +
      ((j - 1) * n_row))
    inds <- c(inds, forced_ind)
  }

  # To calculate the standard deviation when the zscore norm function is
  # selected we keep at least 3 non-NA values from each sample. If the select_n
  # value is greater than (sum(!is.na(abundance_matrix)) - 3 * nsamps) the
  # sample function will throw an error because we are trying to sample more
  # elements than are in the population. To avoid this error we will compare the
  # number of non-NA values to (3 * nsamps + select_n) and subset the
  # abundance_matrix object accordingly.
  if ((3 * nsamps + select_n) > sum(!is.na(abundance_matrix))) {
    # Manually select all non-NA values to calculate the normalization
    # parameters because the number of randomly selected non-NA values is larger
    # than the actual number of non-NA values.
    inds <- which(!is.na(abundance_matrix))
  } else {
    # randomly assign the rest of the indices
    inds <- c(inds, sample(
      setdiff(which(!is.na(abundance_matrix)), inds),
      select_n
    ))
  }

  # randomly select a normalization method
  rand_norm <- sample(norm_fn, 1)

  # get normalization parameters from subsetted matrix. normalize_global_matrix
  # is not intended to be used outside this function, it returns the location
  # and (if applicable) scale parameters for normalization.
  norm_params <- normalize_global_basic(
    edata = cbind(
      omicsData$e_data %>% dplyr::select(dplyr::all_of(edata_cname)),
      replace(abundance_matrix, -inds, NA)
    ),
    norm_fn = rand_norm
  )

  # Apply the normalization to the abundance_matrix. This must be done outside
  # the normalize_global_basic function because a reduced data set (randomly
  # selected non-NA values are changed to NA) is passed to
  # normalize_global_basic and the normalization is carried out on the full
  # data set.
  lapply(1:ncol(abundance_matrix), function(i) {
    abundance_matrix[, i] <<- abundance_matrix[, i] - norm_params$location[i]
    if (!is.null(norm_params$scale)) abundance_matrix[, i] <<- abundance_matrix[, i] / norm_params$scale[i]
  })

  # run Kruskal-Wallis on the normalized dataset
  hs_mpv <- kw_rcpp(abundance_matrix[sig_inds, ], group_vector)
  ns_mpv <- kw_rcpp(abundance_matrix[nonsig_inds, ], group_vector)

  # NA values are from peptides with no observations in 1 group
  hs_mpv <- hs_mpv[!is.na(hs_mpv)]
  ns_mpv <- ns_mpv[!is.na(ns_mpv)]

  # log transformed median p-values
  hs_mpv <- -log10(median(hs_mpv))
  ns_mpv <- log10(median(ns_mpv))

  return(list(hs_mpv, ns_mpv))
}

#' Normalize e_data within SPANS
#'
#' This function is intended to be used in SPANS only. It is a VERY trimmed down
#' version of normalize_global. It is trimmed down because within SPANS we only
#' need the norm_params element from the output of the normalize_global
#' function. All of the other options and output can be ignored.
#'
#' @param edata a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, lipids, or metabolites and \eqn{n} is the number of samples.
#'   Each row corresponds to data for a peptide, protein, lipid, or metabolite,
#'   with a column giving the identifer name.
#'
#' @param norm_fn character string indicating the normalization function to use
#'   for normalization. See details for the current offerings.
#'
#' @return A list containing the location and scale parameters for normalizing
#'   the data.
#'
normalize_global_basic <- function(edata, norm_fn) {
  # Select the function name for the normalizing function that will be used.
  fn_to_use <- switch(norm_fn,
    mean = "mean_center",
    median = "median_center",
    zscore = "zscore_transform",
    mad = "mad_transform"
  )

  # Attach the normalizing function to the name fn_to_use.
  fn_to_use <- get(fn_to_use,
    envir = asNamespace("pmartR"),
    mode = "function"
  )

  # Normalize the data with the specified function.
  norm_results <- fn_to_use(
    e_data = edata,
    edata_id = names(edata)[1],
    subset_fn = "all",
    feature_subset = edata[, 1],
    backtransform = FALSE,
    apply_norm = FALSE
  )

  # Just return the normalizing parameters.
  return(norm_results$norm_params)
}

#' Gets the parameters for the highest ranked methods from spans.
#'
#' @param SPANSRes_obj an object of the class SPANSRes obtained by calling
#'   \code{spans_procedure()}
#' @param sort_by_nmols a logical indicator of whether to sort by number of
#'   molecules used in the normalization (see \code{\link{spans_procedure}} for
#'   info about the 'mols_used_in_norm' column)
#'
#' @return A list of lists, where there are multiple sublists only if there were
#'   ties for the top SPANS score.  Each sublist contains named elements for the
#'   subset and normalization methods, and the parameters used for the subset
#'   method. \cr
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#'
#' # data must be log transformed and grouped
#' myobject <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' myobject <- group_designation(omicsData = myobject, main_effects = "Phenotype")
#'
#' spans_result <- spans_procedure(omicsData = myobject)
#'
#' # a list of the parameters for any normalization procedure with the best SPANS score
#' best_params <- get_spans_params(spans_result)
#'
#' # extract the arguments from the first list element
#' subset_fn = best_params[[1]]$subset_fn
#' norm_fn = best_params[[1]]$norm_fn
#' params = best_params[[1]]$params
#' if (is.null(params[[1]])) {
#'   params = NULL
#' }
#'
#' # pass arguments to normalize global
#' norm_object <- normalize_global(omicsData = myobject, subset_fn = subset_fn, norm_fn = norm_fn, params = params)
#' }
#'
#' @export
#'
get_spans_params <- function(SPANSRes_obj, sort_by_nmols = FALSE) {
  if (!inherits(SPANSRes_obj, "SPANSRes")) {
    # Suffer the wrath of Dread Pirate Roberts!!!!!
    stop("object must be of class 'SPANSRes'")
  }

  if (all(is.na(SPANSRes_obj$SPANS_score))) stop("No methods were selected for scoring, there is no 'best' set of parameters to return.")

  # get rows that are tied for top score
  best_df <- SPANSRes_obj %>%
    dplyr::top_n(1, wt = SPANS_score)

  if (sort_by_nmols) {
    best_df <- best_df %>%
      dplyr::top_n(1, mols_used_in_norm)
  }

  ## populate a list with the subset method, normalization method, and subset parameters.
  params <- vector("list", nrow(best_df))

  for (i in 1:nrow(best_df)) {
    params[[i]][["subset_fn"]] = best_df[i, "subset_method"]
    params[[i]][["norm_fn"]] = best_df[i, "normalization_method"]

    pars_from_df = as.numeric(strsplit(best_df[i, "parameters"], ";")[[1]])

    if (params[[i]][["subset_fn"]] == "ppp_rip") {
      params[[i]][["params"]] = list(ppp_rip = list(ppp = pars_from_df[1], rip = pars_from_df[2]))
    } else if (params[[i]][["subset_fn"]] == "all") {
      params[[i]][["params"]] = list(NULL)
    } else if (params[[i]][["subset_fn"]] %in% c("los", "rip", "ppp")) {
      sublist = list()
      sublist[[params[[i]][["subset_fn"]]]] <- pars_from_df
      params[[i]][["params"]] <- sublist
    }
  }

  return(params)
}
