#' Calculate Normalization Parameters
#'
#' Calculates normalization parameters based on the data using the specified subset and normalization functions with possibility of apply the normalization to the data.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#' @param subset_fn character vector indicating which subset functions to test. See details for the current offerings.
#' @param norm_fn character vector indicating the normalization functions to test. See details for the current offerings.
#' @param params additional arguments passed to the chosen subset functions. See details for parameter specification and default values.
#' @param group character vector that specifies group assignment of the samples.  Defaults to NULL, in which case the grouping structure given in \code{attr(omicsData, 'group_DF')} is used.
#' @param n_iter number of iterations used in calculating the background distribution in step 0 of SPANS.  Defaults to 1000.
#' @param sig_thresh numeric value that specifies the maximum p-value for which a biomolecule can be considered highly significant based on a Kruskal-Wallis test.  Defaults to 0.0001.
#' @param nonsig_thresh numeric value that specifies the minimum p-value for which a biomolecule can be considered non-significant based on a Kruskal-Wallis test.  Defaults to 0.5.
#' @param min_sig integer value specifying the minimum number of highly significant biomolecules identified in in step 0 order to proceed.  sig_thresh will be adjusted to the minimum value that gives this many biomolecules.
#' @param min_nonsig integer value specifying the minimum number of non-significant biomolecules identified in step 0 of SPANS.  nons0
#'
#'@details Below are details for specifying function and parameter options.
#' @section Subset Functions:
#' Specifying a subset function indicates the subset of features (rows of \code{e_data}) that should be used for computing normalization factors. The following are valid options: "all", "los", "ppp", "rip", and "ppp_rip". The option "all" is the subset that includes all features (i.e. no subsetting is done). The option "los" identifies the subset of the features associated with the top \code{L}, where \code{L} is a proportion between 0 and 1, order statistics. Specifically, the features with the top \code{L} proportion of highest absolute abundance are retained for each sample, and the union of these features is taken as the subset identified (Wang et al., 2006). The option "ppp" (orignally stands for percentage of peptides present) identifies the subset of features that are present/non-missing for a minimum \code{proportion} of samples (Karpievitch et al., 2009; Kultima et al., 2009). The option "rip" identifies features with complete data that have a p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) when subjected to a Kruskal-Wallis test based (non-parametric one-way ANOVA) on group membership (Webb-Robertson et al., 2011). The option "ppp_rip" is equivalent to "rip" however rather than requiring features with complete data, features with at least a \code{proportion} of non-missing values are subject to the Kruskal-Wallis test.
#'
#' @section Normalization Functions:
#' Specifying a normalization function indicates how normalization scale and location parameters should be calculated. The following are valid options: "median", "mean", "zscore", and "mad". Parameters for median centering are calculated if "median" is specified. The location estimates are the sample-wise medians of the subset data. There are no scale estimates for median centering. Parameters for mean centering are calculated if "mean" is specified. The location estimates are the sample-wise means of the subset data. There are no scale estimates for median centering. Parameters for z-score transformation are calculated if "zscore" is specified. The location estimates are the subset means for each sample. The scale estimates are the subset standard deviations for each sample. Parameters for median absolute deviation (MAD) transformation are calculated if "mad" is specified.
#'
#' @section Specifying Subset Parameters Using the \code{params} Argument:
#' Parameters for the chosen subset function should be specified in a list 
#'
#' The following subset functions have parameters that can be specified:
#' \tabular{ll}{
#' los \tab a value between 0 and 1 indicating the top proportion of order statistics. Defaults to 0.05 if unspecified. \cr
#' \tab \cr
#' ppp \tab a value between 0 and 1 specifying the proportion of samples that must have non-missing values for a feature to be retained. Defaults to 0.5 if unspecified. \cr
#' \tab \cr
#' rip \tab a value between 0 and 1 specifying the p-value threshold for determining rank invariance. Defaults to 0.2 if unspecified. \cr
#' \tab \cr
#' ppp_rip \tab two values corresponding to the RIP and PPP parameters above. Defaults to 0.5 and 0.2, respectively.
#' \cr
#' }
#'
#' @author Daniel Claborne
#'
#' @export

spans_procedure <- function(omicsData, norm_fn = c("median", "mean", "zscore", "mad"), subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"), 
                            params = NULL, group = NULL, n_iter = 1000, sig_thresh = 0.0001, nonsig_thresh = 0.5, min_nonsig = 20, min_sig = 20){
  
  # error checks
  
  if(!inherits(omicsData, c("pepData", "proData"))) stop("omicsData must be of class 'pepData', or 'proData'")
  
  # params defaults if none specified
  if(is.null(params)){
    params <- list("los" = list(0.05,0.1,0.2,0.3), 
                   "ppp" = list(0.1,0.25,0.50,0.75), 
                   "rip" = list(0.1,0.15,0.2,0.25), 
                   "ppp_rip" = list(c(0.1,0.1), c(0.25, 0.15), c(0.5, 0.2), c(0.75,0.25)))
  }
  
  # simple function to check if an element of params contains all values between 0 and 1
  checkvals <- function(list){
    sapply(list, function(x){
      isTRUE(all(x >= 0) & all(x<=1))
    })
  }
  
  if(!isTRUE(all(names(params) %in% c("los", "ppp", "rip", "ppp_rip")))) stop("params must be a named list with names of the normalization functions to be tested, one or more of: 'los', 'ppp', 'rip', 'ppp_rip'")
  
  # subset function parameter value checks
  if("ppp_rip" %in% subset_fn){
    if(is.null(params$ppp_rip)) stop("params must contain a list element named 'ppp_rip")
    if(!all(sapply(params$ppp_rip, length) == 2) | !all(checkvals(params$ppp_rip))) stop("each element of 'ppp_rip' in params must be a numeric vector of length 2 with entries between 0 and 1")
  }
  if("los" %in% subset_fn){
    if(is.null(params$los)) stop("params must contain a list element named 'los'")
    if(!all(checkvals(params$los))) stop("each element of 'los' in params must be a numeric value between 0 and 1")
  }
  if("ppp" %in% subset_fn){
    if(is.null(params$ppp)) stop("params must contain a list element named 'ppp'")
    if(!all(checkvals(params$ppp))) stop("each element of 'ppp' in params must be a numeric value between 0 and 1")
  }
  if("rip" %in% subset_fn){
    if(is.null(params$rip)) stop("params must contain a list element named 'rip'")
    if(!all(checkvals(params$rip))) stop("each element of 'rip' in params must be a numeric value between 0 and 1")
  }
    
  # misc input checks
  if(any(c(min_nonsig, min_sig) < 1) | any(c(min_nonsig, min_sig) > nrow(omicsData$e_data))) stop("min_nonsig and min_sig must be an integer value no greater than the number of observed biomolecules")
  if(any(c(sig_thresh, nonsig_thresh) > 1) | any(c(sig_thresh, nonsig_thresh) < 0)) stop("sig_thresh and nonsig_thresh must be numeric values between 0 and 1")
  if(isTRUE(attributes(omicsData)$norm_info$norm_type == "global")) stop("a global normalization scheme has already been applied to this data, SPANS should be run on an unnormalized log-transformed pepData or proData object.")
  if(!isTRUE(grepl("log", attr(omicsData, "data_info")$data_scale))) stop("omicsData object must have been transformed to the log scale.  Either specify the attribute attr(omicsData, 'data_info')$data_scale or call edata_transform() on the omicsData object.")
  if(!all(subset_fn %in% c("all", "los", "ppp", "rip", "ppp_rip"))) stop("subset_fn must be a character vector containing more than one of the elements 'all', 'los', 'ppp', 'rip', 'ppp_rip'")
  if(!all(norm_fn %in% c("median", "mean", "zscore", "mad"))) stop("norm_fn must be a character vector containing more than one of the elements 'median', 'mean', 'zscore', 'mad'")
  if(is.null(attr(omicsData, "group_DF"))) stop("omicsData object must have a grouping structure set by calling group_designation()")
  
  ### end error checking ###
  
  edata_cname = get_edata_cname(omicsData)
  group = attr(omicsData, "group_DF")$Group
  nsamps = attributes(omicsData)$data_info$num_samps
  
  # get indices of significant and nonsignificant p-values
  kw_pvals <- kw_rcpp(omicsData$e_data %>% dplyr::select(-edata_cname) %>% as.matrix(), as.character(group))
  
  # initial storage of both vectors
  sig_inds <- (kw_pvals <= sig_thresh & !is.na(kw_pvals))
  nonsig_inds <- (kw_pvals >= nonsig_thresh & !is.na(kw_pvals))
  
  ## these two while loops change the p-value threshold until at least min_sig and min_nonsig indices are selected
  
  # low to high p-values from KW test
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>% unique() %>% sort()
  iter <- 1
  # while we dont have enough indices, add the indices corresponding to the 1:iter lowest p-values
  while(sum(sig_inds) < min_sig){
    sig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }
  
  # same thing but for high p-values
  ordered_vals <- kw_pvals[!is.na(kw_pvals)] %>% unique() %>% sort(decreasing = TRUE)
  iter <- 1
  while(sum(nonsig_inds) < min_nonsig){
    nonsig_inds <- kw_pvals %in% ordered_vals[1:iter]
    iter <- iter + 1
  }
  
  # get a vector of n_iter sample sizes for randomly selecting peptides to determine normalization factors
  scaling_factor <- sum(!is.na(omicsData$e_data %>% dplyr::select(-edata_cname)))/100
  select_n <- round(runif(n_iter, nsamps/scaling_factor, 100)*scaling_factor) - nsamps
  
  ### produce a list with all combinations of subset funcctions, normalization functions, and parameters ###
  all_calls <- list()
  
  # for each normalization/subset method combination...
  for(nf in norm_fn){
    for(sf in subset_fn){
      #...if it is a normalization function that doesn't have parameters specified, just append a list of the normalization and subset function names...
      if(is.null(params[[sf]])){
        all_calls[[length(all_calls)+1]] <- list(norm_fn = nf, subset_fn = sf)
      }
      #...otherwise append the same list with a parameters variable as well
      else{
        for(par in params[[sf]]){
          if(sf == "rip") temp_par <- list("rip" = par)
          if(sf == "los") temp_par <- list("los" = par)
          if(sf == "ppp") temp_par <- list("ppp" = par)
          if(sf == "ppp_rip") temp_par <- list("ppp_rip" = list("ppp" = par[1], "rip" = par[2]))
          all_calls[[length(all_calls)+1]] <- list(norm_fn = nf, subset_fn = sf, params = temp_par)
        }
      }
    }
  }
  
  n_methods <- length(all_calls)
  
  ### STEP 0:  create random distribution ####
  
  # set up parallel backend
  cores<- parallel::detectCores()
  cl<- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  # get a median significant and non-significant p-value for n_iter iterations
  background_distribution <- foreach::foreach(i = 1:n_iter, .packages = "pmartR", .export = c("spans_make_distribution", "kw_rcpp", "normalize_global_matrix")) %dopar% {
    # function defined at bottom of this script
    spans_make_distribution(omicsData, norm_fn, sig_inds, nonsig_inds, select_n[i])
  }
  
  # make empirical cdfs based on vectors of n_iter median p-values
  sig_cdf <- sapply(background_distribution, function(el){el[[1]]}) %>% ecdf()
  nonsig_cdf <- sapply(background_distribution, function(el){el[[2]]}) %>% ecdf()
  
  print("Finished creating background distribution, moving to STEP 1")
  
  #### STEP 1:  
  # determine which methods (subset function + normalization function + parameters combination) will be assessed in step 2.
  # returned list contains information on the method applied and a T/F value for whether it passed to step 2.
  which_spans <- foreach::foreach(el = all_calls, .packages = "pmartR", .export = "kw_rcpp") %dopar% {
    norm_object <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params)
    params <- norm_object$parameters[[1]]
    
    p_location <- kw_rcpp(matrix(params$location, nrow = 1), group = as.character(group))
    
    if(!is.null(params$scale)){
      p_scale <- kw_rcpp(matrix(params$scale, nrow = 1), group = as.character(group))
      if(any(c(p_location, p_scale) < 0.05)) res <- list(passfail = FALSE, step1_pvals = c(p_location, p_scale)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, p_scale))
    }
    else if(p_location < 0.05) res <- list(passfail = FALSE, step1_pvals = c(p_location, NA)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, NA))
    
    res <- c(el, res, list(n_features_calc = norm_object$n_features_calc))
    
    return(res)
    
  }
  
  print("Finished method candidate selection")
  
  # STEP 2: Score each method that passed step 1 by normalizing the full data and getting median Kruskal-Wallis p-values for significant and nonsignificant peptides
  
  scores <- foreach::foreach(el = which_spans, .packages = "pmartR", .export = c("kw_rcpp"), .combine = "cbind") %dopar% {
    if(el$passfail){
      norm_data <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params, apply_norm = TRUE)
      abundance_matrix <- norm_data$e_data %>% dplyr::select(-edata_cname) %>% as.matrix()
      
      sig_score <- -log10(median(kw_rcpp(abundance_matrix[sig_inds,], group = as.character(group))))
      non_sig_score <- log10(median(kw_rcpp(abundance_matrix[nonsig_inds,], group = as.character(group))))
      
      score <- (sig_cdf(sig_score) + nonsig_cdf(non_sig_score))/2
      
      return(score)
    }
    else return(NA)
  }
  
  print("Finished scoring selected methods")
  
  # close the cluster
  parallel::stopCluster(cl)
  
  # create dataframe with selected methods
  results <- data.frame("subset method" = character(n_methods), "normalization method" = character(n_methods), "SPANS_score" = character(n_methods), "Num_norm_peps" = numeric(n_methods),
                        "parameters" = character(n_methods), "location_p_value" = numeric(n_methods), "scale_p_value" = numeric(n_methods), "Passed_Step_1" = logical(n_methods), stringsAsFactors = FALSE)
  
  # populate the dataframe from which_spans
  for(i in 1:n_methods){
    ss <- which_spans[[i]]$subset_fn
    norm <- which_spans[[i]]$norm_fn
    score <- scores[i]
    num_peps <- which_spans[[i]]$n_features_calc
    params <- which_spans[[i]]$params %>% unlist() %>% as.character() %>% paste(collapse = ";")
    p_loc <- which_spans[[i]]$step1_pvals[1]
    p_scale <- which_spans[[i]]$step1_pvals[2]
    pass_fail <- which_spans[[i]]$passfail
    results[i,] <- c(ss, norm, score, num_peps, params, p_loc, p_scale, pass_fail)
  }
  
  results <- dplyr::arrange(results, desc(SPANS_score))
  
  return(results)
}

#' Creates the list of median p-values used to make the background distribution used to compute the SPANS score in step 2.
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'lipidData', or 'metabData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or  \code{\link{as.metabData}}, respectively.
#' @param norm_fn a character vector of normalization methods to choose from.  Current options are 'mean', 'median', 'zscore', and 'mad'.
#' @param sig_inds significant peptide indices (row indices) based on a Kruskal-Wallis test on the un-normalized data
#' @param nonsig_inds non-significant peptide indices (row indices) based on a Kruskal-Wallis test on the un-normalized data
#' @param select_n number of peptide by sample indices in the data to randomly select to determine normalization parameters
#' 
#' @return a list with 2 elements. The median of highly significant p-values, and the median of nonsignificant p-values.
#'  These are obtained from a SINGLE Kruskal-Wallis test on data normalized by scale/location factors determined from a randomly selected subset of peptides and normalization method
spans_make_distribution <- function(omicsData, norm_fn, sig_inds, nonsig_inds, select_n){
  edata_cname = get_edata_cname(omicsData)
  group_vector = as.character(attr(omicsData, "group_DF")$Group)
  nsamps = attributes(omicsData)$data_info$num_samps
  
  # need a matrix to pass to kw_rcpp
  abundance_matrix <- omicsData$e_data %>% dplyr::select(-edata_cname) %>% as.matrix()
  
  # indices vector that will be used for subsetting
  inds <- NULL
  
  # for each sample, randomly select an index to include for that sample, this ensures each sample gets at least 1 observation
  # remember matrices are stored as vectors, we select elements using a single number
  for(j in 1:nsamps){
    forced_ind <- sample(which(!is.na(abundance_matrix[,j])), 1)*j
    inds <- c(inds, forced_ind)
  }
  
  # randomly assign the rest of the indices
  inds <- c(inds, sample(setdiff(1:length(abundance_matrix), inds), select_n - nsamps))
  
  # randomly select a normalization method
  rand_norm <- sample(norm_fn ,1)
  
  # get normalization parameters from subsetted matrix.  
  # normalize_global_matrix is not intended to be used outside this function, it returns the location and (if applicable) scale parameters for normalization.
  norm_params <- normalize_global_matrix(cbind(omicsData$e_data %>% dplyr::select(edata_cname), replace(abundance_matrix, -inds, NA)),
                                                              rand_norm, apply_norm = FALSE)
  
  # apply the normalization
  lapply(1:ncol(abundance_matrix), function(i){
    abundance_matrix[,i] <<- abundance_matrix[,i] - norm_params$location[i]
    if(!is.null(norm_params$scale)) abundance_matrix[,i] <<- abundance_matrix[,i]/norm_params$scale[i]
  })
  
  # run Kruskal-Wallis on the normalized dataset
  hs_mpv <- kw_rcpp(abundance_matrix[sig_inds,], group_vector)
  ns_mpv <- kw_rcpp(abundance_matrix[nonsig_inds,], group_vector)
  
  # NA values are from peptides with no observations in 1 group
  hs_mpv <- hs_mpv[!is.na(hs_mpv)]
  ns_mpv <- ns_mpv[!is.na(ns_mpv)]
  
  # log transformed median p-values
  hs_mpv <- -log10(median(hs_mpv))
  ns_mpv <- log10(median(ns_mpv))

  return(list(hs_mpv, ns_mpv))
}