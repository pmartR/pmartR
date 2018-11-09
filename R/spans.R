#' Calculate SPANS Score for a Number of Normalization Methods
#'
#' Ranks different combinations of subsetting and normalization methods based on a score that captures how much bias a particular normalization procedure introduces into the data.  Higher score implies less bias.
#'
#' @param omicsData an object of the class 'pepData' or 'proData' created by \code{\link{as.pepData}} or \code{\link{as.proData}} respectively. The data must be log transformed (using edata_transform()) and have a grouping structure, usually set by calling group_designation() on the object. 
#' @param subset_fn character vector indicating which subset functions to test. See details for the current offerings.
#' @param norm_fn character vector indicating the normalization functions to test. See details for the current offerings.
#' @param params list of additional arguments passed to the chosen subset functions. See details for parameter specification and default values.
#' @param group character vector that specifies group assignment of the samples.  Defaults to NULL, in which case the grouping structure given in \code{attr(omicsData, 'group_DF')} is used.
#' @param n_iter number of iterations used in calculating the background distribution in step 0 of SPANS.  Defaults to 1000.
#' @param sig_thresh numeric value that specifies the maximum p-value for which a biomolecule can be considered highly significant based on a Kruskal-Wallis test.  Defaults to 0.0001.
#' @param nonsig_thresh numeric value that specifies the minimum p-value for which a biomolecule can be considered non-significant based on a Kruskal-Wallis test.  Defaults to 0.5.
#' @param min_sig integer value specifying the minimum number of highly significant biomolecules identified in in step 0 order to proceed.  sig_thresh will be adjusted to the minimum value that gives this many biomolecules.
#' @param min_nonsig integer value specifying the minimum number of non-significant biomolecules identified in step 0 of SPANS.  nonsig_thresh will be adjusted to the maximum value that gives this many biomolecules.
#'
#' @param ... Additional arguments
#' \tabular{ll}{
#' \code{location_thresh, scale_thresh} The minimum p-value resulting from a Kruskal-Wallis test on the location and scale parameters resulting from a normalization method in order for that method to be considered a candidate for scoring.  
#' \code{verbose} Logical specifying whether to print the completion of SPANS procedure steps to console.  Defaults to TRUE.
#' \code{parrallel} Logical specifying whether to use a parrallel backend.  Depending on the size of your data, setting this FALSE can cause the algorithm to be very slow.  Defaults to TRUE.
#' }
#'@details Below are details for specifying function and parameter options.
#' @section Subset Functions:
#' Specifying a subset function indicates the subset of features (rows of \code{e_data}) that should be used for computing normalization factors. The following are valid options: "all", "los", "ppp", "rip", and "ppp_rip". The option "all" is the subset that includes all features (i.e. no subsetting is done). The option "los" identifies the subset of the features associated with the top \code{L}, where \code{L} is a proportion between 0 and 1, order statistics. Specifically, the features with the top \code{L} proportion of highest absolute abundance are retained for each sample, and the union of these features is taken as the subset identified (Wang et al., 2006). The option "ppp" (orignally stands for percentage of peptides present) identifies the subset of features that are present/non-missing for a minimum \code{proportion} of samples (Karpievitch et al., 2009; Kultima et al., 2009). The option "rip" identifies features with complete data that have a p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) when subjected to a Kruskal-Wallis test based (non-parametric one-way ANOVA) on group membership (Webb-Robertson et al., 2011). The option "ppp_rip" is equivalent to "rip" however rather than requiring features with complete data, features with at least a \code{proportion} of non-missing values are subject to the Kruskal-Wallis test.
#'
#' @section Normalization Functions:
#' Specifying a normalization function indicates how normalization scale and location parameters should be calculated. The following are valid options: "median", "mean", "zscore", and "mad". Parameters for median centering are calculated if "median" is specified. The location estimates are the sample-wise medians of the subset data. There are no scale estimates for median centering. Parameters for mean centering are calculated if "mean" is specified. The location estimates are the sample-wise means of the subset data. There are no scale estimates for median centering. Parameters for z-score transformation are calculated if "zscore" is specified. The location estimates are the subset means for each sample. The scale estimates are the subset standard deviations for each sample. Parameters for median absolute deviation (MAD) transformation are calculated if "mad" is specified.
#'
#' @section Specifying Subset Parameters Using the \code{params} argument:
#' Parameters for the chosen subset function should be specified in a list.  The list elements should have names corresponding to the subset function inputs and contain a \emph{list} of numeric values.  The elements of ppp_rip will be length 2 numeric vectors, corresponding to the parameters for ppp and rip.  See examples.
#'
#' The following subset functions have parameters that can be specified:
#' \tabular{ll}{
#' los \tab list of values between 0 and 1 indicating the top proportion of order statistics. Defaults to list(0.05,0.1,0.2,0.3) if unspecified. \cr
#' \tab \cr
#' ppp \tab list of values between 0 and 1 specifying the proportion of samples that must have non-missing values for a feature to be retained. Defaults to list(0.1,0.25,0.50,0.75) if unspecified. \cr
#' \tab \cr
#' rip \tab list of values between 0 and 1 specifying the p-value threshold for determining rank invariance. Defaults to list(0.1,0.15,0.2,0.25) if unspecified. \cr
#' \tab \cr
#' ppp_rip \tab list of length 2 numeric vectors corresponding to the RIP and PPP parameters above. Defaults list(c(0.1,0.1), c(0.25, 0.15), c(0.5, 0.2), c(0.75,0.25)) if unspecified.
#' \cr
#' }
#' 
#' @examples 
#' 
#' library(pmartR)
#' library(pmartRdata)
#' 
#' pep_object
#' pep_object <- edata_transform(pep_object, data_scale = "log2")
#' pep_object <- group_designation(pep_object, main_effects = "Condition")
#' 
#' ## default parameters
#' # spans_res <- spans_procedure(pep_object)
#' 
#' ## specify only certain subset and normalization functions
#' # spans_res <- spans_procedure(pep_object, norm_fn = c("median", "zscore"), subset_fn = c("all", "los", "ppp"))
#' 
#' ## specify parameters for supplied subset functions, notice ppp_rip takes a vector of two numeric arguments.
#' # spans_res <- spans_procedure(pep_object, subset_fn = c("all", "los", "ppp"), params = list(los = list(0.25, 0.5), ppp = list(0.15, 0.25)))
#' # spans_res <- spans_procedure(pep_object, subset_fn = c("all", "rip", "ppp_rip"), params = list(rip = list(0.3, 0.4), ppp_rip = list(c(0.15, 0.5), c(0.25, 0.5))))
#'
#' @author Daniel Claborne
#'
#' @references Webb-Robertson BJ, Matzke MM, Jacobs JM, Pounds JG, Waters KM. A statistical selection strategy for normalization procedures in LC-MS proteomics experiments through dataset-dependent ranking of normalization scaling factors. Proteomics. 2011;11(24):4736-41. 
#'
#' @export

spans_procedure <- function(omicsData, norm_fn = c("median", "mean", "zscore", "mad"), subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"), 
                            params = NULL, group = NULL, n_iter = 1000, sig_thresh = 0.0001, nonsig_thresh = 0.5, min_nonsig = 20, min_sig = 20, ...){
  .spans_procedure(omicsData, norm_fn = norm_fn, subset_fn = subset_fn, params = params, group = group, 
                   n_iter = n_iter, sig_thresh = sig_thresh, nonsig_thresh = nonsig_thresh, min_nonsig = min_nonsig, min_sig = min_sig, ...)
}

  .spans_procedure <- function(omicsData, norm_fn = c("median", "mean", "zscore", "mad"), subset_fn = c("all", "los", "ppp", "rip", "ppp_rip"), 
                            params = NULL, group = NULL, n_iter = 1000, sig_thresh = 0.0001, nonsig_thresh = 0.5, min_nonsig = 20, min_sig = 20,
                            location_thresh = 0.05, scale_thresh = 0.05, verbose = TRUE, parrallel = TRUE){ 
    
    edata_cname = get_edata_cname(omicsData)
    fdata_cname = get_fdata_cname(omicsData)
    nsamps = attributes(omicsData)$data_info$num_samps
    
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
    
    if(is.null(group)){
      group = attr(omicsData, "group_DF")$Group
    }
    else{
      if(!is.character(group) | !(group %in% colnames(omicsData$f_data[,-fdata_cname]))) stop("group must be a string specifying a column in f_data by which to group by")
      omicsData <- group_designation(omicsData, main_effects = group)
      group = attr(omicsData, "group_DF")$Group
    }
      
    # misc input checks
    if(!inherits(attr(omicsData, "group_DF"), "data.frame")) stop("the omicsData object must have a grouping structure, usually set by calling group_designation() on the object")
    if(any(c(min_nonsig, min_sig) < 1) | any(c(min_nonsig, min_sig) > nrow(omicsData$e_data))) stop("min_nonsig and min_sig must be an integer value no greater than the number of observed biomolecules")
    if(any(c(sig_thresh, nonsig_thresh) > 1) | any(c(sig_thresh, nonsig_thresh) < 0)) stop("sig_thresh and nonsig_thresh must be numeric values between 0 and 1")
    if(isTRUE(attributes(omicsData)$norm_info$norm_type == "global")) stop("a global normalization scheme has already been applied to this data, SPANS should be run on an unnormalized log-transformed pepData or proData object.")
    if(!isTRUE(grepl("log", attr(omicsData, "data_info")$data_scale))) stop("omicsData object must have been transformed to the log scale.  Either specify the attribute attr(omicsData, 'data_info')$data_scale or call edata_transform() on the omicsData object.")
    if(!all(subset_fn %in% c("all", "los", "ppp", "rip", "ppp_rip"))) stop("subset_fn must be a character vector containing more than one of the elements 'all', 'los', 'ppp', 'rip', 'ppp_rip'")
    if(!all(norm_fn %in% c("median", "mean", "zscore", "mad"))) stop("norm_fn must be a character vector containing more than one of the elements 'median', 'mean', 'zscore', 'mad'")
    if(is.null(attr(omicsData, "group_DF"))) stop("omicsData object must have a grouping structure set by calling group_designation()")
    
    ### end main error checking ###
    
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
    select_n <- ceiling(runif(n_iter, nsamps/scaling_factor, 100)*scaling_factor) - nsamps
    
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
    
    if(length(all_calls) < 2) stop("Your input parameters did not result in more than 1 normalization method.")
    
    n_methods <- length(all_calls)
    
    ### STEP 0:  create random distribution ####
    
    # set up parallel backend
    if(parrallel){
      cores<- parallel::detectCores()
      cl<- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      doParallel::registerDoParallel(cl)
    } 
    else foreach::registerDoSEQ()
    
    # get a median significant and non-significant p-value for n_iter iterations
    background_distribution <- foreach::foreach(i = 1:n_iter, .packages = "pmartR", .export = c("spans_make_distribution", "kw_rcpp", "normalize_global_matrix")) %dopar% {
      # function defined at bottom of this script
      spans_make_distribution(omicsData, norm_fn, sig_inds, nonsig_inds, select_n[i])
    }
    
    # make empirical cdfs based on vectors of n_iter median p-values
    sig_cdf <- sapply(background_distribution, function(el){el[[1]]}) %>% ecdf()
    nonsig_cdf <- sapply(background_distribution, function(el){el[[2]]}) %>% ecdf()
    
    if(verbose) print("Finished creating background distribution, moving to method candidate selection")
    
    #### STEP 1:  
    # determine which methods (subset function + normalization function + parameters combination) will be assessed in step 2.
    # returned list contains information on the method applied and a T/F value for whether it passed to step 2.
    which_spans <- foreach::foreach(el = all_calls, .packages = "pmartR", .export = "kw_rcpp") %dopar% {
      norm_object <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params)
      params <- norm_object$parameters[[1]]
      
      p_location <- kw_rcpp(matrix(params$location, nrow = 1), group = as.character(group))
      
      if(!is.null(params$scale)){
        p_scale <- kw_rcpp(matrix(params$scale, nrow = 1), group = as.character(group))
        if(any(c(p_location, p_scale) < c(location_thresh, scale_thresh))) res <- list(passfail = FALSE, step1_pvals = c(p_location, p_scale)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, p_scale))
      }
      else if(p_location < location_thresh) res <- list(passfail = FALSE, step1_pvals = c(p_location, NA)) else res <- list(passfail = TRUE, step1_pvals = c(p_location, NA))
      
      res <- c(el, res, list(n_features_calc = norm_object$n_features_calc))
      
      return(res)
      
    }
    
    if(verbose) print("Finished method candidate selection, proceeding to score selected methods.")
    
    # STEP 2: Score each method that passed step 1 by normalizing the full data and getting median Kruskal-Wallis p-values for significant and nonsignificant peptides
    
    scores <- foreach::foreach(el = which_spans, .packages = "pmartR", .export = c("kw_rcpp"), .combine = "cbind") %dopar% {
      if(el$passfail){
        norm_data <- normalize_global(omicsData, el$subset_fn, el$norm_fn, params = el$params, apply_norm = TRUE)
        abundance_matrix <- norm_data$e_data %>% dplyr::select(-edata_cname) %>% as.matrix()
        
        sig_score <- -log10(median(kw_rcpp(abundance_matrix[sig_inds,], group = as.character(group)), na.rm = TRUE))
        non_sig_score <- log10(median(kw_rcpp(abundance_matrix[nonsig_inds,], group = as.character(group)), na.rm = TRUE))
        
        score <- (sig_cdf(sig_score) + nonsig_cdf(non_sig_score))/2
        
        return(score)
      }
      else return(NA)
    }
    
    if(verbose) print("Finished scoring selected methods")
    
    # create dataframe with selected methods
    spansres_obj <- data.frame("subset_method" = character(n_methods), "normalization_method" = character(n_methods), "SPANS_score" = numeric(n_methods), 
                          "parameters" = character(n_methods), "mols_used_in_norm" = numeric(n_methods), "passed_selection" = logical(n_methods), stringsAsFactors = FALSE)
    
    extra_info <- data.frame("subset method" = character(n_methods), "normalization method" = character(n_methods), "parameters" = character(n_methods), "location_p_value" = numeric(n_methods), "scale_p_value" = numeric(n_methods), stringsAsFactors = FALSE)
    
    # populate the dataframe from which_spans
    for(i in 1:n_methods){
      ss <- which_spans[[i]]$subset_fn
      norm <- which_spans[[i]]$norm_fn
      score <- scores[i]
      num_mols <- which_spans[[i]]$n_features_calc
      params <- which_spans[[i]]$params %>% unlist() %>% as.character() %>% paste(collapse = ";")
      p_loc <- which_spans[[i]]$step1_pvals[1]
      p_scale <- which_spans[[i]]$step1_pvals[2]
      pass_fail <- which_spans[[i]]$passfail
      
      # store into row of df
      spansres_obj[i,] <- c(ss, norm, score, params, num_mols, pass_fail)
      extra_info[i, ] <- c(ss, norm, params, p_loc, p_scale)
    }
    
    spansres_obj <- dplyr::arrange(spansres_obj, desc(SPANS_score))
    extra_info <- dplyr::arrange(extra_info, desc(spansres_obj$SPANS_score))
    
    spansres_obj <- spansres_obj
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
  inds <- c(inds, sample(setdiff(which(!is.na(abundance_matrix)), inds), select_n))
  
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