# This function is intended to be used in SPANS only, it is a very trimmed down version of normalize_global.  It is necessary because the kw_rcpp function only operates on matrices

normalize_global_matrix <- function(edata, norm_fn, params = NULL, apply_norm = FALSE, backtransform = FALSE, min_prop = NULL){
  
  # check for valid normalization function choice #
  if(!(norm_fn %in% c("mean", "median", "zscore", "mad")))stop(paste(norm_fn, " is not a valid subset option", sep = ""))

  # check that apply_norm is T/F #
  if(!inherits(apply_norm, "logical")) stop("apply_norm must be a logical argument")
  
  #check that backtransform is T/F #
  if(!inherits(backtransform, "logical")) stop("backtransform must a logical argument")

  fn_to_use <- switch(norm_fn, mean = "mean_center", median = "median_center", zscore = "zscore_transform", mad = "mad_transform")
  
  fn_to_use <- get(fn_to_use, envir = asNamespace("pmartR"), mode = "function")
  
  norm_results <- fn_to_use(e_data = edata, edata_id = names(edata)[1], feature_subset = edata[,1], backtransform = backtransform, apply_norm = apply_norm, check.names=FALSE)
  
  res <- norm_results$norm_params
  
  return(res)
}






