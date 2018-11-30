#' Calculate Normalization Parameters
#'
#' Calculates normalization parameters based on the data using the specified subset and normalization functions with possibility of apply the normalization to the data.
#'
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#' @param subset_fn character string indicating the subset function to use for normalization. See details for the current offerings.
#' @param norm_fn character string indicating the normalization function to use for normalization. See details for the current offerings.
#' @param params additional arguments passed to the chosen subset functions. See details for parameter specification and default values.
#' @param apply_norm logical argument indicating if the normalization should be applied to the data. Defaults to FALSE. If TRUE, the normalization is applied to the data and an S3 object of the same class as \code{omicsData} (e.g. 'pepData') with normalized values in \code{e_data} is returned.
#' @param backtransform logical argument indicating if parameters for back transforming the data, after normalization, should be calculated.  Defaults to FALSE. If TRUE, the parameters for back transforming the data after normalization will be calculated, and subsequently included in the data normalization if \code{apply_norm} is TRUE or \code{\link{apply.normRes}} is used downstream. See details for an explanation of how these factors are calculated.
#' @param min_prop numeric threshold between 0 and 1 giving the minimum value for the proportion of features subset (rows of \code{e_data})
#'
#'@details Below are details for specifying function and parameter options.
#'@section Subset Functions:
#' Specifying a subset function indicates the subset of features (rows of \code{e_data}) that should be used for computing normalization factors. The following are valid options: "all", "los", "ppp", "rip", and "ppp_rip". The option "all" is the subset that includes all features (i.e. no subsetting is done). The option "los" identifies the subset of the features associated with the top \code{L}, where \code{L} is a proportion between 0 and 1, order statistics. Specifically, the features with the top \code{L} proportion of highest absolute abundance are retained for each sample, and the union of these features is taken as the subset identified (Wang et al., 2006). The option "ppp" (orignally stands for percentage of peptides present) identifies the subset of features that are present/non-missing for a minimum \code{proportion} of samples (Karpievitch et al., 2009; Kultima et al., 2009). The option "rip" identifies features with complete data that have a p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) when subjected to a Kruskal-Wallis test based (non-parametric one-way ANOVA) on group membership (Webb-Robertson et al., 2011). The option "ppp_rip" is equivalent to "rip" however rather than requiring features with complete data, features with at least a \code{proportion} of non-missing values are subject to the Kruskal-Wallis test.
#'
#' @section Normalization Functions:
#' Specifying a normalization function indicates how normalization scale and location parameters should be calculated. The following are valid options: "median", "mean", "zscore", and "mad". Parameters for median centering are calculated if "median" is specified. The location estimates are the sample-wise medians of the subset data. There are no scale estimates for median centering. Parameters for mean centering are calculated if "mean" is specified. The location estimates are the sample-wise means of the subset data. There are no scale estimates for median centering. Parameters for z-score transformation are calculated if "zscore" is specified. The location estimates are the subset means for each sample. The scale estimates are the subset standard deviations for each sample. Parameters for median absolute deviation (MAD) transformation are calculated if "mad" is specified.
#'
#' @section Specifying Subset Parameters Using the \code{params} Argument:
#' Parameters for the chosen subset function should be specified in a list with the function specification followed by an equal sign and the desired parameter value. For example, if LOS with 0.1 is desired, one should use \code{params = list(los = 0.1)}. ppp_rip can be specified in one of two ways: specify the parameters with each separate function or combine using a nested list (e.g. \code{params = list(ppp_rip = list(ppp = 0.5, rip = 0.2))}).
#'
#' The following functions have parameters that can be specified:
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
#'
#' @section Backtransform:
#' The purpose of back transforming data is to ensure values are on a scale similar to their raw values before normaliztion. The following values are calculated and/or applied for backtransformation purposes:
#' \tabular{ll}{
#' \code{median} \tab scale is NULL and location parameter is a global median across all samples \cr
#' \tab \cr
#' \code{mean} \tab scale is NULL and location parameter is a global median across all samples \cr
#' \tab \cr
#' \code{zscore} \tab scale is pooled standard deviation and location is global mean across all samples \cr
#' \tab \cr
#' \code{mad} \tab scale is pooled median absolute deviation and location is global median across all samples
#' \cr
#' }
#'
#'
#'
#' @return If apply_norm is FALSE, an S3 object of type 'normRes' is returned. This object contains a list with: subset method, normalization method, normalization parameters, number of features used in normalization, and proportion of features used in normalization. plot() and summary() methods are available for this object. If apply_norm is TRUE, then the normalized data is returned in an object of the appropriate S3 class (e.g. pepData).
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object <- edata_transform(omicsData = lipid_object, data_scale="log2")
#' lipid_object <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' norm_object <- normalize_global(omicsData = lipid_object, subset_fn = "all", norm_fn = "median")
#' norm_object <- normalize_global(omicsData = lipid_object, subset_fn = "all", norm_fn = "median", apply_norm = FALSE, backtransform = TRUE)
#' norm_data <- normalize_global(omicsData = lipid_object, subset_fn = "all", norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
#'}
#'
#' @author Lisa Bramer
#' @references
#'
#' @export
normalize_global <- function(omicsData, subset_fn, norm_fn, params = NULL, apply_norm = FALSE, backtransform = FALSE, min_prop = NULL){
  
  ## initial checks ##
  
  #Store data class as it will be referred to often
  dat_class <- class(omicsData)
  
  # check for group information #
  grp_pres = length(grep("group_DF", names(attributes(omicsData))))
  
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # data should be log transformed #
  if(!attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")){
    stop("omicsData$e_data should be log transformed prior to calling normalize_data. See documentation for edata_transform function for more information.")
  }
  
  #check that min_prop is greater than 0 but less than or equal to 1
  if(!is.null(min_prop)){
    if(min_prop <= 0 | min_prop > 1) stop("min_prop must be greater than zero but less than or equal to 1")
  }
  
  check_names = getchecknames(omicsData)
  edata_id <- attr(omicsData, "cnames")$edata_cname
  samp_id <- attr(omicsData, "cnames")$fdata_cname
  
  #Use default normalization and subsetting if not specified
  if(missing(subset_fn)) stop("subset_fn wasn't specified")
  if(missing(norm_fn)) stop("norm_fn wasn't specified")
  
  # check for valid subset function choice #
  if(!(subset_fn %in% c("all", "los", "ppp", "rip", "ppp_rip")))stop(paste(subset_fn, " is not a valid subset option", sep = ""))
  
  # check for valid normalization function choice #
  if(!(norm_fn %in% c("mean", "median", "zscore", "mad")))stop(paste(norm_fn, " is not a valid subset option", sep = ""))
  
  # check that group designation was run if "rip" is involved
  if(subset_fn %in% c("rip", "ppp_rip") & grp_pres == 0)stop("group_designation() must be run on the data if subset_fn is 'rip' or 'ppp_rip'")
  
  # check that parameter was specified for subset functions other than all #
  if(subset_fn != "all"){
    
    # ppp_rip #
    if(subset_fn == "ppp_rip"){
      
      # set default parameter values, will be overwritten if user specifies #
      params_ppp = 0.5
      params_rip = 0.2
      
      # check that parameters are either specified with "ppp_rip" or separate options, if both are specified, error out #
      if(sum(c("ppp", "rip") %in% names(params)) > 1 & "ppp_rip" %in% names(params))stop("Too many arguments in 'params'. Specifying only one of 'ppp_rip = ' or 'ppp'/'rip'")
      
      if("ppp" %in% names(params)){
        
        # check that parameters are valid #
        if(params$ppp < 0 | params$ppp > 1)stop("ppp parameter must be between 0 and 1")
        
        params_ppp = params$ppp
      }
      
      if("rip" %in% names(params)){
        if(params$rip < 0 | params$rip > 1)stop("rip parameter must be between 0 and 1")
        
        params_rip = params$rip
      }
      
      if("ppp_rip" %in% names(params)){
        
        # check that at least one of the two parameters is specified #
        if(sum(c("ppp", "rip") %in% names(params[["ppp_rip"]])) < 1)stop("Invalid parameter specification for 'ppp_rip'")
        
        if("ppp" %in% names(params[["ppp_rip"]])){
          # check that parameters are valid #
          if(params[["ppp_rip"]]$ppp < 0 | params[["ppp_rip"]]$ppp > 1)stop("ppp parameter must be between 0 and 1")
          
          params_ppp = params[["ppp_rip"]]$ppp
        }
        if("rip" %in% names(params)){
          if(params[["ppp_rip"]]$rip < 0 | params[["ppp_rip"]]$rip > 1)stop("rip parameter must be between 0 and 1")
          
          params_rip = params[["ppp_rip"]]$rip
        }
        
        
      }
    }else{
      param_val = NULL
      # if something other than "all" or "ppp_rip" is specified #
      if(subset_fn %in% names(params)){
        if(params[[subset_fn]] < 0 | params[[subset_fn]] > 1)stop("Specified parameter must be between 0 and 1")
        param_val = params[[subset_fn]]
      }
    }
    
  }

  # check that apply_norm is T/F #
  if(!inherits(apply_norm, "logical")) stop("apply_norm must be a logical argument")
  
  #check that backtransform is T/F #
  if(!inherits(backtransform, "logical")) stop("backtransform must a logical argument")

  
  # subset data using current subset method #
  if(subset_fn == "all"){
    peps = all_subset(omicsData$e_data, edata_id)
  }
  if(subset_fn == "los"){
    if(!is.null(param_val)){
      peps = los(omicsData$e_data, edata_id, param_val)
    }else{
      peps = los(omicsData$e_data, edata_id)
    }
  }
  if(subset_fn == "rip"){
    group_df = attr(omicsData, "group_DF")
    if(!is.null(param_val)){
      peps = rip(omicsData$e_data, edata_id, samp_id, group_df, param_val)
    }else{
      peps = rip(omicsData$e_data, edata_id, samp_id, group_df)
    }
  }
  if(subset_fn == "ppp"){
    if(!is.null(param_val)){
      peps = ppp(omicsData$e_data, edata_id, param_val)
    }else{
      peps = ppp(omicsData$e_data, edata_id)
    }
  }
  if(subset_fn == "ppp_rip"){
    group_df = attr(omicsData, "group_DF")
    peps = ppp_rip(omicsData$e_data, edata_id, samp_id, group_df, alpha = params_rip, proportion = params_ppp)
  }
  
  fn_to_use <- switch(norm_fn, mean = mean_center, median = median_center, zscore = zscore_transform, mad = mad_transform)
  
  norm_results <- fn_to_use(e_data = omicsData$e_data, edata_id=edata_id, feature_subset = peps, backtransform = backtransform, apply_norm = apply_norm, check.names=check_names)
  
  if(apply_norm == FALSE){
    res = list(subset_fn = subset_fn, norm_fn = norm_fn, parameters = list(normalization = norm_results$norm_params, backtransform = norm_results$backtransform_params), n_features_calc = length(peps), prop_features_calc = length(peps)/nrow(omicsData$e_data))
    
    class(res) = "normRes"
    attributes(res)$omicsData = omicsData
    
    if(!is.null(min_prop)){
      if(res[[5]] < min_prop) message(paste("The minimum proportion of features allowed (min_prop) was specified as", min_prop, "but the actual proportion of features used to calculate the normalization factors using the given subset function (subset_fn) was", round(res[[5]], 3), sep = " "))
    }
    
  }else{
    omicsData$e_data = norm_results$transf_data
    attributes(omicsData)$data_info$norm_info$is_normalized = TRUE
    attributes(omicsData)$data_info$norm_info$norm_type = "global"  # added 12/21/17 by KS
    attributes(omicsData)$data_info$norm_info$subset_fn = subset_fn
    attributes(omicsData)$data_info$norm_info$subset_params = params
    attributes(omicsData)$data_info$norm_info$norm_fn = norm_fn
    attributes(omicsData)$data_info$norm_info$n_features_calc = length(peps)
    attributes(omicsData)$data_info$norm_info$prop_features_calc = length(peps)/nrow(omicsData$e_data)
    attributes(omicsData)$data_info$norm_info$params$norm_scale = norm_results$norm_params$scale
    attributes(omicsData)$data_info$norm_info$params$norm_location = norm_results$norm_params$location
    attributes(omicsData)$data_info$norm_info$params$bt_scale = norm_results$backtransform_params$scale
    attributes(omicsData)$data_info$norm_info$params$bt_location = norm_results$backtransform_params$location
    res = omicsData
    
    if(!is.null(min_prop)){
      prop_features_calc = attributes(omicsData)$data_info$norm_info$prop_features_calc
        if(prop_features_calc < min_prop) stop(paste("The minimum proportion of features allowed (min_prop) was specified as", min_prop, "but the actual proportion of features used to calculate the normalization factors using the given subset function (subset_fn) was", round(prop_features_calc, 3), "hence normalization was not carried out." , sep = " "))
    }
  }
  
  return(res)
}






