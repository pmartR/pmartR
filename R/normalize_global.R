#' Calculate Normalization Parameters
#'
#' Calculates normalization parameters based on the data using the specified subset and normalization functions with possibility of apply the normalization to the data.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively. The function \code{\link{group_designation}} must have been run on omicsData to use several of the subset functions (i.e. rip and ppp_rip).
#' @param subset_fn character string indicating the subset function to use for normalization. See details for the current offerings.
#' @param norm_fn character string indicating the normalization function to use for normalization. See details for the current offerings.
#' @param params additional arguments passed to the chosen subset functions. See details for parameter specification and default values.
#' @param apply_norm logical argument indicating if the normalization should be applied to the data. Defaults to FALSE. If TRUE, the normalization is applied to the data and an S3 object of the same class as \code{omicsData} (e.g. 'pepData') with normalized values in \code{e_data} is returned.
#' @param backtransform logical argument indicating if parameters for back transforming the data, after normalization, should be calculated.  Defaults to FALSE. If TRUE, the parameters for back transforming the data after normalization will be calculated, and subsequently included in the data normalization if \code{apply_norm} is TRUE or \code{\link{apply.normRes}} is used downstream. See details for an explanation of how these factors are calculated.
#' @param min_prop numeric threshold between 0 and 1 giving the minimum value for the proportion of features subset (rows of \code{e_data})
#'
#' @details Below are details for specifying function and parameter options.
#' @section Subset Functions:
#' Specifying a subset function indicates the subset of features (rows of \code{e_data}) that should be used for computing normalization factors. The following are valid options: "all", "los", "ppp", "complete", "rip", and "ppp_rip". The option "all" is the subset that includes all features (i.e. no subsetting is done). The option "los" identifies the subset of the features associated with the top \code{L}, where \code{L} is a proportion between 0 and 1, order statistics. Specifically, the features with the top \code{L} proportion of highest absolute abundance are retained for each sample, and the union of these features is taken as the subset identified (Wang et al., 2006). The option "ppp" (orignally stands for percentage of peptides present) identifies the subset of features that are present/non-missing for a minimum \code{proportion} of samples (Karpievitch et al., 2009; Kultima et al., 2009). The option "complete" retains molecules with no missing data across all samples, equivalent to "ppp" with proportion = 1.  The option "rip" identifies features with complete data that have a p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) when subjected to a Kruskal-Wallis test based (non-parametric one-way ANOVA) on group membership (Webb-Robertson et al., 2011). The option "ppp_rip" is equivalent to "rip" however rather than requiring features with complete data, features with at least a \code{proportion} of non-missing values are subject to the Kruskal-Wallis test.
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
#' @return If apply_norm is FALSE, an S3 object of type 'normRes' is returned. This object contains a list with: subset method, normalization method, normalization parameters, number of features used in normalization, and proportion of features used in normalization. plot() and summary() methods are available for this object. If apply_norm is TRUE, then the normalized data is returned in an object of the appropriate S3 class (e.g. pepData).
#' 
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object <- edata_transform(omicsData = lipid_object, data_scale="log2")
#' lipid_object <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' norm_object <- normalize_global(omicsData = lipid_object, subset_fn = "all", norm_fn = "median")
#' norm_data <- normalize_global(omicsData = lipid_object, subset_fn = "all", norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
#'}
#'
#' @author Lisa Bramer
#' 
#' @references
#'
#' @export
#' 
normalize_global <- function (omicsData, subset_fn, norm_fn, params = NULL,
                              apply_norm = FALSE, backtransform = FALSE,
                              min_prop = NULL) {
  
  ## initial checks ##
  
  # Ensure omicsData is an appropriate class
  if (!inherits(omicsData, c("pepData", "proData", "lipidData", "metabData",
                             "nmrData"))) {
    
    # Throw down an error for blatant abuse of the pmartR standards.
    stop (paste("omicsData must be of class 'pepData', 'proData', 'lipidData',",
                "'metabData' or 'nmrData'",
                sep = " "))
    
  }
  
  # data should be log transformed #
  if (!(get_data_scale(omicsData) %in% c("log2", "log10", "log"))) {
    
    # Tell the user to read the documentation before running functions willy
    # nilly.
    stop (paste("omicsData$e_data should be log transformed prior to calling",
                "normalize_data. See documentation for the edata_transform",
                "function for more information.",
                sep = " "))
    
  }
  
  #check that min_prop is greater than 0 but less than or equal to 1
  if (!is.null(min_prop)) {
    
    if (min_prop < 0 || min_prop >= 1) {
      
      # Stop the user for not understanding what a proportion is.
      stop ("min_prop must be greater than zero but less than or equal to 1")
      
    }
    
  }
  
  check_names = getchecknames(omicsData)
  edata_id <- attr(omicsData, "cnames")$edata_cname
  samp_id <- attr(omicsData, "cnames")$fdata_cname
  
  #Use default normalization and subsetting if not specified
  if(missing(subset_fn)) stop("subset_fn wasn't specified")
  if(missing(norm_fn)) stop("norm_fn wasn't specified")
  
  # check for valid subset function choice #
  if (!(subset_fn %in% c("all", "los", "ppp", "rip", "ppp_rip", "complete"))) {
    
    # There is no hope!!!
    stop (paste(subset_fn, " is not a valid subset option", sep = ""))
    
  }
  
  # Check if subset_fn is "all" or "complete" and params is not NULL.
  if (subset_fn %in% c("all", "complete") && !is.null(params)) {
    
    # Throw an error because "all" and "complete" do not require any parameters
    # to be specified.
    stop (paste("The subset functions 'all' and 'complete' do not require",
                "params to be specified. Either set params = NULL or change",
                "subset_fn to another option.",
                sep = " "))
    
  }
  
  # check for valid normalization function choice #
  if (!(norm_fn %in% c("mean", "median", "zscore", "mad"))) {
    
    # It is really quite simple to select one of the four options.
    stop (paste(norm_fn, " is not a valid normalization option", sep = "")) 
    
  }
  
  # check for group information #
  grp_pres = length(grep("group_DF", names(attributes(omicsData))))
  
  # check that group designation was run if "rip" is involved
  if (subset_fn %in% c("rip", "ppp_rip") && grp_pres == 0) {
    
    stop (paste("group_designation() must be run on the data if subset_fn is",
                "'rip' or 'ppp_rip'",
                sep = " "))
    
  }
  
  # check that parameter was specified for subset functions other than all #
  if(subset_fn != "all"){
    
    # ppp_rip #
    if (subset_fn == "ppp_rip") {
      
      # set default parameter values, will be overwritten if user specifies #
      params_ppp = 0.5
      params_rip = 0.2
      
      # check that parameters are either specified with "ppp_rip" or separate
      # options, if both are specified, error out #
      if (length(names(params)) > 1 && "ppp_rip" %in% names(params)) {
        
        stop (paste("Too many arguments in 'params'. Specify either",
                    "'ppp_rip = ' or 'ppp = ' and 'rip = '",
                    sep = " "))
        
      }
      
      # Check if the parameters were specified with ppp. If it wasn't the
      # default value will be used.
      if ("ppp" %in% names(params)) {
        
        # check that parameters are valid #
        if (params$ppp < 0 || params$ppp > 1) {
          
          # STOP!! JUST STOP!!!
          stop ("The ppp parameter must be between 0 and 1")
          
        }
        
        # Change the ppp value to the one the user specified.
        # IF they got it right.
        params_ppp = params$ppp
        
        
      }
      
      # Check if the parameters were specified with rip. If it wasn't the
      # default value will be used.
      if ("rip" %in% names(params)) {
        
        # Make sure rip is correctly specified. (You can't ever trust the user.)
        if (params$rip < 0 || params$rip > 1) {
          
          # *SIGH*
          stop ("The rip parameter must be between 0 and 1")
          
        }
        
        # Change the rip value to the one the user specified.
        params_rip = params$rip
        
      }
      
      # Check if the parameters were specified using the ppp_rip option.
      if ("ppp_rip" %in% names(params)) {
        
        # check that at least one of the two parameters is specified #
        if (sum(c("ppp", "rip") %in% names(params[["ppp_rip"]])) < 1) {
          
          stop ("Invalid parameter specification for 'ppp_rip'")
          
        }
        
        # Check if ppp was specified. If it wasn't the default value will be
        # used.
        if ("ppp" %in% names(params[["ppp_rip"]])) {
          
          # check that parameters are valid #
          if (params[["ppp_rip"]]$ppp < 0 || params[["ppp_rip"]]$ppp > 1) {
            
            stop ("The ppp parameter must be between 0 and 1")
            
          }
          
          # Change the ppp value to the one the user specified.
          # IF they got it right.
          params_ppp = params[["ppp_rip"]]$ppp
          
        }
        
        # Check if rip was specified.
        if ("rip" %in% names(params)) {
          
          # Check that rip is a valid number.
          if (params[["ppp_rip"]]$rip < 0 || params[["ppp_rip"]]$rip > 1) {
            
            stop ("The rip parameter must be between 0 and 1")
            
          }
          
          # Change the ppp value to the one the user specified.
          params_rip = params[["ppp_rip"]]$rip
          
        }
        
      }
      
      # An option other than "all" or "ppp_rip" was input.
    } else {
      
      # Set to NULL in case the user specified a subset function that requires a
      # value for params but they did NOT specify a value for param.
      param_val = NULL
      
      # Check if the name of the subset function is also the name specified in
      # the params argument.
      if (subset_fn %in% names(params)) {
        
        # Make sure they got the value right.
        if (params[[subset_fn]] < 0 || params[[subset_fn]] > 1) {
          
          stop ("Specified parameter must be between 0 and 1")
          
        }
        
        # Set the value to the one the user input.
        param_val = params[[subset_fn]]
        
      }
      
    }
    
  }

  # check that apply_norm is T/F #
  if (!inherits(apply_norm, "logical")) {
    
    # Stop the illogical user in their tracks!
    stop ("apply_norm must be a logical argument")
    
  }
  
  #check that backtransform is T/F #
  if (!inherits(backtransform, "logical")) {
    
    # Stop the illogical user again!!
    stop ("backtransform must a logical argument")
    
  }

  # Check if the apply_norm and backtransform arguments agree.
  if (apply_norm == FALSE & backtransform == TRUE) {
    
    # STOP! You cannot backtransform if no normalization occurred!!!
    stop (paste("apply_norm is set to FALSE and backtransform is set to TRUE;",
                "a backtransform cannot be applied if the normalization is not",
                "being applied.",
                sep = " "))
    
  }
  
  # Subset the data prior to normalization -------------------------------------
  
  # subset data using current subset method #
  if (subset_fn == "all") {
    
    # keep all biomolecule IDs.
    peps = all_subset(omicsData$e_data, edata_id)
    
  } else if (subset_fn == "los") {
    
    # Check if the user specified a cutoff value.
    if (!is.null(param_val)) {
      
      # Apply the users cutoff value.
      peps = los(omicsData$e_data, edata_id, param_val)
      
      # No cutoff value was specified by the user.
    } else {
      
      # Use the default value to subset.
      peps = los(omicsData$e_data, edata_id)
      
    }
    
  } else if (subset_fn == "rip") {
    
    # Extricate the group information.
    group_df = attr(omicsData, "group_DF")
    
    # Check if the user specified a cutoff value.
    if (!is.null(param_val)) {
      
      # Apply the users cutoff value.
      peps = rip(omicsData$e_data, edata_id, samp_id, group_df, param_val)
      
      # No cutoff value was specified by the user.
    } else {
      
      # Use the default value to subset.
      peps = rip(omicsData$e_data, edata_id, samp_id, group_df)
      
    }
    
  } else if (subset_fn == "complete") {
    
    # Only keep rows with non-missing values for all samples.
    peps = complete_mols(omicsData$e_data, edata_id)
    
  } else if (subset_fn == "ppp") {
    
    # Check if the user specified a cutoff value.
    if (!is.null(param_val)) {
      
      # Apply the users cutoff value.
      peps = ppp(omicsData$e_data, edata_id, param_val)
    } else {
      
      # Use the default value to subset.
      peps = ppp(omicsData$e_data, edata_id)
      
    }
    
  } else if (subset_fn == "ppp_rip") {
    
    # Separate the group attribute from the pack.
    group_df = attr(omicsData, "group_DF")
    
    # Use the parameters defined above (either the default values or the values
    # set by the user) to subset the data.
    peps = ppp_rip(omicsData$e_data,
                   edata_id,
                   samp_id,
                   group_df,
                   alpha = params_rip,
                   proportion = params_ppp)
    
  }
  
  # Normalize the data ---------------------------------------------------------
  
  # Pull out the name of the normalization function from the input.
  fn_to_use <- switch (norm_fn,
                       mean = mean_center,
                       median = median_center,
                       zscore = zscore_transform,
                       mad = mad_transform)
  
  # Normalize the data according to the method selected.
  norm_results <- fn_to_use(e_data = omicsData$e_data,
                            edata_id = edata_id,
                            subset_fn = subset_fn,
                            feature_subset = peps,
                            backtransform = backtransform,
                            apply_norm = apply_norm,
                            check.names = check_names)
  
  # Calculate the proportion of features used to calculate the normalization
  # parameters.
  prop_features_calc <- length(peps) / nrow(omicsData$e_data)
  
  # If min_prop is specified check if it is below the actual proportion of
  # features used for normalizing.
  if (!is.null(min_prop) && (prop_features_calc < min_prop)) {
    
    # Throw an error because the proportion of features used for normalizing is
    # less than the threshold.
    stop (paste("The minimum proportion of features allowed (min_prop) was",
                "specified as", min_prop, "but the actual proportion of",
                "features used to calculate the normalization parameters using",
                "the given subset function (subset_fn) was",
                round(prop_features_calc, 3), "hence normalization was not",
                "carried out." , sep = " "))
    
  }
  
  # Check if the normalization will be applied. If return all of the information
  # needed to normalize (which biomolecules to normalize with, normalization
  # method and parameters, ...).
  if (apply_norm == FALSE) {
    
    # Return all of the information needed to normalize the data. This will be
    # turned into its own class.
    res = list(subset_fn = subset_fn,
               norm_fn = norm_fn,
               parameters = list(
                 normalization = norm_results$norm_params,
                 backtransform = norm_results$backtransform_params
               ), 
               n_features_calc = length(peps),
               feature_subset = peps,
               prop_features_calc = prop_features_calc)
    
    # Create a normRes object.
    class(res) = "normRes"
    
    # Add the input data as an attribute to the normRes object.
    attributes(res)$omicsData = omicsData

    # Return the normRes object.
    return (res)

    # This will run if apply_norm is TRUE.
  } else {

    # Update e_data with the normalized data.
    omicsData$e_data = norm_results$transf_data
    
    # Update the norm_info attribute.
    attributes(omicsData)$data_info$norm_info <- list(
      is_normalized = TRUE,
      norm_type = "global",  # added 12/21/17 by KS
      subset_fn = subset_fn,
      subset_params = params,
      norm_fn = norm_fn,
      n_features_calc = length(peps),
      prop_features_calc = prop_features_calc,
      params = list(norm_scale = norm_results$norm_params$scale,
                    norm_location = norm_results$norm_params$location,
                    bt_scale = norm_results$backtransform_params$scale,
                    bt_location = norm_results$backtransform_params$location)
    )

  }

  # Return the normalized omicsData object with all of its updated attributes!
  return (omicsData)
  
}
