#' Molecule filter object
#'
#' This function returns a moleculeFilt object for use with
#' \code{\link{applyFilt}}
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param min_num an integer value specifying the minimum number of times each
#'   feature must be observed across all samples. Default value is 2.
#' @param use_group logical indicator for whether to utilize group information
#'   from \code{\link{group_designation}} when calculating the molecule filter.
#'   Defaults to FALSE.
#' @param use_batch logical indicator for whether to utilize batch information
#'   from \code{\link{group_designation}} when calculating the molecule filter.
#'   Defaults to FALSE. Necessary to be set to TRUE if running ComBat batch
#'   correction.
#'
#' @details Attribute of molecule_filt object is "total_poss_obs", the number of
#'   total possible observations for each feature (same as the number of
#'   samples)
#'
#' @return Object of class moleculeFilt (also a data.frame) that contains the
#'   molecule identifier and the number of samples for which the molecule was
#'   measured (not NA)
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' to_filter <- molecule_filter(omicsData = pep_object)
#' summary(to_filter, min_num = 2)
#' }
#'
#' @author Kelly Stratton
#'
#' @export
#'
molecule_filter <- function (omicsData,use_groups = FALSE, use_batch = FALSE) {
  ## some initial checks ##
  # test#

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {

    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }

  # Make sure the arguemnts are logical.
  if (!is.logical(use_groups)) stop ("use_groups must be logical.")
  if (!is.logical(use_batch)) stop ("use_batch must be logical.")

  # check that omicsData has batch_id data if specified
  if(is.null(attributes(attr(omicsData,"group_DF"))$batch_id) && use_batch == TRUE){
    stop (paste("omicsData must have batch_id specified if use_batch = TRUE"))
  }

  if(is.null(attr(omicsData,"group_DF")) && use_groups == TRUE){
    stop (paste("omicsData must have groups specified if use_groups = TRUE"))
  }

  # find the column which has the edata cname
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))

  ordering = omicsData$e_data[,id_col]
  

  # SCENARIO 1: use_groups = FALSE, use_batch = FALSE
  # we run the scenario as before
  if((use_batch == FALSE | is.null(attributes(attr(omicsData,"group_DF"))$batch_id)) & (use_groups == FALSE | is.null(attr(omicsData,"group_DF")))){
    # Extricate the column number of the ID column.

    # Compute the number of non-missing values
    num_obs <- rowSums(!is.na(omicsData$e_data[, -id_col]))

    # Create a data frame with the ID column and the number of non-missing values.
    output <- data.frame(omicsData$e_data[, id_col], num_obs)
  }

  # SCENARIO 2: use_groups = FALSE, use_batch = TRUE
  else if((use_batch == TRUE & !is.null(attributes(attr(omicsData,"group_DF"))$batch_id)) & (use_groups == FALSE | is.null(attr(omicsData,"group_DF")))){
    # create a data frame with ID columns and the number of non-missing values per group
    # save the group data frame
    batchDat <- attributes(attr(omicsData,"group_DF"))$batch_id
    colnames(batchDat)[2] <- "Batch"
    # Create a data frame with the ID columns and the minimum number of non-missing values per grouping
    output <- omicsData$e_data %>%
      tidyr::pivot_longer(cols = -tidyselect::all_of(id_col), names_to = names(batchDat)[1], values_to = "value") %>%
      dplyr::left_join(batchDat, by = pmartR::get_fdata_cname(omicsData)) %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col)), Batch) %>%
      dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col))) %>%
      dplyr::summarise(min_num_obs = as.numeric(min(num_obs)),.groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::rename(molecule = tidyselect::all_of(id_col)) %>% 
      dplyr::arrange(match(molecule,ordering)) %>% 
      data.frame()
    colnames(output)[1] <- get_edata_cname(omicsData)
  }

  # SCENARIO 3: use_groups = TRUE, use_batch = FALSE
  else if((use_batch == FALSE| is.null(attributes(attr(omicsData,"group_DF"))$batch_id)) & (use_groups == TRUE & !is.null(attr(omicsData,"group_DF")))){
    # create a data frame with ID columns and the number of non-missing values per group
    # save the group data frame
    groupDat <- attr(omicsData,"group_DF")
    # Create a data frame with the ID columns and the minimum number of non-missing values per grouping
    output <- omicsData$e_data %>%
      tidyr::pivot_longer(cols = -tidyselect::all_of(id_col), names_to = names(groupDat)[1], values_to = "value") %>%
      dplyr::left_join(groupDat, by = pmartR::get_fdata_cname(omicsData)) %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col)), Group) %>%
      dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col))) %>%
      dplyr::summarise(min_num_obs = as.numeric(min(num_obs)),.groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::rename(molecule = tidyselect::all_of(id_col)) %>%
      dplyr::arrange(match(molecule,ordering)) %>%
      data.frame()
    colnames(output)[1] <- get_edata_cname(omicsData)
  }

  # SCENARIO 4: use_groups = TRUE, use_batch = TRUE
  else {
    groupDat <- attr(omicsData,"group_DF")
    batchDat <- attributes(attr(omicsData,"group_DF"))$batch_id
    colnames(batchDat)[2] <- "Batch"

    output <- omicsData$e_data %>%
      tidyr::pivot_longer(cols = -tidyselect::all_of(id_col), names_to = names(groupDat)[1], values_to = "value") %>%
      dplyr::left_join(groupDat, by = pmartR::get_fdata_cname(omicsData)) %>%
      dplyr::left_join(batchDat, by = pmartR::get_fdata_cname(omicsData)) %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col)), Group, Batch) %>%
      dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
      dplyr::group_by(dplyr::across(tidyselect::all_of(id_col))) %>%
      dplyr::summarise(min_num_obs = as.numeric(min(num_obs)),.groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::rename(molecule = tidyselect::all_of(id_col)) %>%
      dplyr::arrange(match(molecule,ordering)) %>%
      data.frame()
    colnames(output)[1] <- get_edata_cname(omicsData)
  }

  # change the names of the data.frame
  names(output) <- c(get_edata_cname(omicsData), "Num_Observations")

  # Extract the 'data.frame' class from the the output data frame.
  orig_class <- class(output)

  # Create the moleculeFilt class and attach the data.frame class to it as well.
  class(output) <- c("moleculeFilt", orig_class)

  # Fabricate an attribute that has the total number of samples (columns in
  # e_data minus the ID column). This will be used to ensure someone doesn't try
  # to filter e_data using a threshold larger than the number of samples.
  attr(output, "num_samps") <- get_data_info(omicsData)$num_samps

  # Add the group designation information to the attributes.
  attr(output, "group_DF") <- attr(omicsData, "group_DF")

  # Fabricate an attribute that states whether or not we have added a batch_id
  attr(output, "use_batch") <- ifelse(use_batch == FALSE,FALSE,TRUE)
  attr(output, "use_groups") <- ifelse(use_groups == FALSE,FALSE,TRUE)

  # Return the completed object!!!
  return(output)
}

#'Filter Based on Pooled Coefficient of Variation (CV) Values
#'
#'A pooled CV is calculated for each biomolecule.
#'
#'@param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'  'lipidData', or 'nmrData' created by \code{\link{as.pepData}},
#'  \code{\link{as.proData}}, \code{\link{as.metabData}},
#'  \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively. Note,
#'  if \code{\link{group_designation}} has not been run, the CV is calculated
#'  based on all samples for each biomolecule.
#'@param use_groups logical indicator for whether to utilize group information
#'  from \code{\link{group_designation}} when calculating the CV. Defaults to
#'  TRUE. If use_groups is set to TRUE but \code{\link{group_designation}} has
#'  not been run on the omicsData object, use_groups will be treated as FALSE.
#'
#'@return  An S3 object of class 'cvFilt' giving the pooled CV for each
#'  biomolecule and additional attributes used for plotting a data.frame with a
#'  column giving the biomolecule name and a column giving the pooled CV value.
#'
#'@details For each biomolecule, the CV of each group is calculated as the
#'  standard deviation divided by the mean, excluding missing values. A pooled
#'  CV estimate is then calculated based on the methods of Ahmed (1995). Any
#'  groups consisting of a single sample are excluded from the CV calculation,
#'  and thus, from the cv_filter result. If group_designation has not been run
#'  on the omicsData object, all samples are considered to belong to the same
#'  group.
#'@references Ahmed, S.E. (1995). \emph{A pooling methodology for coefficient of
#'  variation}. The Indian Journal of Statistics. 57: 57-75.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_object2 <- group_designation(omicsData = pep_object,
#'                                  main_effects = "Condition")
#' to_filter <- cv_filter(omicsData = pep_object2, use_groups = TRUE)
#' summary(to_filter, cv_threshold = 30)
#'}
#'
#'@author Lisa Bramer, Kelly Stratton
#'
#'@export
#'
cv_filter <- function(omicsData, use_groups = TRUE) {

  # Run some preliminary checks ------------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {

    # Follow the instructions foul creature!!!
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))

  }

  # check that use_groups is valid #
  if (!is.logical(use_groups)) {

    # Let the user know with an error that they are not logical.
    stop ("Argument 'use_groups' must be either TRUE or FALSE")

  }

  # Check if use_groups is TRUE but the group_designation function has not been
  # run yet.
  if (use_groups == TRUE && is.null(attr(omicsData, "group_DF"))) {

    # Set use groups to FALSE and continue running. At this point it is clear
    # the user doesn't get it. We will help them out a little by doing some of
    # the work for them.
    use_groups <- FALSE

  }

  # Prepare data for CV calculations -------------------------------------------

  # Extricate the column number of the ID column.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))

  # Check the data scale.
  if (get_data_scale(omicsData) == "abundance") {

    # Create a copy of the original data. This copy is created so we don't have
    # to repeat the code below (with a different name for omicsData$e_data)
    # depending on what scale the input data are in.
    cur_edata <- omicsData$e_data[, -id_col]

    # The following code is run if the data is NOT on the abundance scale.
  } else {

    # Convert the data back to the abundance scale to perform the CV
    # calculations. Remove the column containing the biomolecule IDs.
    cur_edata <- edata_transform(omicsData = omicsData,
                                 data_scale = "abundance")$e_data[, -id_col]

  }

  # Conduct CV calculations ----------------------------------------------------

  # Calculate the unpooled CV if group_designation has not been run or if
  # use_groups = FALSE.
  if (use_groups == FALSE || is.null(attr(omicsData, "group_DF"))) {

    # Calculate the unpooled coefficient of variation.
    cvs <- unpooled_cv_rcpp(as.matrix(cur_edata))

    # For mysterious reasons multiply the unpooled CV by 100.
    cvs <- 100 * cvs

    # Set the is_pooled variable to FALSE (because we didn't pool the CV).
    is_pooled <- FALSE

    # Calculate the pooled CV if use_groups is TRUE and the group_designation
    # function has been run.
  } else if (use_groups == TRUE && !is.null(attr(omicsData, "group_DF"))) {

    # Extract the group_DF attribute to reduce typing below. Curse you S3
    # objects for not having a convenient way of extracting attributes!!
    groupDF <- attr(omicsData, "group_DF")

    # From the group_DF attribute extract the non-singleton group names.
    nonsingletons <- attr(groupDF, "nonsingleton_groups")

    # Check if any groups are singletons.
    # The following if statement removes any singleton samples from cur_edata.
    # These samples will not be part of the CV calculation and they will not be
    # filtered out of the omicsData object when the filter is applied.
    if (!setequal(unique(groupDF$Group), nonsingletons)) {

      # Give a warning that singleton groups will not be used to determine
      # which biomolecules will be filtered.
      warning (paste("Grouping information is being utilized when calculating",
                     "the CV, and there are group(s) consisting of a single",
                     "sample. The singleton group(s) will be ignored by this",
                     "filter.",
                     sep = " "))

      # Keep rows in groupDF that correspond to non-singleton groups.
      groupDF <- groupDF[which(groupDF$Group %in% nonsingletons), ]

      # Fish out the name of the column that contains the sample names.
      sID <- get_fdata_cname(omicsData)

      # Keep columns in cur_edata corresponding to non-singleton groups.
      cur_edata <- cur_edata[, which(names(cur_edata) %in% groupDF[, sID])]

    }

    # Make sure the order of the groups in group_DF matches the order of the
    # sample names in cur_edata.
    groupie <- groupDF$Group[match(names(cur_edata), groupDF$SampleID)]

    # Calculate the pooled CV. The data needs to be converted to a matrix and
    # the group names need to be converted to a character vector for
    # pooled_cv_rcpp to run properly.
    cvs <- pooled_cv_rcpp(as.matrix(cur_edata), as.character(groupie))

    # For mystifying reasons multiply the pooled CV values by 100.
    cvs <- cvs * 100

    # Set the is_pooled variable to TRUE (because we pooled the CV).
    is_pooled <- TRUE

  }

  # Create and add attributes to the cv filter data frame ----------------------

  # Create a data frame with the ID column from e_data and the CV values. This
  # data frame is called pool_cv even though the CV may not be pooled (this
  # makes us mysterious).
  pool_cv <- data.frame(omicsData$e_data[, id_col],
                        CV = cvs)
  names(pool_cv)[1] <- get_edata_cname(omicsData)

  ## determine plotting window cutoff ##
  # calculate percentage of observations with CV <= 200 #
  prct.less200 <- (sum(pool_cv$CV <= 200, na.rm = T) /
                     length(pool_cv$CV[!is.na(pool_cv$CV)]))

  if (prct.less200 > 0.95) {
    x.max = min(200, quantile(pool_cv$CV, 0.99, na.rm = TRUE))
  } else{
    x.max = quantile(pool_cv$CV, 0.95, na.rm = TRUE)
  }

  ## generate some summary stats for CV values, for PMART purposes only ##
  tot.nas <- sum(is.na(pool_cv$CV))

  output <- data.frame(pool_cv, row.names = NULL)

  orig_class <- class(output)

  class(output) <- c("cvFilt", orig_class)

  # Add the group designation information to the attributes.
  attr(output, "group_DF") <- attr(omicsData, "group_DF")

  attr(output, "pooled") <- is_pooled
  attr(output, "max_x_val") <- x.max
  attr(output, "tot_nas") <- tot.nas
  attr(output, "use_groups") <- ifelse(use_groups == FALSE,FALSE,TRUE)


  # Return the completed object. We did it!!!
  return (output)

}

#' RMD Runs
#'
#' The method computes a robust Mahalanobis distance that can be mapped to a
#'p-value and used to identify outlying samples
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'  'lipidData', or 'nmrData' usually created by \code{\link{as.pepData}},
#'  \code{\link{as.proData}}, \code{\link{as.metabData}},
#'  \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @param ignore_singleton_groups logical indicator of whether to remove
#'  singleton groups or not; defaults to TRUE. If TRUE, rmd_filter results are
#'  returned only for samples in groups of size greater than 1. This is used
#'  when calculating the correlation.
#'
#' @param metrics A character vector indicating which metrics should be used when
#'       calculating the robust Mahalanobis distance. This vector must contain
#'       between two and five of the following options: "MAD" (Median Absolute
#'       Deviation), "Kurtosis", "Skewness", "Correlation", and
#'       "Proportion_Missing". The default is NULL. When NULL a combination of
#'       metrics will be chosen depending on the class of omicsData.
#'
#' @return a data.frame containing columns for the Sample ID, log2 robust
#'  Mahalanobis distance, p-values, and robust Mahalanobis distance
#'
#' @details The metrics on which the log2 robust Mahalanobis distance is based
#'  can be specified using the \code{metrics} argument. \tabular{ll}{ pepData
#'  \tab For pepData objects, all five of the metrics "MAD", "Kurtosis",
#'  "Skewness", "Correlation", "Proportion_Missing" may be used (this is the
#'  default). \cr proData \tab For proData objects, all five of the metrics
#'  "MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing" may be
#'  used (this is the default). \cr metabData \tab For metabData objects, the
#'  use of "Proportion_Missing" is discouraged due to the general lack of
#'  missing data in metabolomics datasets (the default behavior omits
#'  "Proportion_Missing" from the metrics). \cr lipidData \tab For lipidData
#'  objects, , the use of "Proportion_Missing" is discouraged due to the general
#'  lack of missing data in metabolomics datasets (the default behavior omits
#'  "Proportion_Missing" from the metrics). \cr }
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' metab_object2 <- edata_transform(omicsData = metab_object, data_scale = "log2")
#' metab_object3 <- group_designation(omicsData = metab_object2, main_effects = "Condition")
#' rmd_results <- rmd_filter(omicsData = metab_object3, metrics=c("MAD", "Skewness", "Correlation"))
#' rmd_results <- rmd_filter(omicsData = metab_object2)
#'
#' data(pep_pepData)
#' pep_pepData2 <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' pep_pepData3 <- group_designation(omicsData = pep_pepData2, main_effects = "Condition")
#' rmd_results <- rmd_filter(omicsData = pep_pepData3)
#' }
#'
#' @references Matzke, M., Waters, K., Metz, T., Jacobs, J., Sims, A., Baric, R.,
#'  Pounds, J., and Webb-Robertson, B.J. (2011), \emph{Improved quality control
#'  processing of peptide-centric LC-MS proteomics data}. Bioinformatics.
#'  27(20): 2866-2872.
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
#'
#' @rdname rmd_filter
#' @name rmd_filter
#'
rmd_filter <- function (omicsData,
                        ignore_singleton_groups = TRUE,
                        metrics = NULL) {

  # Create functions to be used later ------------------------------------------

  # This function will be used in two different places later depending on the
  # metrics used in the input. It is here so changes only need to be made in one
  # location of the code.
  mal.fun <- function(x) (t(as.vector(x) - med.mat) %*%
                            solve(as.matrix(cov.mat)) %*%
                            (as.vector(x) - med.mat))

  # Run some preliminary checks ------------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {

    # Throw down an error like Zeus with his lightning bolts.
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))

  }

  # check that ignore_singleton_groups is logical #
  if (!is.logical(ignore_singleton_groups)) {

    # Stop the illogical user with an error.
    stop ("ignore_singleton_groups must be either TRUE or FALSE")

  }

  # group_DF attribute is required #
  if (is.null(attr(omicsData, "group_DF"))) {

    # Kindly tell the user they are wrong and point them in the right direction
    # to get some much needed help.
    stop (paste("omicsData must contain attribute information for 'group_DF'.",
                "See documentation for group_designation function for more",
                "information.",
                sep = " "))

  }

  # data should be log transformed #
  if (!attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")) {

    # Yet another error message for the user. Poor thing.
    stop (paste("omicsData$e_data should be log transformed prior to calling",
                "rmd_filter. See documentation for edata_transform function",
                "for more information.",
                sep = " "))

  }

  # Check the number of observations in e_data.
  if (attributes(omicsData)$data_info$num_edata < 50) {

    # Warn the user that skimping on the sample size is ALWAYS a bad idea!
    warning (paste("Use the results of the RMD filter with caution due to a",
                   "small number of biomolecules (<50).",
                   sep = " "))

  }

  # Check if metrics is NULL.
  if (is.null(metrics)) {

    # Take a looksie at the omicsData class.
    if (inherits(omicsData, c("pepData", "proData"))) {

      # Set default metrics values for peptide and protein data.
      metrics <- c("MAD", "Kurtosis", "Skewness", "Correlation",
                   "Proportion_Missing")

    } else if (inherits(omicsData, c("lipidData", "metabData", "nmrData"))) {

      # Set default metrics values for lipid and metabolite data.
      metrics <- c("MAD", "Kurtosis", "Skewness", "Correlation")

    }

    # Check they didn't mess up the metrics vector (if it isn't NULL).
  } else {

    # Make sure the entries in metrics are an acceptable type.
    if (!all(metrics %in% c("MAD", "mad", "m", "median absolute deviation",
                            "median_absolute_deviation",
                            "Kurtosis", "kurtosis", "kurt", "k",
                            "Skewness", "skew", "s", "skewness",
                            "Correlation", "corr", "c", "cor", "correlation",
                            "Proportion_Missing", "proportion_missing", "p",
                            "prop_missing", "prop_miss", "proportion_miss",
                            "prop", "proportion", "proportion missing",
                            "prop miss", "proportion miss"))) {

      # Stop the user for using unholy elements in the metrics vector.
      stop (paste("One or more of the elements in metrics does not match",
                  "the acceptable values. See details for metrics in the",
                  "rmd_filter documenation.",
                  sep = " "))

    }

    # Make sure there are not any repeated metrics.
    if (length(unique(metrics)) != length(metrics)) {

      # Stop the user for trying to use a metric more than once.
      stop ("One or more of the elements in metrics occurs more than once.")

    }

    # Make sure the metrics vector as an appropriate number of elements
    if (length(metrics) < 2) {

      # Just being Zeus again with the error throwing.
      stop ("The metrics vector must contain at least two elements.")

    } else if (length(metrics) > 5) {

      # Zeus is nothing compared to pmartR!! BWAHAHAHA!!!
      stop ("The metrics vector cannot contain more than five elements.")

    }

  }

  # Carry out preliminary group calculations -----------------------------------

  # Create an object for the group_DF attribute with a name that is almost as
  # long as using the call to the attr function.
  mintR_groupDF = attr(omicsData, "group_DF")

  ### Aug 2020: deal with any groups that have only a single sample in them ###
  if (any(table(mintR_groupDF$Group) == 1) && ignore_singleton_groups) {

    # which groups have a single sample in them #
    singleton_groups <- names(which(table(mintR_groupDF$Group) == 1))
    singleton_rows <- which(mintR_groupDF$Group %in% singleton_groups)

    # filter out samples corresponding to groups of size 1 #
    to_remove <- as.character(mintR_groupDF[singleton_rows,
                                            get_fdata_cname(omicsData)])
    myfilter <- custom_filter(omicsData, f_data_remove = to_remove)
    omicsData <- applyFilt(myfilter, omicsData)
    mintR_groupDF <- attr(omicsData, "group_DF")

  }

  # Compare the number of samples to the number of metrics used.
  if (dim(omicsData$f_data)[1] < 2*length(metrics)) {

    # Warn the user that little sample sizes are the bane of their existence.
    warning (paste("Use the results of the RMD filter with caution due to a",
                   "small number of samples relative to the number of metrics",
                   "being used. Consider reducing the number of metrics being",
                   "used.",
                   sep = " "))
  }

  # Compute the metrics used to calculate the RMD ------------------------------

  # Extricate the column number from e_data that contains the IDs.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))

  # Convert metrics to lower case for matching purposes.
  metrics <- tolower(metrics)

  # Initialize a vector for counting the number of metrics used.
  metrics_final <- rep(NA, 5)

  # Initialize a data frame with the sample ID column from f_data.
  rmd.vals <- data.frame(Sample.ID = names(omicsData$e_data[, -id_col]))

  # Compute the median absolute deviation across the samples (columns).
  if (any(metrics %in% c("mad", "m", "median absolute deviation",
                         "median_absolute_deviation"))){
    ind = which(metrics %in% c("mad", "m", "median absolute deviation",
                               "median_absolute_deviation"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'MAD'.")}
    metrics = metrics[-ind]
    metrics_final[1] = 1
    rmd.vals$MAD = run_mad(omicsData$e_data[, -id_col])[,2]
  }

  # Compute the kurtosis across the samples (columns).
  if(any(metrics %in% c("kurtosis", "kurt", "k"))){
    ind = which(metrics %in% c("kurtosis", "kurt", "k"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Kurtosis'.")}
    metrics = metrics[-ind]
    metrics_final[2] = 1
    rmd.vals$Kurtosis = run_kurtosis(omicsData$e_data[, -id_col])[,2]
  }

  # Compute the skewness across the samples (columns).
  if(any(metrics %in% c("skew", "s", "skewness"))){
    ind = which(metrics %in% c("skew", "s", "skewness"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Skewness'.")}
    metrics = metrics[-ind]
    metrics_final[3] = 1
    rmd.vals$Skewness = run_skewness(omicsData$e_data[, -id_col])[,2]
  }

  # Compute the correlation across the samples (columns).
  if(any(metrics %in% c("corr", "c", "cor", "correlation"))){
    ind = which(metrics %in% c("corr", "c", "cor", "correlation"))
    if(length(ind)>1){
      stop("More than one of the entries in metrics matches 'Correlation'.")}
    metrics = metrics[-ind]
    metrics_final[4] = 1
    rmd.vals$Corr = suppressWarnings(run_group_meancor(omicsData,
                                                       mintR_groupDF,
                                                       ignore_singleton_groups)[,2])
  }

  # Compute the proportion of missing values across the samples (columns).
  if (any(metrics %in% c("proportion_missing", "p", "prop_missing", "prop_miss",
                         "proportion_miss", "prop", "proportion",
                         "proportion missing", "prop miss", "proportion miss"))) {

    # check to see whether there is any missing data
    if (attr(omicsData, 'data_info')$num_miss_obs == 0) {

      # if no missing data, we flat out refuse to include prop_missing in the
      # metrics, so give a warning to this effect.
      warning (paste("There are no missing values in e_data, therefore",
                     "Proportion_Missing will not be used as one of the",
                     "metrics for RMD-Runs.",
                     sep = " "))

      # If there are missing values compute the proportion of missing values for
      # each sample (column) as usual.
    } else {

      ind = which(metrics %in% c("proportion_missing", "p", "prop_missing",
                                 "prop_miss", "proportion_miss", "prop",
                                 "proportion", "proportion missing", "prop miss",
                                 "proportion miss"))

      if(length(ind)>1){
        stop("More than one of the entries in metrics matches 'Proportion_Missing'.")}

      metrics = metrics[-ind]
      metrics_final[5] = 1
      rmd.vals$Proportion_Missing = run_prop_missing(omicsData$e_data[, -id_col])[,2]

      ## proceed, to check the rank of cov.mat ##

      ## Conduct Robust PCA ##
      robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1],
                                     k = (ncol(rmd.vals)-1),
                                     mcd = FALSE,
                                     scale = FALSE)

      # Check for errors when conducting robust PCA.
      if (inherits(robpca.res, "try-error")) {

        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "prop_missing as one of the metrics. Try again, excluding",
                    "Proportion_Missing from the metrics vector.",
                    sep = " "))
      }

      ## Calculate Covariance Matrix #
      cov.mat = (robpca.res@loadings %*%
                   diag(robpca.res@eigenvalues) %*%
                   t(robpca.res@loadings))

      # Extract the rank of the covariance matrix.
      myrank = qr(cov.mat)$rank

      # Check the rank of the covariance matrix compared to the number of
      # metrics used to calculate it.
      if (myrank != sum(metrics_final, na.rm = TRUE)) {

        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "prop_missing as one of the metrics. Try again, excluding",
                    "Proportion_Missing from the metrics vector.",
                    sep = " "))

      }

      ## Calculate Robust Mahalanobis Distance ##
      med.mat = matrix(apply(rmd.vals[,-1], 2, median),
                       nrow = (ncol(rmd.vals)-1))

      rob.dist.vals = try (apply(rmd.vals[,-1], 1, mal.fun))

      # Check if there was an error when calculating the the RMD.
      if (inherits(rob.dist.vals, "try-error")) {

        # Throw an error because there is not enough missing data. Good problem
        # to have I guess.
        stop (paste("There are not enough missing values in e_data to use",
                    "Proportion_Missing as one of the metrics. Try again,",
                    "excluding Proportion_Missing from the metrics vector.",
                    sep = " "))

      }

    }

  }

  # Check the number of metrics used. This value could be less than two if only
  # two metrics were input and one of them was Proportion_Missing. For example,
  # if Proportion_Missing was removed because there were no missing values then
  # there would only be one metric left and a PCA and RMD cannot be calculated
  # for only one metric.
  if (sum(metrics_final, na.rm = TRUE) < 2) {

    # Throw an error if only two metrics were given, one of them was
    # Proportion_Missing, and there is no missing data. Give the user a helpful
    # suggestion. We are quite nice :)
    stop (paste("Vector of metrics must contain at least two valid entries.",
                "Try including a metric other than Proportion_Missing.",
                sep = " "))

  }

  # Extract the names of the metrics that have been used. We cannot simply report
  # the metrics vector used as the input because it is possible that the input
  # metrics vector has had elements removed. For example, if the input data does
  # not have any missing values, and Proportion_Missing was an input metric, then
  # Proportion_Missing will not be calculated. Therefore, this metric will be
  # removed from the metrics_final_txt vector.
  metrics_final_txt <- names(rmd.vals[,-1])

  # RMD the heck out of the data -----------------------------------------------

  # Check if robust PCA, covariance matrix, and RMD objects have been created
  if (!exists("robpca.res")) {

    ## Conduct Robust PCA ##
    robpca.res = rrcov:::PcaHubert(x = rmd.vals[,-1],
                                   k = (ncol(rmd.vals)-1),
                                   mcd = FALSE,
                                   scale = FALSE)

    ## Calculate Covariance Matrix #
    cov.mat = (robpca.res@loadings %*%
                 diag(robpca.res@eigenvalues) %*%
                 t(robpca.res@loadings))

    ## Calculate Robust Mahalanobis Distance ##
    med.mat = matrix(apply(rmd.vals[,-1], 2, median),
                     nrow = (ncol(rmd.vals)-1))

    rob.dist.vals = apply(rmd.vals[,-1], 1, mal.fun)

  }

  log2.dist.vals = log(rob.dist.vals, base = 2)

  rmd.pvals = 1 - pchisq(rob.dist.vals, df = (ncol(rmd.vals)-1))

  temp.res = data.frame(Sample.ID = rmd.vals[,1],
                        Log2.md = log2.dist.vals,
                        pvalue = rmd.pvals,
                        rmd.vals[,-1])

  names(temp.res)[1] = get_fdata_cname(omicsData)

  output = merge(x = mintR_groupDF,
                 y = temp.res,
                 by = get_fdata_cname(omicsData),
                 all = TRUE)

  # The order of the rows in output can be scrambled from the order of the rows
  # in either mintR_groupDF or temp.res because of the merge function. The
  # following lines will reorder the rows of the output data frame to match the
  # order of the samples from the group_DF attribute. The merge function didn't
  # behave this way until I made up some long and silly sample names when
  # creating a VizSampNames attribute for the RMD filter object. It lead to some
  # plots that didn't make sense. Very frustrating!
  # When matching the sample names I hard coded the sample ID column from output
  # because how it is created (from mintR_groupDF and temp.res) the sample ID
  # column will always be first. ... I should say that in every scenario I have
  # tested the sample ID column is first.
  output <- output[match(names(omicsData$e_data[, -id_col]), output[, 1]), ]

  orig_class <- class(output)

  class(output) <- c("rmdFilt", orig_class)

  # Add the sample names in e_data to attributes.
  attr(output, "sample_names") <- names(omicsData$e_data[, -id_col])

  # If there are custom sample names in f_data add them as an attribute.
  if ("VizSampNames" %in% names(omicsData$f_data)) {

    # Nab the order of the samples in f_data in relation to the order of the
    # sample names in e_data which are the column names of e_data minus the
    # column containing the biomolecule IDs.
    oder <- match(names(omicsData$e_data[, -id_col]),
                  omicsData$f_data[[get_fdata_cname(omicsData)]])

    # Add the custom sample names as its own attribute. This attribute can/will
    # be used in the plot.rmdFilt function. The sample names of this attribute
    # will be in the same order as the sample names in e_data (column names of
    # e_data).
    attr(output, "VizSampNames") <- omicsData$f_data$VizSampNames[oder]

  }

  # Add the group designation information to the attributes.
  attr(output, "group_DF") <- attr(omicsData, "group_DF")

  # Save fdata as an attribute. This is used in the summary.rmdFilt method.
  attr(output, "fdata") <- omicsData$f_data

  # We also need the name of the sample ID column in fdata in the
  # summary.rmdFilt method.
  attr(output, "fdata_cname") <- get_fdata_cname(omicsData)

  # Sum the number of metrics actually calculated.
  attr(output, "df") <- sum(!is.na(metrics_final))

  # Report the metrics actually used.
  attr(output, "metrics") <- metrics_final_txt

  return(output)

}

#' Calculate the Fraction of Missing Data of Sample Runs
#'
#' This function calculates the fraction of missing data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number
#'   of peptides and \eqn{n} is the number of samples.
#'
#' @return data.frame with two elements: Sample, a character vector giving the
#'   sample names; and Prop_missing, a numeric vector giving the fraction of
#'   missing values per run
#'
#' @author Lisa Bramer
#'
run_prop_missing <- function(data_only){

  # calculate number of missing values per run #
  nummiss <- apply(is.na(data_only), 2, sum)

  # calculate the fraction of missing values per run #
  fracmiss <- nummiss/nrow(data_only)

  # store data #
  res.final <- data.frame(Sample = names(data_only), Prop_missing = fracmiss, row.names = NULL)

  return(res.final)
}

#' Calculate the Median Absolute Deviance (MAD) of Sample Runs
#'
#' This function calculates the median absolute deviance across data for each
#' sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number
#'   of peptides and \eqn{n} is the number of samples.
#'
#' @details When calculating the MAD within a sample NA values are ignored. If
#'   all peptide abundance values are missing within a sample, the MAD is
#'   replaced by the overall mean MAD values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the
#'   sample names; and MAD, a numeric vector giving the MAD values
#'
#' @author Lisa Bramer
#'
run_mad <- function(data_only){

  # calculate MAD #
  mad_val = apply(data_only, 2, function(x) median(abs(x - median(x, na.rm = T)), na.rm = T))

  # calculate the number of samples with MAD equal to NA #
  num.miss <- sum(is.na(mad_val))

  # if at least one sample has a MAD of NA, replace it with mean MAD value #
  if(num.miss > 0){
    mad_val[is.na(mad_val)] <- mean(mad_val, na.rm = T)
  }

  # store data #
  res.final <- data.frame(Sample = names(data_only), MAD = mad_val, row.names = NULL)

  return(res.final)
}

#' Calculate the Kurtosis of Sample Runs
#'
#' This function calculates the kurtosis across data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number
#'   of peptides and \eqn{n} is the number of samples.
#'
#' @details Kurtosis is calculated by method 2 in the \code{e1071} package,
#'   which is unbiased under normality. Within a sample NA values are ignorned
#'   in the kurtosis calculation. If all peptide abundance values are missing
#'   within a sample, the kurtosis is replaced by the overall mean of nonmissing
#'   kurtosis values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the
#'   sample names; and Kurtosis, a numeric vector giving the kurtosis
#'
#' @author Lisa Bramer
#'
run_kurtosis <- function(data_only){

  # calculate kurtosis #
  kurt_res <- apply(data_only, 2, e1071::kurtosis, na.rm = TRUE, type = 2)

  # calculate the number of samples with kurtosis equal to NA #
  num.miss <- sum(is.na(kurt_res))

  # if at least one sample has a kurtosis of NA, replace it with mean kurtosis #
  if(num.miss > 0){
    kurt_res[is.na(kurt_res)] <- mean(kurt_res, na.rm = T)
  }

  # store data #
  res.final <- data.frame(Sample = names(data_only), Kurtosis = kurt_res, row.names = NULL)

  return(res.final)
}

#' Calculate the Skewness of Sample Runs
#'
#' This function calculates the skewness across data for each sample run.
#'
#' @param data_only a \eqn{p \times n} data.frame, where \eqn{p} is the number
#'   of peptides and \eqn{n} is the number of samples.
#'
#' @details Skewness is calculated as a bias-corrected calculation given by
#'   method 2 in the \code{e1071} package. Within a sample NA values are
#'   ignorned in the skewness calculation. If all peptide abundance values are
#'   missing within a sample, the skewness is replaced by the overall mean of
#'   nonmissing skewness values for the data.
#'
#' @return data.frame with two elements: Sample, a character vector giving the
#'   sample names; and Skewness, a numeric vector giving the skewness values
#'
#' @author Lisa Bramer
#'
run_skewness <- function(data_only){

  # calculate skewness #
  skew_res <- apply(data_only, 2, e1071::skewness, na.rm = TRUE, type = 2)

  # calculate the number of samples with skewness equal to NA #
  num.miss <- sum(is.na(skew_res))

  # if at least one sample has a skewness of NA, replace it with mean skewness #
  if(num.miss > 0){
    skew_res[is.na(skew_res)] <- mean(skew_res, na.rm = T)
  }

  # store data #
  res.final <- data.frame(Sample = names(data_only), Skewness = skew_res, row.names = NULL)

  return(res.final)
}

#' Calculate the Mean Correlation of a Sample with Respect to Group
#'
#' This function calculates the mean correlation of a sample with all other
#' samples that have the same group membership
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or
#'   'lipidData' usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}}, or
#'   \code{\link{as.lipidData}}, respectively.
#' @param mintR_groupDF data.frame created by \code{\link{group_designation}}
#'   with columns for sample.id and group.
#' @param use_singletons logical indicator of whether to include singleton
#'   groups or not; defaults to FALSE. If FALSE, rmd_filter results are returned
#'   only for samples in groups of size greater than 1. This is a pass-through
#'   argument from \code{\link{rmd_filter}}.
#'
#' @details Correlation calculations use only complete pairwise observations.
#'
#' @return data.frame with two elements: Sample.ID, a character vector giving
#'   the sample names; and Mean_Correlation, a numeric vector giving the mean
#'   correlation values
#'
#' @author Lisa Bramer
#'
run_group_meancor <- function(omicsData, mintR_groupDF, ignore_singleton_groups = TRUE){

  # if group.data has column for TimeCourse, re-compute group.data including TimeCourse as main effect
  # if mintR_groupDF has TimeCourse variable, re-compute group_data including TimeCourse as main effect
  if(!is.null(attr(mintR_groupDF, "time_course"))){
    if(!is.null(attr(mintR_groupDF, "main_effects"))){
      temp_maineff = c(attr(mintR_groupDF, "main_effects"), attr(mintR_groupDF, "time_course"))}else{
        temp_maineff = attr(mintR_groupDF, "time_course")
      }
    group_data = group_designation(omicsData, temp_maineff)
  }

  # create a list of unique groups #
  grps = unique(as.character(mintR_groupDF$Group))

  samp_id = attr(omicsData, "cnames")$fdata_cname

  # compute all pairwise correlations between samples of the same group #
  nonsingleton_groups <- attributes(mintR_groupDF)$nonsingleton_groups
  singleton_groups <- setdiff(grps, nonsingleton_groups)
  if(length(singleton_groups) > 0 & ignore_singleton_groups == FALSE){
    ## there are singelton groups and we don't want to ignore them in creating the rmd filter ##

    # calculate the correlation of the singleton sample to all other samples and
    # average that value to get the correlation for that particular sample/group
    omicsData_singletons <- omicsData
    omicsData_singletons$f_data$Dummy <- "dummy" # create a dummy grouping variable so that all samples belong to same group
    omicsData_singletons <- group_designation(omicsData_singletons, main_effect = "Dummy")

    prwse.grp.cors.all <- cor(omicsData_singletons$e_data[, -which(names(omicsData_singletons$e_data) == get_edata_cname(omicsData_singletons))], use = "pairwise.complete.obs")

    prwse.grp.cors.singletons <- list()
    for(i in 1:length(singleton_groups)){
      prwse.grp.cors.singletons[[i]] <- prwse.grp.cors.all
    }

    # do the "usual" thing for the nonsingleton groups:
    # calculate the mean correlation of a sample with all other samples that have the same group membership
    myfilt <- custom_filter(omicsData, f_data_keep = as.character(mintR_groupDF[which(mintR_groupDF$Group %in% nonsingleton_groups), samp_id]))
    omicsData_nonsingletons <- applyFilt(myfilt, omicsData)

    # make a list of which columns belong to which groups #
    grp.col.ids = list()
    for(i in 1:length(nonsingleton_groups)){

      # pull sample names from group.data in current group #
      nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == nonsingleton_groups[i]), samp_id])

      # pull column numbers corresponding to above names #
      grp.col.ids[[i]] = which(names(omicsData_nonsingletons$e_data) %in% nms)
    }

    prwse.grp.cors.nonsingletons <- lapply(grp.col.ids, function(x){
      cor(omicsData_nonsingletons$e_data[,x], use = "pairwise.complete.obs")
    })


    ## combine the singleton and nonsingleton results...make sure the order is correct
    prwse.grp.cors <- c(prwse.grp.cors.nonsingletons, prwse.grp.cors.singletons)
    prws.grp.cors.grpnames <- c(nonsingleton_groups, singleton_groups)

  }else{

    # make a list of which columns belong to which groups #
    grp.col.ids = list()
    for(i in 1:length(grps)){

      # pull sample names from group.data in current group #
      nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == grps[i]), samp_id])

      # pull column numbers corresponding to above names #
      grp.col.ids[[i]] = which(names(omicsData$e_data) %in% nms)
    }

    prwse.grp.cors = lapply(grp.col.ids, function(x){
      cor(omicsData$e_data[,x], use = "pairwise.complete.obs")
    })
    # structure of prws.grp.cors: list w/number elements equal to number groups
    # each element is correlation matrix
    # [[1]]
    # Mock1     Mock2     Mock3
    # Mock1 1.0000000 0.9670513 0.9703483
    # Mock2 0.9670513 1.0000000 0.9858240
    # Mock3 0.9703483 0.9858240 1.0000000
    #
    # [[2]]
    # Infection1 Infection2 Infection3 Infection4 Infection5
    # Infection1  1.0000000  0.9802595  0.9757083  0.9790413  0.9732911
    # Infection2  0.9802595  1.0000000  0.9837869  0.9804636  0.9824165
    # Infection3  0.9757083  0.9837869  1.0000000  0.9786230  0.9801930
    # Infection4  0.9790413  0.9804636  0.9786230  1.0000000  0.9784125
    # Infection5  0.9732911  0.9824165  0.9801930  0.9784125  1.0000000
  }

  # turn diagonal into NAs, so we don't include a sample's correlation with itself #
  grp.cors = lapply(prwse.grp.cors, function(x) x*(matrix(1, nrow = nrow(x), ncol = ncol(x)) + diag(NA, nrow(x))))

  # compute mean correlation for each sample #
  mean.cor = lapply(grp.cors, function(x) apply(x, 1, mean, na.rm = T))

  ## need to adjust the list elements of mean.cor for any singleton groups,
  ## to just contain the value for the sample in that group
  if(length(singleton_groups) > 0 & ignore_singleton_groups == FALSE){
    mean.cor2 <- mean.cor
    # when I concatenated the pairwise group correlations, I put the singleton groups last #
    for(i in 1:length(singleton_groups)){
      # get the sample name in the current singleton group, and pull that value out of mean.cor #
      cur_singleton <- singleton_groups[i]
      cur_sample <- as.character(omicsData_singletons$f_data[which(mintR_groupDF$Group == cur_singleton), samp_id])

      # get the element in mean.cor list that corresponds to this singleton group #
      j <- which(prws.grp.cors.grpnames == cur_singleton)
      mean.cor2[[j]] <- mean.cor[[j]][cur_sample]
      # pull sample names from group.data in current group #
      # nms = as.character(mintR_groupDF[which(mintR_groupDF$Group == nonsingleton_groups[i]), samp_id])

    }
  }else{
    mean.cor2 <- mean.cor
  }

  # format results #
  # get order to put results in original sample order based on peptide.data #
  temp = match(names(omicsData$e_data)[-1],names(unlist(mean.cor2)))
  res.cor = data.frame(Sample.ID = names(omicsData$e_data)[-1], Mean_Correlation = unlist(mean.cor2)[temp], row.names = NULL)

  return(res.cor)
}

#' Proteomics filter object
#'
#' This function counts the number of peptides that map to each protein and/or
#' the number of proteins to which each individual peptide maps.
#'
#' @param omicsData an object of class "pepData", the a result of
#'   \code{\link{as.pepData}}. The e_meta component of omicsData must be
#'   nonempty.
#'
#' @return A list with two elements. The first element is a data frame of counts
#'   for each unique peptide. The second element is also a data frame. This data
#'   frame contains the counts for the number of peptides that map to each
#'   unique protein.
#'
#' @examples
#' \dontrun{
#' library(pmartR)
#' data("pep_object")
#' my_filter <- proteomics_filter(omicsData = pep_object)
#' summary(my_filter, min_num_peps = 3)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export
#'
proteomics_filter <- function (omicsData) {

  # Preliminary checks and setup -----------------------------------------------

  # check that omicsData is of class 'pepData' #
  if (!inherits(omicsData, "pepData")) {

    # Let the user know that if the data does not contain peptides or proteins
    # you CANNOT APPLY A FILTER USING PEPTIDES AND PROTEINS. What is going on
    # upstairs!?!
    stop ("omicsData must be of class 'pepData'")

  }

  # check that e_meta is not NULL #
  if(is.null(omicsData$e_meta)) stop("e_meta must be non-NULL")

  # get peptide and protein column names #
  pep_id = attr(omicsData, "cnames")$edata_cname
  pro_id = attr(omicsData, "cnames")$emeta_cname

  # check that peptide and protein column names are non-null #
  if(is.null(pep_id)) stop("Peptide column name is NULL")
  if(is.null(pro_id)) stop("Protein column name is NULL")

  # Count peptides and proteins ------------------------------------------------

  # Count the number of rows for each peptide in e_meta.
  pepCount <- omicsData$e_meta %>%
    # Group the e_meta data frame by the peptide identifiers in the ID column.
    dplyr::group_by(.data[[pep_id]]) %>%
    # Count the number of times each peptide identifier occurs (the number of
    # rows an identifier appears in).
    dplyr::tally() %>%
    # Convert pepCount to a data frame from a tibble.
    data.frame()

  # Count the number of peptides associated with each protein.
  proCount <- omicsData$e_meta %>%
    # Group the e_meta data frame by the protein identifiers in the protein ID
    # column.
    dplyr::group_by(.data[[pro_id]]) %>%
    # Count the number of times each protein identifier occurs (the number of
    # rows an identifier appears in).
    dplyr::tally() %>%
    # Convert proCount to a data frame from a tibble.
    data.frame()

  # Generate a list containing the data frames for the two counts.
  output <- list(counts_by_pep = pepCount,
                 counts_by_pro = proCount)

  # Preserve the list class.
  orig_class <- class(output)

  # Incorporate the proteomicsFilt class along with the list class.
  class(output) <- c("proteomicsFilt", orig_class)

  # Return the list of data frames containing peptide and protein counts.
  # We can count!!
  return(output)

}

#' IMD-ANOVA filter object
#'
#' This function returns an imdanovaFilt object for use with
#'\code{\link{applyFilt}}
#'
#' @param omicsData object of one of the classes "pepData", "isobaricpepData",
#'   "proData", "lipidData", "metabData", or "nmrData", usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}},
#'   \code{\link{as.proData}}, \code{\link{as.lipidData}},
#'   \code{\link{as.metabData}}, or \code{\link{as.nmrData}}, respectively.
#'   Groups (more than one group) must have been specified using the
#'   \code{\link{group_designation}} function prior to using the imdanova_filter
#'   function.
#'
#' @details The output from this function can be used in conjunction with
#'   \code{\link{applyFilt}} to filter out molecules that are not present in
#'   enough samples to do statistical comparisons. If any singleton groups are
#'   present in the omicsData object, those groups are not part of the filter
#'   object that is returned.
#'
#' @return Object of class imdanovaFilt (also a data.frame) containing the
#'   molecule identifier and number of samples in each group with non-missing
#'   values for that molecule.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_pepData2 <- group_designation(omicsData = pep_object,
#'                                   main_effects = "Condition")
#' to_filter <- imdanova_filter(omicsData = pep_pepData2)
#' summary(to_filter, min_nonmiss_anova = 2)
#' }
#'
#' @author Kelly Stratton
#'
#' @export
#'
imdanova_filter <- function (omicsData) {

  # Run some preliminary checks ------------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {

    # Follow the instructions foul creature!!!
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))

  }

  # group_DF attribute is required #
  if (is.null(attr(omicsData, "group_DF"))) {

    # Kindly tell the user they are wrong and point them in the right direction
    # to get some much needed help.
    stop (paste("omicsData must contain attribute information for 'group_DF'.",
                "See documentation for group_designation function for more",
                "information.",
                sep = " "))

  }

  # Check the number of groups. There must be more than one group if the data
  # are not paired.
  if (length(names(get_group_table(omicsData))) < 2 &&
      is.null(attr(attr(omicsData, "group_DF"), "pair_id"))) {

    # Let the user know that they cannot compare statistics between groups if
    # there is only one group!!
    stop (paste("There must be more than one group in order to create an",
                "imdanovaFilt object if the data are not paired.",
                sep = " "))

  }

  # Determine number of samples per group --------------------------------------

  # Extract the group_DF attribute.
  groupDF <- attr(omicsData, "group_DF")

  # if groupDF has column for TimeCourse, do the following within time point
  # (do not re-compute groupDF including TimeCourse as main effect--this was our
  # old strategy but Bobbie-Jo directed us to not do this, and instead loop
  # through time points)
  if (any(names(groupDF) == "TimeCourse")) {

    # added by KGS 9/4/2020 since we have disabled TimeCourse functionality
    stop (paste("Option for TimeCourse in group_designation is not currently",
                "supported.",
                sep = " "))
    #     filt.edata <- vector(mode = "list", length =
    #     length(unique(groupDF$TimeCourse))) names(filt.edata) <-
    #     unique(groupDF$TimeCourse)
    #
    #     for(tind in 1:length(unique(groupDF$TimeCourse))){
    #
    #     t = unique(groupDF$TimeCourse)[tind] t.e_data <- cbind(e_data[,1],
    #     e_data[, names(e_data) %in%
    #     as.character(groupDF[,samp_id][groupDF$TimeCourse==t])])
    #     names(t.e_data)[1] <- names(e_data)[1]
    #
    #     t.groupDF <- groupDF[groupDF$TimeCourse==t, ] #all(names(t.e_data)[-1] ==
    #     t.groupDF$sampleID) # just checking, should be TRUE
    #
    #     nonmiss_per_group <- nonmissing_per_group(omicsData=NULL, e_data=t.e_data,
    #     groupDF=t.groupDF, cname_id=edata_id, samp_id=samp_id)
    #     if(filter_method=="anova"){ filt.edata[[tind]] <-
    #     anova_filter(nonmiss_per_group=nonmiss_per_group,
    #     min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_id) }else{
    #     if(filter_method=="gtest"){ filt.edata[[tind]] <-
    #     gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=t.groupDF,
    #     e_data=t.e_data, alpha=alpha, min_nonmiss_gtest=min_nonmiss_gtest,
    #     cname_id = edata_id, samp_id = samp_id) }else{
    #     if(filter_method=="combined"){ filt.edata.gtest <-
    #     gtest_filter(nonmiss_per_group, groupDF=t.groupDF, e_data=t.e_data,
    #     alpha=alpha, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_id) #
    #     min.nonmiss.allowed <- 2 filt.edata.anova <-
    #     anova_filter(nonmiss_per_group, min_nonmiss_anova, cname_id = edata_id)
    #     filt.edata[[tind]] <- intersect(filt.edata.anova, filt.edata.gtest) } } }
    #     }
    #
    #     filter.edata <- Reduce(base::intersect, filt.edata)

  } else { # end of if-statement for the presence of TimeCourse variable

    # Count the number of non-missing values for all groups. For example, if
    # group A has 5 samples and 4 of the samples have missing values then the
    # count for group A will be 1.
    nonmiss_per_group <- nonmissing_per_group(omicsData = omicsData)

    # Extract the data frame that contains a column for the biomolecule IDs and
    # columns for the counts of nonmissing values for each group.
    output <- nonmiss_per_group$nonmiss_totals

  } # end of else-stament for the absence of TimeCourse variable

  # remove columns of output that correspond to any singleton groups present #
  singleton_groups <- setdiff(unique(groupDF$Group),
                              attr(groupDF, "nonsingleton_groups"))

  # Check for the presence of singleton groups.
  if (length(singleton_groups) > 0) {

    # Remove the columns of the output data frame that correspond to singleton
    # groups (there is only one sample for that particular group).
    output <- output[, -which(names(output) %in% singleton_groups)]

  }

  orig_class <- class(output)
  class(output) <- c("imdanovaFilt", orig_class)

  # Save the entire omicsData object as an attribute. This will be used in the
  # summary.imdanovaFilt method when data are paired.
  attr(output, "omicsData") <- omicsData

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes
  # KS added attribute for nonsingleton groups 12/3/2020 #
  attr(output, "nonsingleton_groups") <-
    nonmiss_per_group$group_sizes$Group[which(
      nonmiss_per_group$group_sizes$n_group > 1
    )]

  return (output)

}

#' Custom Filter
#'
#' This function creates a customFilt S3 object based on user-specified items to
#' filter out of the dataset
#'
#' @param omicsData an object of class "pepData", "proData", "metabData",
#'   "lipidData", or "nmrData, created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#'
#' @param e_data_remove character vector specifying the names of the e_data
#'   identifiers to remove from the data. This argument can only be specified
#'   with other 'remove' arguments.
#'
#' @param f_data_remove character vector specifying the names of f_data
#'   identifiers to remove from the data. This argument can only be specified
#'   with other 'remove' arguments.
#'
#' @param e_meta_remove character vector specifying the names of the e_meta
#'   identifiers to remove from the data. This argument can only be specified
#'   with other 'remove' arguments.
#'
#' @param e_data_keep character vector specifying the names of the e_data
#'   identifiers to keep from the data. This argument can only be specified with
#'   other 'keep' arguments.
#'
#' @param f_data_keep character vector specifying the names of f_data
#'   identifiers to keep from the data. This argument can only be specified with
#'   other 'keep' arguments.
#'
#' @param e_meta_keep character vector specifying the names of the e_meta
#'   identifiers to keep from the data. This argument can only be specified with
#'   other 'keep' arguments.
#'
#' @return An S3 object of class 'customFilt', which is a list with 3 elements:
#'   e_data_remove, f_data_remove, and e_meta_remove.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("metab_object")
#' to_filter <- custom_filter(metab_object, e_data_remove = "fumaric acid",
#'                            f_data_remove = "Infection1")
#' summary(to_filter)
#' to_filter2 <- custom_filter(metab_object, e_data_remove = "fumaric acid")
#' summary(to_filter2)
#' }
#'
#' @author Kelly Stratton
#'
#' @export
#'
custom_filter <- function (omicsData,
                           e_data_remove = NULL,
                           f_data_remove = NULL,
                           e_meta_remove = NULL,
                           e_data_keep = NULL,
                           f_data_keep = NULL,
                           e_meta_keep = NULL ) {

  # Run some preliminary checks ------------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {

    # Follow the instructions foul creature!!!
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = " "))

  }

  # check that not all remove and keep arguments are NULL.
  if (is.null(c(e_data_remove, f_data_remove, e_meta_remove,
                e_data_keep, f_data_keep, e_meta_keep))) {

    # Stop the user because they are trying to take the long way to the exact
    # same data set.
    stop ("No items have been identified for filtering.")

  }

  #check that both keep and remove arguments are not non-NULL
  if (!is.null(c(e_data_remove, f_data_remove, e_meta_remove)) &&
      !is.null(c(e_data_keep, f_data_keep, e_meta_keep))) {

    stop (paste("Cannot have both remove arguments and keep arguments",
                "be non-NULL. Create separate filter objects for the",
                "remove arguments and keep arguments.",
                sep = " "))

  }

  # Extricate the names of the different ID columns.
  edata_id = attr(omicsData, "cnames")$edata_cname
  emeta_id = attr(omicsData, "cnames")$emeta_cname
  samp_id = attr(omicsData, "cnames")$fdata_cname

  # Create filter_object for remove arguments ----------------------------------

  if (!is.null(c(e_data_remove, f_data_remove, e_meta_remove))) {

    # checks for e_data_remove #
    if (!is.null(e_data_remove)) {

      # check that e_data_remove are all in omicsData #
      if (!(all(e_data_remove %in% omicsData$e_data[, edata_id]))) {

        # Throw an error because the user tried to filter nonexistent
        # biomolecules.
        stop ("Not all of the items in e_data_remove are found in e_data.")

      }

      # check that e_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_data[, edata_id] %in% e_data_remove)) {

        # Throw an error because the user tried to filter every biomolecule.
        stop("e_data_remove specifies all the items in the data.")

      }

    }

    # checks for f_data_remove #
    if (!is.null(f_data_remove)) {

      # check that f_data_remove are all in omicsData #
      if (!(all(f_data_remove %in% omicsData$f_data[, samp_id]))) {

        # Throw an error because the user tried to filter nonexistent samples.
        stop ("Not all of the items in f_data_remove are found in f_data.")

      }

      # check that f_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$f_data[, samp_id] %in% f_data_remove)) {

        # Throw an error because the user tried to filter every sample.
        stop ("f_data_remove specifies all the items in f_data.")

      }

    }

    # checks for e_meta_remove #
    if (!is.null(e_meta_remove)) {

      # check that e_meta_remove are all in omicsData #
      if (!(all(e_meta_remove %in% omicsData$e_meta[, emeta_id]))) {

        # Throw an error because the user tried to filter nonexistent mapping
        # variables.
        stop ("Not all of the items in e_meta_remove are found in e_meta.")

      }

      # check that e_meta_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_meta[, emeta_id] %in% e_meta_remove)) {

        # Throw an error because the user tried to filter every single mapping
        # variable. They must be tired of their work and they want to be done
        # as soon as possible. Not having any data to analyze is a very
        # effective way of finishing quickly.
        stop ("e_meta_remove specifies all the items in e_meta.")

      }
    }

    # Fashion a filter object with the remove elements.
    filter_object <- list(e_data_remove = e_data_remove,
                          f_data_remove = f_data_remove,
                          e_meta_remove = e_meta_remove)

  }

  # Create filter_object for keep arguments ----------------------------------

  if (!is.null(c(e_data_keep, f_data_keep, e_meta_keep))) {

    # checks for e_data_keep #
    if(!is.null(e_data_keep)){

      # check that e_data_keep are all in omicsData #
      if (!(all(e_data_keep %in% omicsData$e_data[, edata_id]))) {

        # Stop the greedy user from trying to keep things that aren't theirs!!
        stop ("Not all of the items in e_data_keep are found in e_data.")

      }

      # check that e_data_keep doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_data[, edata_id] %in% e_data_keep)) {

        # Let the greedy user know that they are keeping everything. Why run the
        # filter in the first place!?!
        stop ("e_data_keep specifies all the items in e_data.")

      }

    }

    # checks for f_data_keep #
    if (!is.null(f_data_keep)) {

      # check that f_data_keep are all in omicsData #
      if (!(all(f_data_keep %in% omicsData$f_data[, samp_id]))) {

        # Stop the greedy user from trying to keep samples that aren't theirs!
        stop ("Not all of the items in f_data_keep are found in f_data.")

      }

      # check that f_data_remove doesn't specify ALL the items in omicsData #
      if (all(omicsData$f_data[, samp_id] %in% f_data_keep)) {

        # Stop the greedy user from keeping all of the samples. Why filter?!
        stop ("f_data_keep specifies all the items in f_data.")

      }

    }

    # checks for e_meta_remove #
    if (!is.null(e_meta_keep)) {

      # check that e_meta_keep are all in omicsData #
      if (!(all(e_meta_keep %in% omicsData$e_meta[, emeta_id]))) {

        # Stop the greedy user from trying to keep mapping variables that do not
        # belong to them.
        stop ("Not all of the items in e_meta_keep are found in e_meta.")

      }

      # check that e_meta_keep doesn't specify ALL the items in omicsData #
      if (all(omicsData$e_meta[, emeta_id] %in% e_meta_keep)) {

        # Stop the greedy user from keeping ALL of the mapping variables.
        stop ("e_meta_keep specifies all the items in e_meta")

      }

    }

    # Manufacture a filter object with the keep elements.
    filter_object <- list(e_data_keep = e_data_keep,
                          f_data_keep = f_data_keep,
                          e_meta_keep = e_meta_keep)

  }


  # Set filter_object class and attributes -------------------------------------

  class(filter_object) <- c("customFilt", "list")

  # Save the counts of the biomolecules, samples, and mapping variables.
  attr(filter_object,
       "num_samples") <- length(unique(omicsData$f_data[, samp_id]))
  attr(filter_object,
       "num_edata") <-  length(unique(omicsData$e_data[, edata_id]))
  attr(filter_object,
       "num_emeta") <- if (!is.null(emeta_id)) {

         length(unique(omicsData$e_meta[, emeta_id]))

       }

  # Save the ID column names. These attributes are used in the
  # summary.customFilt method.
  attr(filter_object, "cnames")$edata_cname = edata_id
  attr(filter_object, "cnames")$emeta_cname = emeta_id
  attr(filter_object, "cnames")$fdata_cname = samp_id

  # Save the entire omicsData object as an attribute. This is used in the
  # summary.customFilt method.
  attr(filter_object, "omicsData") = omicsData # added 12/5/2017 by KS #

  # Add the group designation information to the attributes.
  attr(filter_object, "group_DF") <- attr(omicsData, "group_DF")

  # Return the customated filter object. Good on us!!!
  return (filter_object)

}
