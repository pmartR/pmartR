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
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_object2 <- group_designation(omicsData = pep_object, main_effects = "Condition")
#' to_filter <- cv_filter(omicsData = pep_object2, use_groups = TRUE)
#' summary(to_filter, cv_threshold = 30)
#'}
#'
#'@author Lisa Bramer, Kelly Stratton
#'
#'@export

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
  if (use_groups==TRUE && is.null(attr(omicsData, "group_DF"))) {
    
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
    # calculations.
    cur_edata <- edata_transform(omicsData = omicsData,
                                 data_scale = get_data_scale(omicsData))$e_data
    
  }
  
  # Conduct CV calculations ----------------------------------------------------
  
  # Calculate the unpooled CV if group_designation has not been run or if
  # use_groups = FALSE.
  if (use_groups == FALSE || is.null(attr(omicsData, "group_DF"))) {
    
    # Calculate the unpooled coefficient of variation.
    cvs <- apply(cur_edata, 1, cv.calc)
    
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
    
    # added by KGS Sept 2020: filters out any samples corresponding to singleton groups before proceeding with CV filter
    # the samples themselves won't be filtered out of the omicsData object upon application of the filter though
    group_sizes <- data.frame(Group = names(table(groupDF$Group)), n_group = as.numeric(table(groupDF$Group)))
    
    if (any(group_sizes$n_group == 1)) {
      # which group(s) #
      singleton_groups <-
        group_sizes$Group[group_sizes$n_group == 1]
      
      # which sample name(s) #
      samps_to_rm <-
        as.character(groupDF[which(groupDF$Group %in% singleton_groups), get_fdata_cname(omicsData)])
      
      warning("Grouping information is being utilized when calculating the CV, and there are group(s) consisting of a single sample. These singleton group(s) will be ignored by this filter.")
      
      # use custom_filter to remove the sample(s) from the omicsData object #
      my_cust_filt <-
        custom_filter(omicsData, f_data_remove = samps_to_rm)
      omicsData <- applyFilt(my_cust_filt, omicsData)
      
      groupDF <-
        groupDF[-which(groupDF[, get_fdata_cname(omicsData)] %in% samps_to_rm),]
    }
    
    ##########
    #    ## format the data ##
    #    temp_dat = data.table(omicsData$e_data)
    #    melt_dat = melt(temp_dat, id.var = edata_id)
    #
    #    setnames(melt_dat, names(melt_dat)[2], samp_id)
    #
    #    merge_dat = merge(x = melt_dat, y = data.table(groupDF), by = samp_id, all.x = TRUE)
    #
    #    # group the data by peptide and group #
    #    dat_grouped = dplyr::group_by_(data.frame(merge_dat), edata_id, "Group")
    ##########
    
    #added 2/2/17 iobani lines 53-58 and 63,68,73,78
    
    # Reorder the columns of cur_edata so the C++ function that calculates the
    # pooled CV can correctly account for group membership.
    cur_edata <- cur_edata[, order(groupDF$Group)]

    ##########
    #    cv_grouped = dplyr::group_by_(grp_cv, edata_id)
    #
    #    # calculate pooled cv #
    #    pool_cv = dplyr::summarise(cv_grouped,
    #    CV_pooled = 100*sum(CV_group[!is.na(CV_group)]*n_cv[!is.na(CV_group)])/
    #    sum(n_cv[!is.na(CV_group)]))
    ##########
    
    #added 2/2/17 iobani lines 87-95
    
    # Reorder the groups so the C++ function that calculates the pooled CV can
    # correctly account for group membership.
    group_dat <- as.character(groupDF$Group[order(groupDF$Group)])
    
    # Calculate the pooled CV.
    cvs <- pooled_cv_rcpp(as.matrix(cur_edata), group_dat)
    
    # For mystifying reasons multiply the pooled CV values by 100.
    cvs <- cvs * 100
    
    # Set the is_pooled variable to TRUE (because we pooled the CV).
    is_pooled <- TRUE
    
  }
  
  # Create and add attributes to the cv filter data frame ----------------------
  
  # Create a data frame with the ID column from e_data and the CV values. This
  # data frame is called pool_cv even though the CV may not be pooled (this
  # makes us mysterious). Call the column that holds the CV values CV_pooled.
  # Again, the CV may not be pooled (this makes us even more mysterious).
  pool_cv <- data.frame(omicsData$e_data[, id_col],
                        CV_pooled = cvs)
  names(pool_cv)[1] <- get_edata_cname(omicsData)
  
  ## determine plotting window cutoff ##
  # calculate percentage of observations with CV <= 200 #
  prct.less200 <- (sum(pool_cv$CV_pooled <= 200, na.rm = T) /
                     length(pool_cv$CV_pooled[!is.na(pool_cv$CV_pooled)]))
  
  if (prct.less200 > 0.95) {
    x.max = min(200, quantile(pool_cv$CV_pooled, 0.99, na.rm = TRUE))
  } else{
    x.max = quantile(pool_cv$CV_pooled, 0.95, na.rm = TRUE)
  }
  
  ## generate some summary stats for CV values, for PMART purposes only ##
  tot.nas <- sum(is.na(pool_cv$CV_pooled))
  
  output <- data.frame(pool_cv, row.names = NULL)
  
  orig_class <- class(output)
  
  class(output) <- c("cvFilt", orig_class)
  
  attr(output, "pooled") <- is_pooled
  attr(output, "max_x_val") <- x.max
  attr(output, "tot_nas") <- tot.nas
  
  # Return the completed object. You did it!!!
  return (output)

  #else{
  #   # group designation has NOT been called on omicsData #
  #   
  #   samp_id <- attr(omicsData, "cnames")$fdata_cname
  #   edata_id <- attr(omicsData, "cnames")$edata_cname
  #   
  #   cvs <-
  #     apply(omicsData$e_data[, -which(names(omicsData$e_data) == edata_id)], 1, cv.calc)
  #   
  #   pool_cv <-
  #     data.frame(omicsData$e_data[, which(names(omicsData$e_data) == edata_id)], CV_pooled = cvs)
  #   names(pool_cv)[1] <- edata_id
  #   
  #   ## determine plotting window cutoff ##
  #   # calculate percentage of observations with CV <= 200 #
  #   prct.less200 <-
  #     sum(pool_cv$CV_pooled <= 200, na.rm = T) / length(pool_cv$CV_pooled[!is.na(pool_cv$CV_pooled)])
  #   
  #   if (prct.less200 > 0.95) {
  #     x.max = min(200, quantile(pool_cv$CV_pooled, 0.99, na.rm = TRUE))
  #   } else{
  #     x.max = quantile(pool_cv$CV_pooled, 0.95, na.rm = TRUE)
  #   }
  #   
  #   ## generate some summary stats for CV values, for PMART purposes only ##
  #   tot.nas <- sum(is.na(pool_cv$CV_pooled))
  #   
  #   output <- data.frame(pool_cv, row.names = NULL)
  #   
  #   orig_class <- class(output)
  #   
  #   class(output) <- c("cvFilt", orig_class)
  #   
  #   attr(output, "sample_names") <-
  #     names(omicsData$e_data)[-which(names(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]
  #   attr(output, "group_DF") <- NULL
  #   attr(output, "max_x_val") <- x.max
  #   attr(output, "tot_nas") <- tot.nas
  #   
  # }
  
}


#' Calculate the Coefficient of Variation (CV)
#'
#' This function calculates the coefficient of variation for a vector of data
#'
#' @param data a vector of values on the "abundance" scale (not log transformed)
#'
#' @return cv value for the vector of values with missing values ignored in the calculation
#'
#' @author Lisa Bramer

cv.calc <- function(data) {
  # calculate mean of observations #
  ybar = mean(data, na.rm = TRUE)
  
  # calculate standard deviation of observations #
  ystd = sd(data, na.rm = TRUE)
  
  # calculate the sample cv value #
  sampcv = ystd / ybar
  
  # return cv values #
  return(sampcv)
}
