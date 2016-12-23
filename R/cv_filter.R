#' Filter Based on Pooled Coefficient of Variation (CV) Values
#'
#' A pooled CV is calculated for each biomolecule.
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', or 'lipidData', created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{as.lipidData}}, respectively. Note, if \code{\link{group.designation}} has not been run, the CV is calculated based on all samples for each biomolecule.
#'
#' @return  An S3 object of class 'cvFilt' giving the pooled CV for each biomolecule and additional attributes used for plottinga data.frame with a column giving the peptide, protein, or gene name and a column giving the pooled CV value
#'
#'@details For each biomolecule, the CV of each group is calculated as the standard deviation divided by the mean, excluding missing values. A pooled CV estimate is then calculated based on the methods of Ahmed (1995).
#'@references Ahmed, S.E. (1995). \emph{A pooling methodology for coefficient of variation}. The Indian Journal of Statistics. 57: 57-75.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_object2 <- group_designation(omicsData = pep_object, main_effects = "Condition")
#' to_filter <- cv_filter(omicsData = pep_object2)
#' summary(to_filter, cv_threshold = 30)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' @export

cv_filter <- function(omicsData){

  ## some intial checks ##

  # check that omicsData is of appropriate class #
  if(!(class(omicsData) %in% c("pepData", "proData", "lipidData", "metabData"))) stop("omicsData must be an object of class 'pepData', 'proData', 'lipidData', or 'metabData'")

  # check that group_designation has been called already #
  if(!is.null(attr(omicsData, "group_DF"))){

    groupDF = attr(omicsData, "group_DF")

    samp_id = attr(omicsData, "cnames")$fdata_cname
    edata_id = attr(omicsData, "cnames")$edata_cname

    ## format the data ##
    temp_dat = data.table::data.table(omicsData$e_data)
    melt_dat = data.table::melt.data.table(temp_dat, id.var = edata_id)

    data.table::setnames(melt_dat, names(melt_dat)[2], samp_id)

    merge_dat = data.table:::merge.data.table(x = melt_dat, y = data.table::data.table(groupDF), by = samp_id, all.x = TRUE)

    # group the data by peptide and group #
    dat_grouped = dplyr::group_by_(data.frame(merge_dat), edata_id, "Group")

    # calculate cv (on original scale) and number of non-missing values, by peptide and group #
    if(attr(omicsData, "data_info")$data_scale == "log2"){
      grp_cv = dplyr::summarise(dat_grouped, CV_group = cv.calc(2^value), n_cv = sum(!is.na(2^value)))
    }

    if(attr(omicsData, "data_info")$data_scale == "log10"){
      grp_cv = dplyr::summarise(dat_grouped, CV_group = cv.calc(10^value), n_cv = sum(!is.na(10^value)))
    }

    if(attr(omicsData, "data_info")$data_scale == "log"){
      grp_cv = dplyr::summarise(dat_grouped, CV_group = cv.calc(exp(value)), n_cv = sum(!is.na(exp(value))))
    }

    if(attr(omicsData, "data_info")$data_scale == "abundance"){
      grp_cv = dplyr::summarise(dat_grouped, CV_group = cv.calc(value), n_cv = sum(!is.na(value)))
    }

    cv_grouped = dplyr::group_by_(grp_cv, edata_id)

    # calculate pooled cv #
    pool_cv = dplyr::summarise(cv_grouped, CV_pooled = 100*sum(CV_group[!is.na(CV_group)]*n_cv[!is.na(CV_group)])/sum(n_cv[!is.na(CV_group)]))

    ## determine plotting window cutoff ##
    # calculate percentage of observations with CV <= 200 #
    prct.less200 = sum(pool_cv$CV_pooled <= 200,na.rm=T)/length(pool_cv$CV_pooled[!is.na(pool_cv$CV_pooled)])

    if(prct.less200 > 0.95){x.max = min(200, quantile(pool_cv$CV_pooled, 0.99, na.rm=TRUE))}else{x.max = quantile(pool_cv$CV_pooled, 0.95, na.rm = TRUE)}

    ## generate some summary stats for CV values, for PMART purposes only ##
    tot.nas = sum(is.na(pool_cv$CV_pooled))

  output <- data.frame(pool_cv, row.names=NULL)

  orig_class <- class(output)

  class(output) <- c("cvFilt", orig_class)

  attr(output, "sample_names") <- names(omicsData$e_data)[-which(names(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]
  attr(output, "group_DF") <- attr(omicsData, "group_DF")
  attr(output, "max_x_val") <- x.max
  attr(output, "tot_nas") <- tot.nas

  }else{

    samp_id = attr(omicsData, "cnames")$fdata_cname
    edata_id = attr(omicsData, "cnames")$edata_cname

    cvs = apply(omicsData$e_data[,-which(names(omicsData$e_data) == edata_id)], 1, cv.calc)

    pool_cv = data.frame(omicsData$e_data[,which(names(omicsData$e_data) == edata_id)], CV_pooled = cvs)
    names(pool_cv)[1] = edata_id

    ## determine plotting window cutoff ##
    # calculate percentage of observations with CV <= 200 #
    prct.less200 = sum(pool_cv$CV_pooled <= 200,na.rm=T)/length(pool_cv$CV_pooled[!is.na(pool_cv$CV_pooled)])

    if(prct.less200 > 0.95){x.max = min(200, quantile(pool_cv$CV_pooled, 0.99, na.rm=TRUE))}else{x.max = quantile(pool_cv$CV_pooled, 0.95, na.rm = TRUE)}

    ## generate some summary stats for CV values, for PMART purposes only ##
    tot.nas = sum(is.na(pool_cv$CV_pooled))

    output <- data.frame(pool_cv, row.names=NULL)

    orig_class <- class(output)

    class(output) <- c("cvFilt", orig_class)

    attr(output, "sample_names") <- names(omicsData$e_data)[-which(names(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]
    attr(output, "group_DF") <- NULL
    attr(output, "max_x_val") <- x.max
    attr(output, "tot_nas") <- tot.nas

}
  return(output)
}


#' Calculate the Coefficient of Variation (CV)
#'
#' This function calculates the coefficient of variation for a vector of data
#'
#' @param data a vector of values
#'
#' @return cv value for the vector of values with missing values ignored in the calculation
#'
#' @author Lisa Bramer

cv.calc <- function(data){
  # calculate mean of observations #
  ybar = mean(data, na.rm = TRUE)

  # calculate standard deviation of observations #
  ystd = sd(data, na.rm = TRUE)

  # calculate the sample cv value #
  sampcv = ystd/ybar

  # return cv values #
  return(sampcv)
}
