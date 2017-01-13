#' IMD-ANOVA filter object
#'
#' This function returns an imdanovaFilt object for use with \code{\link{applyFilt}}
#'
#' @param omicsData object of one of the classes "pepData", "proData", "lipidData", or "metabData", usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, \code{\link{as.metabData}}, respectively.
#'
#' @details The output from this function can be used in conjunction with \code{\link{applyFilt}} to filter out molecules that are not present in enough samples to do statistical comparisons.
#'
#' @return Object of class imdanovaFilt (also a data.frame) containing the molecule identifier and number of samples in each group with non-missing values for that molecule.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
#' to_filter <- imdanova_filter(omicsData = pep_pepData2)
#' summary(to_filter, min_nonmiss_anova = 3)
#'}
#'
#' @author Kelly Stratton
#'
#' @export
imdanova_filter <- function(omicsData){ #}, filter_method, min_nonmiss_gtest=3, min_nonmiss_anova=2, alpha=NULL){

  ## Initial checks ##

  # omicsData must be of the appropriate class
  if(!(class(omicsData) %in% c("pepData", "proData", "lipidData", "metabData"))){
    stop("Invalid class for omicsData. Valid classes are pepData, proData, lipidData, and metabData.")
  }

  # omicsData must include groupDF information #
  # if(!is.null(omicsData)){
  e_data <- omicsData$e_data
  groupDF <- attr(omicsData, "group_DF")


  ## end of initial checks ##

  samp_id <- attr(omicsData, "cnames")$fdata_cname
  edata_id <- attr(omicsData, "cnames")$edata_cname
  emeta_id <- attr(omicsData, "cnames")$emeta_cname

  # if groupDF has column for TimeCourse, do the following within time point (do not re-compute groupDF including TimeCourse as main effect--this was our old strategy but Bobbie-Jo directed us to not do this, and instead loop through time points)
  if(any(names(groupDF)=="TimeCourse")){

#     filt.edata <- vector(mode = "list", length = length(unique(groupDF$TimeCourse)))
#     names(filt.edata) <- unique(groupDF$TimeCourse)
#
#     for(tind in 1:length(unique(groupDF$TimeCourse))){
#
#       t = unique(groupDF$TimeCourse)[tind]
#       t.e_data <- cbind(e_data[,1], e_data[, names(e_data) %in% as.character(groupDF[,samp_id][groupDF$TimeCourse==t])])
#       names(t.e_data)[1] <- names(e_data)[1]
#
#       t.groupDF <- groupDF[groupDF$TimeCourse==t, ]
#       #all(names(t.e_data)[-1] == t.groupDF$sampleID) # just checking, should be TRUE
#
#       nonmiss_per_group <- nonmissing_per_group(omicsData=NULL, e_data=t.e_data, groupDF=t.groupDF, cname_id=edata_id, samp_id=samp_id)
#       if(filter_method=="anova"){
#         filt.edata[[tind]] <- anova_filter(nonmiss_per_group=nonmiss_per_group, min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_id)
#       }else{
#         if(filter_method=="gtest"){
#           filt.edata[[tind]] <- gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=t.groupDF, e_data=t.e_data, alpha=alpha, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_id, samp_id = samp_id)
#         }else{
#           if(filter_method=="combined"){
#             filt.edata.gtest <- gtest_filter(nonmiss_per_group, groupDF=t.groupDF, e_data=t.e_data, alpha=alpha, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_id)
#             #           min.nonmiss.allowed <- 2
#             filt.edata.anova <- anova_filter(nonmiss_per_group, min_nonmiss_anova, cname_id = edata_id)
#             filt.edata[[tind]] <- intersect(filt.edata.anova, filt.edata.gtest)
#           }
#         }
#       }
#     }
#
#     filter.edata <- Reduce(base::intersect, filt.edata)

  }else{ # end of if-statement for the presence of TimeCourse variable
    nonmiss_per_group <- nonmissing_per_group(omicsData=omicsData)
    output <- nonmiss_per_group$nonmiss_totals

    #attr(output, "filter_method") <- filter_method

#     if(filter_method=="anova"){
#       filter.edata <- anova_filter(nonmiss_per_group=nonmiss_per_group, min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_id)
#     }else{
#       if(filter_method=="gtest"){
#         filter.edata <- gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=groupDF, e_data=e_data, alpha=NULL, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_id, samp_id = samp_id)
#       }else{
#         if(filter_method=="combined"){
#           filter.edata.gtest <- gtest_filter(nonmiss_per_group=nonmiss_per_group, groupDF=groupDF, e_data=e_data, alpha=NULL, min_nonmiss_gtest=min_nonmiss_gtest, cname_id = edata_id, samp_id = samp_id)
#           #           min.nonmiss.allowed <- 2
#           filter.edata.anova <- anova_filter(nonmiss_per_group=nonmiss_per_group, min_nonmiss_anova=min_nonmiss_anova, cname_id = edata_id)
#           filter.edata <- intersect(filter.edata.anova, filter.edata.gtest)
#         }
#       }
#     }
  } # end of else-stament for the absence of TimeCourse variable


  orig_class <- class(output)
  class(output) <- c("imdanovaFilt", orig_class)

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes

  return(output)

}




