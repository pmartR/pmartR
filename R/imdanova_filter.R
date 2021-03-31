#'IMD-ANOVA filter object
#'
#'This function returns an imdanovaFilt object for use with
#'\code{\link{applyFilt}}
#'
#'@param omicsData object of one of the classes "pepData", "isobaricpepData",
#'  "proData", "lipidData", "metabData", or "nmrData", usually created by
#'  \code{\link{as.pepData}}, \code{\link{as.isobaricpepData}},
#'  \code{\link{as.proData}}, \code{\link{as.lipidData}},
#'  \code{\link{as.metabData}}, or \code{\link{as.nmrData}}, respectively.
#'  Groups (more than one group) must have been specified using the
#'  \code{\link{group_designation}} function prior to using the imdanova_filter
#'  function.
#'
#'@details The output from this function can be used in conjunction with
#'  \code{\link{applyFilt}} to filter out molecules that are not present in
#'  enough samples to do statistical comparisons. If any singleton groups are 
#'  present in the omicsData object, those groups are not part of the filter 
#'  object that is returned.
#'
#'@return Object of class imdanovaFilt (also a data.frame) containing the
#'  molecule identifier and number of samples in each group with non-missing
#'  values for that molecule.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_pepData2 <- group_designation(omicsData = pep_object,
#'                                   main_effects = "Condition")
#' to_filter <- imdanova_filter(omicsData = pep_pepData2)
#' summary(to_filter, min_nonmiss_anova = 2)
#'}
#'
#'@author Kelly Stratton
#'
#'@export
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
  
  # groupDF must have more than 1 group #
  if (length(names(get_group_table(omicsData))) < 2) {
    
    # Let the user know that they cannot compare statistics between groups if
    # there is only one group!!
    stop (paste("There must be more than one group in order to create an",
                "imdanovaFilt object.",
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

  } else { # end of if-statement for the presence of TimeCourse variable
    
    # Count the number of nonmissing elements per group per biomolecule. For
    # example if group A has 5 samples and 4 of the samples have missing values
    # then the count for group A will be 1.
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

  attr(output, "group_sizes") <- nonmiss_per_group$group_sizes
  # KS added attribute for nonsingleton groups 12/3/2020 #
  attr(output, "nonsingleton_groups") <- nonmiss_per_group$group_sizes$Group[which(nonmiss_per_group$group_sizes$n_group > 1)]
  

  return (output)

}
