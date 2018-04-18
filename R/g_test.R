#' Implement the G-test 
#'
#' Can test for a qualitiative difference between groups for continuous data (see imd_test) or a quantiative difference between groups
#' for count data.  For continuous data, this is often used to test if the number of missing values is the same across all groups (i.e.a
#' qualitative difference).  For count data, this is used to test if the number of counts is the same for all groups.
#'
#' @param omicsData A mintR data object of any class
#' @param return_sizes Logical, should the group sizes be returned as well?  If true, a separate data.frame containing the group sizes is returned
#'
#' @return a data frame with the following columns: group means, global G-test statistic and corresponding p-value
#'
#' @author Bryan Stanfill
#' @references 
#' Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical analyses of peptide intensities and peptide occurrences improves identification of significant peptides from MS-based proteomics data." Journal of proteome research 9.11 (2010): 5748-5756.
#'
#' @examples 
#' library(mintJansson)
#' rRNA_data <- mintR::group_designation(omicsData = rRNA_data, main_effects = c("site","treatment"))
#' group_df <- attr(rRNA_data, "group_DF")
#' #This is will return a lot of errors because the peptides with too few data haven't 
#' #been filtered out yet
#' g_res <- g_test(omicsData = rRNA_data, groupData = group_df)
#' 
#' @export

g_test <- function(omicsData, return_sizes = FALSE){
  # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }
  
  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("At least two groups are necessary to perform a G-test")
  }
  
  #Double check that count data are provided, give warning if not.
  if(attr(omicsData,"data_info")$data_scale!="count"){
    warning("This function is intendend for count data only; consider using 'anova_test(...)' for the data provided.")
  }
  
  #Remove rows in "groupData" that don't have corresponding columns in 'omicsData$edata'
  groupData <- dplyr::filter(groupData,sample_name%in%colnames(omicsData$e_data))
  
  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData$sample_name)]

  gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))
  raw_results <- gtest_cpp(data.matrix(data),gp)
  #The C++ function returns a list so we need to make it into a data.frame
  results <- data.frame(raw_results$Gstat, raw_results$pvalue,raw_results$group_means)
  
  #Rename the columns to match group names
  group_names <- paste("Group",as.character(unique(groupData$Group)))
  colnames(results) <- c("G-Statistic","p-value",group_names)
  
  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  results <- cbind(omicsData$e_data[edatacname],results)
  
  if(return_sizes){
    sizes_df <- raw_results$group_sizes
    colnames(sizes_df) <- group_names
    return(list(Results=results,Sizes=sizes_df))
  }
  
  return(results)
}



