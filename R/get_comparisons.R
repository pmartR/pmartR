#' Returns data.frame with comparisons to be made 
#' 
#' The method creates a data.frame containing the comparisons to be made using differential statistics.  
#'  
#' @param comp_type string taking one of the following values: control, pairwise, custom. Specifying "control" indicates that all other groups are to be compared to this single control group. Specifying "pairwise" indicates that all pairwise comparisons are to be made. 
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param control_group string indicating the group to use for the control group. Only required when comp_type="control".
#'
#' @details This function takes in the omicsData and type of comparison, and returns a data.frame where each row corresponds to a comparison of interest.
#'
#' @return data.frame with columns for Test and Control. Each row corresponds to a comparison of interest. 
#' #'  
#' @author Kelly Stratton
#' 
#' @export
get_comparisons <- function(comp_type, omicsData, control_group=NULL){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    group_data <- attr(omicsData, "group_DF")
  }
  
  # do we even run this if they want custom groups? 
  # need to address the order of the factor levels somewhere...maybe not in this function, but somewhere...
  
  
  ## initial checks ##
  if(!(comp_type %in% c("control", "pairwise", "custom"))){
    stop("Comparison type must be one of the following: control, pairwise, custom")
  }
  
  if(comp_type=="control" & is.null(control_group)){
    stop("A control_group must be specified when comp_type is set to 'control'.")
  }
  
  #if(comp_type=="control" & (table(group_data$Group)[control_group]<2)){
  #  stop("The control_group selected has fewer than two subjects and cannot be used.")
  #}
  
  ## end of initial checks ##
  
  grps <- unique(group_data$Group)
  n_grps <- length(grps)
  
  
  if(comp_type=="control"){
    query_grps <- grps[grps != control_group]
    control_grps <- rep(control_group, length(query_grps))
  }
  
  if(comp_type=="pairwise"){
    query_grps <- c()
    control_grps <- c()
    for(i in 1:n_grps){
      for(j in (i+1):n_grps){
        if(i != j){
          control_grps <- c(control_grps, as.character(grps[i]))
          query_grps <- c(query_grps, as.character(grps[j]))
        }
      }
    }
    #cbind(control_grps, query_grps)
    # remove any rows with an NA
    inds <- c(which(is.na(control_grps)), which(is.na(query_grps)))
    inds <- unique(inds)
    query_grps <- query_grps[-inds]
    control_grps <- control_grps[-inds]
  }
  
  
  comparisons <- data.frame(Test = query_grps, Control = control_grps)
  return(comparisons)
}