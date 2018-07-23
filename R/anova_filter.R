#' Identifies biomolecules to be filtered in preparation for IMD-ANOVA.
#'
#' The method identifies biomolecules to be filtered specifically according data requirements for running an ANOVA.
#'
#' @param nonmiss_per_group a list of length two. The first element giving the total number of possible samples for each group. The second element giving a data.frame with the first column giving the biomolecule identifier and the second through kth columns giving the number of non-missing observations for each of the \code{k} groups. Usually the result of \code{\link{nonmissing_per_group}}
#' @param min_nonmiss_anova the minimum number of nonmissing biomolecule values required, in each group, in order for the biomolecule to not be filtered. Must be greater than or equal to 2; default value is 2.
#' @param cname_id character string specifying the name of the column containing the biomolecule identifiers in \code{e_data} and \code{e_meta} (if applicable).
#'
#' @details This function filters biomolecules that do not have at least \code{min.nonmiss.allowed} values per group, where groups are from \code{group_designation}.
#'
#' @return filter.peps a character vector of the biomolecules to be filtered out prior to ANOVA or IMD-ANOVA
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' pep_object <- group_designation(pep_object, main_effects = "Condition")
#' filter_results <- anova_filter(nonmiss_per_group = nonmissing_per_group(omicsData = pep_object), cname_id = "Mass_Tag_ID")
#' }
#'
#' @seealso \code{\link{nonmissing_per_group}}
#'
#' @author Kelly Stratton
#'
#' @export
anova_filter <- function(nonmiss_per_group, min_nonmiss_anova=2, cname_id){

  # min_nonmiss_anova must be >=2
  if(!is.null(min_nonmiss_anova)){
    if(min_nonmiss_anova<2){
      stop("min_nonmiss_anova must be >=2")
    }
  }

  # check that nonmiss_per_group is a list of length 2 #
  if(!inherits(nonmiss_per_group, "list") | length(nonmiss_per_group)!=2) stop("nonmiss_per_group must be a list of length 2.")

  # column names of nonmiss_per_group$nonmiss_totals
  my.names <- names(nonmiss_per_group$nonmiss_totals)
  inds.names <- my.names %in% c(cname_id, "NA")
  inds.names <- which(inds.names==TRUE)
  my.names <- my.names[-inds.names] # remove cname_id and "NA" (if there's a column named "NA")

  # need at least n=2 per group in at least 2 groups in order to keep the biomolecule
  temp <- (nonmiss_per_group$nonmiss_totals[,which(names(nonmiss_per_group$nonmiss_totals) %in% my.names)] >= min_nonmiss_anova)

  # sum the number of groups that meet the nonmissing per group requirement
  temp2 <- rowSums(temp)

  # append the vector of biomolecules
  temp3 <- data.frame(nonmiss_per_group$nonmiss_totals[, 1], temp2)
  names(temp3)[1] <- names(nonmiss_per_group$nonmiss_totals)[1]

  # create indicator for which rows do not meet the nonmissing requirement for at least 2 groups
  inds.rm <- which(temp2 < 2) # these are the rows to remove since they do not have at least 2 groups not meeting nonmissing requirements

  # get names of biomolecules to be filtered
  filter.ids <- as.vector(temp3[inds.rm, 1])

  # output names of peptides/proteins/genes to be filtered
  return(filter.ids)
}

