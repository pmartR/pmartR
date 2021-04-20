#' Identifies biomolecules to be filtered in preparation for IMD-ANOVA.
#'
#' The method identifies biomolecules to be filtered specifically according data
#' requirements for running an ANOVA.
#'
#' @param nonmiss_per_group a list of length two. The first element giving the
#' total number of possible samples for each group. The second element giving a
#' data.frame with the first column giving the biomolecule identifier and the
#' second through kth columns giving the number of non-missing observations for
#' each of the \code{k} groups. Usually the result of
#' \code{\link{nonmissing_per_group}}
#' 
#' @param min_nonmiss_anova the minimum number of nonmissing biomolecule values
#' required, in each group, in order for the biomolecule to not be filtered.
#' Must be greater than or equal to 2; default value is 2.
#'
#' @details This function filters biomolecules that do not have at least
#' \code{min.nonmiss.allowed} values per group, where groups are from
#' \code{group_designation}.
#'
#' @return filter.peps a character vector of the biomolecules to be filtered out
#' prior to ANOVA or IMD-ANOVA
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' 
#' pep_object2 <- group_designation(pep_object, main_effects = "Condition")
#' 
#' nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
#' 
#' to_filter <- anova_filter(nonmiss_per_group = nonmissing_result,
#'                           min_nonmiss_anova = 2)
#' }
#'
#' @seealso \code{\link{nonmissing_per_group}}
#'
#' @author Kelly Stratton
#'
#' @export
#' 
anova_filter <- function (nonmiss_per_group,
                          min_nonmiss_anova = 2) {
  
  # check that min_nonmiss_anova is of length 1 #
  if (length(min_nonmiss_anova) != 1) {
    
    # Warn the user of their treachery with an error.
    stop ("min_nonmiss_anova must be of length 1")
    
  }

  # min_nonmiss_anova must be >=2
  if (!is.null(min_nonmiss_anova)) {
    
    if (min_nonmiss_anova < 2) {
      
      stop ("min_nonmiss_anova must be >=2")

    }
    
  }

  # check that nonmiss_per_group is a list of length 2 #
  if (!inherits(nonmiss_per_group, "list") || length(nonmiss_per_group) != 2) {
    
    stop ("nonmiss_per_group must be a list of length 2.")
    
  }

  # Remove the column with the biomolecule IDs and any columns that have an NA
  # as the column name. Then sum (by row) the number of groups that meet the
  # non-missing per group requirement. Then extract the rows that do not meet
  # the requirement of at least two groups with non-missing counts higher than
  # the threshold. These are the rows that will be filtered out later.
  inds_rm <- which(rowSums(
    nonmiss_per_group$nonmiss_totals[, -c(
      1,
      which(names(nonmiss_per_group$nonmiss_totals) %in%
              c(NA, "<NA>", "NA.", "NA"))
    )] >= min_nonmiss_anova
  ) < 2)

  # get names of biomolecules to be filtered
  return (as.character(nonmiss_per_group$nonmiss_totals[inds_rm, 1]))
  
}
