#' Identifies peptides to be filtered out in preparation for IMD-ANOVA.
#'
#' The method identifies peptides, proteins, lipids, or metabolites to be
#' filtered specifically according to the G-test.
#'
#' @param nonmiss_per_group list created by \code{\link{nonmissing_per_group}}.
#' The first element giving the total number of possible samples for each group.
#' The second element giving a data.frame with the first column giving the 
#' biomolecule and the second through kth columns giving the number of
#'  non-missing observations for each of the \code{k} groups.
#' 
#' @param min_nonmiss_gtest the minimum number of non-missing peptide values
#' allowed in a minimum of one group. Default value is 3.
#'
#' @details Two methods are available for determining the peptides to be
#' filtered. The naive approach is based on \code{min.nonmiss.allowed}, and
#' looks for peptides that do not have at least \code{min.nonmiss.allowed}
#' values per group. The other approach also looks for peptides that do not
#' have at least a minimum number of values per group, but this minimum number
#' is determined using the G-test and a p-value threshold supplied by the user.
#' The G-test is a test of independence, used here to test the null hypothesis
#' of independence between the number of missing values across groups.
#'
#' @return filter.peps a character vector of the peptides to be filtered out
#' prior to the G-test or IMD-ANOVA
#'
#' @examples
#' \dontrun{
#' library(pmartR)
#' data(pep_object)
#' 
#' pep_object2 <- group_designation(omicsData = pep_object,
#'                                  main_effects = "Condition")
#' 
#' nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
#' 
#' to_filter <- gtest_filter(nonmiss_per_group = nonmissing_result,
#'                           min_nonmiss_gtest = 3)
#'}
#'
#' @author Kelly Stratton
#'
#' @export
#'
gtest_filter <- function (nonmiss_per_group,
                          min_nonmiss_gtest = 3) {
  
  # check that min_nonmiss_gtest is of length 1 #
  if (length(min_nonmiss_gtest) != 1) {
    
    # Warn the user of their treachery with an error.
    stop ("min_nonmiss_gtest must be of length 1")
    
  }
  
  # min_nonmiss_gtest must be >=2
  if (!is.null(min_nonmiss_gtest)) {
    
    if (min_nonmiss_gtest < 2) {
      
      stop ("min_nonmiss_gtest must be >=2")
      
    }
    
  }
  
  # check that nonmiss_per_group is a list of length 2 #
  if (!inherits(nonmiss_per_group, "list") || length(nonmiss_per_group) != 2) {
    
    stop ("nonmiss_per_group must be a list of length 2.")
    
  }
  
  # remove the column with Peptide info and any NA's
  # for each peptide/row, need max(nonmiss_per_group$nonmiss_totals[,-1]) to be
  # >= min_nonmiss_gtest
  trueFalse <- apply(
    nonmiss_per_group$nonmiss_totals[, -c(
      1,
      which(names(nonmiss_per_group$nonmiss_totals) %in%
              c(NA, "<NA>", "NA.", "NA"))
    )],
    1,
    function (x) max(x) >= min_nonmiss_gtest
  )
  
  # Extract the biomolecule IDs corresponding to counts that fall below the
  # minimum. These are the biomolecules that will be filtered out.
  return (as.character(nonmiss_per_group$nonmiss_totals[!trueFalse, 1]))
  
}
