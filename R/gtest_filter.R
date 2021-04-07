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
#' @param omicsData an optional object of one of the classes "pepData",
#' "proData", "lipidData", or "metabData" usually created by
#' \code{\link{as.pepData}}, \code{\link{as.proData}},
#' \code{\link{as.lipidData}}, or \code{\link{as.metabData}}, respectively.
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
#' pep_object2 <- group_designation(omicsData = pep_object,
#'                                  main_effects = "Condition")
#' nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
#' to_filter <- gtest_filter(nonmiss_per_group = nonmissing_result,
#'                           omicsData = pep_object2,
#'                           min_nonmiss_gtest = 3)
#'}
#'
#' @author Kelly Stratton
#'
#' @export
#'
gtest_filter <- function (nonmiss_per_group,
                          omicsData,
                          min_nonmiss_gtest = 3) {
  
  # Ensure omicsData is an appropriate class
  if (!inherits(omicsData, c("pepData", "proData", "lipidData", "metabData",
                             "nmrData"))) {
    
    # Throw down an error for blatant abuse of the pmartR standards.
    stop (paste("omicsData must be of class 'pepData', 'proData', 'lipidData',",
                "'metabData' or 'nmrData'",
                sep = " "))
    
  }
  
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

  # the order of the peptides in groups has been changed--it no longer matches
  # the order in the data (added by KS on 10/13/2015)
  inds <- match(omicsData$e_data[, 1], nonmiss_per_group$nonmiss_totals[, 1])
  temp <- nonmiss_per_group$nonmiss_totals[inds, ]
  nonmiss_per_group$nonmiss_totals <- temp
  
  # to make sure it's an object of class data.frame and not also data.table
  groups <- data.frame(nonmiss_per_group$nonmiss_totals) 
  
  # remove the column with Peptide info and any NA's
  groups2 <- groups[, -c(1, which(names(groups) %in% c(NA, "<NA>", "NA.")))] 
  # for ea peptide/row, need max(nonmiss_per_group$nonmiss_totals[,-1]) to be
  # >= min_nonmiss_gtest
  trueFalse <- apply(groups2, 1, function (x) max(x) >= min_nonmiss_gtest)
  
  peps <- as.character(groups[, 1])
  myresults <- data.frame(Peptide = peps, MeetsMinReq = trueFalse)
  
  myresults2 <- myresults[myresults$MeetsMinReq == FALSE, ]
  peps.filt <- as.character(myresults2$Peptide)

  return (peps.filt)
  
}
