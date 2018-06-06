#' Identifies peptides to be filtered out in preparation for IMD-ANOVA.
#'
#' The method identifies peptides, proteins, lipids, or metabolites to be filtered specifically according to the G-test.
#'
#' @param nonmiss_per_group list created by \code{\link{nonmissing_per_group}}. The first element giving the total number of possible samples for each group. The second element giving a data.frame with the first column giving the biomolecule and the second through kth columns giving the number of non-missing observations for each of the \code{k} groups.
#' @param groupDF data.frame created by \code{group_designation} with columns for sampleID and Group. If omicsData is supplied, this argument is optional; if e_data is supplied, this argument is required. If two main effects are provided the original main effect levels for each sample are returned as the third and fourth columns of the data.frame.
#' @param omicsData an optional object of one of the classes "pepData", "proData", "lipidData", or "metabData" usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.lipidData}}, or \code{\link{as.metabData}}, respectively.
#' @param e_data an optional \eqn{p} by \eqn{n} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples.
#' @param alpha numeric value specifing the p-value cut-off to be used when determining whether to filter out each peptide. Either \code{alpha} or \code{min.nonmiss.allowed}, but not both, must be specified; this specification determines which approach is used for filtering the peptides. See \link{Details} for more information.
#' @param min_nonmiss_gtest the minimum number of non-missing peptide values allowed in a minimum of one group. Default value is 3. Either \code{alpha} or \code{min.nonmiss.allowed} must be specified; this specification determines which approach is used for filtering the peptides. See \link{Details} for more information.
#' @param cname_id name for biomolecule identifier column found in \code{e_data} and \code{e_meta} (if applicable)
#' @param samp_id name for sample identifier column found in \code{f_data}
#'
#' @details Two methods are available for determining the peptides to be filtered. The naive approach is based on \code{min.nonmiss.allowed}, and looks for peptides that do not have at least \code{min.nonmiss.allowed} values per group. The other approach also looks for peptides that do not have at least a minimum number of values per group, but this minimum number is determined using the G-test and a p-value threshold supplied by the user. The G-test is a test of independence, used here to test the null hypothesis of independence between the number of missing values across groups.
#'
#' @return filter.peps a character vector of the peptides to be filtered out prior to the G-test or IMD-ANOVA
#'
#' @examples
#' dontrun{
#' library(pmartR)
#' data(pep_object)
#' pep_object2 <- group_designation(omicsData = pep_object, main_effects = "Condition")
#' nonmissing_result <- pmartR:::nonmissing_per_group(omicsData = pep_object2)
#' to_filter <- gtest_filter(nonmiss_per_group = nonmissing_result, omicsData = pep_object2, min_nonmiss_gtest = 3)
#'}
#'
#' @author Kelly Stratton
#'
#' @export
#'

gtest_filter <- function(nonmiss_per_group, groupDF=NULL, omicsData=NULL, e_data=NULL, alpha=NULL, min_nonmiss_gtest=NULL, cname_id = NULL, samp_id = NULL){

  # one of omicsData and e_data is required (but not both); if both are provided use e_data; if omicsData is provided, define peptide data; if neither is specified then stop
  if(is.null(omicsData) & is.null(e_data)){
    stop("Neither omicsData nor e_data has been specified. One of these two arguments is required.")
  }else{
    if(!is.null(omicsData) & !is.null(e_data)){
      warning("Both omicsData and e_data were specified. Function will proceed using e_data")
    }else{
      if(!is.null(omicsData) & is.null(e_data)){
        e_data <- omicsData$e_data
      }
    }
  }
  
  # check that groupDF is supplied if e_data is supplied #
  if(!is.null(omicsData)){
    e_data <- omicsData$e_data
    groupDF <- attr(omicsData, "group_DF")

    samp_id <- attr(omicsData, "cnames")$fdata_cname
    cname_id <- attr(omicsData, "cnames")$edata_cname
  }else{
    if(is.null(e_data)) stop("Either omicsData or e_data and groupDF must be provided.")
    if(is.null(groupDF)) stop("Either omicsData or e_data and groupDF must be provided.")

    if(is.null(samp_id)) stop("samp_id must be provided if omicsData is not")
    if(is.null(cname_id)) stop("cname_id must be provided if omicsData is not")
  }

  # one of alpha and min_nonmiss_gtest is required (but not both); if both alpha and min_nonmiss_gtest are specified, use the alpha approach and warn the user that this is happening; if neither is specified then stop #
  if(!is.null(alpha) & !is.null(min_nonmiss_gtest)){
    warning("Both alpha and min_nonmiss_gtest have been specified. Function will proceed using the alpha value and more sophisticated filtering approach.")
    min_nonmiss_gtest=NULL
  }else{
    if(is.null(alpha) & is.null(min_nonmiss_gtest)){
      stop("Neither alpha nor min_nonmiss_gtest has been specified. One of these two arguments is required.")
    }
  }

  

  # the order of the peptides in groups has been changed--it no longer matches the order in the data (added by KS on 10/13/2015)
  inds <- match(e_data[,1], nonmiss_per_group$nonmiss_totals[,1])
  temp <- nonmiss_per_group$nonmiss_totals[inds,]
  nonmiss_per_group$nonmiss_totals <- temp
  #all(e_data$Mass_Tag_ID==nonmiss_per_group$nonmiss_totals$Peptide)

  ### Use min non-missing approach ###
  if(!is.null(min_nonmiss_gtest) & is.null(alpha)){
    groups <- data.frame(nonmiss_per_group$nonmiss_totals) # to make sure it's an object of class data.frame and not also data.table

    groups2 <- groups[, -c(1, which(names(groups) %in% c(NA, "<NA>", "NA.")))] # remove the column with Peptide info and any NA's
    # for ea peptide/row, need max(nonmiss_per_group$nonmiss_totals[,-1]) to be >= min_nonmiss_gtest
    trueFalse <- apply(groups2, 1, function(x) max(x)>=min_nonmiss_gtest)

    peps <- as.character(groups[, 1])
    myresults <- data.frame(Peptide=peps, MeetsMinReq=trueFalse)

    myresults2 <- myresults[myresults$MeetsMinReq==FALSE,]
    peps.filt <- as.character(myresults2$Peptide)
  }

  ### use p-value threshold approach ###
  if(!is.null(alpha)){

    # make sure group.vals.unique is character and not factor
    group.vals.unique <- as.character(unique(groupDF$Group))
    # remove any NA or <NA> group values #
    if(any(group.vals.unique %in% c(NA, '<NA>'))){
      inds <- which(group.vals.unique %in% c(NA, '<NA>'))
      group.vals.unique <- group.vals.unique[-inds]
    }

    num.groups <- length(group.vals.unique)

    result <- rep(NA, num.groups)
    for(i in 1:num.groups){
      # check the pval if everything is present in one group and not in the others; if pval>alpha, skip the while loop and set num.nonmiss to 0
      group.size <- length(groupDF$Group[groupDF$Group==group.vals.unique[i]]) - length(groupDF$Group[groupDF$Group=="<NA>"])
      inds <- groupDF$Group==group.vals.unique[i] # which elements of groupDF are in the ith group?
      inds[is.na(inds)] <- FALSE
      cur.data <- rep(NA, length(groupDF$Sample.ID))
      cur.data[inds==TRUE] <- c(rep(1, group.size))
      cur.data <- t(data.frame(cur.data))
      if(gtest_vector(cur.data, groupDF$Group)$pvalues > alpha){
        # if the pvalue is big (i.e. > alpha) when all data is present in 1 group but not the others, then we don't need to go through the while loop because the G-test will never be significant
        num.nonmiss <- 0 #max(table(groupDF$Group)) + 1000
      }else{ # find the smallest number of peptide identifications that must be present in 1 group in order for the G-test to be significant
        cur.data <- rep(NA, length(groupDF$Sample.ID))
        cur.data <- t(data.frame(cur.data))
        pval <- 1
        num.nonmiss <- 1
        while(pval > alpha && num.nonmiss < length(groupDF$Group[inds])){
          # modify cur.data #
          inds <- groupDF$Group==group.vals.unique[i]
          inds[is.na(inds)] <- FALSE # set any missing values to FALSE

          cur.data[inds] <- c(rep(1, num.nonmiss), rep(NA, length(groupDF$Group[inds])-num.nonmiss))
          pval <- gtest_vector(cur.data, groupDF$Group)$pvalues

          num.nonmiss <- num.nonmiss + 1
        }
      }

      result[i] <- max(0, (num.nonmiss - 1))
    }

    # now use same approach as min.missing.allowed to determine which peps to filter out #
    # the indices of 'result' are the same as the indices of 'group.vals.unique'
    groups = nonmiss_per_group$nonmiss_totals
    groups2 = groups[, - which(names(groups) %in% c(cname_id, NA, "NA", "<NA>", "NA."))] # remove the column with Peptide info


    # for ea peptide/row, need min(nonmiss_per_group$nonmiss_totals[,-1]) to be >= result
    trueFalse <- apply(groups2, 1, function(x) x>=result)
    trueFalse <- t(trueFalse)

    #temp = data.frame(nonmiss_per_group$nonmiss_totals)
    peps = as.character(groups[, 1])
    myresults <- data.frame(Peptide=peps, MeetsMinReq=trueFalse)
    myresults$remove0s <- rowSums(trueFalse) # 0's are the ones we'll filter out

    # which ones to filter out? only the rows that have all FALSE values
    myresults2 <- myresults[myresults$remove0s==0,]
    peps.filt <- as.character(myresults2$Peptide)
  }

  return(peps.filt)
}

