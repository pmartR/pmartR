#' Computes the Number of Non-Missing Data Points by Group
#'
#' This function computes the number of non-missing observations for samples, based on a group designation, for every peptide in the dataset
#'
#' @param omicsData an optional object of one of the classes "pepData", "proData", "metabData", or "lipidData", usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{lipidData}}, respectively. Either omicsData or all of e_data, groupDF, cname_id, and samp_id must be provided.
#' @param e_data an optional a \eqn{p \times n + 1} data.frame of expression data, where \eqn{p} is the number of lipids observed and \eqn{n} is the number of samples. Each row corresponds to data for each lipid. One column specifying a unique identifier for each lipid (row) must be present. Not required if omicsData is provided.
#' @param groupDF data.frame created by \code{group_designation} with columns for the sample identifier and the designated group. Not required if omicsData is provided.
#' @param cname_id character string specifying the name of the column containing the biomolecule identifiers in \code{e_data} and \code{e_meta} (if applicable). Not required if omicsData is provided.
#' @param samp_id character string specifying the name of the column containing the sample identifiers in \code{groupDF}. Not required if omicsData is provided.
#' @param check.names Logical, determines whether check.names argument of the data.frame function is TRUE or FALSE. Not required if omicsData is provided. 
#'
#' @return a list of length two. The first element giving the total number of possible samples for each group. The second element giving a data.frame with the first column giving the peptide and the second through kth columns giving the number of non-missing observations for each of the \code{k} groups.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' pep_object2 <- group_designation(omicsData = pep_object, main_effects = "Condition")
#' nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
#' nonmissing_result <- nonmissing_per_group(e_data = pep_object2$e_data, groupDF = attr(pep_object2, "group_DF"), cname_id = attr(pep_object2, "cnames")$edata_cname, samp_id = attr(pep_object2, "cnames")$fdata_cname)
#'}
#'
#' @author Lisa Bramer, Kelly Stratton
#'
#' 
#'
nonmissing_per_group <- function(omicsData = NULL, e_data = NULL, groupDF=NULL, cname_id=NULL, samp_id=NULL, check.names = NULL){

  # NOTE: Keep e_data & groupDF in this function because it gets called from within other functions and on the output of other functions that only return e_data (modifying all of those is too big of a task to handle, at least for the time being. -KS Jan 14, 2016)

  ## initial checks ##

  # check that omicsData is of an appropriate class #
  if(!is.null(omicsData) && !inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")



  if(!is.null(omicsData)){
    check_names = getchecknames(omicsData)
    
    e_data <- omicsData$e_data
    groupDF <- attr(omicsData, "group_DF")

    samp_id <- attr(omicsData, "cnames")$fdata_cname
    cname_id <- attr(omicsData, "cnames")$edata_cname
  }else{
    if(is.null(e_data)) stop("Either omicsData or e_data and groupDF must be provided.")
    if(is.null(groupDF)) stop("Either omicsData or e_data and groupDF must be provided.")

    if(is.null(samp_id)) stop("samp_id must be provided if omicsData is not")
    if(is.null(cname_id)) stop("cname_id must be provided if omicsData is not")
    if(is.null(check.names)) stop("check.names must be provided if omicsData is not")
    
    check_names = check.names
  }

  # check that groupDF is not NULL #
  if(is.null(groupDF)) stop("group_designation() must be run prior to imdanova_filter()")
  # check that the first column gives a unique sample id #
  if(nrow(groupDF) != length(unique(groupDF[, samp_id]))) stop("the first column of groupDF must be a unique sample identifier")


  ## end of initial checks ##


#  # melt the data by peptide and merge with groupDF then group by peptide and group #
#  melt.data = pre_imdanova_melt(e_data = e_data, groupDF = groupDF, samp_id = samp_id)
#  names(melt.data)[2] <- "Peptide" # set this explicitly now, and then before we return the function output, reset "Peptide" to cname_id
#  # class(melt.data) ## "grouped_dt" "tbl_dt"     "tbl"        "data.table" "data.frame" --> this doesn't work with dplyr::summarise
#  class(melt.data) <- c("grouped_df", "tbl_df", "tbl", "data.frame")

  # summarize total sample size per group #
  tot_samps = data.frame(dplyr::summarise(dplyr::group_by(groupDF, Group), n_group = n()))

#  # summarize the number of non-missing samples per group #
#  nonmiss_dat = dplyr::summarise(melt.data, non_miss = sum(!is.na(value)))

#  # put data back in wide format #
#  nonmiss_res = data.table::dcast.data.table(data.table::data.table(nonmiss_dat), Peptide~Group, value.var = "non_miss") # this command, using "Peptide", is the reason that we had to set the column name in melt.data to "Peptide" (using cname_id did not work)

#  #
#  nonmiss_res <- as.data.frame(nonmiss_res)

#  # format Peptide column as a character vector #
#  nonmiss_res[,"Peptide"] = as.character(nonmiss_res[,"Peptide"])

#  # [Kelly added 12/29/14] set the names(nonmiss_res) here, so the spaces and dashes don't get changed to periods (this messes things up in gtest_filter.R)
#  mynames <- names(nonmiss_res)
#  mynames[1] <- cname_id # re-set "Peptide" to whatever the cname_id is
#  nonmiss_totals <- data.frame(nonmiss_res)
#  names(nonmiss_totals) <- mynames

  
  group_dat<- as.character(groupDF$Group[order(groupDF$Group)])
  
  Mass_Tag_ID<- as.character(e_data[,cname_id])
  
  temp_data<- e_data[ , -which(names(e_data) %in% cname_id)]
  temp_data2<- temp_data[, match(groupDF[ ,samp_id], names(temp_data))]
  temp_data3<- temp_data2[,order(groupDF$Group)]
  
  
  nonmissing<- nonmissing_per_grp(as.matrix(temp_data3),group_dat)
  
  nonmissing<- data.frame(nonmissing, check.names = check_names)
  
  colnames(nonmissing)<- unique(group_dat)
  
  nonmiss_totals<- data.frame(Mass_Tag_ID,nonmissing,stringsAsFactors = FALSE, check.names = check_names)
  
  names(nonmiss_totals)[1] <- cname_id
  

  return(list(group_sizes = tot_samps, nonmiss_totals = nonmiss_totals))
}
