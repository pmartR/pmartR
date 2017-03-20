

#' Identifies the 'subset' of All Features
#'
#' Selects features for normalization via choosing all features currently in the data
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details This function returns the subset of all features. These will be used for normalization.
#'
#' @return Character vector containing all features.
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(pep_edata)
#' keep <- none_subset(e_data = pep_edata, edata_id = "Mass_Tag_ID")
#'}
#'
#' @author Kelly Stratton
#'



all_subset <- function(e_data, edata_id){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  # subset to all cases
  peps <- as.character(e_data[, edata_id_ind])

  if(length(peps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(peps)
}


#' Identify Features from the Top L Order Statistics
#'
#' Select features for normalization via the method of the top L order statistics (LOS)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the column giving the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param L numeric value between 0 and 1, indicating the top proportion of features to be retained (default value 0.05)
#'
#' @details The feature abundances of the top \code{L} order statistics are identified and returned. Specifically, for each sample, the features with the top \code{L} proportion of highest absolute abundance are retained, and the union of these features is taken as the subset identified.
#'
#' @return Character vector containing the features belonging to the subset.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_pepData)
#' pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
#' pep_subset <- los(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, mintR_groupDF = attr(pep_pepData2, "group_DF"))
#' }
#'
#' @author Kelly Stratton, Lisa Bramer
#'


los <- function(e_data, edata_id, L =.05){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  n_samps <- ncol(e_data[, -edata_id_ind])
  n_peps <- nrow(e_data[, -edata_id_ind])

  peps <- e_data[, edata_id_ind]
  row.names(e_data) <- peps
  mydata <- e_data[, -edata_id_ind]

  # get number of peps to keep from L that was provided
  num_kp <- round(L * n_peps)

  kp <- c()
  # for each sample, ... #
  for(i in 1:n_samps){
    # rank the features from highest to lowest (NAs are low) #
    cur_data <- mydata[, i]
    names(cur_data) <- peps

    # keep the top L percent of the features (excluding NAs, if NAs are included)
    ordered_data <- sort(abs(cur_data), decreasing = TRUE) # absolute value!

    # take the union of those features - that is the subset to return #
    kp <- c(kp, names(ordered_data)[1:num_kp])

  }

  kp <- unique(kp)

  if(length(kp)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(kp)
}


#' Identify Proportion of the Peptides Present (PPP) Biomolecules
#'
#' Selects features for normalization via the method of percentage of the peptides (or proteins, metabolites, etc.) present (PPP)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param proportion numeric value between 0 and 1, indicating the proportion at or above which a feature must be present across all samples in order to be retained (default value 0.5)
#'
#' @details Features present across \code{proportion} samples are designated as PPP.
#'
#' @return Character vector containing the features belonging to the PPP subset.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_pepData)
#' pep_subset <- ppp(e_data = pep_pepData$e_data, edata_id = attr(pep_pepData, "cnames")$edata_cname)
#'
#' pep_subset <- ppp(e_data = pep_pepData$e_data, edata_id = attr(pep_pepData, "cnames")$edata_cname, percent = 0.6)
#'}
#'
#' @author Kelly Stratton
#'



ppp <- function(e_data, edata_id, proportion=0.5){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  # subset to features present in at least proportion samples
  peps <- e_data[, edata_id_ind]
  row.names(e_data) <- peps
  mydata <- e_data[, -edata_id_ind]

  # get matrix of !is.na's #
  mydata_present <- !is.na(mydata)

  # get column of proportion present #
  mydata_pct_present <- rowSums(mydata_present)/ncol(mydata_present)

  # which features have proportion present above the value of "proportion" #
  inds <- which(mydata_pct_present >= proportion)
  mydata <- mydata[inds,]

  subset_peps <- as.character(row.names(mydata))

  if(length(subset_peps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(subset_peps)
}


#' Identify Proportion of Peptides Present (PPP) and Rank Invariant Peptides (RIP)
#'
#' Selects features for normalization via the method of proportion of peptides present and rank invariant peptides (ppp_rip)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with a column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param fdata_id character string indicating the name of the sample column name in f_data.
#' @param groupDF data.frame created by \code{group_designation} with columns for sample.id and group. If two main effects are provided the original main effect levels for each sample are returned as the third and fourth columns of the data.frame.
#' @param alpha numeric p-value threshold, above which the features are retained as rank invariant (default value 0.25)
#' @param proportion numeric value between 0 and 1, indicating the percentage at or above which a feature must be present across all samples in order to be retained (default value 0.5)
#'
#' @details Features present across \code{proportion} samples are subjected to a Kruskal-Wallis test (non-parametric one-way ANOVA, where NAs are ignored) on group membership, and those biomolecules with p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) are retained as rank-invariant biomolecules.
#'
#' @return Character vector containing the features belonging to the ppp_rip subset.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_pepData)
#' pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
#' pep_subset <- ppp_rip(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, samp_id = attr(pep_pepData2, "cnames")$fdata_cname, mintR_groupDF = attr(pep_pepData2, "group_DF"))
#'}
#'
#' @author Kelly Stratton
#'



ppp_rip <- function(e_data, edata_id, fdata_id, groupDF, alpha=0.2, proportion=0.5){

  samp_id = fdata_id
  
  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)
  
  # subset to features present in at least proportion samples
  peps <- e_data[, edata_id_ind]
  #row.names(e_data) <- peps
  mydata <- e_data#[, -edata_id_ind]
  
  # get matrix of !is.na's #
  mydata_present <- !is.na(mydata)
  
  # get column of proportion present #
  mydata_pct_present <- rowSums(mydata_present)/ncol(mydata_present)
  
  # which features have proportion present above the value of "proportion" #
  inds <- which(mydata_pct_present >= proportion)
  mydata <- mydata[inds,]
  peps <- peps[inds]
##########
#  pvals <- data.frame(rep(NA, nrow(mydata)))
#  for(i in 1:nrow(mydata)){
#    # check to see whether all observations are in the same group
#    nonmiss <- nonmissing_per_group(e_data = mydata[i,], groupDF= groupDF, cname_id=edata_id, samp_id=samp_id)$nonmiss_totals
#    if(sum(nonmiss[-1]==0) > (length(nonmiss[-1])-2)){
#      # if so, return 0 (we don't want to keep this one)
#      pvals[i,1] <- 0
#    }else{
#      # otherwise, return the p-value
#      pvals[i,1] <- kruskal.test(as.numeric(mydata[i,-1])~as.factor(groupDF$Group), na.action="na.omit")$p.value
#    }
#
#  }
#  #row.names(pvals) <- row.names(mydata)
##########
  
  #added 2/6/17 lines 223-231 iobani
  mydata = mydata[,-which(names(mydata) %in% edata_id)]
  group_dat = as.character(groupDF$Group[order(groupDF$Group)])
  rtemp = mydata[,match(names(mydata),groupDF[,samp_id])]
  rtemp2 = rtemp[,order(groupDF$Group)]
  
  # conduct K-W test using kw_rcpp function 
  pvals = kw_rcpp(as.matrix(rtemp2),group_dat)
  pvals= data.frame(pvals)
  
  RIPeps <- as.character(peps[as.numeric(pvals[,1]) > alpha])

  if(length(RIPeps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")


  return(RIPeps)
}



#' Identify Rank-Invariant Peptides
#'
#' Selects features for normalization via the method of rank-invariant peptides (RIP)
#'
#' @param e_data a \eqn{p \times n} data.frame, where \eqn{p} is the number of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of samples. Each row corresponds to data for a peptide, protein, lipid, or metabolite, with the first column giving the identifer name.
#' @param edata_id character string indicating the name of the peptide, protein, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param fdata_id character string indicating the name of the sample column name in f_data.
#' @param groupDF data.frame created by \code{group_designation} with columns for sample.id and group. If two main effects are provided the original main effect levels for each sample are returned as the third and fourth columns of the data.frame.
#' @param alpha numeric p-value threshold, above which the features are retained as rank invariant (default value 0.25)
#'
#' @details Features with complete data are subjected to a Kruskal-Wallis test (non-parametric one-way ANOVA) on group membership, and those features with p-value greater than a defined threshold \code{alpha} (common values include 0.1 or 0.25) are retained as rank-invariant features.
#'
#' @return Character vector containing the features belonging to the RIP subset.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_pepData)
#' pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
#' pep_subset <- rip(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, mintR_groupDF = attr(pep_pepData2, "group_DF"))
#'}
#'
#' @author Kelly Stratton
#'



rip <- function(e_data, edata_id, fdata_id, groupDF, alpha=.2){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  # subset to complete cases (features with complete data)
  peps <- e_data[, edata_id_ind]
  row.names(e_data) <- as.character(peps)
  mydata <- e_data[, -edata_id_ind]
  inds <- which(complete.cases(mydata) == TRUE)
  mydata <- mydata[inds,]
 
##########
#  pvals <- data.frame(rep(NA, nrow(mydata)))
#  for(i in 1:nrow(mydata)){
#    pvals[i,1] <- kruskal.test(as.numeric(mydata[i,])~as.factor(groupDF$Group))$p.value
#  }
##########
  
  #added 2/6/17 iobani
  samp_id = fdata_id
  group_dat = as.character(groupDF$Group[order(groupDF$Group)])
    
  rtemp =  mydata[,match(names(mydata),groupDF[,samp_id])]
  rtemp2 =  rtemp[,order(groupDF$Group)]
    
  # conduct K-W test on un-normalized data # #used kw_rcpp function 
  pvals = kw_rcpp(as.matrix(rtemp2),group_dat)#here we use the kw_rcpp function to calculate pvals#
  pvals = data.frame(pvals)
 
  
  row.names(pvals) <- row.names(mydata)

  RIPeps <- as.character(row.names(pvals)[as.numeric(pvals[,1]) > alpha])

  if(length(RIPeps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(RIPeps)
}
