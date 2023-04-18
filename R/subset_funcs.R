#' Identify All Biomolecules for Use in Normalization
#'
#' Selects biomolecules for normalization via choosing all biomolecules
#' currently in the data
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details This function returns the subset of all biomolecules. These will be
#'   used for normalization.
#'
#' @return Character vector containing all biomolecules.
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

#' Identify Biomolecules from the Top L Order Statistics for Use in
#' Normalization
#'
#' Select biomolecules for normalization via the method of the top L order
#' statistics (LOS)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the column giving the
#'   peptide, protein, lipid, or metabolite identifier. Usually obtained by
#'   calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param L numeric value between 0 and 1, indicating the top proportion of
#'   biomolecules to be retained (default value 0.05)
#'
#' @details The biomolecule abundances of the top \code{L} order statistics are
#'   identified and returned. Specifically, for each sample, the biomolecules with
#'   the top \code{L} proportion of highest absolute abundance are retained, and
#'   the union of these biomolecules is taken as the subset identified.
#'
#' @return Character vector containing the biomolecules belonging to the subset.
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

#' Identify Biomolecules from the Proportion Present (PPP) for Use in
#' Normalization
#'
#' Selects biomolecules for normalization via the method of percentage of the
#' peptides (or proteins, metabolites, etc.) present (PPP)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#' @param proportion numeric value between 0 and 1, indicating the proportion at
#'   or above which a biomolecule must be present across all samples in order to be
#'   retained (default value 0.5)
#'
#' @details Biomolecules present across \code{proportion} samples are designated as
#'   PPP.
#'
#' @return Character vector containing the biomolecules belonging to the PPP
#'   subset.
#'   
#' @author Kelly Stratton
#'
ppp <- function(e_data, edata_id, proportion=0.5){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  # get column of proportion present #
  mydata_pct_present <- rowSums(!is.na(e_data[, -edata_id_ind])) /
    ncol(e_data[, -edata_id_ind])

  # which features have proportion present above the value of "proportion" #
  inds <- which(mydata_pct_present >= proportion)

  subset_peps <- as.character(e_data[inds, edata_id_ind])

  if(length(subset_peps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(subset_peps)
}

#' Identify Proportion of Peptides Present (PPP) and Rank Invariant Peptides
#' (RIP) for Use in Normalization
#'
#' Selects biomolecules for normalization via the method of proportion of
#' biomolecules present and rank invariant biomolecules (ppp_rip)
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#' @param fdata_id character string indicating the name of the sample column
#'   name in f_data.
#' @param groupDF data.frame created by \code{group_designation} with columns
#'   for sample.id and group. If two main effects are provided the original main
#'   effect levels for each sample are returned as the third and fourth columns
#'   of the data.frame.
#' @param alpha numeric p-value threshold, above which the biomolecules are retained
#'   as rank invariant (default value 0.25)
#' @param proportion numeric value between 0 and 1, indicating the percentage at
#'   or above which a biomolecule must be present across all samples in order to be
#'   retained (default value 0.5)
#'
#' @details Biomolecules present across \code{proportion} samples are subjected to a
#'   Kruskal-Wallis test (non-parametric one-way ANOVA, where NAs are ignored)
#'   on group membership, and those biomolecules with p-value greater than a
#'   defined threshold \code{alpha} (common values include 0.1 or 0.25) are
#'   retained as rank-invariant biomolecules.
#'
#' @return Character vector containing the biomolecules belonging to the ppp_rip subset.
#'
#' @author Kelly Stratton
#'
ppp_rip <- function(e_data, edata_id, fdata_id, groupDF, alpha=0.2, proportion=0.5){

  samp_id = fdata_id
  
  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)
  
  # subset to features present in at least proportion samples
  peps <- e_data[, edata_id_ind]
  
  # get column of proportion present #
  mydata_pct_present <- rowSums(!is.na(e_data[, -edata_id_ind])) /
    ncol(e_data[, -edata_id_ind])
  
  # which features have proportion present above the value of "proportion" #
  inds <- which(mydata_pct_present >= proportion)
  mydata <- e_data[inds, -edata_id_ind]
  peps <- peps[inds]
  
  # make sure groups are aligned with the columns of e_data
  reorder = match(colnames(mydata), as.character(groupDF[,fdata_id]))
  group_dat = as.character(groupDF[reorder,]$Group)
  
  # conduct K-W test using kw_rcpp function 
  pvals = kw_rcpp(as.matrix(mydata), group_dat)
  pvals = data.frame(pvals)
  
  RIPeps <- as.character(peps[as.numeric(pvals[, 1]) > alpha &
                                !is.na(as.numeric(pvals[, 1]))])

  if(length(RIPeps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")


  return(RIPeps)
  
}

#' Identify Rank-Invariant Biomolcules for Use in Normalization
#'
#' Selects biomolecules for normalization via the method of rank-invariant
#' biomolcules (RIP)
#'
#' @param e_data a \eqn{p \times n} data.frame, where \eqn{p} is the number of
#'   peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#' @param fdata_id character string indicating the name of the sample column
#'   name in f_data.
#' @param groupDF data.frame created by \code{group_designation} with columns
#'   for sample.id and group. If two main effects are provided the original main
#'   effect levels for each sample are returned as the third and fourth columns
#'   of the data.frame.
#' @param alpha numeric p-value threshold, above which the biomolecules are
#'   retained as rank invariant (default value 0.25)
#'
#' @details Biomolecules with complete data are subjected to a Kruskal-Wallis
#'   test (non-parametric one-way ANOVA) on group membership, and those
#'   biomolecules with p-value greater than a defined threshold \code{alpha}
#'   (common values include 0.1 or 0.25) are retained as rank-invariant
#'   biomolecules.
#'
#' @return Character vector containing the biomolecules belonging to the RIP
#'   subset.
#'
#' @author Kelly Stratton
#'
rip <- function(e_data, edata_id, fdata_id, groupDF, alpha=.2){

  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)

  # subset to complete cases (features with complete data)
  peps <- e_data[, edata_id_ind]
  inds <- which(complete.cases(e_data[, -edata_id_ind]) == TRUE)
  mydata <- e_data[inds, -edata_id_ind]
  peps <- peps[inds]
  
  #added 2/6/17 iobani
  group_dat = as.character(groupDF$Group[order(groupDF$Group)])
  mydata <- mydata[, order(groupDF$Group)]
    
  # conduct K-W test on un-normalized data, used kw_rcpp function 
  pvals = kw_rcpp(as.matrix(mydata), group_dat)
  pvals = data.frame(pvals)

  RIPeps <- as.character(peps[as.numeric(pvals[, 1]) > alpha &
                                !is.na(as.numeric(pvals[, 1]))])

  if(length(RIPeps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")

  return(RIPeps)
}

#' Identify biomolecules with no missing values across samples
#'
#' Selects biomolecules that have complete rows in e_data, equivalent to 'ppp'
#' with proportion = 1.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame, where \eqn{p} is the number
#'   of peptides, proteins, lipids, or metabolites and \eqn{n} is the number of
#'   samples. Each row corresponds to data for a peptide, protein, lipid, or
#'   metabolite, with one column giving the biomolecule identifier name.
#' @param edata_id character string indicating the name of the peptide, protein,
#'   lipid, or metabolite identifier. Usually obtained by calling
#'   \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @return Character vector containing the biomolecules with no missing values
#'   across all samples.
#'   
complete_mols <- function(e_data, edata_id){
  
  # pull off the column for edata_id
  edata_id_ind <- which(colnames(e_data)==edata_id)
  
  # rows with no missing values
  complete_inds <- which(complete.cases(e_data[, -edata_id_ind]) == TRUE)
  
  # retain only character peptide names for complete rows
  peps <- as.character(e_data[complete_inds, edata_id_ind])
  
  if(length(peps)<2) stop("There are <2 biomolecules in the subset; cannot proceed.")
  
  return(peps)
}
