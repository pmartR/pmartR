#' Tests for the independence of missing data across groups (aka factors, aka main effects)
#'
#' Tests the null hypothesis that the number of missing observations is independent of the groups.  A g-test is used to test this null hypothese
#' against the alternative that the groups and missing data are related.  This is usually performed in conjuction with an ANOVA which tests if
#' the mean response (which varies with data type) is the same across groups; this combination is called IMD_ANOVA.  It's probably a good idea
#' to first filter the data with `imd_anova_filter` to see if there is enough infomration to even do this test.  See Webb-Robertson et al. (2010) for more.
#'
#' @param omicsData A pmartR data object of any class
#' @param comparisons `data.frame` with columns for "Control" and "Test" containing the different comparisons of interest. Comparisons will be made between the Test and the corresponding Control  If left NULL, then all pairwise comparisons are executed.
#' @param pval_adjust character vector specifying the type of multiple comparisons adjustment to implement. The default setting is to not apply an adjustment. Valid options include: holm, bonferonni and none. See \code{\link{p.adjust}} for some details.
#' @param pval_thresh numeric p-value threshold, below or equal to which peptides are considered differentially expressed. Defaults to 0.05
#'
#' @return  a list of `data.frame`s
#' \tabular{ll}{
#' Results  \tab e_data cname, Count of non-missing data for each group, Global G-test statistic and p-value\cr
#'  \tab \cr
#' Gstats  \tab Value of the g statistics for each of the pairwise comparisons specified by the `comparisons` argument \cr
#'  \tab \cr
#' Pvalues  \tab p-values for each of the pairwise comparisons specified by `comparisons` argument \cr
#'  \tab \cr
#' Flags  \tab Indicator of statistical significance where the sign of the flag reflects the difference in the ratio of non-missing observations (0/+-2 to if adjusted p-value>=pval_thresh or p-value<pval_thresh) \cr
#'  }
#'
#' @author Bryan Stanfill
#'
#' @references
#' Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical analyses of peptide intensities and peptide occurrences improves identification of significant peptides from MS-based proteomics data." Journal of proteome research 9.11 (2010): 5748-5756.
#'
imd_test <- function(omicsData, comparisons=NULL, pval_adjust = 'none', pval_thresh = 0.05){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData", "nmrData"))) stop("omicsData is not an object of appropriate class")

  # Check for group_DF attribute #
  if(!("group_DF" %in% names(attributes(omicsData)))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attributes(omicsData)$group_DF
  }

  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("This test cannot be performed with less than two groups.")
  }

  ############
  #Check for gtest filter
  if(is.null(attr(omicsData,"imdanova"))){
    warning("These data haven't been filtered, see `?imdanova_filter` for details.")
  }else{
    cnames <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
    filterrows <- which(cnames%in%attr(omicsData,"imdanova")$test_with_gtest)
    if(length(filterrows)>0)
      omicsData$e_data <- omicsData$e_data[filterrows,]
  }

  #Remove rows in "groupData" that don't have corresponding columns in 'omicsData$edata'
  samp_cname <- attr(omicsData,"cnames")$samp_cname
  if(is.null(samp_cname)){ #Added becuase "samp_cname" has been replaced with "fdata_cname"
    samp_cname <- attr(omicsData,"cnames")$fdata_cname
  }
  group_samp_cname <- groupData[,samp_cname] #The sampleIDs in groupData (possibly more than is in the e_data)
  groupData <- groupData[group_samp_cname%in%colnames(omicsData$e_data),] #Only keep the rows of groupData that have columns in e_data
  #groupData <- dplyr::filter(groupData,SampleID%in%colnames(omicsData$e_data)) #This is faster but too specific to a cname of SampleID

  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData[,samp_cname])]

  #Find the number of times each peptide was absent for each of the experimental groups, matrix of CAk values
  gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))

  #Number of samples associated with each group
  nK <- as.numeric(table(gp))
  label_map <- data.frame(Num_label=1:k,Group=unique(groupData$Group),nK=nK)

  #----- Global IMD test - any pattern in missing across all groups? -----#

  #Rcpp::sourceCpp('src/count_missing_cpp.cpp') #Run for debugging
  absent <- count_missing_cpp(data.matrix(data),gp)
  colnames(absent) <- label_map$Group

  #Number of times each peptide was observed for each of the experimental groups, matrix of COk values
  observed <- matrix(label_map$nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)
  observed <- observed-absent
  colnames(observed) <- label_map$Group

  #compute total number of samples from which each peptide is absent (mA) or observed (mO)
  mA <- rowSums(absent)
  mO <- rowSums(observed)
  N <- mA+mO

  #Expected number of absences (EAk) and observations (EOk) for each group
  EAk <- (mA/N)*matrix(label_map$nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)
  EOk <- (mO/N)*matrix(label_map$nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)

  #Sum is only taken over non-zero counts; zero counts are turned into NaN so they need to be zeroed out
  Gks <- observed*log(observed/EOk)
  Gks[which(is.nan(Gks))] <- 0
  GksA <- absent*log(absent/EAk)
  GksA[which(is.nan(GksA))] <- 0
  Gks <- Gks+GksA

  #Compute test statistic but summing over all nonzero counts (na.rm=TRUE removes the zero counts)
  Global_Gk <- 2*rowSums(Gks)

  #Put the statistic, degrees of freedom and p-value into a data frame
  Global_results <- data.frame(Statistic=Global_Gk, p_value=pchisq(Global_Gk,df=(k-1),lower.tail = FALSE))

  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  colnames(observed) <- paste("Count",colnames(observed),sep='_')
  Global_results <- cbind(omicsData$e_data[edatacname],observed,Global_results)


  #----- Pairwise IMD test - any qualitative difference between groups? -----#
  #Compute test statistic but summing over all nonzero counts (na.rm=TRUE removes the zero counts)
  #source('~/Documents/MinT/MSomics/MSomicsSTAT/R/group_comparison.R') #Run to debug
  if(k==2){
    #If group_comparison_imd doesn't return Gstats, give them global Gstats
    #pairwise_stats <- cbind(omicsData$e_data[edatacname],Global_Gk)
    #pairwise_pvals <- cbind(omicsData$e_data[edatacname],Global_results$p_value)
    pairwise_stats <- data.frame(Global_Gk)
    pairwise_pvals <- data.frame(Global_results$p_value)
    if(is.null(comparisons)){
      colnames(pairwise_stats)[1] <- colnames(pairwise_pvals)[1] <- paste(label_map$Group[1],label_map$Group[2],sep="_vs_")
    }else{
      colnames(pairwise_stats)[1] <- colnames(pairwise_pvals)[1] <- paste(comparisons$Test[1],comparisons$Control[1],sep="_vs_")
    }

    #Signs for group comparisons
    ratio_observed <- observed/(observed+absent)
    cmat <- create_c_matrix(group_df = groupData,to_compare_df = comparisons)
    ratio_diff <- ratio_observed%*%t(cmat$cmat)
    colnames(ratio_diff) <- cmat$names

  }else{
    #source('~/Documents/MinT/MSomics/MSomicsSTAT/R/group_comparison.R') #Run for debugging
    pairwise_gk <- group_comparison_imd(groupData,comparisons,observed,absent)

    pairwise_stats <- data.frame(pairwise_gk$GStats)
    pairwise_pvals <- data.frame(pairwise_gk$pvalues)
    colnames(pairwise_pvals) <- colnames(pairwise_gk$pvalues)
    ratio_diff <- pairwise_gk$signs
    ##Add edata_cname
    #pairwise_stats <- cbind(omicsData$e_data[edatacname],pairwise_gk$GStats)
    #pairwise_pvals <- cbind(omicsData$e_data[edatacname],pairwise_gk$pvalues)
  }

  #----- P-value adjustment -----#
  #Match provided 'pval_adjust' to available options
  pval_adjust <- try(match.arg(tolower(pval_adjust),c("holm","bonferroni","none","tukey","dunnett")),silent=TRUE)

  if(class(pval_adjust)=='try-error')
    stop("Provided 'pval_adjust' argument is invalid, please select 'holm', 'bonferroni', 'tukey', 'dunnett' or 'none'.")

  if(pval_adjust%in%c("tukey","dunnett")){
    warning("Tukey and Dunnett corrections aren't available for the IMD test, replaced by Holm")
    pval_adjust <- 'holm'
  }

  #Implement Bonferonni correction by multiplying by number of tests if requested, otherwise do nothing

  if(ncol(pairwise_pvals)==1)
    pval_adjust <- 'none'

  if(pval_adjust=="bonferroni"){
    adjusted_pvals <- pmin(data.matrix(pairwise_pvals*(ncol(pairwise_pvals)-1)),1)
  }else if(pval_adjust=="holm"){
    adjusted_pvals <- t(apply(pairwise_pvals,1,ranked_holm_cpp))
    colnames(adjusted_pvals) <- colnames(pairwise_pvals)
  }else{
    adjusted_pvals <- pairwise_pvals
  }

  #----- Create significance flag matrix -----#
  sig_indicators <- data.matrix(adjusted_pvals)
  sigs <- which(sig_indicators<pval_thresh)
  if(length(sigs)>0){
    sig_indicators[sigs] <- 2*sign(ratio_diff[sigs])
    sig_indicators[-sigs] <- 0
  }else{
    sig_indicators[which(sig_indicators>=pval_thresh)] <- 0
  }

  sig_indicators <- data.frame(sig_indicators)
  #sig_indicators <- data.frame(adjusted_pvals[,1],sig_indicators)
  #colnames(sig_indicators)[1] <- colnames(adjusted_pvals)[1]

  return(list(Results=Global_results,Gstats=pairwise_stats,Pvalues=adjusted_pvals,Flags=sig_indicators))
}






