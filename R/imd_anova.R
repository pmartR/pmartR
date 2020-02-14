#' Tests for a qualitative and quantiative difference between groups using IMD and ANOVA, respectively
#'
#'  This is IMD-ANOVA test proposed in Webb-Robertson et al. (2010).
#'
#' @param omicsData A pmartR data object of any class, which has a `group_df` attribute that is usually created by the `group_designation()` function
#' @param comparisons data.frame with columns for "Control" and "Test" containing the different comparisons of interest. Comparisons will be made between the Test and the corresponding Control  If left NULL, then all pairwise comparisons are executed.
#' @param test_method character string specifying the filter method to use: "combined", "gtest", or "anova". "combined" implements both the gtest and anova filters.
#' @param pval_adjust character vector specifying the type of multiple comparisons adjustment to implement via \code{\link{p.adjust}}. A NULL value corresponds to no adjustment. Valid options for ANOVA include: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none. Valid options for g-test include: holm, bonferonni and none. See \code{\link{p.adjust}} for some details.
#' @param pval_thresh numeric p-value threshold, below or equal to which peptides are considered differentially expressed. Defaults to 0.05
#' @param covariates data.frame similar to \code{groupData} consisting of two columns: the sample ID variable (with names matching column names in \code{omicsData$e_data}) and a column containing the numeric or group data for each sample
#' @param paired logical - should the data be paired or not; if TRUE the attribute "pairing" cannot be `NULL`
#' 
#' @return  a list of \code{data.frame} objects 
#' \tabular{ll}{
#' Full_Results \tab Columns: e_data cname, group counts, group means, ANOVA p-values, IMD p-values, fold change estimates, fold change significance flags\cr
#'  \tab \cr
#'  Flags (signatures)  \tab Indicator of statistical significance for one (+/-1 for ANOVA, +/-2 for g-test) or neither (0) test \cr
#'  \tab \cr
#' p-values  \tab p-values for pairwise comparisons adjusted for within, e.g., peptide multiple comparisons only (not across peptides yet) \cr
#' \tab \cr
#' groupData  \tab Mapping between sample ID and group \cr
#'  }
#'
#' @author Bryan Stanfill, Kelly Stratton
#' @references 
#' Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical analyses of peptide intensities and peptide occurrences improves identification of significant peptides from MS-based proteomics data." Journal of proteome research 9.11 (2010): 5748-5756.
#'
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' #Transform the data
#' mypepData <- edata_transform(omicsData = pep_object, data_scale = "log2")
#' 
#' #Group the data by condition
#' mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))
#' 
#' #Apply the IMD ANOVA filter
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
#' 
#' #Implement the IMD ANOVA method and compute all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
#' anova_res <- imd_anova(omicsData = mypepData, test_method = 'anova')
#' imd_res <- imd_anova(omicsData = mypepData, test_method = 'gtest')
#' imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='bon')
#' 
#' #Test with really big dataset
#' library(OvarianPepdataBP)
#' tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp$e_data, f_data = tcga_ovarian_pepdata_bp$f_data, e_meta = tcga_ovarian_pepdata_bp$e_meta)
#' tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("race"))
#' ovarian_res <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, test_method = 'anova')
#' #Tukey adjustment is super slow right now because "ptukey" is super slow, not sure how to fix that
#' ovarian_res_holm <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'holm', test_method='gtest')
#' #Dunnett adjustment, should give an error because dunnett correction shouldn't be applied for all pairwise comparisons
#' ovarian_res_dunnett <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'dunnett', test_method='combined')
#' 
#' #Test really big dataset, two factors
#' tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("vital_status","neoplasm_histologic_grade"))
#' ovarian_res_twofac <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, test_method='comb')
#' 
#' #Same but only test main effects (Dead vs Alive, G2 vs G3)
#' comp_df <- data.frame(Control=c("Alive","G2"), Test=c("Dead","G3"))
#' ovarian_res_twofac_main_effects <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, comparisons = comp_df, test_method='comb')
#' }
#' @export
#' 
imd_anova <- function(omicsData, comparisons = NULL, test_method, pval_adjust = 'none', pval_thresh = 0.05, covariates = NULL, paired = FALSE, equal_var = TRUE){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData"))) stop("omicsData is not an object of appropriate class")
  
  # Check for group_DF attribute #
  if(!("group_DF" %in% names(attributes(omicsData)))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attributes(omicsData)$group_DF
  }
  
  #Match provided 'test_method' to available options
  test_method <- try(match.arg(tolower(test_method),c("combined","gtest","anova")),silent=TRUE)
  
  if(!(test_method%in%c("combined","gtest","anova"))){
    #If test-method isn't valid, stop and tell them
    stop("Provided 'test_method' argument is invalid, please select 'anova','gtest' or 'combined'.")
  }

  # Check for log transform #
  if(!(attr(omicsData,"data_info")$data_scale%in%c("log2","log","log10"))&(test_method%in%c("combined","anova"))){
    stop("Data must be log transformed in order to implement ANOVA.")
  }
  
  ############
  #Check for anova filter - give warning if not present then let `imd_test` and `anova_test` do the actual filtering
  if(is.null(attr(omicsData,"imdanova"))){
    warning("These data haven't been filtered, see `?imdanova_filter` for details.")
    #Add attribute so imd_test and anova_test don't return same warning
    attr(omicsData,"imdanova")$test_with_anova <- "No IMD ANOVA Attribute"
  }#else{
  #  cnames <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
  #  filterrows <- which(cnames%in%attr(omicsData,"imdanova")$test_with_anova)
  #  if(length(filterrows)>0) #Remove rows that need to be filtered
  #    omicsData$e_data <- omicsData$e_data[filterrows,]
  #}
  
  ############
  #Use imd_test() to test for independence of missing data (qualitative difference between groups)
  if(test_method=='anova'){
    #If they don't want the g-test done, save some time by removing comparisons and the pval_adjust arguments
    #Also make the gtest_pvalues NULL so nothing's returned
    imd_results_full <- imd_test(omicsData, comparisons = NULL, pval_adjust = 'none', pval_thresh = pval_thresh)
  }else{
    imd_results_full <- imd_test(omicsData, comparisons = comparisons, pval_adjust = pval_adjust, pval_thresh = pval_thresh)
    gtest_pvalues <- imd_results_full$Pvalues
    colnames(gtest_pvalues) <- paste0("P_value_G_",colnames(gtest_pvalues))
    gtest_flags <- imd_results_full$Flags
    colnames(gtest_flags) <- paste0("Flag_",colnames(imd_results_full$Pvalues))
  }
  
  #Get cname and counts by group from imd_results_full
  imd_counts <- imd_results_full[[1]][,c(1,grep("^Count",colnames(imd_results_full[[1]])))] 
  
  
  ############
  #Use anova_test() to test for a global equality of means, test requested comparisons and 
  #compute foldchanges for those comaprisons (quantitative difference between groups)
  if(test_method=='gtest'){
    #If they only want the g-test, used anova_test to compute means, fold changes but 
    #don't return flags or p-values
    anova_results_full <- anova_test(omicsData, comparisons = comparisons, pval_adjust = 'none', pval_thresh = pval_thresh, covariates = covariates, paired = paired, equal_var = equal_var)
  }else{
    anova_results_full <- anova_test(omicsData, comparisons = comparisons, pval_adjust = pval_adjust, pval_thresh = pval_thresh, covariates = covariates, paired = paired, equal_var = equal_var)
    anova_fold_flags <- anova_results_full$Flags
    colnames(anova_fold_flags) <- paste0("Flag_",colnames(anova_fold_flags))
    anova_pvalues <- anova_results_full$Fold_change_pvalues
    colnames(anova_pvalues) <- paste0("P_value_T_",colnames(anova_pvalues))
  }
  ###########
  #Group means and fold changes should be returned even if they only ask for g-test
  anova_results <- anova_results_full[[1]][,c(1,grep("Mean",colnames(anova_results_full[[1]])))]
  anova_fold_change <- anova_results_full$Fold_changes
  colnames(anova_fold_change) <- paste0("Fold_change_",colnames(anova_fold_change))
  
  
  #########
  #Return appropriate results for each test_method
  
  if(test_method=='anova'){

    all_anova <- cbind(anova_results,anova_pvalues,anova_fold_change,anova_fold_flags)
    Full_results <- right_join(imd_counts,all_anova)
    if(nrow(Full_results)>max(nrow(imd_counts),nrow(all_anova))){
      stop("Combining g-test and ANOVA results failed.")
    }
    
    imd_out <- list(Full_results=Full_results, Flags = NULL,P_values = NULL)
    #Remove p-value and flag indicators from column names
    colnames(anova_pvalues) <- gsub("^P_value_T_","",colnames(anova_pvalues))
    imd_out$P_values <- anova_pvalues
    colnames(anova_fold_flags) <- gsub("^Flag_","",colnames(anova_fold_flags))
    imd_out$Flags <- anova_fold_flags

    final_out <- statRes_output(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh) 
    attr(final_out, "cnames") = attr(omicsData, "cnames")
    attr(final_out, "data_class") = attr(omicsData, "class")
    return(final_out)
    
  }else if(test_method=='gtest'){
    
    all_gtest <- cbind(imd_counts,gtest_pvalues, gtest_flags)
    all_anova <- cbind(anova_results,anova_fold_change)
    Full_results <- left_join(all_gtest, all_anova)
    if(nrow(Full_results)>max(nrow(all_gtest),nrow(all_anova))){
      stop("Combining g-test and ANOVA results failed.")
    }

    imd_out <- list(Full_results = Full_results,Flags = NULL, P_values = NULL)
    #Remove p-value and flag indicators from column names
    colnames(gtest_pvalues) <- gsub("^P_value_G_","",colnames(gtest_pvalues))
    imd_out$P_values <- gtest_pvalues
    colnames(gtest_flags) <- gsub("^Flag_","",colnames(gtest_flags))
    imd_out$Flags <- gtest_flags

    final_out <- statRes_output(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh) 
    attr(final_out, "cnames") = attr(omicsData, "cnames")
    attr(final_out, "data_class") = attr(omicsData, "class")
    return(final_out)
  }
  
  #####-----Determine which p-values/flags to return--------######
  
  #Create the list that is returned
  all_gtest <- cbind(imd_counts,gtest_pvalues)
  all_anova <- cbind(anova_results,anova_pvalues,anova_fold_change,anova_fold_flags)
  Full_results <- dplyr::full_join(all_gtest,all_anova)
  
  #Get counts for rows in "anova_results" but not in "imd_counts"
  final_cnts <- Full_results[,grep("Count",colnames(Full_results))]
  msng_cnts <- which(is.na(rowSums(final_cnts)))
  if(length(msng_cnts)>0){
    to_fix <- Full_results[msng_cnts,get_edata_cname(omicsData)]
    omicsData2 <- omicsData
    omicsData2$e_data <- omicsData$e_data%>%dplyr::filter(!!rlang::sym(get_edata_cname(omicsData))%in%as.character(to_fix))
    new_cnts <- imd_test(omicsData = omicsData2, comparisons = NULL, pval_adjust = 'none', pval_thresh = pval_thresh)
    rm(omicsData2)
    #Replace the NA counts with the correct counts
    Full_results[msng_cnts,grep("Count",colnames(Full_results))] <- new_cnts$Results[,grep("Count",colnames(new_cnts$Results))]
  }
  
  ##I added this catch for some reason, but I don't know why and it stops some perfectly fine code
  #if(nrow(Full_results)>max(nrow(all_gtest),nrow(all_anova))){
  #  stop("Combining g-test and ANOVA results failed.")
  #}
  
  #If both IMD and ANOVA are done do the following (from email from Kelly on 7/28 based on convo. with BJ same day)
  #For reporting p-values when the user does both the g-test and ANOVA,
  #(1)   If ANOVA p-value but no g-test p-value, report ANOVA p-value (+/-1 for flag)
  #(2)   If g-test p-value but no ANOVA p-value, report g-test p-value (+/-2 for flag)
  #(3)   If neither, report NA (NA for flag)
  #(4)   If both but neither is significant, report ANOVA (+/-1 for flag)
  #(5)   If both but ANOVA is significant, report ANOVA (+/-1 for flag)
  #(6)   If both but g-test is significant, report g-test (+/-2 for flag)
  
  #Start with the anova p-values (extracted from combined results) and replace if missing [see (2)] or g-test is significant [see (6)]
  imd_pvals <- data.matrix(Full_results[grep("^P_value_T_",colnames(Full_results))])
  imd_flags <- data.matrix(Full_results[grep("^Flag_",colnames(Full_results))])
  
  #Redefine gtest_pvalues after it's been merged
  gtest_pvalues <- data.matrix(Full_results[grep("^P_value_G_",colnames(Full_results))])
  
  #Same for 'gtest_flags', but a little more complicated
  gtest_flags <- imd_results_full$Flags
  colnames(gtest_flags) <- paste0("Flag_",colnames(imd_results_full$Pvalues))
  gtest_flags <- cbind(imd_counts[,1],gtest_flags)
  colnames(gtest_flags)[1] <- colnames(imd_counts)[1]
  gtest_flags <- dplyr::full_join(anova_results,gtest_flags)
  
  # Everything -but- gtest flags is aligned with the ID column from Full_results, reorder it here
  reorder = match(Full_results[,get_edata_cname(omicsData)], gtest_flags[,get_edata_cname(omicsData)])
  gtest_flags <- gtest_flags[reorder,]
  gtest_flags <- data.matrix(gtest_flags[grep("^Flag_",colnames(gtest_flags))])
  
  #Replance missing ANOVA p-values with g-test p-values
  if(any(is.na(imd_pvals))){
    anova_NAs <- which(is.na(imd_pvals))
    imd_pvals[anova_NAs] <- gtest_pvalues[anova_NAs]
    imd_flags[anova_NAs] <- gtest_flags[anova_NAs]
  }
  
  #Replace insignificant ANOVA p-values with significant g-test
  insig_anova <- which(imd_pvals>pval_thresh) #Insignificant anova_pvalues
  if(any(insig_anova)){
    sig_gtest <- which(gtest_pvalues<=pval_thresh) #Significant gtest_pvalues
    overlap <- sig_gtest[sig_gtest%in%insig_anova] #Overlap between the two groups
    if(any(overlap)){
      imd_pvals[overlap] <- gtest_pvalues[overlap] #Replace insignificant p-values & flags with significant gtest ones
      imd_flags[overlap] <- gtest_flags[overlap]
    }
  }
  
  #Remove p-value and flag indicators from column names
  colnames(imd_pvals) <- gsub("^P_value_T_","",colnames(imd_pvals))
  colnames(imd_flags) <- gsub("^Flag_","",colnames(imd_flags))
  
  # R is fidgety about storing a single column dataframe to a column position
  Full_results[,grep("^Flag_",colnames(Full_results))] <- if(ncol(imd_flags) == 1) as.vector(imd_flags) else imd_flags
  
  ##-------- Construct statRes object --------##
  imd_out <- list(Full_results=Full_results,
                  Flags = data.frame(imd_flags, check.names = FALSE),
                  P_values = data.frame(imd_pvals, check.names = FALSE))
  
  final_out <- statRes_output(imd_out,omicsData,comparisons,test_method,pval_adjust,pval_thresh)
  
  attr(final_out, "cnames") = attr(omicsData, "cnames")
  attr(final_out, "data_class") = attr(omicsData, "class")
  return(final_out)
}

