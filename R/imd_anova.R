#' Tests for a qualitative and quantitative difference between groups using IMD
#' and ANOVA, respectively
#'
#' This is IMD-ANOVA test proposed in Webb-Robertson et al. (2010).
#'
#' @param omicsData A pmartR data object of any class, which has a `group_df`
#'   attribute that is usually created by the `group_designation()` function
#' @param comparisons data.frame with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control. If left NULL, then all
#'   pairwise comparisons are executed.
#' @param test_method character string specifying the filter method to use:
#'   "combined", "gtest", or "anova". "combined" implements both the gtest and
#'   anova filters.
#' @param pval_adjust_a A character string specifying the type of multiple
#'   comparison adjustment to implement for ANOVA tests. Valid options include:
#'   "bonferroni", "holm", "tukey", and "dunnett". The default is "none" which
#'   corresponds to no p-value adjustment.
#' @param pval_adjust_g A character string specifying the type of multiple
#'   comparison adjustment to implement for G-test tests. Valid options include:
#'   "bonverroni" and "holm". The default is "none" which corresponds to no
#'   p-value adjustment.
#' @param pval_thresh numeric p-value threshold, below or equal to which
#'   peptides are considered differentially expressed. Defaults to 0.05
#' @param covariates data.frame similar to \code{groupData} consisting of two
#'   columns: the sample ID variable (with names matching column names in
#'   \code{omicsData$e_data}) and a column containing the numeric or group data
#'   for each sample
#' @param paired logical - should the data be paired or not; if TRUE the
#'   attribute "pairing" cannot be `NULL`
#'
#' @return  a list of \code{data.frame} objects \tabular{ll}{ Full_Results \tab
#'   Columns: e_data cname, group counts, group means, ANOVA p-values, IMD
#'   p-values, fold change estimates, fold change significance flags\cr \tab \cr
#'   Flags (signatures)  \tab Indicator of statistical significance for one
#'   (+/-1 for ANOVA, +/-2 for g-test) or neither (0) test \cr \tab \cr p-values
#'   \tab p-values for pairwise comparisons adjusted for within, e.g., peptide
#'   multiple comparisons only (not across peptides yet) \cr \tab \cr groupData
#'   \tab Mapping between sample ID and group \cr }
#'
#' @author Bryan Stanfill, Kelly Stratton
#'
#' @references Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical
#' analyses of peptide intensities and peptide occurrences improves
#' identification of significant peptides from MS-based proteomics data."
#' Journal of proteome research 9.11 (2010): 5748-5756.
#'
#' @examples
#' \dontrun{
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
#'
#' @export
#'
imd_anova <- function (omicsData,
                       comparisons = NULL,
                       test_method,
                       pval_adjust_a = 'none',
                       pval_adjust_g = 'none',
                       pval_thresh = 0.05,
                       covariates = NULL,
                       paired = FALSE,
                       equal_var = TRUE) {

  # Preliminaries --------------------------------------------------------------

  # check that omicsData is of the appropriate class
  if (!inherits(omicsData,
                c("proData", "pepData", "lipidData", "metabData", "nmrData")))
    stop ("omicsData is not an object of appropriate class")

  # Check for group_DF attribute #
  if(!("group_DF" %in% names(attributes(omicsData)))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attributes(omicsData)$group_DF
    # check for any groups of size 1 #

    # if any, filter out the corresponding sample(s) #
    # we don't want to filter them out, we just want to ignore them...

  }

  # Match provided 'test_method' to available options
  test_method <- try(match.arg(tolower(test_method),
                               c("combined", "gtest", "anova")),
                     silent = TRUE)

  if(!(test_method%in%c("combined","gtest","anova"))){
    #If test-method isn't valid, stop and tell them
    stop("Provided 'test_method' argument is invalid, please select 'anova','gtest' or 'combined'.")
  }

  # Check for log transform #
  if(!(attr(omicsData,"data_info")$data_scale%in%c("log2","log","log10"))&(test_method%in%c("combined","anova"))){
    stop("Data must be log transformed in order to implement ANOVA.")
  }

  # Check for anova filter - give warning if not present then let `imd_test` and `anova_test` do the actual filtering
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

  # Check if combined results was selected. If it is make sure the ANOVA and
  # G-test p-value adjustment arguments are both the same.
  if (test_method == "combined" && pval_adjust_a != pval_adjust_g) {

    # Dear pmartR user,
    #
    #   Please stop making my life difficult. I have enough problems as it is. I
    # do not need you complicating things further. I would really appreciate it
    # if you would read the instructions carefully, consider your options
    # thoroughly, and act wisely based on the outcome of the previous two steps.
    #
    #                                                 Sincerely,
    #                                                 pmartR programmer
    message(paste("The p-value adjustment method selected for ANOVA and G-test",
                  "are different. Check the input to pval_adjust_a and",
                  "pval_adjust_g.",
                  sep = " "))

  }

  # Statisticalness!!! ---------------------------------------------------------

  # Use imd_test() to test for independence of missing data (qualitative difference between groups)
  if(test_method=='anova'){
    #If they don't want the g-test done, save some time by removing comparisons and the pval_adjust arguments
    #Also make the gtest_pvalues NULL so nothing's returned.
    # NOTE: The code below doesn't do what the comments above say they want done.
    imd_results_full <- imd_test(omicsData,
                                 comparisons = NULL, # This actually performs all pairwise comparisons.
                                 pval_adjust = 'none',
                                 pval_thresh = pval_thresh)
  }else{
    imd_results_full <- imd_test(omicsData,
                                 comparisons = comparisons,
                                 pval_adjust = pval_adjust_g,
                                 pval_thresh = pval_thresh)
    gtest_pvalues <- imd_results_full$Pvalues
    colnames(gtest_pvalues) <- paste0("P_value_G_",colnames(gtest_pvalues))
    gtest_flags <- imd_results_full$Flags
    colnames(gtest_flags) <- paste0("Flag_G_",colnames(imd_results_full$Pvalues))
  }

  #Get cname and counts by group from imd_results_full
  imd_counts <- imd_results_full[[1]][,c(1,grep("^Count",colnames(imd_results_full[[1]])))]

  #Use anova_test() to test for a global equality of means, test requested comparisons and
  #compute foldchanges for those comaprisons (quantitative difference between groups)
  if(test_method=='gtest'){
    #If they only want the g-test, used anova_test to compute means, fold changes but
    #don't return flags or p-values
    anova_results_full <- anova_test(omicsData,
                                     comparisons = comparisons,
                                     pval_adjust = 'none',
                                     pval_thresh = pval_thresh,
                                     covariates = covariates,
                                     paired = paired,
                                     equal_var = equal_var)
  }else{
    anova_results_full <- anova_test(omicsData,
                                     comparisons = comparisons,
                                     pval_adjust = pval_adjust_a,
                                     pval_thresh = pval_thresh,
                                     covariates = covariates,
                                     paired = paired,
                                     equal_var = equal_var)
    anova_fold_flags <- anova_results_full$Flags
    colnames(anova_fold_flags) <- paste0("Flag_A_",colnames(anova_fold_flags))
    anova_pvalues <- anova_results_full$Fold_change_pvalues
    colnames(anova_pvalues) <- paste0("P_value_A_",colnames(anova_pvalues))
  }

  #Group means and fold changes should be returned even if they only ask for g-test
  anova_results <- anova_results_full[[1]][,c(1,grep("Mean",colnames(anova_results_full[[1]])))]
  anova_fold_change <- anova_results_full$Fold_changes
  colnames(anova_fold_change) <- paste0("Fold_change_",colnames(anova_fold_change))


  #Return appropriate results for each test_method

  if (test_method == 'anova') {

    all_anova <- cbind(anova_results,anova_pvalues,anova_fold_change,anova_fold_flags)
    Full_results <- dplyr::right_join(imd_counts,all_anova)

    if(nrow(Full_results)>max(nrow(imd_counts),nrow(all_anova))){
      stop("Combining g-test and ANOVA results failed.")
    }

    # Nab the comparisons actually performed. These will be used in the
    # statRes_output function when creating the number_significant attribute.
    actual_comparisons <- names(anova_results_full$Flags)

    # Create a variable for the ANOVA flags. This will be added as an attribute
    # of the statRes object and used in the bpquant function to create the
    # signatures variable.
    the_flag <- cbind(Full_results[[get_edata_cname(omicsData)]],
                      anova_results_full$Flags)
    colnames(the_flag)[1] <- get_edata_cname(omicsData)

  } else if (test_method == 'gtest') {

    all_gtest <- cbind(imd_counts,gtest_pvalues, gtest_flags)
    all_anova <- cbind(anova_results,anova_fold_change)
    Full_results <- dplyr::left_join(all_gtest, all_anova)

    if(nrow(Full_results)>max(nrow(all_gtest),nrow(all_anova))){
      stop("Combining g-test and ANOVA results failed.")
    }

    # Nab the comparisons actually performed. These will be used in the
    # statRes_output function when creating the number_significant attribute.
    actual_comparisons <- names(imd_results_full$Flags)

    # Set the_flag variable to NULL because we don't want this attribute to
    # exist when only the G-test was run.
    the_flag <- NULL

  } else if (test_method == "combined") {

    # Combine the G-test and ANOVA results!
    all_gtest <- cbind(imd_counts, gtest_pvalues, gtest_flags)
    all_anova <- cbind(anova_results, anova_pvalues, anova_fold_change,
                       anova_fold_flags)
    Full_results <- dplyr::full_join(all_gtest, all_anova)

    # Nab the comparisons actually performed. These will be used in the
    # statRes_output function when creating the number_significant attribute.
    actual_comparisons <- names(anova_results_full$Flags)

    # Combine ANOVA and G-test flags ---------------

    # Criteria for reporting p-values when combined test is selected:
    # (1)   If ANOVA flag but no G-test flag, report ANOVA flag
    # (2)   If G-test flag but no ANOVA flag, report G-test flag
    # (3)   If neither are present, report NA
    # (4)   If both are present but corresponding p-values are not significant,
    #       report ANOVA flag
    # (5)   If both are present but ANOVA is p-value is significant, report
    #       ANOVA flag
    # (6)   If both are present but G-test p-value is significant, report G-test
    #       flag

    # Start with the anova p-values (extracted from combined results) and
    # replace if missing [see (2)] or G-test is significant [see (6)]
    imd_pvals <- data.matrix(Full_results[grep("^P_value_A_",
                                               colnames(Full_results))])
    imd_flags <- data.matrix(Full_results[grep("^Flag_A_",
                                               colnames(Full_results))])

    # Extract G-test p-values and flags.
    g_pvals <- data.matrix(Full_results[grep("^P_value_G_",
                                             colnames(Full_results))])
    g_flags <- data.matrix(Full_results[grep("^Flag_G_",
                                             colnames(Full_results))])

    # Replance missing ANOVA p-values with g-test p-values
    if (any(is.na(imd_pvals))) {

      anova_NAs <- which(is.na(imd_pvals))
      imd_pvals[anova_NAs] <- g_pvals[anova_NAs]
      imd_flags[anova_NAs] <- g_flags[anova_NAs]

    }

    # Replace insignificant ANOVA p-values with significant g-test:
    # Insignificant ANOVA p-value indices.
    insig_anova <- which(imd_pvals > pval_thresh)

    if (any(insig_anova)) {

      # Nab significant G-test p-value indices.
      sig_gtest <- which(g_pvals <= pval_thresh)

      # Find overlap between the significant G-test p-values and the
      # insignificant ANOVA p-values.
      overlap <- sig_gtest[sig_gtest %in% insig_anova]

      if (any(overlap)) {

        # Replace ANOVA flags (with insignificant p-values) with G-test flags
        # (that have significant p-values).
        imd_flags[overlap] <- g_flags[overlap]

      }

    }

    # Remove "Flag_A_" from column names.
    colnames(imd_flags) <- gsub("^Flag_A_","",colnames(imd_flags))

    # Create a variable for the combined flags. This will be added as an
    # attribute of the statRes object and used in the bpquant function to create
    # the signatures variable.
    the_flag <- data.frame(Full_results[[get_edata_cname(omicsData)]],
                           imd_flags)
    colnames(the_flag)[1] <- get_edata_cname(omicsData)

  }

  # Reorder the columns of Full_results to: biomolecule ID, counts, means, fold
  # changes, p-values, and flags.
  Full_results <- Full_results %>%
    dplyr::relocate(
      dplyr::starts_with("Count_", vars = colnames(Full_results)),
      .after = !!rlang::sym(get_edata_cname(omicsData))
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("Mean_", vars = colnames(.data)),
      .after = dplyr::last_col()
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("Fold_change_", vars = colnames(.data)),
      .after = dplyr::last_col()
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("P_value_A", vars = colnames(.data)),
      .after = dplyr::last_col()
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("P_value_G", vars = colnames(.data)),
      .after = dplyr::last_col()
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("Flag_A", vars = colnames(.data)),
      .after = dplyr::last_col()
    ) %>%
    dplyr::relocate(
      dplyr::starts_with("Flag_G", vars = colnames(.data)),
      .after = dplyr::last_col()
    )

  final_out <- statRes_output(Full_results,
                              omicsData,
                              actual_comparisons,
                              test_method,
                              pval_adjust_a,
                              pval_adjust_g,
                              pval_thresh)

  attr(final_out, "bpFlags") <- the_flag
  attr(final_out, "cnames") = attr(omicsData, "cnames")
  attr(final_out, "data_class") = attr(omicsData, "class")

  return (final_out)

}
