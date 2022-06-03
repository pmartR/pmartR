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
#'   "bonferroni" and "holm". The default is "none" which corresponds to no
#'   p-value adjustment.
#' @param pval_thresh numeric p-value threshold, below or equal to which
#'   peptides are considered differentially expressed. Defaults to 0.05
#' @param use_parallel A logical value indicating whether or not to use a
#'   "doParallel" loop when running the G-Test with covariates. The default is
#'   TRUE.
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
                       equal_var = TRUE,
                       use_parallel = TRUE) {

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

  # Throw an error if there are no main effects and gtest or combined is
  # selected.
  if (
    "no_main_effect" %in% attr(groupData, "main_effects") &&
    test_method %in% c("gtest", "combined")
  ) {

    stop (
      paste(
        "If there are no main effects test_method cannot be gtest or combined.",
        sep = " "
      )
    )

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

  # Check if combined results was selected.
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

  # Farm boy, check the inputs. Farm boy, fix all the problems. Farm boy, hurry
  # up. Farm boy, why is the user an idiot? Farm boy, make everything uniform.
  # Farm boy, do all the tedious crap. AS YOU WISH!
  pval_adjust_a <- try(
    match.arg(tolower(pval_adjust_a),
              c("bonferroni", "tukey", "dunnett", "holm", "none")),
    silent = TRUE
  )
  if (class(pval_adjust_a) == 'try-error') {

    # I clearly told you what you could choose from in the documentation.
    # Obviously it is too much to ask for you to read the instructions!
    stop (
      paste("The available options for pval_adjust_a are: 'bonferroni',",
            "'holm', 'tukey', 'dunnett', and 'none'.",
            sep = " ")
    )

  }
  pval_adjust_g <- try(
    match.arg(tolower(pval_adjust_g),
              c("bonferroni", "holm", "none")),
    silent = TRUE
  )
  if (class(pval_adjust_g) == 'try-error') {

    # I cannot make this any easier for you. Please seek the help you
    # desperately need.
    stop (
      paste("The available options for pval_adjust_g are: 'bonferroni',",
            "'holm', and 'none'.",
            sep = " ")
    )

  }

  # Have a looksie at the covariates attribute.
  if (is.null(attr(attr(omicsData, "group_DF"), "covariates"))) {

    covariates <- NULL

  } else {

    # Grab the names of the covariates from the group_DF attribute. This is
    # necessary because covariates used to be an argument to the imd_anova
    # function. It is easier to create an object that contains the covariate
    # names than to change every instance where covariates show up. The code
    # would be much better if. ... I am going to stop there and not say what I
    # think about how things were done in the past.
    covariates <- names(attr(attr(omicsData, "group_DF"), "covariates"))[-1]

  }

  # Check if the data are paired. The paired object is used by multiple
  # functions within imd_anova.
  if (is.null(attr(attr(omicsData, "group_DF"), "pair_id"))) {
    paired <- FALSE
  } else {
    paired <- TRUE
  }

  # Check the input to the pval_thresh argument.
  if (!is.numeric(pval_thresh)) {

    # Why the heck would you input the p-value threshold as a character
    # vector?!?!?! In what universe does that make any sense?! I am still
    # baffled that we have run into this situation. I can feel my eye starting
    # to twitch.
    stop ("pval_thresh must be numeric.")

  }

  # Make sure pval_thresh is between 0 and 1.
  if (pval_thresh < 0 || pval_thresh > 1) {

    # If Lisa comes to me asking me to find out why the user's plots are
    # obviously not correct and it is because pval_thresh is not between 0 and 1
    # I will lose it. The Hulk will have nothing on me at that point.
    stop ("pval_thresh must be between 0 and 1.")

  }

  # Statisticalness!!! ---------------------------------------------------------

  # Use imd_test() to test for independence of missing data (qualitative difference between groups)
  if(test_method=='anova'){

    # If there are no main effects we just need the biomolecule ID column and
    # the counts.
    if ("no_main_effect" %in% attr(groupData, "main_effects")) {

      # Snag the index in e_data of the biomolecule ID column.
      edata_idx <- which(
        names(omicsData$e_data) == get_edata_cname(omicsData)
      )

      # Create a very watered down list similar to the output from imd_test.
      # This list will only contain the Results data frame with the biomolecule
      # ID and counts across all samples.
      imd_results_full <- list(
        Results = data.frame(
          ids = omicsData$e_data[, edata_idx],
          Count_paired_diff = rowSums(!is.na(omicsData$e_data[, -edata_idx]))
        )
      )

      # Rename the ID column with the appropriate name.
      names(imd_results_full$Results)[1] <- get_edata_cname(omicsData)

    } else {

      #If they don't want the g-test done, save some time by removing
      #comparisons and the pval_adjust arguments Also make the gtest_pvalues
      #NULL so nothing's returned.
      # NOTE: The code below doesn't do what the comments above say they want
      # done.
      imd_results_full <- imd_test(omicsData,
                                   groupData = groupData,
                                   comparisons = NULL, # This actually performs all pairwise comparisons.
                                   pval_adjust = 'none',
                                   pval_thresh = pval_thresh,
                                   covariates = NULL,
                                   paired = paired,
                                   use_parallel = use_parallel)
      # NOTE: covariates = NULL in the above call to imd_test because the
      # covariate portion of the code is the slowest part of the imd_test
      # function. In addition, the covariates portion of the results are not
      # used when test_method = anova.

    }

  }else{
    imd_results_full <- imd_test(omicsData,
                                 groupData = groupData,
                                 comparisons = comparisons,
                                 pval_adjust = pval_adjust_g,
                                 pval_thresh = pval_thresh,
                                 covariates = covariates,
                                 paired = paired,
                                 use_parallel = use_parallel)
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
                                     groupData = groupData,
                                     comparisons = comparisons,
                                     pval_adjust = 'none',
                                     pval_thresh = pval_thresh,
                                     covariates = NULL,
                                     paired = paired,
                                     equal_var = equal_var,
                                     use_parallel = use_parallel)
  }else{
    anova_results_full <- anova_test(omicsData,
                                     groupData = groupData,
                                     comparisons = comparisons,
                                     pval_adjust = pval_adjust_a,
                                     pval_thresh = pval_thresh,
                                     covariates = covariates,
                                     paired = paired,
                                     equal_var = equal_var,
                                     use_parallel = use_parallel)
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

    # Replace missing ANOVA p-values with g-test p-values
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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ANOVA FUNCTIONS --------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Tests for a quantiative difference between groups (aka factors, aka main
#' effects)
#'
#' This is the ANOVA part of the IMD-ANOVA test proposed in Webb-Robertson et
#' al. (2010).
#'
#' The order in which different scenarios are handeled: \enumerate{ \item If the
#' data are paired, then the pairing is accounted for first then each of the
#' next steps is carried out on the new variable that is the difference in the
#' paired individuals.<br> \item If covariates are provided, their effect is
#' removed before testing for group differences though mathematically covariates
#' and grouping effects are accounted for simultaneously \item ANOVA is executed
#' to assess the effect of each main effects, results in a vector of group means
#' for each biomolecule and variance estimate \item Group comparisons defined by
#' `comaprison` argument are implemented use parameter vector and variance
#' estimates in ANOVA step }
#'
#' @param omicsData A pmartR data object of any class
#' @param comparisons data.frame with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param pval_adjust A character string specifying the type of multiple
#'   comparisons adjustment to implement. The default, "none", corresponds to no
#'   adjustment. Valid options include: "bonferroni", "holm", "tukey", and
#'   "dunnett".
#' @param pval_thresh numeric p-value threshold, below or equal to which
#'   peptides are considered differentially expressed. Defaults to 0.05
#' @param covariates A character vector with no more than two variable names
#'   that will be used as covariates in the IMD-ANOVA analysis.
#' @param paired logical; should the data be paired or not? if TRUE then the
#'   `f_data` element of `omicsData` is checked for a "Pair" column, an error is
#'   returned if none is found
#' @param equal_var logical; should the variance across groups be assumed equal?
#' @param use_parallel A logical value indicating if the t test should be run in
#'   parallel.
#'
#'
#' @return  a list of `data.frame`s
#'   \tabular{ll}{
#'   Results  \tab Edata cname,
#'   Variance Estimate, ANOVA F-Statistic, ANOVA p-value, Group means\cr \tab
#'   \cr Fold_changes  \tab Estimated fold-changes for each comparison \cr \tab
#'   \cr Fold_changes_pvalues  \tab P-values corresponding to the fold-changes
#'   for each comparison \cr \tab \cr Fold_change_flags  \tab Indicator of
#'   statistical significance (0/+-2 to if adjusted p-value>=pval_thresh or
#'   p-value<pval_thresh) \cr
#'   }
#'
#' @author Bryan Stanfill
#'
#' @references Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical
#'   analyses of peptide intensities and peptide occurrences improves
#'   identification of significant peptides from MS-based proteomics data."
#'   Journal of proteome research 9.11 (2010): 5748-5756.
#'
anova_test <- function (omicsData, groupData, comparisons, pval_adjust,
                        pval_thresh, covariates, paired, equal_var,
                        use_parallel) {

  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if (k < 2 && !"no_main_effect" %in% attr(groupData, "main_effects")) {
    stop ("At least two groups are necessary to perform an ANOVA if the data are not paired.")
  }

  ###--------Check for anova filter-------------###
  if(is.null(attr(omicsData,"imdanova"))){
    warning("These data haven't been filtered, see `?imdanova_filter` for details.")
  }else{
    cnames <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
    filterrows <- which(cnames%in%attr(omicsData,"imdanova")$test_with_anova)
    if(length(filterrows)>0)
      omicsData$e_data <- omicsData$e_data[filterrows,]
  }

  #----Check to see if appropriate p-value adjustment technique is chosen----#
  #all pairwise -> Tukey
  #Case versus control -> Dunnett
  if(length(unique(comparisons$Control))==1 & tolower(pval_adjust)=="tukey" & length(comparisons$Control)>1){
    stop("Tukey adjustment for multiple comparisions should only be used when all pairwise comparisons are being made; try a 'Dunnett' adjustment.")
  }

  # if((length(unique(comparisons$Control))>1 | (is.null(comparisons) & length(unique(groupData$Group))>2)) & tolower(pval_adjust)=="dunnett"){
  #   stop("Dunnett adjustment for multiple comparisions should only be used when all case-vs-control comparisons are being made; try a 'Tukey' adjustment.")
  # }

  samp_cname <- attr(omicsData,"cnames")$fdata_cname

  # 10/27/2021 Not sure when/if the columns in e_data could be different from
  # the rows in group_DF. We are leaving the following code how it is to prevent
  # any unforeseen and unfortunate cases where groupData (the group_DF
  # attribute) has more columns than e_data (yikes!).

  # Remove rows in "groupData" that don't have corresponding columns in
  # 'omicsData$edata'
  # The sampleIDs in groupData (possibly more than is in the e_data)
  group_samp_cname <- groupData[,samp_cname]
  # Only keep the rows of groupData that have columns in e_data
  groupData <- groupData[group_samp_cname%in%colnames(omicsData$e_data),]

  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData[,samp_cname])]

  # Create a NULL object for covariate data. This will be used to determine if
  # this object has been modified in the paired section.
  cov_data <- NULL

  # Paired stuffs --------------------------------------------------------------

  # If paired==TRUE then use the pairing to create pair adjusted abundances.
  if (paired) {

    # Create a catchy name for the paired variable in f_data. The trick is to
    # make it readable but not too long. For example, using da instead of the
    # creates a name with one character less. These are the types of tricks you
    # learn when obtaining an advanced degree.
    da_pair_name <- attr(attr(omicsData, "group_DF"), "pair_id")

    # Make sure the paired column wasn't deleted between calling
    # group_designation and imd_anova. It is essential that you carry out all of
    # these tedious checks. The user is generally unstable when it comes to
    # following instructions. If you don't create all these safeguards they can
    # really do a number on your hard work. Then they come complaining to you as
    # if it was your fault they couldn't take the time to read and follow your
    # very clear (and amazing) instructions.
    if (!da_pair_name %in% names(omicsData$f_data)) {

      stop ("The pair ID variable is not present in f_data.")


    }

    # New steps for creating pair data:
    # 1. call function to take differences
    # 2. nab sample names according to type of data in pairing column
    # 3. update column (sample) names in data
    # 4. update groupData according to reduced number of samples in paired data
    # 5. update cov_data according to reduced number of samples in paired data

    # 1-3.
    # Take the difference between sample pairs. The columns of this data matrix
    # have been renamed according to the values in the pair ID column of f_data.
    data <- take_diff(omicsData)

    # 4.
    # The data has previously been combined into one value (the difference was
    # taken between each pair of samples). This means groupData also needs to be
    # collapsed. Otherwise the dimensions of the data matrix and groupData will
    # not match and the remaining pipeline will explode because the column names
    # of e_data will not match the sample names in f_data. There will be twice
    # as many rows in groupData as columns in the data matrix. Here we will
    # select the first row in groupData while grouping by the column with the
    # pairing information.
    groupData <- merge(
      groupData,
      omicsData$f_data[, c(samp_cname, da_pair_name)],
      sort = FALSE
    ) %>%
      dplyr::group_by(!!rlang::sym(da_pair_name)) %>%
      dplyr::slice(1) %>%
      data.frame()

    # Update the sample names in the groupData object according to the column
    # names of data (the paired difference data matrix).
    groupData[, samp_cname] <- colnames(data)

    # Add the original main_effects, pairs, and nonsingleton_groups attributes
    # to the updated groupData object. These attributes will be used in various
    # locations throughout the remainder of anova_test().
    attr(groupData, "main_effects") <- attr(
      attr(omicsData, "group_DF"), "main_effects"
    )
    attr(groupData, "pair_id") <- attr(attr(omicsData, "group_DF"), "pair_id")
    attr(groupData, "nonsingleton_groups") <- attr(
      attr(omicsData, "group_DF"), "nonsingleton_groups"
    )

    # Update the number of groups. This should be the same as before but we are
    # proceeding with an "overabundance of caution".
    k <- dplyr::n_distinct(groupData$Group)

    #If provided, covariate data needs to be updated too
    if (!is.null(covariates)) {

      # Add the covariate data to the sample ID and pair ID columns of f_data.
      # The sample ID column will always be the first element of the cols vector
      # and pair ID will always be the second element of that vector.
      cov_data <- dplyr::inner_join(
        x = omicsData$f_data[, c(samp_cname, da_pair_name)],
        y = attr(attr(omicsData, "group_DF"), "covariates")
      )

      # 5.
      # The data has previously been combined into one value (the difference was
      # taken between each pair). This means the covariate data also needs to be
      # collapsed. Otherwise the dimensions of the data matrix and the covariate
      # data will not match. The covariate data frame will have two values for
      # each difference. Here we will select the first row in cov_data while
      # grouping by the column with the pairing information.
      cov_data <- cov_data %>%
        dplyr::group_by(!!rlang::sym(da_pair_name)) %>%
        dplyr::slice(1) %>%
        # Remove the grouping structure. If we don't then the pair ID variable
        # will be added to the beginning of the data frame after removing it
        # with select.
        dplyr::ungroup() %>%
        # The second column will always be the pair ID column because we created
        # the cov_data object to be that way above. Remove this column because
        # it is not needed. We don't need it messing up our hard work in
        # unpredictable ways in other functions.
        dplyr::select(-2) %>%
        data.frame()

      # Update the sample names in the cov_data object according to the column
      # names of data (the paired difference data matrix).
      cov_data[, samp_cname] <- colnames(data)

    }

  }

  # Covariate stuffs -----------------------------------------------------------

  # If covariates are provided, remove their effect before computing group
  # means/variances.
  if(!is.null(covariates)){

    # Check if cov_data was updated/modified in the paired section. If it was
    # then it needs to remain in the format created in the paired section.
    if (is.null(cov_data)) {

      # Extract the covariates data frame from the group_DF attribute. This will
      # be used to combine with the main effects data to create the correct x
      # matrix (when removing the effect of covariates).
      cov_data <- attr(attr(omicsData, "group_DF"), "covariates")

    }

    # The -1 is hard coded because the sample ID name is ALWAYS the first column
    # in the covariates data frame.
    if (any(is.na(cov_data[, -1]))) {

      warning("Missing values were detected in the provided covariates thus covariates will be ignored")
      red_df <- matrix(rep(0,nrow(data)),ncol=1)

    } else {

      # Combine main effect and covariate data. We only want two columns from
      # the main effect data frame (groupData). These columns are the sample ID
      # column and the Group column. The Group column contains all the
      # information we need for creating the X matrix. We also need the sample
      # ID column from the covariate data frame along with the corresponding
      # columns for any covariates in the input.
      covariates <- merge(groupData[c(samp_cname, "Group")],
                          cov_data,
                          sort = FALSE)
      cov_samp_col <- which(colnames(covariates)==samp_cname)
      #source('~/pmartR/R/covariate_supp_funs.R') #Run to debug
      xmatrix <- build_x_mat(covariates[,-cov_samp_col]) #Build the appropriate X matrix based on the covariates data frame
      xmatrix <- reduce_xmatrix(xmatrix,length(unique(groupData$Group)))

      #Remove effect of covariates by projecting into the null space generated by the X matrix
      allproj <- project_to_null_cpp(data.matrix(data),xmatrix,length(unique(groupData$Group)))
      data <- allproj$data_no_x
      red_df <- allproj$lost_df
      #print(red_df)

    }

    # If covariates are adjusted for, then force equal variance assumption even
    # if there are only two levels to main effect
    equal_var <- TRUE

  } else {

    # Set the value of reduced degrees of freedom (the degrees of freedom used
    # up for covariates) to zero because covariates were not used.
    red_df <- matrix(rep(0,nrow(data)),ncol=1)

  }

  # paired t test stuffs -------------------------------------------------------

  if ("no_main_effect" %in% attr(groupData, "main_effects")) {

    return (
      paired_test(data = data,
                  bio_ids = omicsData$e_data[get_edata_cname(omicsData)],
                  cutoff = pval_thresh,
                  use_parallel = use_parallel)
    )

  }

  # ANOVA stuffs ---------------------------------------------------------------

  if (length(attr(get_group_DF(omicsData), "main_effects")) == 1) {
    ##---- One factor ANOVA ----##

    ##Translate the groups into numeric labels for anova_cpp() function
    gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))

    #Rcpp::sourceCpp('~/pmartR/src/anova_helper_funs.cpp') #Run if debugging code
    raw_results <- anova_cpp(data.matrix(data),gp,1-equal_var,red_df)
    group_names <- paste("Mean",as.character(unique(groupData$Group)),sep="_")

  }else{
    ##---- Two factor ANOVA ----##
    #source('~/pmartR/R/anova_helper_fun.R') #Run if debugging
    #source('~/pmartR/R/covariate_supp_funs.R')
    #Rcpp::sourceCpp('~/pmartR/src/two_factor_anova.cpp')
    raw_results <- run_twofactor_cpp(data=data.matrix(data),gpData=groupData,red_df)
    group_names <- paste("Mean",colnames(raw_results$group_means),sep="_")
  }

  #The C++ function returns a list so we need to make it into a data.frame
  results <- data.frame(raw_results$Sigma2, raw_results$group_means, raw_results$Fstats, raw_results$pvalue)

  #Rename the columns to match group names
  colnames(results) <- c("Variance",group_names,"F_Statistic","p_value")

  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  results <- cbind(omicsData$e_data[edatacname],results)

  ###-------Use group_comparison_anova() to compare the groups that were requested--------------##
  #source('~/pmartR/R/group_comparison.R') # Run if debugging
  #Rcpp::sourceCpp('~/pmartR/src/group_comparisons.cpp') #Run if debugging
  group_comp <- group_comparison_anova(groupData=groupData,comparisons=comparisons,anova_results_full=list(Results=results,Sizes=raw_results$group_sizes))

  #If there are only two levels, replace group_comp p-values with those taken from ANOVA function
  if(k==2){
    group_comp$p_values <- raw_results$pvalue
  }

  # Fold change stuffs ---------------------------------------------------------

  #Different way to compute fold change: diff if log scale, ratio otherwise
  if(attr(omicsData,"data_info")$data_scale%in%c("log2","log","log10")){
    fold_change <- group_comp$diff_mat
  }else{
    # CONSIDER MAKING GROUP_COMPARISON RETURN THE CMAT ALONG WITH EVERYTHING
    # ELSE INSTEAD OF RECREATING IT HERE
    Cmat_res <- create_c_matrix(groupData, comparisons)
    Cmat <- Cmat_res$cmat
    fold_change <- fold_change_ratio(raw_results$group_means, Cmat)
    fold_change <- data.frame(fold_change)
    colnames(fold_change) <- colnames(group_comp$diff_mat)
  }


  # p-value correction stuffs --------------------------------------------------

  p_values_adjust <- p_adjustment_anova(p_values = group_comp$p_values,
                                        diff_mean = group_comp$diff_mat,
                                        t_stats = group_comp$t_tests,
                                        sizes = raw_results$group_sizes,
                                        pval_adjust = pval_adjust)
  colnames(p_values_adjust) <- colnames(fold_change)


  # Flag stuffs ----------------------------------------------------------------

  # Create object to be returned
  sig_indicators <- data.matrix(p_values_adjust)

  sigs <- which(sig_indicators<pval_thresh)
  if(length(sigs)>0){
    sig_indicators[sigs] <- sign(data.matrix(group_comp$diff_mat)[sigs])
    sig_indicators[-sigs] <- 0
  }else{
    sig_indicators[which(sig_indicators>=pval_thresh)] <- 0
  }
  sig_indicators <- data.frame(sig_indicators)
  colnames(sig_indicators) <- colnames(fold_change)

  # Return stuffs --------------------------------------------------------------

  # (estimated variance, F-statistic, p-value and group means), estimated fold
  # changes, and adjusted p-values for the requested comparisons
  return (list(Results = results,
               Fold_changes = fold_change,
               Fold_change_pvalues = p_values_adjust,
               Flags = sig_indicators))

}

#Wrapper function for the two factor ANOVA function
run_twofactor_cpp <- function(data,gpData,red_df){
  #Create design matrix for reduced model, i.e., no interaction effect
  colnames(gpData)[-c(1,2)] <- c("Factor1","Factor2")
  gpData <- cbind(gpData,y=1:nrow(gpData))

  #Create design matrix for the full model, i.e., all first order and
  #interaction effects
  Xred <- unname(model.matrix(lm(y~Factor1+Factor2-1,data=gpData)))
  Xfull <- unname(model.matrix(lm(y~Factor1*Factor2-1,data=gpData)))

  #Run the two factor ANOVA model
  #Edit 2/16: add "group_ids" so C++ knows which groups belong to which columns,
  #helps with handling NaNs
  res <- two_factor_anova_cpp(
    data,
    Xfull,
    Xred,
    red_df,
    group_ids=as.numeric(factor(gpData$Group,levels=unique(gpData$Group)))
  )

  #Get the unique group levels to translate ANOVA parameters into group means
  red_gpData <- dplyr::distinct(gpData,Group,.keep_all=TRUE)
  red_gpData <- dplyr::arrange(red_gpData,y)

  means <- res$par_estimates[,1:length(red_gpData$Group)]
  colnames(means) <- red_gpData$Group

  res$group_means <- means
  return(res)
}

#' Group comparisons for the anova test
#'
#' Takes the results of anova_test() and returns group comparison p-values
#'
#' @param groupData data frame that assigns sample names to groups
#' @param comparisons dataframe that defiens the comparsions of interest
#' @param anova_results results of the pmartR::anova_test() function
#'
#' @return A data.frame containing the p-values from the group comparisons.
#'
#' @author Bryan Stanfill
#'
group_comparison_anova <- function(groupData,comparisons,anova_results_full){

  #The group means include the word "Group" in them so that is what will be passed
  #to the group_comparison(...) and fold_change(...) functions along with estimated variance
  anova_results <- anova_results_full$Results
  means <- data.matrix(anova_results[,grep("Mean",colnames(anova_results))])
  sigma2 <- data.matrix(anova_results$Variance)
  sizes <- data.matrix(anova_results_full$Sizes)
  #rm(anova_results_full)

  #If "comparisons" is null, then do all pairwise comparisons, else use the "comparisons" df to
  #create the appropriate C matrix

  ##---- I used to deal with one vs two factor differently, that code is at the end of this script -----##
  Cmat_res <- create_c_matrix(group_df = groupData, to_compare_df = comparisons)

  Cmat <- Cmat_res$cmat
  group_comp <- group_comparison_anova_cpp(means, sizes, sigma2, Cmat)

  group_comp$diff_mat <- data.frame(group_comp$diff_mat)
  colnames(group_comp$diff_mat) <- Cmat_res$names
  colnames(group_comp$diff_ses) <- Cmat_res$names
  colnames(group_comp$t_tests) <- paste(Cmat_res$names,"Test-Stat")
  group_comp$p_values <- data.frame(group_comp$p_values)
  colnames(group_comp$p_values) <- paste(Cmat_res$names,"p-value")
  return(group_comp)
}

# Runs a t test on the difference between paired data when there are no main
# effects.

# @author Evan A Martin
paired_test <- function (data, bio_ids, cutoff, use_parallel) {

  # Create the results data frame. With no main effects and paired data this
  # will consist of two columns. The first contains the biomolecule IDs and the
  # second is the mean of differences.
  results <- cbind(bio_ids, rowMeans(data, na.rm = TRUE))
  names(results)[2] <- "Mean_paired_diff"

  # Set up parallel backend.
  if (use_parallel) {
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores - 1)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  t_pval <- foreach::foreach(v = 1:nrow(data), .combine = "c") %dopar% {

    # Run the t test on each row of the difference data frame. Only the p-value
    # will be kept. This will be used to determine the flags and returned in the
    # final statRes object.
    t.test(data[v, ])$p.value

  }

  # Create a vector of zeros for the flag. If a p-value is significant the flag
  # will be changed to either -1 or 1 depending on the sign of the fold change.
  banner <- rep(0, nrow(data))

  # Determine which p-values fall below the threshold.
  siggy <- which(t_pval < cutoff)

  # Use the fold change (mean across columns) and the p-value to determine the
  # flag.
  banner[siggy] <- sign(results$Mean_paired_diff)[siggy]

  return (
    list(Results = results,
         Fold_changes = data.frame(paired_diff = results$Mean_paired_diff),
         Fold_change_pvalues = data.frame(paired_diff = t_pval),
         Flags = data.frame(paired_diff = banner))
  )

}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# IMD FUNCTIONS ----------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Tests for the independence of missing data across groups (aka factors, aka
#' main effects)
#'
#' Tests the null hypothesis that the number of missing observations is
#' independent of the groups.  A g-test is used to test this null hypothese
#' against the alternative that the groups and missing data are related.  This
#' is usually performed in conjuction with an ANOVA which tests if the mean
#' response (which varies with data type) is the same across groups; this
#' combination is called IMD_ANOVA.  It's probably a good idea to first filter
#' the data with `imd_anova_filter` to see if there is enough infomration to
#' even do this test.  See Webb-Robertson et al. (2010) for more.
#'
#' @param omicsData A pmartR data object of any class
#' @param comparisons `data.frame` with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param pval_adjust A character string specifying the type of multiple
#'   comparisons adjustment to implement. The default setting, "none", is to not
#'   apply an adjustment. Valid options include: "bonferroni" and "holm".
#' @param pval_thresh numeric p-value threshold, below or equal to which
#'   peptides are considered differentially expressed. Defaults to 0.05
#' @param covariates A character vector with no more than two variable names
#'   that will be used as covariates in the IMD-ANOVA analysis.
#' @param use_parallel A logical value indicating whether or not to use a
#'   "doParallel" loop when running the G-Test with covariates. The default is
#'   TRUE.
#'
#' @return  a list of `data.frame`s
#'   \tabular{ll}{
#'   Results  \tab e_data cname,
#'   Count of non-missing data for each group, Global G-test statistic and
#'   p-value\cr \tab \cr Gstats  \tab Value of the g statistics for each of the
#'   pairwise comparisons specified by the `comparisons` argument \cr \tab \cr
#'   Pvalues  \tab p-values for each of the pairwise comparisons specified by
#'   `comparisons` argument \cr \tab \cr Flags  \tab Indicator of statistical
#'   significance where the sign of the flag reflects the difference in the
#'   ratio of non-missing observations (0/+-2 to if adjusted
#'   p-value>=pval_thresh or p-value<pval_thresh) \cr
#'   }
#'
#' @author Bryan Stanfill
#'
#' @references Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical
#' analyses of peptide intensities and peptide occurrences improves
#' identification of significant peptides from MS-based proteomics data."
#' Journal of proteome research 9.11 (2010): 5748-5756.
#'
imd_test <- function (omicsData, groupData, comparisons, pval_adjust,
                      pval_thresh, covariates, paired, use_parallel) {

  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("This test cannot be performed with less than two groups.")
  }

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

  # Create a data matrix that only includes the samples in "groupData" and put
  # them in the order specified by "groupData"
  data <- omicsData$e_data[,as.character(groupData[,samp_cname])]

  # Find the number of times each peptide was absent for each of the
  # experimental groups, matrix of CAk values
  gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))

  #Number of samples associated with each group
  nK <- as.numeric(table(gp))
  label_map <- data.frame(Num_label=1:k,Group=unique(groupData$Group),nK=nK)

  # Global IMD test - any pattern in missing across all groups? ----------------

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
  edataIdx <- which(colnames(omicsData$e_data) == edatacname)
  colnames(observed) <- paste("Count",colnames(observed),sep='_')
  Global_results <- cbind(omicsData$e_data[edatacname],observed,Global_results)

  # Pairwise IMD test - any qualitative difference between groups? -------------

  # Create the c matrix. This matrix will be used to determine the number of
  # comparisons being performed, what the column names of the test statistic
  # and p-value data frames should be, and to determine the sign of the imd
  # flags.
  cmat <- create_c_matrix(group_df = groupData,
                          to_compare_df = comparisons)

  # Compute test statistic but summing over all nonzero counts (na.rm=TRUE
  # removes the zero counts)
  if(k==2){

    # Determine if covariates/paired data should be accounted for.
    if (is.null(covariates) || paired) {

      pairwise_stats <- data.frame(Global_Gk)
      pairwise_pvals <- data.frame(Global_results$p_value)

    } else {

      interim <- imd_cov(
        data = omicsData$e_data[, -edataIdx],
        groupData = groupData,
        fdata = omicsData$f_data,
        cmat = cmat,
        covariates = covariates,
        paired = paired,
        use_parallel = use_parallel
      )

      pairwise_stats <- interim$tstats
      pairwise_pvals <- interim$pvals

    }

    if(is.null(comparisons)){
      colnames(pairwise_stats)[1] <- colnames(pairwise_pvals)[1] <- paste(label_map$Group[1],label_map$Group[2],sep="_vs_")
    }else{
      colnames(pairwise_stats)[1] <- colnames(pairwise_pvals)[1] <- paste(comparisons$Test[1],comparisons$Control[1],sep="_vs_")
    }

    # Signs for group comparisons
    ratio_diff <- sign(
      (observed / (observed + absent)) %*% t(cmat$cmat)
    )
    colnames(ratio_diff) <- cmat$names

  }else{

    # Determine if covariates/paired data should be accounted for.
    if (is.null(covariates) && !paired) {

      pairwise_gk <- group_comparison_imd(groupData,comparisons,observed,absent)
      pairwise_stats <- data.frame(pairwise_gk$GStats)
      pairwise_pvals <- data.frame(pairwise_gk$pvalues)
      colnames(pairwise_pvals) <- colnames(pairwise_gk$pvalues)
      ratio_diff <- pairwise_gk$signs

    } else {

      interim <- imd_cov(
        data = omicsData$e_data[, -edataIdx],
        groupData = groupData,
        fdata = omicsData$f_data,
        cmat = cmat,
        covariates = covariates,
        paired = paired,
        use_parallel = use_parallel
      )

      pairwise_stats <- interim$tstats
      pairwise_pvals <- interim$pvals

      ratio_diff <- sign(
        (observed / (observed + absent)) %*% t(cmat$cmat)
      )
      colnames(ratio_diff) <- cmat$names

    }

  }

  # P-value adjustment ---------------------------------------------------------

  #Implement Bonferonni correction by multiplying by number of tests if requested, otherwise do nothing

  if(ncol(pairwise_pvals)==1)
    pval_adjust <- 'none'

  if(pval_adjust=="bonferroni"){
    adjusted_pvals <- pmin(data.matrix(pairwise_pvals*(ncol(pairwise_pvals))),1)
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
    sig_indicators[sigs] <- sign(ratio_diff[sigs])
    sig_indicators[-sigs] <- 0
  }else{
    sig_indicators[which(sig_indicators>=pval_thresh)] <- 0
  }

  sig_indicators <- data.frame(sig_indicators)
  #sig_indicators <- data.frame(adjusted_pvals[,1],sig_indicators)
  #colnames(sig_indicators)[1] <- colnames(adjusted_pvals)[1]

  return (list(Results=Global_results,
               Gstats=pairwise_stats,
               Pvalues=adjusted_pvals,
               Flags=sig_indicators))

}

#' Group comparisons for the g-test
#'
#' Takes the results of imd_test() and returns group comparison p-values
#'
#' @param groupData data frame that assigns sample names to groups
#' @param comparisons dataframe that defiens the comparsions of interest
#' @param observed matrix of number of observed counts
#' @param absent matrix of number of observed counts
#'
#' @return A data.frame containing the p-values from the group comparisons.
#'
#' @author Bryan Stanfill
#'
group_comparison_imd <- function(groupData,comparisons,observed,absent){

  #If "comparisons" is null, then do all pairwise comparisons, else use the "comparisons" df to
  #create the appropriate C matrix
  Cmat_res <- create_c_matrix(groupData, comparisons)
  Cmat <- abs(Cmat_res$cmat)

  #Number of samples per group
  nK <- observed[1,]+absent[1,]

  moMat <- observed%*%t(Cmat)
  maMat <- absent%*%t(Cmat)
  N <- moMat[1,]+maMat[1,]

  #CokMat <- CakMat <- EokMat <- EakMat <- NULL
  Gstats <- matrix(0,ncol=nrow(Cmat),nrow=nrow(moMat))
  for(i in 1:ncol(moMat)){
    nk_N_ratio <- nK/N[i]*Cmat[i,]

    EokMat <- matrix(moMat[,i],ncol=1)%*%matrix(nk_N_ratio,nrow=1)
    EakMat <- matrix(maMat[,i],ncol=1)%*%matrix(nk_N_ratio,nrow=1)
    CokMat <- observed*matrix(Cmat[i,],nrow=nrow(observed),ncol=ncol(Cmat),byrow=TRUE)
    CakMat <- absent*matrix(Cmat[i,],nrow=nrow(observed),ncol=ncol(Cmat),byrow=TRUE)
    GmatObs <- CokMat*log(CokMat/EokMat)
    GmatAbs <- CakMat*log(CakMat/EakMat)
    GmatObs[is.nan(GmatObs)] <- GmatAbs[is.nan(GmatAbs)] <- 0
    Gstats[,i] <- 2*rowSums(GmatObs+GmatAbs)
  }

  colnames(Gstats) <- Cmat_res$names
  pvals <- pchisq(Gstats, 1, lower.tail = FALSE)

  #Return signs for imd_test flags, i.e., diff in proportion of observed values
  signs <- sign((observed/(observed+absent))%*%t(Cmat_res$cmat))

  return(list(GStats=Gstats,pvalues=pvals,signs=signs))
}

# imd_cov carries out pairwise G-tests on all possible pairs of groups
# accounting for covariates and (when present) paired samples.
#
# @param data The e_data object without the biomolecule ID column.
#
# @param groupData The group_DF attribute.
#
# @param fdata The f_data object. This will only be used if the data are paired.
#
# @param cmat A list output by the create_c_matrix function.
#
# @param covariates A character vector specifying the covariates.
#
# @param paired A logical value indicating whether the data is paired.
#
# @param use_parallel A logical value indicating whether or not to use a
#   "doParallel" loop when running the G-Test with covariates. The default is
#   TRUE.
#
# @return A list of two data frames. The first data frame contains the test
#   statistics for each pairwise comparison across all biomolecules. The second
#   contains the corresponding p-values.
#
# @Author Evan A Martin
imd_cov <- function (data, groupData, fdata, cmat,
                     covariates, paired, use_parallel) {

  # Create an object that will correctly subset the anova(glm()) output given
  # the number of covariates present.
  if (paired && is.null(covariates)) {
    gem <- 3
  } else if (paired && length(covariates) == 1) {
    gem <- 4
  } else if (paired && length(covariates) == 2) {
    gem <- 5
  } else if (!paired && length(covariates) == 1) {
    gem <- 3
  } else {
    gem <- 4
  }

  # Lay hold on the appropriate indices of covariates attribute of the groupData
  # object using the names of the covariates in the input.
  covIdx <- which(colnames(attr(groupData, "covariates")) %in% covariates)

  # Snag the index of the pair ID variable in f_data.
  pairIdx <- which(names(fdata) == attr(groupData, "pair_id"))

  # Create a list for the test statistics and p-values from anova for each
  # pairwise test.
  list_anova <- vector(mode = "list", length = length(cmat$names))


  # Set up parallel backend.
  if (use_parallel) {
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores - 1)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  # Loop through each comparison and compute the test statistic and p-value for
  # each row in e_data.
  for (e in 1:length(cmat$names)) {

    # Determine which groups will be compared by extracting the group names from
    # cmat$names.
    groupNames <- stringr::str_replace(string = cmat$names[[e]],
                                       pattern = "_vs_",
                                       replacement = "|")

    # Nab the indices corresponding to the two groups that are being compared.
    # These indicies will be used to subset the columns of e_data and the rows
    # of both the main effect and covariate data frames.
    groupIdx <- stringr::str_which(groupData$Group, groupNames)

    nova <- foreach::foreach(v = 1:nrow(data)) %dopar% {

      # Account for paired data when covariates are NOT present.
      if (paired && is.null(covariates)) {

        # Run the G-test with the pairing information and no covariates.
        suppressWarnings(
          ninja <- anova(
            glm(
              as.numeric(!is.na(data[v, groupIdx])) ~
                fdata[groupIdx, pairIdx] +
                groupData$Group[groupIdx],
              family = binomial
            ),
            test = "Chisq"
          )
        )

        # Account for paired data and one covariate.
      } else if (paired && length(covariates) == 1) {

        # Run the G-test with the pairing information and the covariate.
        suppressWarnings(
          ninja <- anova(
            glm(
              as.numeric(!is.na(data[v, groupIdx])) ~
                fdata[groupIdx, pairIdx] +
                attr(groupData, "covariates")[groupIdx, covIdx] +
                groupData$Group[groupIdx],
              family = binomial
            ),
            test = "Chisq"
          )
        )

        # Account for paired data and two covariates.
      } else if (paired && length(covariates) == 2) {

        # Run the G-test with the pairing information and the covariates.
        suppressWarnings(
          ninja <- anova(
            glm(
              as.numeric(!is.na(data[v, groupIdx])) ~
                fdata[groupIdx, pairIdx] +
                attr(groupData, "covariates")[groupIdx, covIdx[[1]]] +
                attr(groupData, "covariates")[groupIdx, covIdx[[2]]] +
                groupData$Group[groupIdx],
              family = binomial
            ),
            test = "Chisq"
          )
        )

        # Account for one covariate when data are NOT paired.
      } else if (!paired && length(covariates) == 1) {

        # Run the G-test with one covariate.
        suppressWarnings(
          ninja <- anova(
            glm(
              as.numeric(!is.na(data[v, groupIdx])) ~
                attr(groupData, "covariates")[groupIdx, covIdx] +
                groupData$Group[groupIdx],
              family = binomial
            ),
            test = "Chisq"
          )
        )

        # Account for two covariates when data are NOT paired.
      } else {

        # Run the G-test with two covariates.
        suppressWarnings(
          ninja <- anova(
            glm(
              as.numeric(!is.na(data[v, groupIdx])) ~
                attr(groupData, "covariates")[groupIdx, covIdx[[1]]] +
                attr(groupData, "covariates")[groupIdx, covIdx[[2]]] +
                groupData$Group[groupIdx],
              family = binomial
            ),
            test = "Chisq"
          )
        )

      }

      # Snag the test statistic and p-value from the anova(glm()) output.
      dragon <- data.frame(daStats = ninja$Deviance[[gem]],
                           daPvals = ninja$`Pr(>Chi)`[[gem]])

      dragon

    }

    # Combine each row into one data frame.
    list_anova[[e]] <- data.table::rbindlist(nova)

  }

  # Create separate data frames for the test statistics and the p-values.
  devi <- cbindList(lapply(list_anova, `[[`, "daStats"))
  pchi <- cbindList(lapply(list_anova, `[[`, "daPvals"))

  # Name the columns of the test statistic and p-value data frames with the
  # names of the corresponding groups being compared.
  colnames(devi) <- cmat$names
  colnames(pchi) <- cmat$names

  return (list(tstats = devi,
               pvals = pchi))

}

# cbindList takes a list and combines each element of the list in a data frame
# column-wise.
#
# @param aList A list whose elements are vectors all with the same length.
#
# @return A data frame with each element of the list as a column in the data
#   frame.
#
# @Author Evan A Martin
cbindList <- function (aList) {

  # Create a matrix that will have one column for each element of aList.
  aMatrix <- matrix(0,
                    nrow = length(aList[[1]]),
                    ncol = length(aList))

  # Loop through each element in aList and add it as a column to a data frame.
  for (e in 1:length(aList)) {

    aMatrix[, e] <- aList[[e]]

  }

  # Convert to a data frame.
  aDF <- data.frame(aMatrix)

  return (aDF)

}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMPARISON FUNCTIONS ---------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

create_c_matrix <- function(group_df,to_compare_df=NULL){
  #This function will create the C-matrix that defines the comparisons that are requested
  #and returns a list of names that defines which two groups are being compared for each
  #row of the C matrix

  #group_df - defines which columns (`sample_name`) belong to which groups (`Group`)
  groups_two_fac <- groups <- as.character(unique(group_df$Group))
  ngroups <- length(groups)
  compi <- NULL

  if(is.null(to_compare_df)){

    #Create the Cmatrix that corresponds to all pairwise comparisons
    #With two factors, all pairwise includes the higher order terms
    #THIS COULD PROBABLY BE MADE MUCH FASTER IN C++
    n_comparisons <- choose(ngroups,2)
    Cmat <- matrix(0,n_comparisons,ngroups)
    row <- 1
    for(i in 1:(ngroups-1)){
      for(j in (i+1):ngroups){
        Cmat[row,i] <- 1
        Cmat[row,j] <- -1
        row <- row+1
        compi <- c(compi,paste0(groups[i],"_vs_",groups[j]))
      }
    }
  }else{
    n_comparisons <- nrow(to_compare_df)
    Cmat <- matrix(0,n_comparisons,ngroups)
    all_comp_levels <- as.character(unique(unlist(to_compare_df)))

    if(ncol(group_df)>2){#Two factor case, allow for main effect tests so grow groups to include those
      group_df[,-1] <- apply(group_df[,-1],2,as.factor)
      groups_two_fac <- as.character(unique(unlist(group_df[,-1])))
    }

    #Check that all of the groups in "to_compare_df" are present in "group_df"
    if(!all(all_comp_levels%in%groups_two_fac) & !all(all_comp_levels%in%groups)){
      stop("Some groups listed in the 'comparisons' argument do not appear in the 'groupData' argument.")
    }

    for(i in 1:n_comparisons){
      control_i <- grep(to_compare_df$Control[i],groups)
      test_i <- grep(to_compare_df$Test[i],groups)
      Cmat[i,control_i] <- (-1/length(control_i))
      Cmat[i,test_i] <- 1/length(test_i)
      compi <- c(compi,paste0(to_compare_df$Test[i],"_vs_",to_compare_df$Control[i]))
    }

  }

  return(list(cmat=Cmat, names=compi))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# P-VALUE FUNCTIONS ------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Adjust p-values for multiple comparisons
#'
#' Depending upon the \code{pval_adjust} method selected, the supplied p_values are compared against an adjusted \code{pval_thresh} value or the provided
#' means are used to compute new statistics, p-values are computed and compared against the provided \code{pval_thresh}.  A \code{data.frame} that indicates which
#' of the tests are significant, 1 if significant or 0 if insignificant.  If \code{means} is also provided and the p-value is signficant then the direction
#' of the change is indicated by the sign on 1, i.e., means<0 and p_value<pval_thresh will return -1, similarly for means>0.
#'
#' @param p_values A matrix (or \code{data.frame}) of p-values to be adjusted.
#' @param diff_mean A matrix (or \code{data.frame}) of groups means that are to be compared
#' @param t_stats A matrix (or \code{data.frame}) of t-test statistics resulting from from standard procedures
#' @param sizes A matrix (or \code{data.frame}) of group sizes
#' @param pval_adjust character vector specifying the type of multiple comparisons adjustment to implement. A NULL value corresponds to no adjustment. Valid options include: holm, bonferroni, dunnett, tukey or none.
#'
#' @return a data frame with the following columns: group means, global G-test statistic and corresponding p-value
#'
#' @author Bryan Stanfill
#'
p_adjustment_anova <- function (p_values, diff_mean, t_stats,
                                sizes, pval_adjust) {

  if(pval_adjust=="tukey" | pval_adjust=="dunnett"){

    #Tukey-Kramer statistics are t-statistic/sqrt(2)
    if(is.null(t_stats) | is.null(sizes)){
      stop("The standard t-tests and group sizes need to be supplied in order to apply Tukey-Kramer adjustment.")
    }

    n_compares <- ncol(t_stats)

    if(n_compares==1){
      #Tukey adjustment is not done if only one comparison is provided
      pval_adjust <- 'none'
    }else{

      if(pval_adjust=="tukey"){
        tukey_stats <- data.matrix(t_stats*sqrt(2))
        #Tukey only needs sizes for the rows, not individual group sizes
        if(is.data.frame(sizes)){
          sizes <- rowSums(data.matrix(sizes))
        }else if(is.matrix(sizes)){
          sizes <- rowSums(sizes)
        }

        #Rcpp::sourceCpp('src/tukey_helper.cpp') #Use for debugging
        adjusted_pvals <- ptukey_speed(tukey_stats,sizes)

        # Need to find all rows where only one test was performed and replace
        # the adjusted p-value with the original p-value.
        just_one <- which(rowSums(!is.na(t_stats)) == 1)

        # Find the non-missing values from the original p-value data frame.
        # These are the p-values that need to replace the incorrectly adjusted
        # p-values (because only one test was performed). The numbers returned
        # from this code will be the indices of the data frame as a vector (not
        # a row/column index). R's convention is to number down the rows then
        # across the columns. For example,
        # matrix(1:25, nrow = 5, byrow = TRUE)[c(3, 5)]
        # [1] 11 21
        not_na <- which(!is.na(p_values[just_one, ]))

        # Replace all unjustly adjusted p-values with the original p-value.
        # These are the rows where only one p-value was calculated.
        adjusted_pvals[just_one, ][not_na] <- data.matrix(
          p_values[just_one, ]
        )[not_na]

      }else{#Dunnett adjustment - Needs to be sped up
        adjusted_pvals <- matrix(NA,nrow(t_stats), ncol(t_stats))
        #colnames(adjusted_pvals) <- colnames(t_stats)

        for(i in 1:nrow(adjusted_pvals)){
          k <- length(which(!is.na(p_values[i,])))
          if(k>0){
            dfi <- sum(sizes[i,])-k
            rm_nas <- which(is.na(p_values[i,]))
            #Modified version of an example from "Additional multcomp Examples" viggnette of the multcomp package
            #Until we add covariates the correlation matrix is assumed to be the identity
            if(length(rm_nas)>0){
              adjusted_pvals[i,-rm_nas] <- as.numeric(sapply(abs(t_stats[i,-rm_nas]), function(x,k,df) 1 - mvtnorm::pmvt(-rep(x, k), rep(x, k), df = df),k=k,df=dfi))
            }else{
              adjusted_pvals[i,] <- as.numeric(sapply(abs(t_stats[i,]), function(x,k,df) 1 - mvtnorm::pmvt(-rep(x, k), rep(x, k), df = df),k=k,df=dfi))
            }
          }
        }
      }

      colnames(adjusted_pvals) <- colnames(t_stats)
      colnames(adjusted_pvals) <- gsub("Test-Stat","p-value",colnames(adjusted_pvals)) #This is a band-aid right now, may need to make more general later
    }
  }

  #Don't make this an "else if" because pval_adjust can be changed to 'none' if n_compares==1
  if(pval_adjust%in%c('none',"bonferroni","holm")){

    #p_values needs to be supplied to apply bonferroni adjustment, if it's NULL tell them that's a problem
    if(is.null(p_values)){
      stop("The `p_values` argument must be supplied to perform the selected `pval_adjust` method")
    }

    if(is.null(dim(p_values)) || ncol(p_values)==1){ # leaving this as "||" because only evaluating the 1st argument is fine in this case, as it returns TRUE when p_values is a vector (and thus ncol(p_values) does not compute and using a "|" gives an error)
      pval_adjust='none'
    }

    if(pval_adjust !='holm'){
      #For bonferroni adjustment, multiply p_values by number of comparisons (columns in p_values) else do no adjustment
      multiplier <- ifelse(pval_adjust=="bonferroni",ncol(p_values),1)
      adjusted_pvals <- pmin(data.matrix(multiplier*p_values), 1)
    }else{
      #Rcpp::sourceCpp('src/holm_adjustment.cpp') #For debugging
      #source('~/pmartR/R/support_p_adjustment.R') #For debugging
      adjusted_pvals <- t(apply(p_values,1,ranked_holm_cpp))

      #NaN p-values should stay NaN
      nan_pvals <- lapply(p_values,is.nan)
      for(i in 1:ncol(adjusted_pvals)){
        adjusted_pvals[,i][nan_pvals[[i]]] <- NaN
      }
    }

  }

  #Return the adjusted p-values

  return(pvalues=data.frame(adjusted_pvals))
}

# I haven't found a function that quickly gets the rank of the original p-values
# so do it in R for now, still much faster than a pure R version
ranked_holm_cpp <- function(ps){
  return(holm_cpp(ps)[rank(ps)])
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COVARIATE FUNCTIONS ----------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Function to create an X matrix based on a covariate data frame
build_x_mat <- function(cov_df){

  if(is.null(ncol(cov_df))){
    cov_df <- matrix(cov_df,ncol=1)
  }

  #If all covariates are numeric, simply return the same matrix back
  if(is.numeric(cov_df)){
    return(data.matrix(cov_df))
  }

  #If the covariates are a mix of numeric, factors and characters, return matrix
  #of group identifiers
  Xmatrix <- NULL
  for(i in 1:ncol(cov_df)){

    if(is.numeric(cov_df[,i])){
      #If column i is numeric, append it to X
      Xmatrix <- cbind(Xmatrix,cov_df[,i])
    }else{
      coli_levels <- unique(cov_df[,i])
      n_levels <- length(coli_levels)
      if(n_levels!=length(cov_df[,i])){
        Xcoli <- matrix(0,nrow(cov_df),n_levels)

        for(j in 1:n_levels){
          Xcoli[cov_df[,i]==coli_levels[j],j] <- 1
        }
        Xmatrix <- cbind(Xmatrix,Xcoli)
      }
    }
  }
  return(Xmatrix)
}

#---Function to remove linearly dependent columns in x corresponding to group
#levels----####
reduce_xmatrix <- function(x,ngroups){
  p <- ncol(x)
  orig_rank <- qr(x)$rank
  if(orig_rank<p){
    #Remove constant columns
    col_sds <- apply(x,2,sd)
    if(any(col_sds==0)){
      x <- x[,-which(col_sds==0)]
    }else{
      stdx <- apply(x,2,function(mat) (mat-mean(mat))/sd(mat))
      dmat <- as.matrix(dist(t(stdx)))
      dmat <- dmat[(1:ngroups),-(1:ngroups)]
      if(any(dmat==0)){
        ind_mat <- matrix(1:ngroups,ngroups,(p-ngroups))
        cols_to_remove <- ind_mat[which(dmat==0)]
        x[,cols_to_remove] <- 0
      }
    }
  }
  return(x)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PAIRED FUNCTIONS -------------------------------------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Function to translate paired data into differences
build_pair_mat <- function(pair_df){
  #pair_df - data.frame with two columns: sampleIDs and PairIDs (in that order)

  #Creat matrix of 1, -1 that will turn single observations into paired differences
  coli_levels <- unique(pair_df[,2])
  n_levels <- length(coli_levels)
  Xmatrix <- matrix(0,nrow(pair_df),n_levels)

  for(j in 1:n_levels){
    if(length(which(pair_df[,2]==coli_levels[j]))>2){
      stop(paste("Only two samples can be associated with the same 'pair'.  More than two samples are associated with pair",coli_levels[j],"."))
    }
    Xmatrix[pair_df[,2]==coli_levels[j],j] <- c(1,-1)
  }

  return(Xmatrix)
}

#' Compute pairwise differences
#'
#' Computes the differences for paired data according to the information in the
#' pairing column of f_data. This variable name is also an attribute of the
#' group_DF attribute.
#'
#' @param omicsData Any one of the omicsData objects (pepData, metabData, ...).
#'
#' @return A data.frame containing the differences between paired samples.
#'
#' @author Evan A Martin
#'
take_diff <- function (omicsData) {

  # Apprehend the index in f_data for the sample ID and pair ID columns.
  samp_idx <- which(names(omicsData$f_data) == get_fdata_cname(omicsData))
  pair_idx <- which(
    names(omicsData$f_data) == attr(attr(omicsData, "group_DF"), "pair_id")
  )

  # Create a matrix with -1 and 1 in the corresponding rows to calculate the
  # difference between the paried samples. This matrix is created according to
  # the column in f_data containing the paired information.
  pid_matrix <- build_pair_mat(
    pair_df = omicsData$f_data[, c(samp_idx, pair_idx)]
  )

  # Find the biomolecule ID index in e_data.
  bio_idx <- which(
    names(omicsData$e_data) == get_edata_cname(omicsData)
  )

  # Compute the differences between paired samples according to the pid matrix.
  diff_data <- fold_change_diff_na_okay(
    data = data.matrix(omicsData$e_data[, -bio_idx]),
    C = t(pid_matrix)
  )

  # Extract the first occurance of each pair variable in f_data. These will be
  # used later to rename the columns of diff_data.
  moniker <- omicsData$f_data %>%
    dplyr::group_by(dplyr::across(pair_idx)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(pair_idx)

  # Check if the data from the pair ID column is numeric. If it is the word
  # "Pair" will be added before each number. That way data.frame() won't
  # complain because the column names of the difference data will not be
  # numbers. Everyone wins!!!
  if (is.numeric(moniker)) {

    moniker <- paste("Pair", moniker,
                     sep = "_")

  }

  # Rename the columns of e_data to reflect that the paired differences were
  # calculated.
  colnames(diff_data) <- moniker

  return (diff_data)

}
