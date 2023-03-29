#' Differential Expression for seqData
#' 
#' Performs statistical analysis for differential expression of seqData objects,
#' using methods from one of: edgeR, DESeq2, or limma-voom
#'
#' @param omicsData object of type 'seqData' created by
#'   \code{\link{as.seqData}}
#' @param method character string of length one specifying which wrapper to
#'   use. Can be 'edgeR', 'DESeq2', or 'voom'
#' @param p_adjust character string for p-value correction method, refer to
#'   ?p.adjust() for valid options. Defaults to "BH" (Benjamini & Hochberg).
#' @param comparisons data frame with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param p_cutoff numeric value between 0 and 1 for setting p-value
#'   significance threshold
#' @param ... additional arguments passed to methods functions. Note, formatting
#'   option changes will interfere with wrapping functionality.
#'   
#'
#' @return object of class statRes 
#'
#' @details Runs default differential expression workflows.
#' 
#' Flags (signatures) -	Indicator of statistical significance.
#' Zeroes indicate no significance, while +/- 1 indicates direction of
#' significance.
#'
#' Method "edgeR" - Runs default edgeR workflow with empirical Bayes
#' quasi-likelihood F-tests. Additional arguments can be passed for use in the
#' function. Refer to calcNormFactors() and glmQLFit() in edgeR package.
#'
#' Method "DESeq2" - Runs default DESeq workflow. Defaults to Wald test, no
#' independent filtering, and running in parallel. Additional arguments can be
#' passed for use in the function. Refer to DESeq() and results() in DESeq2
#' package. Requires 'survival' package to run.
#'
#' Method "voom" - Runs default limma-voom workflow using empirical Bayes
#' moderated t-statistics. Additional arguments can be passed for use in the
#' function. Refer to calcNormFactors() in edgeR package.
#' 
#' @references 
#'  Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package 
#'  for differential expression analysis of digital gene expression data.” 
#'  Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.
#'  
#'  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and 
#'  dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#'  
#'  Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, 
#'  G.K. (2015). limma powers differential expression analyses for 
#'  RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#' 
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' myseqData <- group_designation(omicsData = rnaseq_object, main_effects = "Virus")
#' edger_results <- diffexp_seq(omicsData = myseqData, method = "edgeR")
#' deseq_results <- diffexp_seq(omicsData = myseqData, method = "DESeq2")
#' voom_results <- diffexp_seq(omicsData = myseqData, method = "voom")
#' }
#' @export
#' 
diffexp_seq <- function(omicsData, method = "edgeR", p_adjust = "BH", 
                        comparisons = NULL, p_cutoff = 0.05, ...){
  
  ## check inputs
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("seqData"))) {
    # Throw an error that the input for omicsData is not the appropriate class.
    stop("omicsData must be of class 'seqData'")
  }
  
  # check p_adjust #
  if (length(p_adjust) != 1 || !(p_adjust %in% p.adjust.methods)) {
    
    stop(paste("p_adjust must a single character string of length 1 in one of the following: ", p.adjust.methods))
    
  }
  
  # check p_cutoff #
  if (length(p_cutoff) != 1 || !is.numeric(p_cutoff) || p_cutoff > 1 || p_cutoff < 0) {
    
    stop("p_cutoff must numeric of length 1 between 0 and 1.")
    
  }
  
  # check method #
  if(length(method) != 1 || !(method %in% c('edgeR', 'DESeq2','voom'))){
    stop("method must a single character string of length 1 in 'edgeR', 'DESeq2', or 'voom'")
    
  }
  
  # must have group_desig #
  if(is.null(get_group_DF(omicsData))){
    stop("Running group_designation is required before statistical analysis")
    
  }
  
  
  if(method == 'edgeR'){
    
    edgeR_wrapper(omicsData = omicsData, p_adjust = p_adjust, 
                  comparisons = comparisons, p_cutoff = p_cutoff, ...)
    
  } else if(method == "DESeq2"){
    
    DESeq2_wrapper(omicsData = omicsData, p_adjust = p_adjust, 
                   comparisons = comparisons, p_cutoff = p_cutoff, ...)
    
  } else if(method == "voom"){
    
    voom_wrapper(omicsData = omicsData, p_adjust = p_adjust, 
                 comparisons = comparisons, p_cutoff = p_cutoff, ...)
    
  }
}



#' Wrapper for DESeq2 workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param test either "Wald" or "LRT", which will then use either Wald 
#' significance tests, or the likelihood ratio test on the difference in 
#' deviance between a full and reduced model formula 
#' @param method a character string of length one specifying which wrapper to use. 
#' Can be 'edgeR', 'DESeq2', or 'voom' 
#' @param p_adjust Character string for p-value correction method, refer to 
#' ?p.adjust() for valid options. Defaults to "BH" (Benjamini & Hochberg)
#' @param comparisons `data.frame` with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param p_cutoff Numeric value between 0 and 1 for setting p-value
#'   significance threshold
#' @param ... additional arguments passed to methods functions. Note, formatting
#'   option changes will interfere with wrapping functionality.
#' @param ... additional arguments passed to function
#'
#' @details Runs default DESeq workflow. Defaults to Wald test, no independent
#'   filtering, and running in parallel. Additional arguments can be passed for
#'   use in the function, refer to DESeq() and results() in DESeq2 package.
#'   Requires 'survival' package to run.
#'
#'   Flags (signatures) -	Indicator of statistical significance for computed
#'   test. Zeros indicate no significance, while +/- 1 indicates direction of
#'   significance.
#'
#' @references 
#'  
#'  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and 
#'  dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#'  
#'
#' @return statRes object
#' 
#' 
DESeq2_wrapper <- function(
    omicsData, test = "Wald",  p_adjust = "BH", comparisons = NULL,
    p_cutoff = 0.05, ...
){
  
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("package 'survival' required for DESeq2 processing")
  }
  
  l <- list(...)
  
  DESeq_args <- l[names(l) %in% names(formals(DESeq2::DESeq))]
  results_args <- l[names(l) %in% names(formals(DESeq2::results))]
  
  ## Normal p_value adjustments, test can either be "Wald" or liklihood ratio "LRT"
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  group_res <- get_group_formula(omicsData)
  grouping_info <- group_res[[1]]
  grouping_formula <- group_res[[2]]
  
  all_cols <- stringr::str_trim(
    stringr::str_split(
      stringr::str_remove(
        grouping_formula, "~0 \\+"), " \\+ ")[[1]])
  grouping_contrasts <- grouping_info[all_cols]
  grouping_contrasts <- grouping_contrasts[!apply(grouping_contrasts, 2, is.numeric)]
  contrast_levels <- as.character(unique(unlist(grouping_contrasts)))
  cob_list <- combn(contrast_levels, 2)
  all_contrasts <- apply(cob_list, 2, paste, collapse = "-")
  
  ## Make sure comparisons are doing the thing ##
  if(is.null(comparisons) && 
     !identical(attr(grouping_info, "main_effects"), "no_main_effect")
     ){
    
    comparisons <- c()
    comparisons <- unique(grouping_info[["Group"]])
    cob_list_group <- combn(comparisons, 2)
    all_contrasts_group <- apply(cob_list_group, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_group)
    
  } else if (is.null(comparisons)){
    
    comparisons <- c()
    comparisons <- c(unique(as.character(grouping_info[[attr(grouping_info, "pair_group")]])))
    cob_list_pair <- combn(comparisons, 2)
    all_contrasts_pair <- apply(cob_list_pair, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_pair)
    
  } else {
    
    if(ncol(cob_list) == 1){
      all_contrasts <- unique(c(all_contrasts, paste(cob_list[nrow(cob_list):1,], collapse = "-")))
      
    } else {
      all_contrasts <- unique(c(all_contrasts,
                                apply(cob_list[nrow(cob_list):1,], 2, paste, collapse = "-")))
    }
    
    cob_list <- cbind(cob_list, cob_list[nrow(cob_list):1,])
    all_contrasts_comp <- paste(comparisons$Control, comparisons$Test, sep = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_comp)
    if(length(interesting_comparisons) == 0) stop("Invalid comparisons given")
    
  }
  
  ## DESeq workflow
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  grouping_info <- grouping_info[colnames(grouping_info) != fdata_cname]
  
    edata_deseq <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
      e_data_counts, 
      colData = grouping_info,
      design = as.formula(stringr::str_remove(
        grouping_formula, "0 \\+")) ## The zero does wonky things with covariates
    ))
    
    dds <- DESeq2::estimateSizeFactors(edata_deseq) # run this or DESeq() first
    norm_factors_apply <- DESeq2::counts(dds, normalized=TRUE)
    
    group_means <- purrr::map_dfc(1:length(unique(grouping_info$Group)), function(group_n){
      group <- unique(grouping_info$Group)[group_n]
      group_cols <- as.data.frame(norm_factors_apply)[grouping_info$Group == group]
      means_group <- apply(group_cols, 1, mean, na.rm = T)
      df <- data.frame(means_group)
      colnames(df) <- group
      df
    })
    
    group_means <- cbind(omicsData$e_data[edata_cname], group_means)
    
  ## All possible arguments ugh
    
    list_defaults <- list(
      object = edata_deseq,
      test = test, ## reasonable to allow change in pmartR
      fitType = "parametric", #fitting of dispersions to the mean intensity
      quiet = TRUE,
      minReplicatesForReplace = Inf, ## cook's distance used to flag outliers, require at least n replicates to replace flagged outliers-- defined as .99 quantile of the F(p, m - p) distribution, where p is the number of parameters and m is the number of samples. Replacement is values predicted by the trimmed mean over all samples (and adjusted by size factor or normalization factor). Inf disables replacement
      modelMatrixType = "standard", ## If we need anything fancier than what is defined by "Group" we should consider this argument for the glm use case -- betapriors req for expanded
      parallel = TRUE,
      
      ## Waldy
      betaPrior = F, 
      #whether or not to put a zero-mean normal prior on the non-intercept coefficients See nbinomWaldTest for description of the calculation of the beta prior.In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
      useT = F ## only for waldtest, uses a t-distribution instead of a normal distribution)
    )
    
    run_deseq <- c(DESeq_args, list_defaults)
    run_deseq <- run_deseq[!duplicated(names(run_deseq))]
    
  run_stats_deseq <- do.call(DESeq2::DESeq, args = run_deseq)
 
  ## Args
  
  if(stringr::str_detect(grouping_formula, "Group")){
    contrast_front <- "Group"
  } else {
    contrast_front <- attr(get_group_DF(omicsData), "pair_group")
  }
  
  list_defaults <- list(
    object = run_stats_deseq, 
    lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    cooksCutoff = FALSE, 
    ## Fiters by distance, default it filters .99 quantile of the F(p, m-p) 
    ## distribution, where p is the number of coefficients being fitted 
    ## and m is the number of samples. Excludes groups w/ only 2 samples
    independentFiltering = FALSE,
    pAdjustMethod = p_adjust,
    tidy = TRUE,
    parallel = T
  )
  
  run_results <- c(results_args, list_defaults)
  run_results <- run_results[!duplicated(names(run_results))]
  
  all_res <- purrr::map(interesting_comparisons, function(combo_n){
    
    combo <- cob_list[,combo_n]
    
    add_list <- list(contrast = c(contrast_front, combo[1], combo[2]))
    run_results <- c(run_results, add_list)
    run_results <- run_results[!duplicated(names(run_results))]
    

    ## Run tests
    res <- do.call(DESeq2::results, args = run_results)
    
    res[["row"]] <- NULL
    # Flag stuffs ----------------------------------------------------------------
    sigs <- which(res[["padj"]] < p_cutoff)
    res[["Flag"]] <- 0
    if(length(sigs)>0){
      res[["Flag"]][sigs] <- sign(res$log2FoldChange[sigs])
    }
    
    res[["row"]] <- NULL
    colnames(res) <- paste0(
      colnames(res), 
      paste0("_", combo[1], "_vs_", combo[2])
    )
    
    ## Non-zero counts ##
    cmb1 <- e_data_counts[grouping_info$Group == combo[1]]
    res[[paste0("NonZero_Count_", combo[1])]] <- rowSums(cmb1 != 0)
    
    cmb2 <- e_data_counts[grouping_info$Group == combo[2]]
    res[[paste0("NonZero_Count_", combo[2])]] <- rowSums(cmb2 != 0)
    row.names(res) <- NULL
    res[[edata_cname]] <- as.character(omicsData$e_data[[edata_cname]])
    res[c(ncol(res), 1:(ncol(res) - 1))]
  })
  
  
  #merge all data frames in list
  all_cont <- purrr::reduce(all_res, dplyr::full_join)
  
  count_cols <- grep("^NonZero_Count_", colnames(all_cont))

  lfc_cols <- grep("^log2FoldChange", colnames(all_cont))
  padj_cols <- grep("^padj", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  
  # colnames(all_cont)[-1] <- gsub("^baseMean", "Mean", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^log2FoldChange", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^padj", "P_value", colnames(all_cont)[-1])
  
  all_cont[[1]] <- as.character(all_cont[[1]])

  results <- cbind(all_cont[c(1, count_cols)],
                   all_cont[c(lfc_cols, padj_cols, flag_cols)])
  
  # attr_list <- c("cnames", "data_info", "filters", "group_DF")
  # keep_attr <- attributes(omicsData)[names(attributes(results)) %in% attr_list]
  # attributes(results) <- c(attributes(results), keep_attr)
  
  attr(results, "cnames") <- attr(omicsData, "cnames")
  attr(results, "data_info") <- attr(omicsData, "data_info")
  attr(results, "filters") <- attr(omicsData, "filters")
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  
  ## sig totes
  flag_df <- reshape2::melt(
    all_cont[flag_cols], 
    variable.name = "Comparison", 
    value.name = "Flags"
  )
  flag_df$Comparison <- gsub("Flag_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = TRUE),
      Down_total = sum(Flags < 0, na.rm = TRUE),
      row.names = NULL
    )
  
  attr(results, "comparisons") <- purrr::map_chr(interesting_comparisons, function(n){
    combo <- cob_list[,n]
    paste0(combo[1], "_vs_", combo[2])
  })
  attr(results, "MA_means") <- group_means
  attr(results, "statistical_test") <- paste0("DESeq_", test)
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  return(results)
  
}

#' Wrapper for edgeR workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by
#'   \code{\link{as.seqData}}
#' @param method a character string of length one specifying which wrapper to
#'   use. Can be 'edgeR', 'DESeq2', or 'limma-voom'
#' @param p_adjust Character string for p-value correction method, refer to
#'   ?p.adjust() for valid options. Defaults to "BH" (Benjamini & Hochberg).
#' @param comparisons `data.frame` with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param p_cutoff Numeric value between 0 and 1 for setting p-value
#'   significance threshold
#' @param ... additional arguments passed to methods functions. Note, formatting
#'   option changes will interfere with wrapping functionality.
#'
#' @return statRes object
#'
#' @details Runs default edgeR workflow with empirical Bayes quasi-likelihood
#'   F-tests. Additional arguments can be passed for use in the function, refer
#'   to calcNormFactors() and glmQLFit() in edgeR package.
#' 
#' Flags (signatures) -	Indicator of statistical significance for computed test. 
#' Zeroes indicate no significance, while +/- 1 indicates direction of significance.
#' 
#' @references 
#'  Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package 
#'  for differential expression analysis of digital gene expression data.” 
#'  Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.
#' 
edgeR_wrapper <- function(
    omicsData, p_adjust = "BH", comparisons = NULL, p_cutoff = 0.05,
    ...
){
  
  l <- list(...)
  
  NF_args <- l[names(l) %in% names(formals(edgeR::calcNormFactors.DGEList))]
  QLFit_args <- l[names(l) %in% names(formals(edgeR::glmQLFit.DGEList))]
  
  ## Get useful variables ##
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  
  ## Snag appropriate formula ##
  group_res <- get_group_formula(omicsData)
  grouping_info <- group_res[[1]]
  grouping_formula <- group_res[[2]]
  design_matrix_edgeR <- model.matrix(as.formula(grouping_formula), grouping_info)
  
  ## Warning for levels that are not "correct" in R
  if(!is.numeric(grouping_info$Group) && 
     !identical(grouping_info$Group, make.names(grouping_info$Group))){
    stop("Main effect levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  ## Warning for levels that are not "correct" in R
  if(!is.null(attr(grouping_info, "covariates"))){
    
    covar <- attr(grouping_info, "covariates")[-1]
    for(i in 1:ncol(covar)){
      if(!is.numeric(covar[[i]]) && !identical(covar[[i]], make.names(covar[[i]]))){
        stop("Covariate levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
      }
    }
  }
  
  ## Warning for levels that are not "correct" in R
  if(identical(attr(grouping_info, "main_effects"), "no_main_effect") && 
     !identical(as.character(grouping_info[[attr(grouping_info, "pair_group")]]), make.names(
       grouping_info[[attr(grouping_info, "pair_group")]]))){
    stop("Pair group column levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  ## Warning for levels that are not "correct" in R
  if(identical(attr(grouping_info, "main_effects"), "no_main_effect") && 
     !identical(as.character(grouping_info[[attr(grouping_info, "pair_id")]]), make.names(
       grouping_info[[attr(grouping_info, "pair_id")]]))){
    stop("Pair id column levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  all_cols <- stringr::str_trim(
    stringr::str_split(
      stringr::str_remove(
        grouping_formula, "~0 \\+"), " \\+ ")[[1]])
  grouping_contrasts <- grouping_info[all_cols]
  grouping_contrasts <- grouping_contrasts[!apply(grouping_contrasts, 2, is.numeric)]
  contrast_levels <- as.character(unique(unlist(grouping_contrasts)))
  cob_list <- combn(contrast_levels, 2)
  all_contrasts <- apply(cob_list, 2, paste, collapse = "-")
  
  ## Make sure comparisons are doing the thing ##
  if(is.null(comparisons) && 
     !identical(attr(grouping_info, "main_effects"), "no_main_effect")){
    
    comparisons <- c()
    comparisons <- unique(grouping_info[["Group"]])
    cob_list_group <- combn(comparisons, 2)
    all_contrasts_group <- apply(cob_list_group, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_group)
    
  } else if (is.null(comparisons)){
    
    comparisons <- c()
    comparisons <- c(unique(as.character(grouping_info[[attr(grouping_info, "pair_group")]])))
    cob_list_pair <- combn(comparisons, 2)
    all_contrasts_pair <- apply(cob_list_pair, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_pair)
    
  } else {
    
    if(ncol(cob_list) == 1){
      all_contrasts <- unique(c(all_contrasts, paste(cob_list[nrow(cob_list):1,], collapse = "-")))
      
    } else {
      all_contrasts <- unique(c(all_contrasts,
                                apply(cob_list[nrow(cob_list):1,], 2, paste, collapse = "-")))
    }

    cob_list <- cbind(cob_list, cob_list[nrow(cob_list):1,])
    all_contrasts_comp <- paste(comparisons$Control, comparisons$Test, sep = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_comp)
    if(length(interesting_comparisons) == 0) stop("Invalid comparisons given")
  }
  
  ## Slice frames
  edata_egdeR <- edgeR::DGEList(omicsData$e_data[-1]) 
  
  ## EdgeR workflow ##
  list_defaults <- list(
    object = edata_egdeR
  )
  
  run_NF <- c(NF_args, list_defaults)
  run_NF <- run_NF[!duplicated(names(run_NF))]
  
  norm_factors_edgeR <- do.call(edgeR::calcNormFactors, run_NF)
  
  ## Stick in attributes for MA plots
  norm_factors_find <- norm_factors_edgeR
  norm_factors_apply <- as.data.frame(sweep(as.data.frame(norm_factors_find$counts), MARGIN=2,
        FUN="/",STATS= norm_factors_find$samples$norm.factors))
  group_means <- suppressMessages(purrr::map_dfc(1:ncol(design_matrix_edgeR), function(col){
    group_cols <- norm_factors_apply[as.logical(design_matrix_edgeR[,col])]
    means_group <- apply(group_cols, 1, mean, na.rm = T)
    data.frame(means_group)
  }))
  
  colnames(group_means) <- unique(grouping_info$Group)
  
  group_means <- cbind(omicsData$e_data[get_edata_cname(omicsData)], group_means)
  

  D_edgeR <- edgeR::estimateDisp(norm_factors_edgeR,
                      design_matrix_edgeR)
  
  list_defaults <- list(
    y = D_edgeR,
    design = design_matrix_edgeR
  )
  
  run_QLFit <- c(QLFit_args, list_defaults)
  run_QLFit <- run_QLFit[!duplicated(names(run_QLFit))]
  
  fit_edgeR <- do.call(edgeR::glmQLFit, run_QLFit)
  
  ## Get all the contrasts
  res_contrasts <- purrr::map(interesting_comparisons, function(n){
    
    combo <- cob_list[,n]
    checker <- colnames(fit_edgeR$coefficients)
    
    text_remove <- c("Group", 
                     attr(grouping_info, "pair_id"), 
                     attr(grouping_info, "pair_group"),
                     colnames(attr(grouping_info, "covariates"))[-1])
    
    text_remove <- text_remove[!(text_remove %in% checker)]
    text_remove <- rev(text_remove[order(nchar(text_remove), text_remove)])
    
    
    for(el in text_remove){
      checker <- sub(el, "", checker)
    }
    
    checkin <- purrr::map_lgl(purrr::map(combo, stringr::str_detect, string = checker), any)

    
    if(all(checkin)){
      
      CONTRASTS <- limma::makeContrasts(
        contrasts = all_contrasts[n], 
        levels = checker)
      
      res_stats <- edgeR::glmQLFTest(fit_edgeR, contrast = CONTRASTS)
      
    } else {
      
      get_coef <- which(checker %in% combo[checkin])
      res_stats <- edgeR::glmQLFTest(fit_edgeR, coef = get_coef)
      
    }
    
    res <- edgeR::topTags(res_stats, n = Inf, adjust.method = p_adjust, 
                          sort.by = "none")
    res <- as.data.frame(res$table)
    
    sig_col <- if("FDR" %in% colnames(res)) "FDR" else "FWER"
    
    # Flag stuffs ----------------------------------------------------------------
    sigs <- which(res[[sig_col]] < p_cutoff)
    res[["Flag"]] <- 0
    if(length(sigs) > 0){
      res[["Flag"]][sigs] <- sign(res[["logFC"]][sigs])
    }
    
    colnames(res) <- paste0(colnames(res), 
                            paste0("_", combo[1], "_vs_", combo[2]))
    
    ## Non-zero counts and means##
    cmb1 <- e_data_counts[grouping_info$Group == combo[1]]
    if(length(cmb1) == 0) cmb1 <- e_data_counts[grouping_info[[attr(grouping_info, "pair_group")]] == combo[1]]
    
    res[[paste0("NonZero_Count_", combo[1])]] <- rowSums(cmb1 != 0)
    # res[[paste0("Mean_", combo[1])]] <- apply(cmb1, 1, mean, na.rm = TRUE)
    
    cmb2 <- e_data_counts[grouping_info$Group == combo[2]]
    if(length(cmb2) == 0) cmb2 <- e_data_counts[grouping_info[[attr(grouping_info, "pair_group")]] == combo[2]]
    
    res[[paste0("NonZero_Count_", combo[2])]] <- rowSums(cmb2 != 0)
    # res[[paste0("Mean_", combo[2])]] <- apply(cmb2, 1, mean, na.rm = TRUE)
    
    res[[get_edata_cname(omicsData)]] <- as.character(omicsData$e_data[[get_edata_cname(omicsData)]])
    row.names(res) <- NULL
    res[c(ncol(res), 1:(ncol(res) - 1))]
  })
  
  all_cont <- res_contrasts[purrr::map_int(res_contrasts, nrow) != 0]
  all_cont <-  purrr::reduce(all_cont, dplyr::full_join)
  
  count_cols <- grep("^NonZero_Count_", colnames(all_cont))
  # mean_cols <- grep("^Mean", colnames(all_cont))
  lfc_cols <- grep("^logFC", colnames(all_cont))
  # pval_cols <- grep(colnames(all_cont), "_pvalue")
  padj_cols <- grep("^(FDR|FWER)", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  
  colnames(all_cont)[-1] <- gsub("^logFC", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^(FDR|FWER)", "P_value", colnames(all_cont)[-1])
  
  ## Ordering
  results <- all_cont[c(1, count_cols, #mean_cols, 
                        lfc_cols, padj_cols, flag_cols)]
  
  # attr_list <- c("cnames", "data_info", "filters", "group_DF")
  # keep_attr <- attributes(omicsData)[names(attributes(results)) %in% attr_list]
  # attributes(results) <- c(attributes(results), keep_attr)
  
  attr(results, "cnames") <- attr(omicsData, "cnames")
  attr(results, "data_info") <- attr(omicsData, "data_info")
  attr(results, "filters") <- attr(omicsData, "filters")
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  
  ## sig totes
  flag_df <- reshape2::melt(all_cont[flag_cols], 
                            variable.name = "Comparison", 
                            value.name = "Flags")
  flag_df$Comparison <- gsub("Flag_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = TRUE),
      Down_total = sum(Flags < 0, na.rm = TRUE),
      row.names = NULL
    )
  
  attr(results, "comparisons") <- purrr::map_chr(interesting_comparisons, function(n){
    combo <- cob_list[,n]
    paste0(combo[1], "_vs_", combo[2])
  })
  attr(results, "MA_means") <- group_means
  attr(results, "statistical_test") <- "EdgeR_F"
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  return(results)
  
}


#' Wrapper for limma-voom workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by
#'   \code{\link{as.seqData}}
#' @param method a character string of length one specifying which wrapper to
#'   use. Can be 'edgeR', 'DESeq2', or 'limma-voom'
#' @param p_adjust Character string for p-value correction method, refer to
#'   ?p.adjust() for valid options
#' @param comparisons `data.frame` with columns for "Control" and "Test"
#'   containing the different comparisons of interest. Comparisons will be made
#'   between the Test and the corresponding Control  If left NULL, then all
#'   pairwise comparisons are executed.
#' @param p_cutoff Numeric value between 0 and 1 for setting p-value
#'   significance threshold
#' @param ... additional arguments passed to methods functions. Note, formatting
#'   option changes will interfere with wrapping functionality.
#'
#' @details Runs default limma-voom workflow using empirical Bayes moderated
#'   t-statistics. Additional arguments can be passed for use in the function,
#'   refer to calcNormFactors() in edgeR package.
#'   
#' @return statRes object
#' 
#' @references 
#'  Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, 
#'  G.K. (2015). limma powers differential expression analyses for 
#'  RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#' 
#' 
voom_wrapper <- function(
    omicsData,  p_adjust = "BH", comparisons = NULL, p_cutoff = 0.05, ...
){
  
  l <- list(...)
  
  NF_args <- l[names(l) %in% names(formals(edgeR::calcNormFactors.DGEList))]
  
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  
  ## Snag appropriate formula ##
  group_res <- get_group_formula(omicsData)
  grouping_info <- group_res[[1]]
  grouping_formula <- group_res[[2]]
  
  ## Warning for levels that are not "correct" in R
  if(!identical(as.character(grouping_info$Group), make.names(grouping_info$Group))){
    stop("Main effect levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  ## Warning for levels that are not "correct" in R
  if(identical(attr(grouping_info, "main_effects"), "no_main_effect") && 
     !identical(as.character(grouping_info[[attr(grouping_info, "pair_group")]]), make.names(
       grouping_info[[attr(grouping_info, "pair_group")]]))){
    stop("Pair group column levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  ## Warning for levels that are not "correct" in R
  if(identical(attr(grouping_info, "main_effects"), "no_main_effect") && 
     !identical(as.character(grouping_info[[attr(grouping_info, "pair_id")]]), make.names(
       grouping_info[[attr(grouping_info, "pair_id")]]))){
    stop("Pair id column levels are not in R-acceptable format (A syntactically valid name consists of letters, numbers and the dot or underline characters and starts with a letter or the dot not followed by a number). Limma-voom processing will not be available unless all main effects meet this condition.")
  }
  
  design_matrix_limma <- model.matrix(as.formula(grouping_formula), grouping_info)
  
  # ## Get comparisons

  all_cols <- stringr::str_trim(
    stringr::str_split(
      stringr::str_remove(
        grouping_formula, "~0 \\+"), " \\+ ")[[1]])
  grouping_contrasts <- grouping_info[all_cols]
  grouping_contrasts <- grouping_contrasts[!apply(grouping_contrasts, 2, is.numeric)]
  contrast_levels <- as.character(unique(unlist(grouping_contrasts)))
  cob_list <- combn(contrast_levels, 2)
  all_contrasts <- apply(cob_list, 2, paste, collapse = "-")
  
  ## Make sure comparisons are doing the thing ##
  if(is.null(comparisons) && 
     !identical(attr(grouping_info, "main_effects"), "no_main_effect")){
    
    comparisons <- c()
    comparisons <- unique(grouping_info[["Group"]])
    cob_list_group <- combn(comparisons, 2)
    all_contrasts_group <- apply(cob_list_group, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_group)
    
  } else if (is.null(comparisons)){
    
    comparisons <- c()
    comparisons <- c(unique(as.character(grouping_info[[attr(grouping_info, "pair_group")]])))
    cob_list_pair <- combn(comparisons, 2)
    all_contrasts_pair <- apply(cob_list_pair, 2, paste, collapse = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_pair)
    
  } else {
    
    if(ncol(cob_list) == 1){
      all_contrasts <- unique(c(all_contrasts, paste(cob_list[nrow(cob_list):1,], collapse = "-")))
                              
    } else {
      all_contrasts <- unique(c(all_contrasts,
                                apply(cob_list[nrow(cob_list):1,], 2, paste, collapse = "-")))
    }
    
    cob_list <- cbind(cob_list, cob_list[nrow(cob_list):1,])
    all_contrasts_comp <- paste(comparisons$Control, comparisons$Test, sep = "-")
    interesting_comparisons <- which(all_contrasts %in% all_contrasts_comp)
    if(length(interesting_comparisons) == 0) stop("Invalid comparisons given")
    
  }
  
  
  ## Limma-voom pipeline
  edata_limma <- edgeR::DGEList(e_data_counts) 
  
  list_defaults <- list(
    object = edata_limma
  )
  
  run_NF <- c(NF_args, list_defaults)
  run_NF <- run_NF[!duplicated(names(run_NF))]
  
  norm_factors_limma <- do.call(edgeR::calcNormFactors, run_NF)
  
  norm_factors_find <- norm_factors_limma
  norm_factors_apply <- as.data.frame(sweep(as.data.frame(norm_factors_find$counts), MARGIN=2,
                                            FUN="/",STATS= norm_factors_find$samples$norm.factors))
  norm_factors_apply <- as.data.frame(sweep(as.data.frame(norm_factors_find$counts), MARGIN=2,
                                            FUN="/",STATS= norm_factors_find$samples$norm.factors))
  group_means <- suppressMessages(purrr::map_dfc(1:ncol(design_matrix_limma), function(col){
    group_cols <- norm_factors_apply[as.logical(design_matrix_limma[,col])]
    means_group <- apply(group_cols, 1, mean, na.rm = T)
    data.frame(means_group)
  }))
  
  colnames(group_means) <- unique(grouping_info$Group)
  
  group_means <- cbind(omicsData$e_data[get_edata_cname(omicsData)], group_means)
  
  limma_voom <- limma::voom(norm_factors_limma, design_matrix_limma)
  limma_vfit <- limma::lmFit(limma_voom, design_matrix_limma)
  
  res_contrasts <- purrr::map(interesting_comparisons, function(n){
    
    combo <- cob_list[,n]

    checker <- colnames(limma_vfit$coefficients)
    text_remove <- c("Group", attr(grouping_info, "pair_id"), 
                     attr(grouping_info, "pair_group"),
                     colnames(attr(grouping_info, "covariates"))[-1])
    text_remove <- text_remove[!(text_remove %in% checker)]
    text_remove <- rev(text_remove[order(nchar(text_remove), text_remove)])
    
    for(el in text_remove){
      checker <- sub(el, "", checker)
    }
    
    checkin <- purrr::map_lgl(purrr::map(combo, stringr::str_detect, string = checker), any)

    
    if(all(checkin)){
      
      combo <- cob_list[,n]
      CONTRASTS <- limma::makeContrasts(
        contrasts = all_contrasts[n], 
        levels = checker)
      
      ## Suppress name difference warnings
      tmp <- suppressWarnings(limma::contrasts.fit(limma_vfit, CONTRASTS))
      
    } else {
      
      get_coef <- which(checker %in% combo[checkin])
      ## Suppress name difference warnings
      tmp <- suppressWarnings(limma::contrasts.fit(limma_vfit, coefficients = get_coef))
      
    }

    fit <- limma::eBayes(tmp)
    
    # topTable uses adjust method  = "BH"
    res <- limma::topTable(fit, n = Inf, adjust.method = p_adjust, sort.by = "none")
    if(nrow(res) > 1){
      # Flag stuffs ----------------------------------------------------------------
      sigs <- which(res[["adj.P.Val"]] < p_cutoff)
      res[["Flag"]] <- 0
      if(length(sigs) >0 ){
        res[["Flag"]][sigs] <- sign(res[["logFC"]][sigs])
      }
      
      colnames(res) <- paste0(colnames(res), paste0("_", combo[1], "_vs_", combo[2]))
      
      ## Non-zero counts and means##
      cmb1 <- e_data_counts[grouping_info$Group == combo[1]]
      if(length(cmb1) == 0) cmb1 <- e_data_counts[grouping_info[[attr(grouping_info, "pair_group")]] == combo[1]]
      res[[paste0("NonZero_Count_", combo[1])]] <- rowSums(cmb1 != 0)
      # res[[paste0("Mean_", combo[1])]] <- apply(cmb1, 1, mean, na.rm = TRUE)
      
      cmb2 <- e_data_counts[grouping_info$Group == combo[2]]
      if(length(cmb2) == 0) cmb2 <- e_data_counts[grouping_info[[attr(grouping_info, "pair_group")]] == combo[2]]
      res[[paste0("NonZero_Count_", combo[2])]] <- rowSums(cmb2 != 0)
      # res[[paste0("Mean_", combo[2])]] <- apply(cmb2, 1, mean, na.rm = TRUE)

      res[[get_edata_cname(omicsData)]] <- as.character(omicsData$e_data[[get_edata_cname(omicsData)]])
      row.names(res) <- NULL
      return(res[c(ncol(res), 1:(ncol(res) - 1))])
    } else return()
  })
  
  all_cont <- res_contrasts[purrr::map_int(res_contrasts, nrow) != 0]
  all_cont <-  purrr::reduce(all_cont, dplyr::full_join)
  
  results <- list(Full_results = all_cont)
  
  count_cols <- grep("^NonZero", colnames(all_cont))
  # mean_cols <- grep("^Mean", colnames(all_cont))
  lfc_cols <- grep("^logFC", colnames(all_cont))
  # pval_cols <- grep(colnames(all_cont), "_pvalue")
  padj_cols <- grep("^adj.P.Val", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  colnames(all_cont)[-1] <- gsub("^logFC", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^adj.P.Val", "P_value", colnames(all_cont)[-1])
  
  results <- all_cont[c(1, count_cols, #mean_cols, 
                        lfc_cols, padj_cols, flag_cols)]
  
  # attr_list <- c("cnames", "data_info", "filters", "group_DF")
  # keep_attr <- attributes(omicsData)[names(attributes(results)) %in% attr_list]
  # attributes(results) <- c(attributes(results), keep_attr)
  
  attr(results, "cnames") <- attr(omicsData, "cnames")
  attr(results, "data_info") <- attr(omicsData, "data_info")
  attr(results, "filters") <- attr(omicsData, "filters")
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  
  flag_df <- reshape2::melt(all_cont[flag_cols], 
                            variable.name = "Comparison", 
                            value.name = "Flags")
  flag_df$Comparison <- gsub("Flag_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = TRUE),
      Down_total = sum(Flags < 0, na.rm = TRUE),
      row.names = NULL
    )
  
  attr(results, "comparisons") <- purrr::map_chr(interesting_comparisons, function(n){
    combo <- cob_list[,n]
    paste0(combo[1], "_vs_", combo[2])
  })
  attr(results, "MA_means") <- group_means
  attr(results, "statistical_test") <- "Voom_T"
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  return(results)
  
}


#' Get formula for group design
#' 
#' For generating group design formulas and correctly ordered group data for 
#' seqData statistical functions
#' 
#' @param omicsData an object of type 'seqData', created by \code{\link{as.seqData}}
#' 
#' @examples
#' \dontrun{
#' 
#' library(pmartRdata)
#' myseqData <- rnaseq_object
#' myseqData <- group_designation(omicsData = myseqData, main_effects = "Virus")
#' grouping_info <- pmartR:::get_group_formula(omicsData = myseqData)
#' 
#' }
#' 
#' @rdname get_group_formula
#' @name get_group_formula
#' 
get_group_formula <- function(omicsData){
  
  grouping_info <- get_group_DF(omicsData)
  if(is.null(grouping_info)) stop("group_designation has not been run")
  e_data_index <- which(colnames(omicsData$e_data) %in% get_edata_cname(omicsData))
  collist <- colnames(omicsData$e_data[-e_data_index])
  collist_group <- grouping_info[[get_fdata_cname(omicsData)]]
  grouping_info <- grouping_info[match(collist, collist_group),]

  
  ## If pairs, add to group_df
  pairs <- attr(grouping_info, "pair_id")
  pair_group <- attr(grouping_info, "pair_group")
  pair_denom <- attr(grouping_info, "pair_denom")
  covariates <- attr(grouping_info, "covariates")
  main_effects <- attr(grouping_info, "main_effects")
  
  if(!is.null(pairs)){
    keepcols <- which(colnames(omicsData$f_data) %in% c(pairs, pair_group, 
                                                        get_fdata_cname(omicsData)))
    grouping_info <- dplyr::left_join(grouping_info, omicsData$f_data[keepcols])
    
    grouping_info[[pairs]] <- as.factor(grouping_info[[pairs]] )
    grouping_info[[pair_group]] <- factor(
      grouping_info[[pair_group]], 
      levels = unique(c(pair_denom, grouping_info[[pair_group]]))
      )
    
    if(identical(attr(grouping_info, "main_effects"), "no_main_effect")) {
      design_matrix_add_pairs <- paste(pairs, pair_group, sep = " + ")
    } else {
      
      grouping_info <- dplyr::arrange(grouping_info, 
                               Group, 
                               !!rlang::sym(pair_group), !!rlang::sym(pairs))
      
      all_mult_levels <- table(apply(grouping_info[c("Group", pair_group)], 
                                     1, paste, collapse = ""))
      
      grouping_info[[pairs]] <- unlist(purrr::map(all_mult_levels, function(el){
        paste0("new_pair_name", 1:el)
      }))
      
      design_matrix_add_pairs <- paste(paste(c(pairs, pair_group), " + "), 
                                       collapse = "")
    }
    
    if(!is.null(covariates)){
      
      redun <- purrr::map_lgl(2:ncol(covariates), function(x){
        identical(as.numeric(as.factor(covariates[[x]])), 
                  as.numeric(as.factor(grouping_info[[pairs]])))
      })
   
      if(any(redun)){
        warning("At least 1 detected covariate is confounded with pair_ID specification and will not be used in final model.")
        covariates <- covariates[c(TRUE, !redun)]
      }
    } else design_matrix_add_covariates <- NULL
    
  } else {
    
    design_matrix_add_pairs <- NULL
    
  }
  
  if(!identical(main_effects, "no_main_effect")){
    design_matrix_add_group <- "Group"
    
    if(!is.null(covariates)){
      
      redun <- purrr::map_lgl(2:ncol(covariates), function(x){
        identical(as.numeric(as.factor(covariates[[x]])), 
                  as.numeric(as.factor(grouping_info[["Group"]])))
      })
      
      if(any(redun)){
        warning("At least 1 detected covariate is confounded with Group and will not be used in final model.")
        covariates <- covariates[c(TRUE, !redun)]
      }
    }
    
  } else {
    design_matrix_add_group <- NULL
  }
  
  if(!is.null(covariates) && ncol(covariates) > 1){
    grouping_info <- dplyr::left_join(grouping_info, covariates)
    design_matrix_add_covariates <- paste(paste(colnames(covariates)[-1], "+ "), collapse = "")
  } else design_matrix_add_covariates <- NULL
  
  
  formula_string <- paste0("~0 + ", design_matrix_add_pairs, 
                           design_matrix_add_covariates, 
                           design_matrix_add_group)
  
  return(list(grouping_info, formula_string))
  
}

#' Diagnostic plot for seqData
#'
#' For generating statistics for 'seqData' objects
#'
#' @param omicsData seqData object used to terst dispersions
# @param comparisons Comparisons used in dispersion estimates
#' @param method either "DESeq2", "edgeR", or "voom" for testing dispersion
#' @param interactive Logical. If TRUE produces an interactive plot.
#' @param x_lab A character string specifying the x-axis label when the metric
#'   argument is NULL. The default is NULL in which case the x-axis label will
#'   be "count".
#' @param y_lab A character string specifying the y-axis label. The default is
#'   NULL in which case the y-axis label will be the metric selected for the
#'   \code{metric} argument.
#' @param x_lab_size An integer value indicating the font size for the x-axis.
#'   The default is 11.
#' @param y_lab_size An integer value indicating the font size for the y-axis.
#'   The default is 11.
#' @param x_lab_angle An integer value indicating the angle of x-axis labels.
#' @param title_lab A character string specifying the plot title when the
#'   \code{metric} argument is NULL.
#' @param title_lab_size An integer value indicating the font size of the plot
#'   title. The default is 14.
#' @param legend_lab A character string specifying the legend title.
#' @param legend_position A character string specifying the position of the
#'   legend. Can be one of "right", "left", "top", or "bottom". The default is
#'   "right".
#' @param point_size An integer specifying the size of the points. The default
#'   is 0.2.
#' @param bw_theme Logical. If TRUE uses the ggplot2 black and white theme.
#' @param palette A character string indicating the name of the RColorBrewer
#'   palette to use. For a list of available options see the details section in
#'   \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param custom_theme a ggplot `theme` object to be applied to non-interactive
#'   plots, or those converted by plotly::ggplotly().
#'
#' @details DESeq2 option requires package "survival" to be available.
#'
#' @return plot result
#' 
#' @examples 
#' \dontrun{
#' library(pmartRdata)
#' myseqData <- group_designation(omicsData = rnaseq_object, main_effects = "Virus")
#' dispersion_est(omicsData = myseqData, method = "edgeR")
#' dispersion_est(omicsData = myseqData, method = "DESeq2")
#' dispersion_est(omicsData = myseqData, method = "voom")
#' }
#' 
#' @references 
#'  Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package 
#'  for differential expression analysis of digital gene expression data.” 
#'  Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.
#'  
#'  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and 
#'  dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)
#'  
#'  Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, 
#'  G.K. (2015). limma powers differential expression analyses for 
#'  RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
#'
#' @export
#' @rdname dispersion_est
#' @name dispersion_est
#'
dispersion_est <- function(omicsData, method,
                           interactive = FALSE,
                           x_lab = NULL,
                           x_lab_size = 11,
                           x_lab_angle = NULL,
                           y_lab = NULL,
                           y_lab_size = 11,
                           title_lab = NULL,
                           title_lab_size = 14,
                           legend_lab = NULL,
                           legend_position = "right",
                           bw_theme = TRUE,
                           palette = NULL,
                           point_size = 0.2,
                           custom_theme = NULL){
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("seqData"))) {
    # Throw an error that the input for omicsData is not the appropriate class.
    stop("omicsData must be of class 'seqData'")
  }
  
  # check method #
  if(length(method) != 1 || !(method %in% c('edgeR', 'DESeq2','voom'))){
    stop("method must a single character string of length 1 in 'edgeR', 'DESeq2', or 'voom'")
  }

  ## Set theme ##

  if(!is.null(custom_theme)){
    if(bw_theme)
      warning(paste("Setting both bw_theme to TRUE and specifying a custom",
                    "theme may cause undesirable results",
                    sep = " "))
    if(!inherits(custom_theme, c("theme", "gg")))
      stop("custom_theme must be a valid 'theme' object as used in ggplot")
    mytheme <-  custom_theme
  } else mytheme <-  ggplot2::theme(
    plot.title = ggplot2::element_text(size = title_lab_size),
    axis.title.x = ggplot2::element_text(size = x_lab_size),
    axis.title.y = ggplot2::element_text(size = y_lab_size),
    axis.text.x = ggplot2::element_text(angle = x_lab_angle),
    legend.position = legend_position
  )
  
  ## Get nice things
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  
  if(is.null(get_group_DF(omicsData))) stop(
    "OmicsData requires group_designation for statistical analysis"
  )
  
  ## Snag appropriate formula ##
  group_res <- get_group_formula(omicsData)
  grouping_info <- group_res[[1]]
  grouping_formula <- group_res[[2]]
  design_matrix <- model.matrix(as.formula(grouping_formula), grouping_info)
  
  # Both packages will do the job when it comes to detecting DE genes in a routine analysis. 
  # limma (+ voom) is faster and has access to more methodology (e.g., duplicateCorrelation) 
  # compared to edgeR, courtesy of lots of things being easier when you assume normality. 
  # However, voom does rely on the presence of a well-fitted mean-variance trend to estimate 
  # the precision weights. For some applications (though not usually DE analyses), the spread 
  # of abundances is too low to stably fit the trend, and in such cases edgeR will give better 
  # performance. edgeR also handles low counts better, which is worth considering if you want 
  # to focus on lowly-expressed genes or are dealing with very low-coverage data (e.g., single-cell stuff).
  
  if(method == "DESeq2"){
    
    if(!requireNamespace("survival", quietly = TRUE)){
      stop("package 'survival' required for DESeq2 processing")
    }
    
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab)) "Mean of Normalized Counts" else x_lab
    the_y_label <- if (is.null(y_lab)) "Dispersion" else y_lab
    the_title_label <- if (is.null(title_lab)) "DESeq2 dispersion fit" else title_lab
    the_legend_label <- if (is.null(legend_lab)) "" else legend_lab
    cols <- if (is.null(palette)) c("black","blue", "red") else RColorBrewer::brewer.pal(3, palette)
    
    ## DESeq workflow
    edata_deseq <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(
      e_data_counts, 
      colData = grouping_info,
      design = as.formula(grouping_formula)
    ))
    
    dds <- DESeq2::estimateSizeFactors(edata_deseq)
    dds <- DESeq2::estimateDispersions(dds)
    
    ## only plots for those above 0
    df1 <- as.data.frame(S4Vectors::mcols(dds))
    
    ## Plot
    p <- ggplot2::ggplot(data = df1, 
                         ggplot2::aes(x = baseMean, y = dispGeneEst)) +
      ggplot2::geom_point(color = "black", size = point_size) +
      ggplot2::geom_point(ggplot2::aes(y = dispersion, color = "blue"), alpha = 0.25, size = point_size) + 
      ggplot2::geom_point(ggplot2::aes(y = dispFit, color = "red"), size = point_size) + 
      ggplot2::scale_y_continuous(trans='log10') + 
      ggplot2::scale_x_continuous(trans='log10') + 
      ggplot2::scale_colour_manual(name = legend_lab, 
                                   values =c('black'=cols[1],'red'=cols[3], "blue" = cols[2]), 
                                   labels = c('Gene-est','Fitted', 'Final')) +
      ggplot2::labs(x = the_x_label, y = the_y_label, 
                    title = the_title_label, color = the_legend_label)
    
  } else if (method == "edgeR"){
    
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab)) "Average Log2 CPM" else x_lab
    the_y_label <- if (is.null(y_lab)) "Quarter-Root Mean Deviance" else y_lab
    the_title_label <- if (is.null(title_lab)) "EdgeR dispersion fit" else title_lab
    the_legend_label <- if (is.null(legend_lab)) "" else legend_lab
    cols <- if (is.null(palette)) c("black","blue", "red") else RColorBrewer::brewer.pal(3, palette)
    
    ## EdgeR workflow
    edata_egdeR <- edgeR::DGEList(e_data_counts) 
    norm_factors_edgeR <- edgeR::calcNormFactors(edata_egdeR)
    
    D_edgeR <- edgeR::estimateDisp(norm_factors_edgeR,
                                   design_matrix)
    
    fit_edgeR <- edgeR::glmQLFit(D_edgeR, design_matrix)

    df2 <- data.frame(
      CD = D_edgeR$common.dispersion,
      TD = D_edgeR$trended.dispersion,
      TagD = D_edgeR$tagwise.dispersion,
      fitD1 = fit_edgeR$var.prior,
      fitD2 = fit_edgeR$var.post,
      AveLogCPM = fit_edgeR$AveLogCPM
    )
    
    ## Remake of plotBCV(GTagD_edgeR, log = "y")
    df2_melt <- reshape2::melt(df2[c("AveLogCPM", "TagD", "TD", "CD")], id.var = "AveLogCPM")
    p <- ggplot2::ggplot(data = df2_melt,
                         ggplot2::aes(x = AveLogCPM, y = sqrt(value), color = variable)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::scale_colour_manual(name = legend_lab, 
                                   values = c( TagD = cols[1], 
                                               TD = cols[2], 
                                               CD = cols[3]),
                                   labels = c(TagD = 'Tagwise', CD = 'Common',TD = 'Trend')) +
      ggplot2::scale_y_continuous(trans='log10') + 
      ggplot2::labs(x = the_x_label, y = the_y_label, 
                    title = the_title_label, color = the_legend_label)
    
    # p <- ggplot2::ggplot(data = df2, 
    #                      ggplot2::aes(x = AveLogCPM, y = (fitD2)^(1/4), color = "black")) +
    #   ggplot2::geom_point( size = point_size) +
    #   ggplot2::geom_line(ggplot2::aes(y = (fitD1)^(1/4), color = "red")) +
    #   ggplot2::scale_colour_manual(name = legend_lab, 
    #                       values =c('black'=cols[1],'red'=cols[2]), 
    #                       labels = c('Squeezed Dispersion','Trend')) +
    #   ggplot2::labs(x = the_x_label, y = the_y_label, 
    #                 title = the_title_label, color = the_legend_label)
    
    
  } else if (method == "voom"){
    
    # Farm boy, do all the tedious label crap. As you wish.
    the_x_label <- if (is.null(x_lab)) "Sqrt (Standard Deviation)" else x_lab
    the_y_label <- if (is.null(y_lab)) "Log (count size + 0.5)" else y_lab
    the_title_label <- if (is.null(title_lab)) "Limma-Voom dispersion Fit" else title_lab
    the_legend_label <- if (is.null(legend_lab)) "" else legend_lab
    cols <- if (is.null(palette)) c("black", "red") else RColorBrewer::brewer.pal(2, palette)
    
    edata_limma <- edgeR::DGEList(e_data_counts) 
    norm_factors_limma <- edgeR::calcNormFactors(edata_limma)
    mean_norm_counts <- norm_factors_limma
    
    limma_voom <- limma::voom(norm_factors_limma, design_matrix, save.plot = TRUE)
    limma_vfit <- limma::lmFit(limma_voom, design_matrix)
    
    efit <- limma::eBayes(limma_vfit)
    
    df3 <- data.frame(
      x_disp = limma_voom$voom.xy$x,
      y_disp = limma_voom$voom.xy$y,
      x_fit = limma_voom$voom.line$x,
      y_fit = limma_voom$voom.line$y
    )
    
    y_fitted <- efit$sigma
    p <- ggplot2::ggplot(data = df3, 
                         ggplot2::aes(x = x_disp, y = y_disp, color = "black")) +
      ggplot2::geom_point( size = point_size) +
      ggplot2::geom_point(ggplot2::aes(x = x_fit, y = y_fit, color = "red"), size = point_size) +
      ggplot2::scale_y_continuous(trans='log10') +
      ggplot2::scale_colour_manual(name = legend_lab, 
                                   values =c('black'=cols[1],'red'=cols[2]), 
                                   labels = c('Mean-Varience','Trend')) +
      ggplot2::labs(x = the_x_label, y = the_y_label, 
                    title = the_title_label, color = the_legend_label)
    
  }
  
  if(bw_theme) p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"))
  
  p <- p + mytheme
  
  if(interactive)
    return(plotly::ggplotly(p, tooltip = c("text"))) else
      return(p)
  
}

