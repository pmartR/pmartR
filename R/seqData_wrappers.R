## File contains wrappers for DESeq2, EdgeR, and limma-voom as well as plotting function for limma-voom


#' Wrapper for DESeq2 workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param test an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param p_adjust an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param comparisons an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param fitType an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param sfType an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param full an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param reduced an object of type 'seqData', created by \code{\link{as.seqData}} 
#' @param quiet an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param minReplicatesForReplace an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param modelMatrixType an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param useT an object of type 'seqData', created by \code{\link{as.seqData}} 
#' @param parallel an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param BPPARAM an object of type 'seqData', created by \code{\link{as.seqData}} 
#' @param lfcThreshold an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param altHypothesis an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param independentFiltering an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param alpha an object of type 'seqData', created by \code{\link{as.seqData}} 
#' @param format an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param addMLE an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param plotMA an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param plotDispEsts an object of type 'seqData', created by \code{\link{as.seqData}} 
#' @param p_cutoff an object of type 'seqData', created by \code{\link{as.seqData}}  
#' 
#' @return data.frame object
#' 
#' 
#' @examples
#' \dontrun{
#' }
#' 
#' @export
#' @rdname Deseq2_wrapper
#' @name Deseq2_wrapper
#' 
Deseq2_wrapper <- function(
    omicsData, test = "Wald",  p_adjust = "BH", comparisons = NULL,
    fitType = NULL, sfType = NULL, betaPrior = FALSE,
    full = NULL, reduced = NULL, quiet = TRUE,
    minReplicatesForReplace = Inf, modelMatrixType = NULL,
    useT = FALSE, minmu = 0.5,
    parallel = FALSE, BPPARAM = NULL, lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    # cooksCutoff, filter, theta, filterFun, 
    independentFiltering = TRUE, alpha = 0.1, 
    format = "DataFrame",
    addMLE = FALSE,
    plotMA = F, plotDispEsts = F,
    p_cutoff = 0.05
){
  
  ## Normal p_value adjustments, test can either be "Wald" or liklihood ratio "LRT"
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  grouping_info <- get_group_DF(omicsData)[colnames(omicsData$f_data) != fdata_cname]
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  
  if(is.null(comparisons)) comparisons <- unique(grouping_info[["Group"]])
  
  edata_deseq <- DESeq2::DESeqDataSetFromMatrix(
    e_data_counts, 
    colData = grouping_info,
    design = ~Group ## Figure out design designation
  )
  
  run_stats_deseq <- DESeq2::DESeq(
    edata_deseq,
    test = test, ## reasonable to allow change in pmartR
    
    ## maybe not all the options
    fitType = "parametric", #fitting of dispersions to the mean intensity
    
    #    parametric - fit a dispersion-mean relation of the form:
    # 
    #         dispersion = asymptDisp + extraPois / mean
    # 
    #       via a robust gamma-family GLM. The coefficients asymptDisp and extraPois are given in the attribute      coefficients of the dispersionFunction of the object.
    # 
    #   local - use the locfit package to fit a local regression of log dispersions over log base mean (normal scale means and dispersions are input and output for dispersionFunction). The points are weighted by normalized mean count in the local regression.
    # 
    # mean - use the mean of gene-wise dispersion estimates.
    # 
    # glmGamPoi - use the glmGamPoi package to fit the gene-wise dispersion, its trend and calculate the MAP based on the quasi-likelihood framework. The trend is calculated using a local median regression.
    
    # sfType = c("ratio", "poscounts", "iterate"), # type of size factor estimation,
    #Method for estimation: either "ratio", "poscounts", or "iterate". "ratio" uses the standard median ratio method introduced in DESeq. The size factor is the median ratio of the sample over a "pseudosample": for each gene, the geometric mean of all samples. "poscounts" and "iterate" offer alternative estimators, which can be used even when all genes contain a sample with a zero (a problem for the default method, as the geometric mean becomes zero, and the ratio undefined). The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts. This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples. The "iterate" estimator iterates between estimating the dispersion with a design of ~1, and finding a size factor vector by numerically optimizing the likelihood of the ~1 model.
    
    quiet = T,
    
    ## Yeah probs set to inf
    minReplicatesForReplace = Inf, ## cook's distance used to flag outliers, require at least n replicates to replace flagged outliers-- defined as .99 quantile of the F(p, m - p) distribution, where p is the number of parameters and m is the number of samples. Replacement is values predicted by the trimmed mean over all samples (and adjusted by size factor or normalization factor). Inf disables replacement
    
    modelMatrixType = "standard", ## If we need anything fancier than what is defined by "Group" we should consider this argument for the glm use case -- betapriors req for expanded
    
    parallel = T,
    
    ## Waldy
    betaPrior = F, 
    #whether or not to put a zero-mean normal prior on the non-intercept coefficients See nbinomWaldTest for description of the calculation of the beta prior.In versions >=1.16, the default is set to FALSE, and shrunken LFCs are obtained afterwards using lfcShrink.
    
    useT = F ## only for waldtest, uses a t-distribution instead of a normal distribution
    
    ## Full vs reduced design models to use for testingg
    
    # minmu = NULL ##  lower bound on the estimated count for fitting gene-wise dispersion and for use with nbinomWaldTest and nbinomLRT. If fitType="glmGamPoi", then 1e-6 will be used (as this fitType is optimized for single cell data, where a lower minmu is recommended), otherwise the default value as evaluated on bulk datasets is 0.5
  )
  
  cob_list <- utils::combn(comparisons, 2)
  
  all_res <- purrr::map(1:ncol(cob_list), function(combo_n){
    
    combo <- cob_list[,combo_n]
    
    ## Run tests
    res <- results(
      run_stats_deseq, 
      lfcThreshold = 0,
      altHypothesis = "greaterAbs",
      contrast=c("Group", combo[1], combo[2]),
      cooksCutoff = FALSE, 
      ## Fiters by distance, default it filters .99 quantile of the F(p, m-p) 
      ## distribution, where p is the number of coefficients being fitted 
      ## and m is the number of samples. Excludes groups w/ only 2 samples
      independentFiltering = FALSE,
      pAdjustMethod = p_adjust,
      # format = "DataFrame",
      tidy = T,
      parallel = T
      
      # alpha = 0.1, used for independent filtering
      # filterFun = NULL, ## You can specify a custom filtering method,
      # addMLE = F ## backwards compatibility argument, specifies if the "unshrunken" maximum likelihood estimates (MLE) of log2 fold change should be added as a column to the results table 
      # minmu = NULL ## Lower bound on estimated count (used when calulating contrasts)
    )
    
    # Flag stuffs ----------------------------------------------------------------
    sigs <- which(res[["padj"]] < p_cutoff)
    res[[paste0("Flag_", test)]] <- 0
    if(length(sigs)>0){
      res[[paste0("Flag_", test)]][sigs] <- sign(res$log2FoldChange[sigs])
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
    res[[edata_cname]] <- omicsData$e_data[[edata_cname]]
    res[c(ncol(res), 1:(ncol(res) - 1))]
  })
  
  #merge all data frames in list
  all_cont <- all_res %>% purrr::reduce(full_join)
  
  count_cols <- grep("^NonZero_Count_", colnames(all_cont))
  mean_cols <- grep("^Mean", colnames(all_cont))
  lfc_cols <- grep("^log2FoldChange", colnames(all_cont))
  # pval_cols <- grep(colnames(all_cont), "_pvalue")
  padj_cols <- grep("^padj", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  
  colnames(all_cont)[-1] <- gsub("^baseMean", "Mean", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^log2FoldChange", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^padj", paste0("P_value_", test), colnames(all_cont)[-1])
  
  results <- all_cont[c(1, count_cols, mean_cols, 
                        lfc_cols, padj_cols, flag_cols)]
  
  attr_list <- c("cnames", "data_info", "filters", "group_DF")
  keep_attr <- attributes(results)[names(attributes(results)) %in% attr_list]
  attributes(results) <- c(attributes(results), keep_attr)
  
  ## sig totes
  flag_df <- reshape2::melt(
    all_cont[flag_cols], 
    variable.name = "Comparison", 
    value.name = "Flags"
    )
  flag_df$Comparison <- gsub("Flag_(Wald|LRT)_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = T),
      Down_total = sum(Flags < 0, na.rm = T),
      row.names = NULL
  )
  
  attr(results, "comparisons") <- comparisons
  attr(results, "statistical_test") <- paste0("DESeq_", test)
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  if(plotDispEsts){
    attr(results, "plotDispEsts") <- run_stats_deseq
    DESeq2::plotDispEsts(run_stats_deseq)
  }
  if(plotMA){
    attr(results, "plotMA") <- run_stats_deseq
    DESeq2::plotMA(run_stats_deseq)
  }
  
  return(results)
  
}

#' Wrapper for EdgeR workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param p_adjust an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param comparisons an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param p_cutoff an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param plotMDS an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param plotBCV an object of type 'seqData', created by \code{\link{as.seqData}} 
#' 
#'
#' @return data.frame object
#' 
#' 
#' @examples
#' \dontrun{
#' }
#' 
#' @export
#' @rdname EdgeR_wrapper
#' @name EdgeR_wrapper
#' 
EdgeR_wrapper <- function(
    omicsData, p_adjust = "BH", comparisons = NULL, p_cutoff = 0.05,
    plotMDS = F, plotBCV = F
){
  
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  grouping_info <- get_group_DF(omicsData)[colnames(omicsData$f_data) != fdata_cname]
  
  ### Combinations ###
  if(is.null(comparisons)) comparisons <- unique(grouping_info$Group)
  cob_list <- combn(comparisons, 2)
  all_contrasts <- purrr::map_chr(
    1:ncol(cob_list), 
    function(combo_n) paste0(cob_list[,combo_n], collapse = "-")
  )
  
  edata_egdeR <- edgeR::DGEList(e_data_counts) 
  norm_factors_edgeR <- edgeR::calcNormFactors(edata_egdeR)
  
  design_matrix_edgeR <- model.matrix(~Group, grouping_info)
  
  GCD_edgeR <- edgeR::estimateGLMCommonDisp(norm_factors_edgeR,
                                     design_matrix_edgeR)
  
  GTD_edgeR <- edgeR::estimateGLMTrendedDisp(GCD_edgeR,
                                      design_matrix_edgeR,
                                      method="power")
  
  GTagD_edgeR <- edgeR::estimateGLMTagwiseDisp(GTD_edgeR, design_matrix_edgeR)
  
  fit_edgeR <- edgeR::glmQLFit(GTagD_edgeR, design_matrix_edgeR)
  

  #### Make p-value adjustments across # of comparisons too ??
  ## Row numbers in excess of edata nrow
  res_contrasts <- purrr::map(1:length(all_contrasts), function(n){
    combo <- cob_list[,n]
    CONTRASTS <- limma::makeContrasts(
      contrasts = all_contrasts[n], 
      levels = comparisons)
    res_stats <- edgeR::glmQLFTest(fit_edgeR, contrast = CONTRASTS)
    res <- edgeR::topTags(res_stats, n = Inf, adjust.method = p_adjust, 
                   sort.by = "none")
    res <- as.data.frame(res$table)
    
    sig_col <- if("FDR" %in% colnames(res)) "FDR" else "FWER"

    # Flag stuffs ----------------------------------------------------------------
    sigs <- which(res[[sig_col]] < p_cutoff)
    res[["Flag_LRT"]] <- 0
    if(length(sigs) > 0){
      res[["Flag_LRT"]][sigs] <- sign(res[["logFC"]][sigs])
    }
    
    colnames(res) <- paste0(colnames(res), 
                            paste0("_", combo[1], "_vs_", combo[2]))
    
    ## Non-zero counts and means##
    cmb1 <- e_data_counts[grouping_info$Group == combo[1]]
    res[[paste0("NonZero_Count_", combo[1])]] <- rowSums(cmb1 != 0)
    res[[paste0("Mean_", combo[1])]] <- apply(cmb1, 1, mean, na.rm = T)
    
    cmb2 <- e_data_counts[grouping_info$Group == combo[2]]
    res[[paste0("NonZero_Count_", combo[2])]] <- rowSums(cmb2 != 0)
    res[[paste0("Mean_", combo[2])]] <- apply(cmb2, 1, mean, na.rm = T)
    
    res[[get_edata_cname(omicsData)]] <- row.names(res)
    row.names(res) <- NULL
    res[c(ncol(res), 1:(ncol(res) - 1))]
  })
  
  all_cont <- res_contrasts[map_int(res_contrasts, nrow) != 0] %>% 
    reduce(full_join)
  
  # mean_cols <- str_detect(colnames(all_cont), "^logCPM")
  # f_cols <- str_detect(colnames(all_cont), "^F_")
  
  count_cols <- grep("^NonZero_Count_", colnames(all_cont))
  mean_cols <- grep("^Mean", colnames(all_cont))
  lfc_cols <- grep("^logFC", colnames(all_cont))
  # pval_cols <- grep(colnames(all_cont), "_pvalue")
  padj_cols <- grep("^(FDR|FWER)", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  
  colnames(all_cont)[-1] <- gsub("^logFC", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^(FDR|FWER)", "P_value_LRT", colnames(all_cont)[-1])
  
  ## Ordering
  results <- all_cont[c(1, count_cols, mean_cols, 
                        lfc_cols, padj_cols, flag_cols)]
  
  attr_list <- c("cnames", "data_info", "filters", "group_DF")
  keep_attr <- attributes(results)[names(attributes(results)) %in% attr_list]
  attributes(results) <- c(attributes(results), keep_attr)
  
  ## sig totes
  flag_df <- reshape2::melt(all_cont[flag_cols], 
                            variable.name = "Comparison", 
                            value.name = "Flags")
  flag_df$Comparison <- gsub("Flag_LRT_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = T),
      Down_total = sum(Flags < 0, na.rm = T),
      row.names = NULL
    )
  
  attr(results, "comparisons") <- comparisons
  attr(results, "statistical_test") <- "EdgeR_LRT"
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  # viz
  if(plotMDS){
    attr(results, "plotMDS") <- norm_factors_edgeR
    plotMDS(norm_factors_edgeR)
  }
  if(plotBCV){
    attr(results, "plotBCV") <- GTagD_edgeR
    plotBCV(GTagD_edgeR)
  }

  return(results)
  
}


#' Wrapper for limma-voom workflow
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param omicsData an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param p_adjust an object of type 'seqData', created by \code{\link{as.seqData}}
#' @param comparisons an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param p_cutoff an object of type 'seqData', created by \code{\link{as.seqData}}  
#' @param plotVoom an object of type 'seqData', created by \code{\link{as.seqData}}
#' 
#'
#' @return data.frame object
#' 
#' 
#' @examples
#' \dontrun{
#' }
#' 
#' @export
#' @rdname Voom_wrapper
#' @name Voom_wrapper
#' 
Voom_wrapper <- function(
    omicsData,  p_adjust = "BH", comparisons = NULL, p_cutoff = 0.05,
    plotVoom = F
){
  
  edata_cname <- get_edata_cname(omicsData)
  fdata_cname <- get_fdata_cname(omicsData)
  e_data_counts <- omicsData$e_data[colnames(omicsData$e_data) != edata_cname]
  grouping_info <- get_group_DF(omicsData)[colnames(omicsData$f_data) != fdata_cname]
  
  if(is.null(comparisons)) comparisons <- unique(grouping_info$Group)
  cob_list <- combn(comparisons, 2)
  
  edata_limma <- edgeR::DGEList(e_data_counts) 
  norm_factors_limma <- edgeR::calcNormFactors(edata_limma)
  design_matrix_limma <- model.matrix(~0 + Group, grouping_info)
  limma_voom <- limma::voom(norm_factors_limma, design_matrix_limma, save.plot = plotVoom)
  limma_vfit <- limma::lmFit(limma_voom, design_matrix_limma)
  
  all_contrasts <- purrr::map_chr(1:ncol(cob_list), 
                           function(combo_n) paste0(paste0("Group", 
                                                           cob_list[,combo_n]), 
                                                    collapse = "-"))
  
  res_contrasts <- purrr::map(1:length(all_contrasts), function(n){
    combo <- cob_list[,n]
    CONTRASTS <- limma::makeContrasts(
      contrasts = all_contrasts[n], 
      levels = paste0("Group", comparisons))
    tmp <- limma::contrasts.fit(limma_vfit, CONTRASTS)
    fit <- limma::eBayes(tmp)
    
    # topTable uses adjust method  = "BH"
    res <- limma::topTable(fit, n = Inf, adjust.method = p_adjust, sort.by = "none")
    
    if(nrow(res) > 1){
      
      # Flag stuffs ----------------------------------------------------------------
      sigs <- which(res[["adj.P.Val"]] < p_cutoff)
      res[["Flag_T"]] <- 0
      if(length(sigs) >0 ){
        res[["Flag_T"]][sigs] <- sign(res[["logFC"]][sigs])
      }
      
      colnames(res) <- paste0(colnames(res), paste0("_", combo[1], "_vs_", combo[2]))
      
      ## Non-zero counts and means##
      cmb1 <- e_data_counts[grouping_info$Group == combo[1]]
      res[[paste0("NonZero_Count_", combo[1])]] <- rowSums(cmb1 != 0)
      res[[paste0("Mean_", combo[1])]] <- apply(cmb1, 1, mean, na.rm = T)
      
      cmb2 <- e_data_counts[grouping_info$Group == combo[2]]
      res[[paste0("NonZero_Count_", combo[2])]] <- rowSums(cmb2 != 0)
      res[[paste0("Mean_", combo[2])]] <- apply(cmb2, 1, mean, na.rm = T)
      
      res[[get_edata_cname(omicsData)]] <- row.names(res)
      row.names(res) <- NULL
      return(res[c(ncol(res), 1:(ncol(res) - 1))])
    } else return()
  })
  
  all_cont <- res_contrasts[map_int(res_contrasts, nrow) != 0] %>% 
    reduce(full_join)
  
  count_cols <- grep("^NonZero", colnames(all_cont))
  mean_cols <- grep("^Mean", colnames(all_cont))
  lfc_cols <- grep("^logFC", colnames(all_cont))
  # pval_cols <- grep(colnames(all_cont), "_pvalue")
  padj_cols <- grep("^adj.P.Val", colnames(all_cont))
  flag_cols <- grep("^Flag", colnames(all_cont))
  
  colnames(all_cont)[-1] <- gsub("^logFC", "Fold_change", colnames(all_cont)[-1])
  colnames(all_cont)[-1] <- gsub("^adj.P.Val", "P_value_T", colnames(all_cont)[-1])
  
  
  results <- all_cont[c(1, count_cols, mean_cols, 
                        lfc_cols, padj_cols, flag_cols)]
  
  attr_list <- c("cnames", "data_info", "filters", "group_DF")
  keep_attr <- attributes(results)[names(attributes(results)) %in% attr_list]
  attributes(results) <- c(attributes(results), keep_attr)
  
  ## sig totes
  flag_df <- reshape2::melt(all_cont[flag_cols], 
                            variable.name = "Comparison", 
                            value.name = "Flags")
  flag_df$Comparison <- gsub("Flag_T_", "", flag_df$Comparison )
  attr(results, "number_significant") <- flag_df %>%
    dplyr::group_by(Comparison) %>%
    dplyr::summarise(
      Up_total = sum(Flags > 0, na.rm = T),
      Down_total = sum(Flags < 0, na.rm = T),
      row.names = NULL
    )
  
  attr(results, "comparisons") <- comparisons
  attr(results, "statistical_test") <- "Voom T"
  attr(results, "adjustment_method") <- p_adjust
  attr(results, "pval_thresh") <- p_cutoff
  attr(results, "data_class") <- attr(omicsData, "class")
  class(results) <- c("statRes", class(results))
  
  # viz
  if(plotVoom){
    attr(results, "plotVoom") <- limma_voom
    plotVoom(limma_voom)
  }
  
  return(results)
  
}


#' Plot for voom wrapper
#' 
#' For generating statistics for 'seqData' objects
#' 
#' @param limma_voom attribute passed in statRes object from seqData where LimmaVoom was run
#'
#' @return plot result
#' 
#' @examples
#' \dontrun{
#' }
#' 
#' @export
#' @rdname plotVoom
#' @name plotVoom
#' 
plotVoom <- function(limma_voom){
  df_points <- data.frame(
    x =  limma_voom$voom.xy$x,
    y = limma_voom$voom.xy$y
  )
  
  df_line <- data.frame(
    x_line =  limma_voom$voom.line$x,
    y_line = limma_voom$voom.line$y
  )
  
  ggplot2::ggplot(
    data = df_points, 
    mapping = ggplot2::aes(x = x, y = y)) + 
    ggplot2::geom_point() +
    ggplot2::geom_line(inherit.aes = F, data = df_line, 
                       mapping = ggplot2::aes(x = x_line, y = y_line), color = "red") +
    ggplot2::labs(x =  limma_voom$voom.xy$xlab, y = limma_voom$voom.xy$ylab) +
    ggplot2::theme_bw()
}
