#' Function to take raw output of `imd_anova` and create output for `statRes` object
#'
#' @param imd_anova_out data frame containing the results of the
#'   \code{imd_anova} call.
#' @param omicsData pmartR data object of any class, which has a `group_df`
#'   attribute that is usually created by the `group_designation()` function
#' @param comparisons character vector of comparison names, e.g. c("A_vs_B", "B_vs_C", ...)
#' @param test_method test method used ("anova", "gtest", or "combined")
#' @param pval_adjust_a_multcomp character string specifying which type of multiple
#'   comparison adjustment was implemented for ANOVA tests. Valid options include:
#'   "bonferroni", "holm", "tukey", and "dunnett".
#' @param pval_adjust_g_multcomp character string specifying which type of multiple
#'   comparison adjustment was implemented for G-tests. Valid options include:
#'   "bonferroni" and "holm".
#' @param pval_adjust_a_fdr character string specifying which type of FDR
#'   adjustment was implemented for ANOVA tests. Valid options include:
#'   "bonferroni", "BH", "BY", and "fdr".
#' @param pval_adjust_g_fdr character string specifying which type of FDR
#'   adjustment was implemented for G-tests. Valid options include:
#'   "bonferroni", "BH", "BY", and "fdr".
#' @param pval_thresh numeric p-value threshold value
#'
#' @return object of class statRes
#'
statRes_output <- function (imd_anova_out,
                            omicsData,
                            comparisons,
                            test_method,
                            pval_adjust_a_multcomp,
                            pval_adjust_g_multcomp,
                            pval_adjust_a_fdr,
                            pval_adjust_g_fdr,
                            pval_thresh) {

  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData", "nmrData"))) stop("omicsData is not an object of appropriate class")

   # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }

  # Create data matrix with just G-test flags to determine number of significant
  # fold changes.
  imd_flags <- imd_anova_out[, grepl("^Flag_G", names(imd_anova_out))] %>%
    data.matrix()

  # Check if imd_flags has zero columns. If it does this means there are no
  # columns for G-test flags and it needs to be set to a zero matrix with one
  # row and column. This will allow it to be used in the apply function later.
  if (dim(imd_flags)[[2]] == 0) {

    imd_flags <- matrix(0, nrow = 1, ncol = 1)

  }

  anova_flags <- imd_anova_out[, grepl("^Flag_A", names(imd_anova_out))] %>%
    data.matrix()

  # Check if anova_flags has zero columns. If it does this means there are no
  # columns for ANOVA flags and it needs to be set to a zero matrix with one
  # row and column. This will allow it to be used in the apply function later.
  if (dim(anova_flags)[[2]] == 0) {

    anova_flags <- matrix(0, nrow = 1, ncol = 1)

  }

  # Add attributes -------------------------------------------------------------

  attr(imd_anova_out, "group_DF") <- groupData
  attr(imd_anova_out, "comparisons") <- comparisons

  # All this flag crap is really more trouble than it is worth!!
  up_a <- apply(anova_flags, 2, function(x){length(which(x == 1))})
  down_a <- apply(anova_flags, 2, function(x){length(which(x == -1))})
  up_g <- apply(imd_flags, 2, function(x){length(which(x == 1))})
  down_g <- apply(imd_flags, 2, function(x){length(which(x == -1))})
  attr(imd_anova_out, "number_significant") <-  data.frame(
    Comparison = comparisons,
    Up_total = up_a + up_g,
    Down_total = down_a + down_g,
    Up_anova = up_a,
    Down_anova = down_a,
    Up_gtest = up_g,
    Down_gtest = down_g,
    row.names = NULL
  )
  attr(imd_anova_out,"statistical_test") <- test_method
  attr(imd_anova_out, "adjustment_method_a_multcomp") <- pval_adjust_a_multcomp
  attr(imd_anova_out, "adjustment_method_g_multcomp") <- pval_adjust_g_multcomp
  attr(imd_anova_out, "adjustment_method_a_fdr") <- pval_adjust_a_fdr
  attr(imd_anova_out, "adjustment_method_g_fdr") <- pval_adjust_g_fdr
  attr(imd_anova_out, "pval_thresh") <- pval_thresh
  attr(imd_anova_out, "data_info") <- attr(omicsData, "data_info")

  class(imd_anova_out) <- c("statRes", "data.frame")

  return (imd_anova_out)

}

#' Summary of statRes Object
#'
#' Provide summary information about statRes objects
#'
#' @name statRes-class
#' @seealso See \code{\link{imd_anova}}
#'
#' @exportClass statRes
setOldClass("statRes")

#' @export
#' @method summary statRes
summary.statRes <- function(x, ...) {
  
  # add if-statement for seqData class since the ANOVA and g-test language is 
  # not applicable for seqData
  if (attributes(x)$data_class == "seqData") {
    cat("Type of test:", attr(x, "statistical_test"), "\n\n")
    cat("Multiple comparison adjustment:",
        attr(x, "adjustment_method"),
        "\n\n")
    cat("p-value threshold:", attr(x, "pval_thresh"), "\n\n")
    cat(
      "Number of significant biomolecules by comparison. Columns specify fold change direction and type of test:\n\n"
    )
    
    table <- attr(x, "number_significant")
    names(table) <- lapply(names(table), function(name) {
      switch(name,
             Up_total = "Total:Positive",
             Down_total = "Total:Negative",
             name)
    })
    
    rownames(table) <- NULL
    
    print(table)
    
    return(invisible(
      list(
        test_type = attr(x, "statistical_test"),
        adjustment = attr(x, "adjustment_method"),
        pval_thresh = attr(x, "pval_thresh"),
        sig_table = table,
        comparisons = attr(x, "comparisons"),
        group_DF = attr(x, "group_DF")
      )
    ))
    
    
  }else{
  # not statRes object from seqData
  cat("Type of test:", attr(x, "statistical_test"), "\n\n")
  cat("Multiple comparison adjustment ANOVA:",
      attr(x, "adjustment_method_a"),
      "\n\n")
  cat("Multiple comparison adjustment G-test:",
      attr(x, "adjustment_method_g"),
      "\n\n")
  cat("p-value threshold:", attr(x, "pval_thresh"), "\n\n")
  cat(
    "Number of significant biomolecules by comparison. Columns specify fold change direction and type of test:\n\n"
  )
  
  table <- attr(x, "number_significant")
  names(table) <- lapply(names(table), function(name) {
    switch(
      name,
      Up_total = "Total:Positive",
      Down_total = "Total:Negative",
      Up_anova = "Positive:ANOVA",
      Down_anova = "Negative:ANOVA",
      Up_gtest = "Positive:G-test",
      Down_gtest = "Negative:G-test",
      name
    )
  })
  
  rownames(table) <- NULL
  
  print(table)
  
  return(invisible(
    list(
      test_type = attr(x, "statistical_test"),
      adjustment_a = attr(x, "adjustment_method_a"),
      adjustment_g = attr(x, "adjustment_method_g"),
      pval_thresh = attr(x, "pval_thresh"),
      sig_table = table,
      comparisons = attr(x, "comparisons"),
      group_DF = attr(x, "group_DF")
    )
  ))
  
  
}

}

#' @export
#' @method print statRes
print.statRes <- function (x) {
  x <- format.data.frame(x)
  
  blank_row = rep("----", 5)
  
  if (nrow(x) >= 9) {
    # Nab the first and last four rows of the statRes data frame.
    statman_head <- head(x, 4)[, 1:min(ncol(x), 5)]
    statman_tail <- tail(x, 4)[, 1:min(ncol(x), 5)]
    
    # Combine the first and last four rows of the statRes data frame with some
    # pretty dashes for a pleasing look when the statRes object is printed in
    # the console.
    statman <- rbind(statman_head, blank_row, statman_tail)
    
  } else {
    # Just print the entire data frame because it is only a little guy.
    statman <- x
    
  }
  
  if (ncol(x) > 5)
    message("only first 5 columns are shown")
  cat("statRes\n")
  cat(capture.output(statman), sep = "\n")
  cat("\n")
  
}
