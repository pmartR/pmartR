#' Panomics Marketplace - Quality Control and Statistical Analysis for Panomics Data
#'
#' @description Provides functionality for quality control processing and
#'  statistical analysis of mass spectrometry (MS) omics data, in particular 
#'  proteomic (either at the peptide or the protein level), lipidomic, and 
#'  metabolomic data. This includes data transformation, specification of 
#'  groups that are to be compared against each other, filtering of features 
#'  and/or samples, data normalization, data summarization (correlation, PCA), 
#'  and statistical comparisons between defined groups.
#'
#' @docType package
#' @importFrom Rcpp evalCpp
#'
#' @importFrom stats anova aov as.formula binomial complete.cases cor dist ecdf formula glm hclust lm mad median model.matrix na.omit p.adjust.methods pbinom pchisq qchisq quantile runif sd setNames t.test var
#' @importFrom utils capture.output combn data head stack tail
#'
#' @useDynLib pmartR
#' @name pmartR
NULL

utils::globalVariables(
  c(".", ":=", ":::", "::")
)
