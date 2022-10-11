#' Panomics Marketplace - Quality Control and Statistical Analysis for Panomics Data
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
