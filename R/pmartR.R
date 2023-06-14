#' Panomics Marketplace - Quality Control and Statistical Analysis for Panomics Data
#'
#' @description Provides functionality for quality control processing and
#'  statistical analysis of mass spectrometry (MS) omics data, in particular 
#'  proteomic (either at the peptide or the protein level), lipidomic, and 
#'  metabolomic data, as well as RNA-seq based count data and nuclear magnetic
#'  resonance (NMR) data. This includes data transformation, specification of 
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
  c(
    ".", ":=", ":::", "::",
    "Absent", "Absent Count", "Absent Proportion", "Abundance", 
    "AveLogCPM", "baseMean", "Batch", "bins", "Biomolecule", "both_cov", 
    "both_me", "color", "Color", "Comparison", "Count", "Count_biomolecules", 
    "Count_First_Group", "Count_Second_Group", "counts", "CV", "Data Type", 
    "dcast", "desc", "Direction", "dispersion", "dispFit", "dispGeneEst", 
    "el", "Flags", "fold_change", "Fold_change", "Fold_change_flag", 
    "Freq", "frequency_counts", "Group", "Group.x", "Group.y", "i", 
    "insufficient", "lcpm", "LibrarySize", "Log2.md", "melt", 
    "metabolite_cname", "Metric", "Min_obs", "molecule", "mols_used_in_norm",
    "n", "N", "n_combine", "n_control", "n_cov", "n_groups", "n_me",
    "n_peps_used", "n_test", "name", "Nested_DF", "NonZero",
    "normalization_method", "nuff", "num_NA", "num_observations", "P_value",
    "p_value_anova", "paired_diff", "Panel By Choice", "parameters", "PC1",
    "PC2", "peps_per_pro", "Peptide", "posneg", "Present", "Present Count", 
    "Present Proportion", "ProportionNonZero", "Protein_Isoform", 
    "pval", "res", "Sample", "SampleID", "setnames", "sig", "Significance",
    "skew", "SPANS_score", "ss_par", "Statistic", "subset_method", "survfit", 
    "Total", "Total_Counts", "Type", "v", "value", "values", "var1", 
    "var2", "variable", "whichtest", "x_disp", "x_fit", "y", "y_disp", 
    "y_fit"
  )
)
