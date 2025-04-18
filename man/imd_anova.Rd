% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imd_anova.R
\name{imd_anova}
\alias{imd_anova}
\title{Test for a qualitative and quantitative difference between groups using IMD
and ANOVA, respectively}
\usage{
imd_anova(
  omicsData,
  comparisons = NULL,
  test_method,
  pval_adjust_a_multcomp = "none",
  pval_adjust_g_multcomp = "none",
  pval_adjust_a_fdr = "none",
  pval_adjust_g_fdr = "none",
  pval_thresh = 0.05,
  equal_var = TRUE,
  model_selection = "auto",
  parallel = TRUE
)
}
\arguments{
\item{omicsData}{pmartR data object of any class, which has a `group_df`
attribute created by the `group_designation()` function}

\item{comparisons}{data frame with columns for "Control" and "Test"
containing the different comparisons of interest. Comparisons will be made
between the Test and the corresponding Control. If left NULL, then all
pairwise comparisons are executed.}

\item{test_method}{character string specifying the filter method to use:
"combined", "gtest", or "anova". Specifying "combined" implements both the gtest and
anova filters.}

\item{pval_adjust_a_multcomp}{character string specifying the type of multiple
comparison adjustment to implement for ANOVA tests. Valid options include:
"bonferroni", "holm", "tukey", and "dunnett". The default is "none" which
corresponds to no p-value adjustment.}

\item{pval_adjust_g_multcomp}{character string specifying the type of multiple
comparison adjustment to implement for G-test tests. Valid options include:
"bonferroni" and "holm". The default is "none" which corresponds to no
p-value adjustment.}

\item{pval_adjust_a_fdr}{character string specifying the type of FDR
adjustment to implement for ANOVA tests. Valid options include:
"bonferroni", "BH", "BY", and "fdr". The default is "none" which corresponds
to no p-value adjustment.}

\item{pval_adjust_g_fdr}{character string specifying the type of FDR
adjustment to implement for G-test tests. Valid options include:
"bonferroni", "BH", "BY", and "fdr". The default is "none" which corresponds
to no p-value adjustment.}

\item{pval_thresh}{numeric p-value threshold, below or equal to which
biomolecules are considered differentially expressed. Defaults to 0.05}

\item{equal_var}{logical; should the variance across groups be assumed equal?}

\item{model_selection}{Character, one of 'full', 'reduced', or 'auto'
indicating the model to be used in the ANOVA analysis. The default 'auto' performs an F-test
to determine if the full model is necessary. If the F-test is significant,
the full model is used, otherwise the reduced model is used.  'full' uses all main
effects, covariates, and interactions between main effects. 'reduced' does 
not consider interactions between main effects (only covariates and marginal main effects).}

\item{parallel}{logical value indicating whether or not to use a
"doParallel" loop when running the G-Test with covariates. Defaults to
TRUE.}
}
\value{
An object of class 'statRes', which is a data frame containing
  columns (when relevant based on the test(s) performed) for: e_data cname,
  group counts, group means, ANOVA p-values, IMD p-values, fold change
  estimates on the same scale as the data (e.g. log2, log10, etc.), and fold
  change significance flags (0 = not significant; +1 = significant and
  positive fold change (ANOVA) or more observations in test group relative to
  reference group (IMD); -1 = significant and negative fold change (ANOVA) or
  fewer observations in test group relative to reference group (IMD))
}
\description{
This is the IMD-ANOVA test defined in Webb-Robertson et al. (2010).
}
\examples{
\dontshow{if (identical(tolower(Sys.getenv("NOT_CRAN")), "true") & requireNamespace("pmartRdata", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
library(pmartRdata)
# Transform the data
mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")

# Group the data by condition
mymetab <- group_designation(omicsData = mymetab, main_effects = c("Phenotype"))

# Apply the IMD ANOVA filter
imdanova_Filt <- imdanova_filter(omicsData = mymetab)
mymetab <- applyFilt(filter_object = imdanova_Filt, omicsData = mymetab, min_nonmiss_anova = 2)

# Implement IMD ANOVA and compute all pairwise comparisons 
# (i.e. leave the comparisons argument NULL), with FDR adjustment
anova_res <- imd_anova(omicsData = mymetab, test_method = "anova",
                       pval_adjust_a_multcomp = "Holm", pval_adjust_a_fdr = "BY")
imd_res <- imd_anova(omicsData = mymetab, test_method = "gtest",
                     pval_adjust_g_multcomp = "bon", pval_adjust_g_fdr = "BY")
imd_anova_res <- imd_anova(omicsData = mymetab, test_method = "combined",
                           pval_adjust_a_fdr = "BY", pval_adjust_g_fdr = "BY")
imd_anova_res <- imd_anova(omicsData = mymetab, test_method = "combined",
                           pval_adjust_a_multcomp = "bon", pval_adjust_g_multcomp = "bon",
                           pval_adjust_a_fdr = "BY", pval_adjust_g_fdr = "BY")

# Two main effects and a covariate
mymetab <- group_designation(omicsData = mymetab, main_effects = c("Phenotype", "SecondPhenotype"),
                             covariates = "Characteristic")
imd_anova_res <- imd_anova(omicsData = mymetab, test_method = 'comb')

# Same but with custom comparisons
comp_df <- data.frame(
  check.names = FALSE, 
  Control = c("Phenotype1", "A"), 
  Test = c("Phenotype2", "B")
)
custom_comps_res <- imd_anova(omicsData = mymetab, comparisons = comp_df, test_method = "combined")
\dontshow{\}) # examplesIf}
}
\references{
Webb-Robertson, Bobbie-Jo M., et al. "Combined statistical
analyses of peptide intensities and peptide occurrences improves
identification of significant peptides from MS-based proteomics data."
Journal of proteome research 9.11 (2010): 5748-5756.
}
\author{
Bryan Stanfill, Kelly Stratton
}
