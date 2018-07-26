library(pmartR)
library(pmartRdata)
library(testthat)
#Trasnform the data
attr(pep_object,"cnames")$fdata_cname <- "SampleID"
mypepData <- edata_transform(omicsData = pep_object, data_scale = "log2")

#Group the data by condition
mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))

#Apply the IMD ANOVA filter
imdanova_Filt <- imdanova_filter(omicsData = mypepData)
mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)

##---- Check that warnings and errors are produced in various situations --------##
context("Check p-value adjustments")
expect_warning(imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='dunnett'))
expect_warning(imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='tukey'))

#Test with really big dataset
library(OvarianPepdataBPsubset)
suppressWarnings(tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp_subset$e_data, f_data = tcga_ovarian_pepdata_bp_subset$f_data, e_meta = tcga_ovarian_pepdata_bp_subset$e_meta,
                                      edata_cname = 'Peptide', fdata_cname = "sampleID", emeta_cname = "Protein", data_scale = "log2"))
suppressWarnings(tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("race")))
expect_warning(expect_error(ovarian_res_dunnett <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'dunnett', test_method='combined')))
