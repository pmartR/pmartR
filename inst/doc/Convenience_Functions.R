## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(pmartR)
library(pmartRdata)
mypro <- pro_object
mymetab <- metab_object
mynmr <- nmr_identified_object

## -----------------------------------------------------------------------------
# get a statRes object using the proData object
mypro <- group_designation(omicsData = mypro, main_effects = "Phenotype")
myfilter <- imdanova_filter(omicsData = mypro)
mypro <- applyFilt(filter_object = myfilter, omicsData = mypro, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)
mystats <- imd_anova(omicsData = mypro, test_method = "combined")

get_data_class(mystats)

## -----------------------------------------------------------------------------
get_data_info(omicsData = mypro)

## -----------------------------------------------------------------------------
get_data_scale(omicsObject = mypro)
get_data_scale(omicsObject = mymetab)
get_data_scale(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_data_scale_orig(omicsObject = mypro)
get_data_scale_orig(omicsObject = mymetab)
get_data_scale_orig(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_edata_cname(omicsObject = mypro)
get_edata_cname(omicsObject = mymetab)
get_edata_cname(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_emeta_cname(omicsObject = mypro)
get_emeta_cname(omicsObject = mymetab)
get_emeta_cname(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_fdata_cname(omicsObject = mypro)
get_fdata_cname(omicsObject = mymetab)
get_fdata_cname(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_group_DF(omicsData = mypro)
get_group_DF(omicsData = mymetab)
get_group_DF(omicsData = mynmr)

## -----------------------------------------------------------------------------
get_group_table(omicsObject = mypro)

## -----------------------------------------------------------------------------
myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
myiso_normalized <- normalize_isobaric(
  omicsData = myiso,
  exp_cname = "Plex",
  apply_norm = TRUE,
  refpool_cname = "Virus",
  refpool_notation = "Pool"
)
get_isobaric_info(omicsData = myiso_normalized)

## -----------------------------------------------------------------------------
get_meta_info(omicsData = mynmr)

## -----------------------------------------------------------------------------
get_nmr_info(omicsData = mynmr)

mynmr <- edata_transform(
  omicsData = nmr_identified_object,
  data_scale = "log2"
)

## -----------------------------------------------------------------------------
# Normalization using a "spiked in" metabolite
nmr_norm <- normalize_nmr(
  omicsData = mynmr, apply_norm = TRUE,
  metabolite_name = "unkm1.53",
  backtransform = TRUE
)
get_nmr_info(omicsData = nmr_norm)

## -----------------------------------------------------------------------------
# Normalization using a sample property
nmr_norm <- normalize_nmr(
  omicsData = mynmr, apply_norm = TRUE,
  sample_property_cname = "Concentration",
  backtransform = TRUE
)
get_nmr_info(omicsData = nmr_norm)

## -----------------------------------------------------------------------------
get_filters(omicsData = mypro)

## -----------------------------------------------------------------------------
get_data_norm(omicsObject = mypro)
get_data_norm(omicsObject = mymetab)
get_data_norm(omicsObject = mynmr)

## -----------------------------------------------------------------------------
get_isobaric_norm(myiso)
get_isobaric_norm(myiso_normalized)

## -----------------------------------------------------------------------------
get_nmr_norm(omicsData = mynmr)
get_nmr_norm(omicsData = nmr_norm)

## -----------------------------------------------------------------------------
imd_anova_res <- imd_anova(
  omicsData = mypro,
  test_method = "comb",
  pval_adjust_a_multcomp = "bon",
  pval_adjust_g_multcomp = "bon"
)

get_comparisons(compObj = imd_anova_res)

