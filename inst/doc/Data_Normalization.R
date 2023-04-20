## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)
library(pmartR)
library(pmartRdata)

## ----load_example_isobaricpepData---------------------------------------------
data(isobaric_object)
class(isobaric_object)

isobaric_object <- edata_transform(isobaric_object,
  data_scale = "log2"
)

## ----specify_refpool----------------------------------------------------------
# do not apply the normalization until we have looked at the distribution of peptides within the reference pool samples
iso_norm <- normalize_isobaric(isobaric_object,
  exp_cname = "Plex",
  apply_norm = FALSE,
  refpool_cname = "Virus",
  refpool_notation = "Pool"
)

# look at the distribution of peptides within the reference pool samples
plot(iso_norm)

# now apply the normalization
isobaric_object_normalized <- normalize_isobaric(isobaric_object,
  exp_cname = "Plex",
  apply_norm = TRUE,
  refpool_cname = "Virus",
  refpool_notation = "Pool"
)

# look at boxplots of the data after normalization to the reference pool samples
plot(isobaric_object_normalized)

## ----specify_channels, eval = FALSE-------------------------------------------
#  # not run because this data does not actually contain channel_cname or refpool_channel information
#  iso_norm <- normalize_isobaric(isobaric_object,
#    exp_cname = "Plex",
#    apply_norm = FALSE,
#    channel_cname = "SampleChannel", # this column in f_data would specify a number corresponding to the channel for each sample
#    refpool_channel = 4
#  ) # this value in the SampleChannel column would always correspond to a reference pool sample

## -----------------------------------------------------------------------------
summary(nmr_identified_object)

# log2 transform the values
nmr_object <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")

## -----------------------------------------------------------------------------
# don't apply the normalization yet
normRes_property <-
  normalize_nmr(nmr_object,
    apply_norm = FALSE,
    sample_property_cname = "Concentration"
  )
plot(normRes_property)

## -----------------------------------------------------------------------------
# now apply the normalization
nmr_norm_property <-
  normalize_nmr(nmr_object,
    apply_norm = TRUE,
    sample_property_cname = "Concentration",
    backtransform = TRUE
  )
plot(nmr_norm_property)

## -----------------------------------------------------------------------------
# don't apply the normalization yet
normRes_reference <-
  normalize_nmr(nmr_object,
    apply_norm = FALSE,
    metabolite_name = "unkm1.53"
  )
plot(normRes_reference)

## -----------------------------------------------------------------------------
# now apply the normalization
nmr_norm_reference <- normalize_nmr(nmr_object,
  apply_norm = TRUE,
  metabolite_name = "unkm1.53"
)
plot(nmr_norm_reference)

## -----------------------------------------------------------------------------
# global median centering - don't apply the norm, show info in the normRes object
mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")

mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "all",
  norm_fn = "median",
  apply_norm = FALSE
)
class(mynorm)
# plot(mynorm)

## -----------------------------------------------------------------------------
# global median centering - apply the norm
mypep_norm <- normalize_global(
  omicsData = mypep,
  subset_fn = "all",
  norm_fn = "median",
  apply_norm = TRUE,
  backtransform = TRUE
)
class(mypep_norm)
plot(mypep_norm)

## -----------------------------------------------------------------------------
# run group_designation
mypep <- group_designation(omicsData = mypep, main_effects = "Phenotype")

# RIP mean centering - don't apply the norm (so we can look at the number of molecules in the subset)
mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "rip",
  params = list(rip = (0.2)),
  norm_fn = "median",
  apply_norm = FALSE
)

# we can see how many biomolecules are included in this subset
mynorm$n_features_calc
mynorm$prop_features_calc

## -----------------------------------------------------------------------------
# LOS MAD
mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "los",
  params = list(los = (0.1)),
  norm_fn = "mad",
  apply_norm = FALSE
)
mynorm$prop_features_calc

## -----------------------------------------------------------------------------
# PPP z-score
mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "ppp",
  params = list(ppp = (0.5)),
  norm_fn = "zscore",
  apply_norm = FALSE,
  min_prop = 0.2
)
mynorm$prop_features_calc

## -----------------------------------------------------------------------------
# PPP-RIP median
mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "ppp_rip",
  params = list(ppp = (0.5), rip = (0.2)),
  norm_fn = "median",
  apply_norm = FALSE
)
mynorm$prop_features_calc

## -----------------------------------------------------------------------------
# Complete mean
mynorm <- normalize_global(
  omicsData = mypep,
  subset_fn = "complete",
  norm_fn = "mean",
  apply_norm = FALSE
)
mynorm$prop_features_calc

## ----spans, eval = FALSE------------------------------------------------------
#  # returns a data frame arranged by descending SPANS score
#  # not run due to long runtime; SPANS plot generated and saved to include in vignette
#  spans_result <- spans_procedure(mypep)
#  
#  plot(spans_result)

## ----spans plot, out.width = "600px", echo=FALSE------------------------------
knitr::include_graphics("SPANS_pep.png")

## -----------------------------------------------------------------------------
mypep_norm <- normalize_global(
  omicsData = mypep,
  subset_fn = "all",
  norm_fn = "median",
  apply_norm = TRUE,
  backtransform = TRUE
)

## -----------------------------------------------------------------------------
plot(mypep, order_by = "Phenotype", color_by = "Phenotype")
mypep_loess <- normalize_loess(omicsData = mypep)
plot(mypep_loess, order_by = "Phenotype", color_by = "Phenotype")

## -----------------------------------------------------------------------------
mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")

# restrict to complete samples using the molecule filter
myfilter <- molecule_filter(omicsData = mymetab)
nsamps <- get_data_info(mymetab)$num_samps
mymetab <- applyFilt(filter_object = myfilter, omicsData = mymetab, min_num = nsamps)

plot(mymetab, order_by = "Phenotype", color_by = "Phenotype")
mymetab_quantile <- normalize_quantile(omicsData = mymetab)
plot(mymetab_quantile, order_by = "Phenotype", color_by = "Phenotype")

