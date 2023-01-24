## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

library(pmartR)
library(pmartRdata)

## ----load_example_isobaricpepData---------------------------------------------
data(isobaric_object)
class(isobaric_object)

isobaric_object <- edata_transform(isobaric_object, 
                                   data_scale = "log2")

## ----specify_refpool----------------------------------------------------------
# do not apply the normalization until we have looked at the distribution of peptides within the reference pool samples
iso_norm <- normalize_isobaric(isobaric_object, 
                               exp_cname = "Plex",
                               apply_norm = FALSE,
                               refpool_cname = "Virus",
                               refpool_notation = "Pool")

# look at the distribution of peptides within the reference pool samples
plot(iso_norm)

# now apply the normalization
isobaric_object_normalized <- normalize_isobaric(isobaric_object, 
                               exp_cname = "Plex",
                               apply_norm = TRUE,
                               refpool_cname = "Virus",
                               refpool_notation = "Pool")

# look at boxplots of the data after normalization to the reference pool samples
plot(isobaric_object_normalized)

## ----specify_channels, eval = FALSE-------------------------------------------
#  # not run because this data does not actually contain channel_cname or refpool_channel information
#  iso_norm <- normalize_isobaric(isobaric_object,
#                                 exp_cname = "Plex",
#                                 apply_norm = FALSE,
#                                 channel_cname = "SampleChannel", # this column in f_data would specify a number corresponding to the channel for each sample
#                                 refpool_channel = 4) # this value in the SampleChannel column would always correspond to a reference pool sample

## -----------------------------------------------------------------------------
summary(nmr_identified_object)

# log2 transform the values
nmr_object <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")

## -----------------------------------------------------------------------------
# don't apply the normalization yet
normRes_property <-
  normalize_nmr(nmr_object,
                apply_norm = FALSE,
                sample_property_cname = "Concentration")
plot(normRes_property)

## -----------------------------------------------------------------------------
# now apply the normalization
nmr_norm_property <-
  normalize_nmr(nmr_object,
                apply_norm = TRUE,
                sample_property_cname = "Concentration",
                backtransform = TRUE)
plot(nmr_norm_property)

## -----------------------------------------------------------------------------
# don't apply the normalization yet
normRes_reference <-
  normalize_nmr(nmr_object,
                apply_norm = FALSE,
                metabolite_name = "unkm1.53")
plot(normRes_reference)

## -----------------------------------------------------------------------------
# now apply the normalization
nmr_norm_reference <- normalize_nmr(nmr_object,
                                    apply_norm = TRUE,
                                    metabolite_name = "unkm1.53")
plot(nmr_norm_reference)

