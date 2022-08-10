## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

library(pmartR)
library(pmartRdata)

## ----specify_channels---------------------------------------------------------


## ----specify_refpool----------------------------------------------------------


## -----------------------------------------------------------------------------
data(nmr_object_identified)
summary(nmr_object_identified)

# log2 transform the values
nmr_object <- edata_transform(nmr_object_identified, "log2")

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

## -----------------------------------------------------------------------------
# global median centering - don't apply the norm, show info in the normRes object



## -----------------------------------------------------------------------------
# global median centering - apply the norm



## -----------------------------------------------------------------------------
# RIP mean centering - don't apply the norm (so we can look at the number of molecules in the subset)



## -----------------------------------------------------------------------------
# LOS MAD



## -----------------------------------------------------------------------------
# PPP z-score



## -----------------------------------------------------------------------------
# PPP-RIP median



## -----------------------------------------------------------------------------
# Complete mean



## -----------------------------------------------------------------------------
# 



## -----------------------------------------------------------------------------
# 



## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------


