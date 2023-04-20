## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(warning = FALSE, fig.width = 8, fig.height = 6)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----logo, out.width = "100px", echo=FALSE------------------------------------
knitr::include_graphics("pmartR_logo_final.jpg")

library(pmartR)
library(pmartRdata)
library(ggplot2)
library(reshape2)

## -----------------------------------------------------------------------------
mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
plot(mymetab, color_by = "Phenotype", order_by = "Phenotype")

## -----------------------------------------------------------------------------
head(mymetab$f_data$SampleID)

# specify new names using delim and components arguments
mymetab_shorter_names <- custom_sampnames(omicsData = mymetab, delim = "_", components = c(1, 2))
plot(mymetab, use_VizSampNames = TRUE, color_by = "Phenotype", order_by = "Phenotype")

## -----------------------------------------------------------------------------
plot(mymetab, color_by = "Phenotype", order_by = "Phenotype")
plot(mymetab, use_VizSampNames = FALSE, color_by = "Phenotype", order_by = "Phenotype")

## -----------------------------------------------------------------------------
plot(mymetab, order_by = "Phenotype", color_by = "Phenotype") +
  theme_dark()

