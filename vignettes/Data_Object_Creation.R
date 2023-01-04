## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(fig.width=8)
knitr::opts_chunk$set(fig.height=6)

## -----------------------------------------------------------------------------
library(pmartR)

## -----------------------------------------------------------------------------
# install the pmartRdata package, if needed
# devtools::install_github("pmartR/pmartRdata")

# load the pmartRdata package
library(pmartRdata)

# load example e_data, f_data, e_meta for the lipid negative ionization mode dataset
edata <- lipid_neg_edata
fdata <- lipid_neg_fdata
emeta <- lipid_neg_emeta

## -----------------------------------------------------------------------------
head(edata)

## -----------------------------------------------------------------------------
head(fdata)

## -----------------------------------------------------------------------------
all(fdata$SampleID %in% names(edata)[-which(names(edata) == "Lipid")])
all(names(edata)[-which(names(edata) == "Lipid")] %in% fdata$SampleID)

## -----------------------------------------------------------------------------
head(emeta)

## -----------------------------------------------------------------------------
all(edata$Lipid %in% emeta$Lipid)

## -----------------------------------------------------------------------------
mylipid <- as.lipidData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = "Lipid", 
                        fdata_cname = "SampleID",
                        emeta_cname = "Lipid", 
                        data_scale = "abundance",
                        data_types = "Negative Ion",
                        check.names = FALSE)

## -----------------------------------------------------------------------------
class(mylipid)

summary(mylipid)

plot(edata_transform(mylipid, data_scale = "log2"))


