## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(pmartR)
library(pmartRdata)

## ----data0, include=FALSE-----------------------------------------------------
# load example e_data, f_data, e_meta for the RNAseq dataset
edata <- rnaseq_edata
fdata <- rnaseq_fdata
emeta <- rnaseq_emeta

## -----------------------------------------------------------------------------
head(edata)

## -----------------------------------------------------------------------------
head(fdata)

## -----------------------------------------------------------------------------
head(emeta)

## -----------------------------------------------------------------------------
myseq <- as.seqData(e_data = rnaseq_edata,
                        e_meta = rnaseq_emeta,
                        f_data = rnaseq_fdata,
                        edata_cname = "Transcript",
                        fdata_cname = "SampleName",
                        emeta_cname = "Transcript")

## ----data2--------------------------------------------------------------------
summary(myseq)

## -----------------------------------------------------------------------------
plot(myseq, transformation = "lcpm")

## ----format3------------------------------------------------------------------
myseq <- group_designation(omicsData = myseq, main_effects = "Virus")
summary(myseq)
plot(myseq, order_by = "Virus", color_by = "Virus", transformation = "lcpm")

## ----format4------------------------------------------------------------------
attributes(myseq)$group_DF

## -----------------------------------------------------------------------------
mytotcountfilt <- total_count_filter(omicsData = myseq)
summary(mytotcountfilt, min_count = 10)
plot(mytotcountfilt)
plot(mytotcountfilt, min_count = 10)

## -----------------------------------------------------------------------------
myseq <- applyFilt(filter_object = mytotcountfilt, omicsData = myseq, min_count = 10)

## -----------------------------------------------------------------------------
# create RNA filter object
myrnafilt <- RNA_filter(omicsData = myseq)

# filter by library size
plot(myrnafilt, plot_type = "library")
plot(myrnafilt, plot_type = "library", size_library = 1000000)


## -----------------------------------------------------------------------------
# filter based on number or proportion of non-zero counts
plot(myrnafilt, plot_type = "biomolecule", min_nonzero = 5000)

## -----------------------------------------------------------------------------
mycor <- cor_result(omicsData = myseq)
plot(mycor, interactive = TRUE)

## -----------------------------------------------------------------------------
mypca <- dim_reduction(omicsData = myseq, k=2)
plot(mypca, interactive = TRUE)

mypca_df <- data.frame(SampleID = mypca$SampleID, PC1 = mypca$PC1, PC2 = mypca$PC2)

## -----------------------------------------------------------------------------
dispersion_est(omicsData = myseq, method = "edgeR")

## -----------------------------------------------------------------------------
stats_edger <- diffexp_seq(omicsData = myseq, method = "edgeR")

## ----warning=FALSE------------------------------------------------------------
plot(stats_edger, plot_type = "volcano")
plot(stats_edger, plot_type = "bar")
plot(stats_edger, plot_type = "ma") 

## -----------------------------------------------------------------------------
dispersion_est(omicsData = myseq, method = "DESeq2")

## -----------------------------------------------------------------------------
stats_deseq <- diffexp_seq(omicsData = myseq, method = "DESeq2")

## ----warning=FALSE------------------------------------------------------------
plot(stats_deseq, plot_type = "volcano")
plot(stats_deseq, plot_type = "bar")
plot(stats_deseq, plot_type = "ma") 

## -----------------------------------------------------------------------------
dispersion_est(omicsData = myseq, method = "voom")

## -----------------------------------------------------------------------------
stats_voom <- diffexp_seq(omicsData = myseq, method = "voom")

## ----warning=FALSE------------------------------------------------------------
plot(stats_voom, plot_type = "volcano")
plot(stats_voom, plot_type = "bar")
plot(stats_voom, plot_type = "ma") 

