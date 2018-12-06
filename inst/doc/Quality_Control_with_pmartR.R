## ----global_options, include = FALSE-------------------------------------
knitr::opts_chunk$set(warning=FALSE, fig.width = 8, fig.height = 6)


## ----logo, out.width = "100px", echo=FALSE-------------------------------
knitr::include_graphics("pmartR_logo_final.jpg")

library(pmartR)
library(pmartRdata)
library(ggplot2)
library(reshape2)

## ----data0, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=8)
knitr::opts_chunk$set(fig.height=6)

data("pep_edata")
data("pep_fdata")

## ----data1---------------------------------------------------------------
data("pep_edata")
dim(pep_edata)
pep_edata[1:6,1:5]
data("pep_fdata")
dim(pep_fdata)
head(pep_fdata)
data("pep_emeta")
dim(pep_emeta)
head(pep_emeta)

## ----data2---------------------------------------------------------------
mypepData <- as.pepData(e_data = pep_edata, f_data = pep_fdata, e_meta = pep_emeta, edata_cname = "Mass_Tag_ID", emeta_cname = "Protein", fdata_cname = "SampleID", data_scale = "abundance")
class(mypepData)
summary(mypepData)

## ----data3---------------------------------------------------------------
data("pep_object")
class(pep_object)

rm(pep_object)

## ----data4---------------------------------------------------------------
plot(mypepData)

## ----checknamesFALSE-----------------------------------------------------
edata <- read.csv(system.file("extdata", "metab_edata_sample_names.csv", package="pmartRdata"), header=TRUE)
names(edata)

fdata <- read.csv(system.file("extdata", "metab_fdata_sample_names.csv", package="pmartRdata"), header=TRUE)

# Do the sample names match each other? #
all(names(edata)[-1] == fdata$SampleID)


## ----checknamesTRUE------------------------------------------------------
edata <- read.csv(system.file("extdata", "metab_edata_sample_names.csv", package="pmartRdata"), header=TRUE, check.names=FALSE)
names(edata)

# Do the sample names match each other? #
all(names(edata)[-1] == fdata$SampleID)

## ----workflowFig, out.width = "600px", echo=FALSE, fig.cap = "Figure 1. Quality control and processing workflow in pmartR package."----
knitr::include_graphics("pmartR_graphic-final.jpg")

## ----format1-------------------------------------------------------------
mypepData <- edata_replace(mypepData, x = 0, y = NA)

## ----format2-------------------------------------------------------------
mypepData <- edata_transform(mypepData, data_scale = "log10")
attributes(mypepData)$data_info$data_scale

mypepData <- edata_transform(mypepData, data_scale = "abundance")
attributes(mypepData)$data_info$data_scale

mypepData <- edata_transform(mypepData, data_scale = "log2")
attributes(mypepData)$data_info$data_scale

plot(mypepData)

## ----format3-------------------------------------------------------------
mypepData <- group_designation(mypepData, main_effects = "Condition", covariates = NULL)

## ----format4-------------------------------------------------------------
attributes(mypepData)$group_DF

plot(mypepData, color_by = "Condition", bw_theme=TRUE)

## ----proteomics_filter---------------------------------------------------
myfilter <- proteomics_filter(mypepData)
summary(myfilter)
plot(myfilter, bw_theme=TRUE)

## ----molecule_filter-----------------------------------------------------
myfilter <- molecule_filter(mypepData)
plot(myfilter, bw_them = TRUE)
summary(myfilter, min_num = 3)

## ----molecule_filter2----------------------------------------------------
summary(myfilter, min_num = 2)
#plot(myfilter, min_num = 2)

mypepData <- applyFilt(filter_object = myfilter, mypepData, min_num = 2)
summary(mypepData)

## ----cv_filter-----------------------------------------------------------
myfilter <- cv_filter(mypepData)
summary(myfilter, cv_threshold = 150)
plot(myfilter, cv_threshold = 150, title_size = 15)

mypepData <- applyFilt(filter_object = myfilter, mypepData, cv_threshold = 150)
summary(mypepData)


## ----custom_filter-------------------------------------------------------
myfilter <- custom_filter(mypepData, e_data_remove = 1047, e_meta_remove = NULL, f_data_remove = "Infection1")
summary(myfilter)

mypepData_temp <- applyFilt(filter_object = myfilter, mypepData)
summary(mypepData_temp)


rm(mypepData_temp)

## ----custom_filter2------------------------------------------------------
myfilter2<- custom_filter(mypepData, e_data_keep = NULL, e_meta_keep = NULL, f_data_keep = c("Infection1", "Infection2", "Infection3", "Infection4", "Infection5"))
summary(myfilter2)

mypepData_temp2<- applyFilt(filter_object = myfilter2, mypepData)
summary(mypepData_temp2)

rm(mypepData_temp2)

## ----imdanova_filter1----------------------------------------------------
myfilter <- imdanova_filter(mypepData)

## ----imdanova_filter2----------------------------------------------------
summary(myfilter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

#plot(myfilter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

## ----imdanova_filter3----------------------------------------------------
summary(myfilter, min_nonmiss_anova = 2)
#plot(myfilter, min_nonmiss_anova = 2)

## ----imdanova_filter4----------------------------------------------------
summary(myfilter, min_nonmiss_gtest = 3)
#plot(myfilter, min_nonmiss_gtest = 3)

## ------------------------------------------------------------------------
mypepData <- applyFilt(filter_object = myfilter, mypepData, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

summary(mypepData)

## ----rmd1----------------------------------------------------------------
myfilter <- rmd_filter(mypepData, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"))
plot(myfilter, bw_theme=TRUE)

summary(myfilter, pvalue_threshold = 0.001)
plot(myfilter, pvalue_threshold = 0.001, bw_theme=TRUE)



## ----rmd2----------------------------------------------------------------
plot(myfilter, sampleID = "Infection5", bw_theme=TRUE)
plot(myfilter, sampleID = "Infection8", bw_theme=TRUE)

## ----norm2---------------------------------------------------------------

# Normalize using all peptides and median centering #
norm_object <- normalize_global(omicsData = mypepData, subset_fn = "all", norm_fn = "median", apply_norm = FALSE, backtransform = TRUE)
norm_data <- normalize_global(omicsData = mypepData, subset_fn = "all", norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
plot(norm_data, title_plot="Normalized Data: Median Centering Using all Peptides", title_size=12)

# Normalize using RIP 0.2 and median centering #
norm_object <- normalize_global(omicsData = mypepData, subset_fn = "rip", params=list(rip=0.2), norm_fn = "median", apply_norm = FALSE, backtransform = TRUE)
norm_data <- normalize_global(omicsData = mypepData, subset_fn = "rip", params=list(rip=0.2), norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
plot(norm_data, title_plot="Normalized Data: Median Centering Using RIP 0.2 Peptides", title_size=12)


## ----spans---------------------------------------------------------------
# returns a dataframe arranged by descending SPANS score
spans_result <- spans_procedure(mypepData)


## ----spans_plot_summarize------------------------------------------------
summary(spans_result)
plot(spans_result)

## ----spans_normalize, eval = FALSE---------------------------------------
#  # a list of the parameters for any normalization procedure with the best SPANS score
#  best_params <- get_spans_params(spans_result)
#  
#  # there are a few ties, all using ppp_rip, well select the method that uses median normalization and parameters ppp = 0.1 and rip = 0.1
#  subset_fn = best_params[[1]]$subset_fn
#  norm_fn = best_params[[1]]$norm_fn
#  params = best_params[[1]]$params
#  
#  # we now pass the extracted subset function, normalization function, and parameters from SPANS to normalize_global()
#  norm_object <- normalize_global(omicsData = mypepData, subset_fn = subset_fn, norm_fn = norm_fn, params = params)

## ----edata_summary-------------------------------------------------------
edata_summary(omicsData = norm_data, by = "sample", groupvar = NULL)

#edata_summary(omicsData = norm_data, by = "molecule", groupvar = "Condition")
#edata_summary(omicsData = norm_data, by = "molecule", groupvar = NULL)


## ----sppPCA--------------------------------------------------------------
pca_results <- dim_reduction(norm_data, k = 2)

summary(pca_results) 
plot(pca_results, bw_theme=TRUE)

## ----corr----------------------------------------------------------------
correlation_results <- cor_result(norm_data)

summary(correlation_results)
plot(correlation_results)



## ----missingval----------------------------------------------------------

results <- missingval_result(mypepData)


plot(results, type = "bySample")
plot(results, type = "byMolecule")


## ----missingval_plots----------------------------------------------------

missingval_scatterplot(mypepData, palette = "Set1")

missingval_heatmap(mypepData)


## ----protein_quant, eval = FALSE-----------------------------------------
#  print(mypepData)
#  
#  mypepdata_median<- protein_quant(pepData = mypepData, method = "rollup", combine_fn = "median")
#  print(mypepdata_median)
#  rm(mypepdata_median)
#  
#  mypepdata_rrollup<- protein_quant(pepData = mypepData, method = "rrollup")
#  print(mypepdata_rrollup)
#  rm(mypepdata_rrollup)
#  

