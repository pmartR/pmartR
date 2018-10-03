## ----logo, out.width = "100px", echo=FALSE-------------------------------
knitr::include_graphics("pmartR_logo_final.jpg")

library(pmartR)
library(pmartRdata)
library(ggplot2)

## ----data0, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=5)
knitr::opts_chunk$set(fig.height=5)

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

plot(mypepData)
plot(mypepData, color_by = "Condition", bw_theme=TRUE)
plot(mypepData, color_by = "Condition", title_size=10, legend_position="bottom")
plot(mypepData, facet_by = "Condition", color_by = "Condition")

## ----proteomics_filter---------------------------------------------------
myfilter <- proteomics_filter(mypepData)
summary(myfilter)
plot(myfilter)
plot(myfilter, bw_theme=TRUE)
plot(myfilter, title_size=10, x_lab_size=9, y_lab_size=9)

## ----molecule_filter-----------------------------------------------------
myfilter <- molecule_filter(mypepData)
plot(myfilter)
summary(myfilter, min_num = 3)
plot(myfilter, min_num = 3)
plot(myfilter, min_num = 3, bw_theme = TRUE)
plot(myfilter, min_num = 3, title_size = 10, x_lab_size = 6, y_lab_size = 8)

## ----molecule_filter2----------------------------------------------------
summary(myfilter, min_num = 2)
plot(myfilter, min_num = 2)

mypepData <- applyFilt(filter_object = myfilter, mypepData, min_num = 2)
summary(mypepData)

## ----cv_filter-----------------------------------------------------------
myfilter <- cv_filter(mypepData)
plot(myfilter)
plot(myfilter, bw_theme = TRUE)
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

plot(myfilter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3, x_lab_size = 9, y_lab_size = 9)

## ----imdanova_filter3----------------------------------------------------
summary(myfilter, min_nonmiss_anova = 2)
plot(myfilter, min_nonmiss_anova = 2)

## ----imdanova_filter4----------------------------------------------------
summary(myfilter, min_nonmiss_gtest = 3)
plot(myfilter, min_nonmiss_gtest = 3)

## ------------------------------------------------------------------------
mypepData <- applyFilt(filter_object = myfilter, mypepData, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

summary(mypepData)

## ----rmd1----------------------------------------------------------------
myfilter <- rmd_filter(mypepData, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"))

plot(myfilter)
plot(myfilter, title_size=7, bw_theme=TRUE)

summary(myfilter, pvalue_threshold = 0.001)
plot(myfilter, pvalue_threshold = 0.001, bw_theme=TRUE)
plot(myfilter, pvalue_threshold = 0.001, bw_theme=TRUE, legend_position="bottom")
plot(myfilter, pvalue_threshold = 0.001, legend_position="top")
plot(myfilter, pvalue_threshold = 0.001, bw_theme=TRUE, legend_position="left")
plot(myfilter, pvalue_threshold = 0.001, x_lab_size = 8, y_lab_size = 8)
plot(myfilter, pvalue_threshold = 0.001, x_lab_size = 8, y_lab_size = 8, point_size = 6)
plot(myfilter, pvalue_threshold = 0.001, x_lab_size = 8, y_lab_size = 8, point_size = 7)


## ----rmd2----------------------------------------------------------------
plot(myfilter, sampleID = "Infection5")
plot(myfilter, sampleID = "Infection5", bw_theme=TRUE)
plot(myfilter, sampleID = "Infection8", title_size = 10, y_lab_size=8)

## ----norm2---------------------------------------------------------------

# Normalize using all peptides and median centering #
norm_object <- normalize_global(omicsData = mypepData, subset_fn = "all", norm_fn = "median", apply_norm = FALSE, backtransform = TRUE)
norm_data <- normalize_global(omicsData = mypepData, subset_fn = "all", norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
plot(norm_data, title_plot="Normalized Data: Median Centering Using all Peptides", title_size=12)

# Normalize using RIP 0.2 and median centering #
norm_object <- normalize_global(omicsData = mypepData, subset_fn = "rip", params=list(rip=0.2), norm_fn = "median", apply_norm = FALSE, backtransform = TRUE)
norm_data <- normalize_global(omicsData = mypepData, subset_fn = "rip", params=list(rip=0.2), norm_fn = "median", apply_norm = TRUE, backtransform = TRUE)
plot(norm_data, title_plot="Normalized Data: Median Centering Using RIP 0.2 Peptides", title_size=12)


## ----norm1---------------------------------------------------------------
# data_norm <- normalize_data(mypepData, subset_fn="ppp.5", norm_fun="median", backtransform=TRUE)




## ----edata_summary-------------------------------------------------------
edata_summary(omicsData = norm_data, by = "sample", groupvar = NULL)
#edata_summary(omicsData = norm_data, by = "sample", groupvar = "Condition")

#edata_summary(omicsData = norm_data, by = "molecule", groupvar = NULL)


## ----sppPCA--------------------------------------------------------------
pca_results <- dim_reduction(norm_data, k = 2)

summary(pca_results) 
plot(pca_results)
plot(pca_results, title_size = 20, bw_theme=TRUE, legend_position="left")

## ----corr----------------------------------------------------------------
correlation_results <- cor_result(norm_data)

summary(correlation_results)
plot(correlation_results)
plot(correlation_results, title_size = 20)
plot(correlation_results, x_lab=FALSE, title_plot="My Correlation Heatmap")
plot(correlation_results, x_lab=FALSE, y_lab=FALSE)
plot(correlation_results, colorbar_lim=c(0,1))


## ----missingval----------------------------------------------------------

results <- missingval_result(mypepData)


plot(results, type = "bySample")
plot(results, type = "bySample", title_plot="Hey look at this!")
plot(results, type = "byMolecule", x_lab_angle = 0)


## ----missingval_plots----------------------------------------------------

missingval_scatterplot(mypepData, palette = "Set1")
missingval_scatterplot(mypepData, title_plot="Hey look at this!", palette = "Set2")

missingval_heatmap(mypepData)
missingval_heatmap(mypepData, title_size = 20, title_plot = "Hey look at this!")


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

