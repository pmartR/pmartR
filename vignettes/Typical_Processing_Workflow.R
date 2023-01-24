## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(warning=FALSE, fig.width = 8, fig.height = 6)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----logo, out.width = "100px", echo=FALSE------------------------------------
knitr::include_graphics("pmartR_logo_final.jpg")

knitr::opts_chunk$set(fig.width=8)
knitr::opts_chunk$set(fig.height=6)

library(pmartR)
library(pmartRdata)
library(ggplot2)
library(reshape2)

## ----data0, include=FALSE-----------------------------------------------------
# load example e_data, f_data, e_meta for the labeled peptide dataset
edata <- isobaric_edata
fdata <- isobaric_fdata
emeta <- isobaric_emeta

## -----------------------------------------------------------------------------
head(edata)

## -----------------------------------------------------------------------------
head(fdata)

## -----------------------------------------------------------------------------
head(emeta)

## ----data1--------------------------------------------------------------------
mypep <- as.isobaricpepData(e_data = edata,
                            f_data = fdata,
                            e_meta = emeta,
                            edata_cname = "Peptide",
                            fdata_cname = "SampleID",
                            emeta_cname = "Protein",
                            data_scale = "abundance")

# log2 transform the data
mypep <- edata_transform(omicsData = mypep, data_scale = "log2")


## ----data2--------------------------------------------------------------------
summary(mypep)

## ----data4--------------------------------------------------------------------
# Don't apply the normalization quite yet; can use summary() and plot() to view reference pool samples
mypep_refpools <- normalize_isobaric(omicsData = mypep, 
                                     exp_cname = "Plex",
                                     apply_norm = FALSE,
                                     refpool_cname = "Virus",
                                     refpool_notation = "Pool")

## -----------------------------------------------------------------------------
summary(mypep_refpools)                     
  
plot(mypep_refpools)

## -----------------------------------------------------------------------------
# Now apply the normalization; can use plot() to view the study samples after reference pool normalization
mypep <- normalize_isobaric(omicsData = mypep, 
                            exp_cname = "Plex",
                            apply_norm = TRUE,
                            refpool_cname = "Virus",
                            refpool_notation = "Pool")

## -----------------------------------------------------------------------------
plot(mypep)

## ----workflowFig, out.width = "600px", echo=FALSE, fig.cap = "Figure 1. Quality control and processing workflow in pmartR package."----
knitr::include_graphics("pmartR_graphic-final.jpg")

## ----format3------------------------------------------------------------------
mypep <- group_designation(omicsData = mypep, main_effects = "Virus")
summary(mypep)
plot(mypep, order_by = "Virus", color_by = "Virus")

## ----format4------------------------------------------------------------------
attributes(mypep)$group_DF

## ----proteomics_filter--------------------------------------------------------
myfilter <- proteomics_filter(mypep)
summary(myfilter, degen_peps = TRUE)
plot(myfilter, bw_theme=TRUE)

## ----proteomics_filter2-------------------------------------------------------
summary(myfilter, min_num_peps = 2)

# apply the filter, requiring proteins to have at least 2 peptides that map to them
mypep <- applyFilt(filter_object = myfilter, omicsData = mypep, min_num_peps = 2)

## ----molecule_filter----------------------------------------------------------
myfilter <- molecule_filter(mypep)
summary(myfilter, min_num = 3)
plot(myfilter, bw_them = TRUE)

## ----molecule_filter2---------------------------------------------------------
summary(myfilter, min_num = 2)

mypep <- applyFilt(filter_object = myfilter, mypep, min_num = 2)
summary(mypep)

## ----cv_filter----------------------------------------------------------------
myfilter <- cv_filter(mypep)
summary(myfilter, cv_threshold = 140)
plot(myfilter, cv_threshold = 140)

mypep <- applyFilt(filter_object = myfilter, omicsData = mypep, cv_threshold = 140)
summary(mypep)

## ----custom_filter------------------------------------------------------------
remove_prots <- mypep$e_meta$Protein[grep("XXX", mypep$e_meta$Protein)]
myfilter <- custom_filter(mypep, e_meta_remove = remove_prots)
summary(myfilter)

mypep <- applyFilt(filter_object = myfilter, omicsData = mypep)

## ----imdanova_filter1---------------------------------------------------------
myfilter <- imdanova_filter(mypep)

## ----imdanova_filter2---------------------------------------------------------
summary(myfilter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

## ----imdanova_filter3---------------------------------------------------------
summary(myfilter, min_nonmiss_anova = 2)

## ----imdanova_filter4---------------------------------------------------------
summary(myfilter, min_nonmiss_gtest = 3)

## -----------------------------------------------------------------------------
mypep <- applyFilt(filter_object = myfilter, omicsData = mypep, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

summary(mypep)

## ----rmd1---------------------------------------------------------------------
myfilter <- rmd_filter(omicsData = mypep, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness"))
plot(myfilter)

summary(myfilter, pvalue_threshold = 0.0001)
plot(myfilter, pvalue_threshold = 0.0001, bw_theme=TRUE)

## ----rmd2---------------------------------------------------------------------
# get vector of potential outliers
potential_outliers <- summary(myfilter, pvalue_threshold = 0.0001)$filtered_samples

# loop over potential outliers to generate plots
if(length(potential_outliers) > 0){
  for(i in 1:length(potential_outliers)){
    print(plot(myfilter, sampleID = potential_outliers[i]))
  }
}



## ----correlation heatmap------------------------------------------------------
mycor <- cor_result(omicsData = mypep)
plot(mycor, interactive = TRUE)

## ----pca----------------------------------------------------------------------
mypca <- dim_reduction(omicsData = mypep)

plot(mypca, interactive = TRUE)

## -----------------------------------------------------------------------------
mypep <- applyFilt(filter_object = myfilter, omicsData = mypep, pvalue_threshold = 0.0001)
summary(mypep)

## ----edata_summary------------------------------------------------------------
edata_summary(omicsData = mypep, by = "sample", groupvar = NULL)

#edata_summary(omicsData = norm_data, by = "molecule", groupvar = "Condition")
#edata_summary(omicsData = norm_data, by = "molecule", groupvar = NULL)


## ----missingval---------------------------------------------------------------
results <- missingval_result(omicsData = mypep)

plot(naRes_obj = results, omicsData = mypep, plot_type = "bar")
plot(naRes_obj = results, omicsData = mypep, plot_type = "scatter")

## ----spans, eval = FALSE------------------------------------------------------
#  # returns a data frame arranged by descending SPANS score
#  # not run due to long runtime; SPANS plot generated and saved to include in vignette
#  spans_result <- spans_procedure(mypep)
#  
#  plot(spans_result)

## ----spans plot, out.width = "600px", echo=FALSE------------------------------
knitr::include_graphics("SPANS_isobaricpep.png")

## -----------------------------------------------------------------------------
#readRDS("/vignettes/spans_result.RDS") 

mynorm <- normalize_global(omicsData = mypep, 
                           subset_fn = "rip", params = list(rip = 0.25), 
                           norm_fn = "median", 
                           apply_norm = FALSE)

# plot(mynorm) # currently errors out

## ----spans_normalize----------------------------------------------------------
# apply the normalization
mynorm <- normalize_global(omicsData = mypep, 
                           subset_fn = "rip", params = list(rip = 0.25), 
                           norm_fn = "median", 
                           apply_norm = FALSE)


## ----protein_quant------------------------------------------------------------
mypro <- protein_quant(pepData = mypep, 
                       method = "rrollup", 
                       combine_fn = "median")

plot(mypro, order_by = "Virus", color_by = "Virus")

## -----------------------------------------------------------------------------
mystats <- imd_anova(omicsData = mypro, test_method = "combined")

summary(mystats)

## -----------------------------------------------------------------------------
plot(mystats, plot_type = "volcano")
plot(mystats, plot_type = "bar")
plot(mystats, plot_type = "gheatmap") # exclude this one?
plot(mystats, plot_type = "fcheatmap")

