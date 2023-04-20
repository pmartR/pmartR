## ----global_options, include = FALSE------------------------------------------
knitr::opts_chunk$set(warning = FALSE, fig.width = 8, fig.height = 6, echo = TRUE, include = TRUE, eval = TRUE)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----logo, out.width = "100px", echo=FALSE------------------------------------
knitr::include_graphics("pmartR_logo_final.jpg")

library(pmartR)
library(pmartRdata)
library(ggplot2)
library(reshape2)

## ----molecule filter----------------------------------------------------------
mypep <- pep_object

# create a molecule filter object
mymolfilt <- molecule_filter(omicsData = mypep)
class(mymolfilt)
summary(mymolfilt)

# summary of molecule filter with a specific threshold
summary(mymolfilt, min_num = 2)

mypep_molfilt <- applyFilt(mymolfilt, omicsData = mypep, min_num = 2)

# now our data object contains the reduced number of biomolecules
summary(mypep_molfilt)

## ----molecule filter groups---------------------------------------------------
# define groups based on a "main effect"
mypep_groups <- group_designation(omicsData = mypep, main_effects = "Phenotype")

# create a molecule filter object
mymolfilt_groups <- molecule_filter(omicsData = mypep_groups, use_groups = TRUE)
class(mymolfilt_groups)
summary(mymolfilt_groups)

# summary of molecule filter with a specific threshold
summary(mymolfilt_groups, min_num = 2)

mypep_molfilt_groups <- applyFilt(mymolfilt_groups, omicsData = mypep_groups, min_num = 2)

# now our data object contains the reduced number of biomolecules
summary(mypep_molfilt_groups)

## ----molecule filter batches--------------------------------------------------
# define groups based on a "main effect" and "batch id" -- here we pretend that "SecondPhenotype" is describing the batches of samples
mypep_batches <- group_designation(omicsData = mypep, main_effects = "Phenotype", batch_id = "SecondPhenotype")

# create a molecule filter object using the batch information -- we can use both groups and batch or just batch
mymolfilt_batches <- molecule_filter(omicsData = mypep_batches, use_groups = TRUE, use_batch = TRUE)
mymolfilt_batches <- molecule_filter(omicsData = mypep_batches, use_batch = TRUE)
class(mymolfilt_batches)
summary(mymolfilt_batches)

# summary of molecule filter with a specific threshold
summary(mymolfilt_batches, min_num = 2)

mypep_molfilt_batches <- applyFilt(mymolfilt_batches, omicsData = mypep_batches, min_num = 2)

# now our data object contains the reduced number of biomolecules
summary(mypep_molfilt_batches)

## ----cv filter----------------------------------------------------------------
# use the example peptide data with Phenotype as main effect
mycvfilt <- cv_filter(omicsData = mypep_groups)
plot(mycvfilt, cv_threshold = 97)

# we get the same graph if we log2 transform our data
mypep_groups_log2 <- edata_transform(mypep_groups, data_scale = "log2")
mycvfilt_log2 <- cv_filter(omicsData = mypep_groups_log2)
plot(mycvfilt_log2, cv_threshold = 97)

summary(mycvfilt_log2, cv_threshold = 97)

# apply the filter
mypep_cvfilt <- applyFilt(filter_object = mycvfilt_log2, omicsData = mypep_groups_log2, cv_threshold = 97)
summary(mypep_cvfilt)

## -----------------------------------------------------------------------------
myproteomicsfilt <- proteomics_filter(omicsData = mypep_groups_log2)
summary(myproteomicsfilt, degen_peps = TRUE)

# This results in an error message since there are no redundant peptides
# mypep_redundant <- applyFilt(filter_object = myredundantfilt, omicsData = mypep_groups_log2)

## -----------------------------------------------------------------------------
# use the same filter object as above
# setting min_num_peps = 2 removes any proteins with just a single observed peptide mapping to them; the associated peptides are also removed
summary(myproteomicsfilt, min_num_peps = 2)

mypep_proteomicsfilt <- applyFilt(filter_object = myproteomicsfilt, omicsData = mypep_groups_log2, min_num_peps = 2)

## ----combined test------------------------------------------------------------
myimdanovafilt <- imdanova_filter(omicsData = mypep_groups_log2)
summary(myimdanovafilt, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

mypep_imdanovafilt <- applyFilt(filter_object = myimdanovafilt, omicsData = mypep_groups_log2, min_nonmiss_anova = 2, min_nonmiss_gtest = 3)

## ----anova only---------------------------------------------------------------
# start with the same filter object as above
summary(myimdanovafilt, min_nonmiss_anova = 2)

mypep_anovafilt <- applyFilt(filter_object = myimdanovafilt, omicsData = mypep_groups_log2, min_nonmiss_anova = 2)

## ----gtest only---------------------------------------------------------------
# start with the same filter object as above
summary(myimdanovafilt, min_nonmiss_gtest = 3)

mypep_imdfilt <- applyFilt(filter_object = myimdanovafilt, omicsData = mypep_groups_log2, min_nonmiss_gtest = 3)

## ----remove from e_data-------------------------------------------------------
remove_peps <- mypep$e_data$Peptide[1:10]

mycustomfilt <- custom_filter(omicsData = mypep_groups_log2, e_data_remove = remove_peps)
summary(mycustomfilt)

mypep_custom <- applyFilt(filter_object = mycustomfilt, omicsData = mypep_groups_log2)

## ----remove from e_meta-------------------------------------------------------
remove_prots <- mypep$e_meta$RazorProtein[1:10]

mycustomfilt <- custom_filter(omicsData = mypep_groups_log2, e_meta_remove = remove_prots)
summary(mycustomfilt)

mypep_custom <- applyFilt(filter_object = mycustomfilt, omicsData = mypep_groups_log2)

## ----create rmd and plot------------------------------------------------------
myrmdfilt <- rmd_filter(mypep_groups_log2)
plot(myrmdfilt, pvalue_threshold = 0.0001)

## ----rmd plot individ outliers------------------------------------------------
plot(myrmdfilt, sampleID = "Sample_35_Phenotype2_A")

## -----------------------------------------------------------------------------
mycustomfilt <- custom_filter(omicsData = mypep_groups_log2, f_data_remove = "Sample_35_Phenotype2_A")

## -----------------------------------------------------------------------------
myseq <- rnaseq_object
summary(myseq)

## -----------------------------------------------------------------------------
mytotcountfilt <- total_count_filter(omicsData = myseq)
summary(mytotcountfilt, min_count = 10)

myseq_totcount <- applyFilt(filter_object = mytotcountfilt, omicsData = myseq, min_count = 10)
summary(myseq_totcount)

## -----------------------------------------------------------------------------
# create RNA filter object
myrnafilt <- RNA_filter(omicsData = myseq)

# filter by library size
plot(myrnafilt, plot_type = "library")
plot(myrnafilt, plot_type = "library", size_library = 10000)
summary(myrnafilt, size_library = 10000)

myseq_librarysize <- applyFilt(filter_object = myrnafilt, omicsData = myseq, size_library = 10000)
summary(myseq_librarysize)

# filter based on number or proportion of non-zero counts
plot(myrnafilt, plot_type = "biomolecule")
summary(myrnafilt, min_nonzero = 5000)
summary(myrnafilt, min_nonzero = 0.2)

myseq_nonzero <- applyFilt(filter_object = myrnafilt, omicsData = myseq, min_nonzero = 0.2)
summary(myseq_nonzero)

