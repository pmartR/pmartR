pkgname <- "pmartR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('pmartR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MSnSet2pepData")
### * MSnSet2pepData

flush(stderr()); flush(stdout())

### Name: MSnSet2pepData
### Title: Convert Data of class MSnSet to pmartR pepData Class
### Aliases: MSnSet2pepData

### ** Examples

## Not run: 
##D library(pmartR)
##D library(MSnbase)
##D data("msnset")
##D result = MSnSet2pepData(msnset,
##D                         data_scale = "log2",
##D                         edata_cname = "UniqueID",
##D                         fdata_cname = "SampleID",
##D                         emeta_cname = "UniqueID",
##D                         check.names = TRUE)
## End(Not run)




cleanEx()
nameEx("RNA_filter")
### * RNA_filter

flush(stderr()); flush(stdout())

### Name: RNA_filter
### Title: RNA filter object
### Aliases: RNA_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("seq_object")
##D to_filter <- RNA_filter(omicsData = seq_object)
##D plot(to_filter)
##D summary(to_filter, size_library = 10000)
##D summary(to_filter, min_nonzero = 5000)
##D summary(to_filter, min_nonzero = .2)
## End(Not run)




cleanEx()
nameEx("all_subset")
### * all_subset

flush(stderr()); flush(stdout())

### Name: all_subset
### Title: Identifies the 'subset' of All Features
### Aliases: all_subset

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_edata)
##D keep <- none_subset(e_data = pep_edata, edata_id = "Mass_Tag_ID")
## End(Not run)




cleanEx()
nameEx("analysis_log")
### * analysis_log

flush(stderr()); flush(stdout())

### Name: analysis_log
### Title: Creates R markdown document report
### Aliases: analysis_log

### ** Examples

dontrun{
library(pmartRdata)
data(metab_object)
result = analysis_log(metab_object)
}




cleanEx()
nameEx("anova_filter")
### * anova_filter

flush(stderr()); flush(stdout())

### Name: anova_filter
### Title: Identifies biomolecules to be filtered in preparation for
###   IMD-ANOVA.
### Aliases: anova_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D 
##D pep_object2 <- group_designation(pep_object, main_effects = "Condition")
##D 
##D nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
##D 
##D to_filter <- anova_filter(nonmiss_per_group = nonmissing_result,
##D                           min_nonmiss_anova = 2)
## End(Not run)




cleanEx()
nameEx("applyFilt")
### * applyFilt

flush(stderr()); flush(stdout())

### Name: applyFilt
### Title: Apply a S3 filter object to a pmartR S3 object
### Aliases: applyFilt applyFilt.moleculeFilt applyFilt.totalCountFilt
###   applyFilt.RNAFilt applyFilt.cvFilt applyFilt.rmdFilt
###   applyFilt.proteomicsFilt applyFilt.imdanovaFilt applyFilt.customFilt

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D to_filter <- molecule_filter(omicsData = pep_object)
##D pep_object2 <- applyFilt(filter_object = to_filter,
##D                          omicsData = pep_object, min_num = 2)
##D print(str(attributes(pep_object2)$filters))
##D pep_object2 <- group_designation(pep_object2, main_effects = "Condition")
##D to_filter2 <- imdanova_filter(omicsData = pep_object2)
##D pep_object3 <- applyFilt(filter_object = to_filter2,
##D                          omicsData = pep_object2,
##D                          min_nonmiss_anova = 3)
##D print(str(attributes(pep_object3)$filters))
## End(Not run)




cleanEx()
nameEx("as.isobaricpepData")
### * as.isobaricpepData

flush(stderr()); flush(stdout())

### Name: as.isobaricpepData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.isobaricpepData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("isobaric_edata")
##D data("isobaric_fdata")
##D data("isobaric_emeta")
##D mypepData <- as.isobaricpepData(e_data = isobaric_edata,
##D                                 e_meta = isobaric_emeta,
##D                                 f_data = isobaric_fdata,
##D                                 edata_cname = "Peptide",
##D                                 fdata_cname = "Sample",
##D                                 emeta_cname = "Protein")
## End(Not run)




cleanEx()
nameEx("as.lipidData")
### * as.lipidData

flush(stderr()); flush(stdout())

### Name: as.lipidData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.lipidData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("lipid_edata")
##D data("lipid_fdata")
##D mylipidData <- as.lipidData(e_data = lipid_edata,
##D                             f_data = lipid_fdata,
##D                             edata_cname = "LipidCommonName",
##D                             fdata_cname = "Sample_Name")
## End(Not run)




cleanEx()
nameEx("as.metabData")
### * as.metabData

flush(stderr()); flush(stdout())

### Name: as.metabData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.metabData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("metab_edata")
##D data("metab_fdata")
##D mymetabData <- as.metabData(e_data = metab_edata,
##D                             f_data = metab_fdata,
##D                             edata_cname = "Metabolite",
##D                             fdata_cname = "SampleID")
## End(Not run)




cleanEx()
nameEx("as.multiData")
### * as.multiData

flush(stderr()); flush(stdout())

### Name: as.multiData
### Title: Create a 'multiData' object from multiple omicsData objects
### Aliases: as.multiData

### ** Examples


## Not run: 
##D library(pmartRdata)
##D library(pmartR)
##D 
##D # Combine lipid and protein object into multidata, both must be log2 + normalized.
##D mylipid_object <- edata_transform(lipid_object, "log2")
##D mylipid_object <- normalize_global(mylipid_object, "all", "median", apply_norm = TRUE)
##D 
##D # Combine without specifically supplying f_meta, either directly, or as one
##D # of the f_datas in any object.
##D mymultidata <- as.multiData(pro_object, mylipid_object, auto_fmeta = TRUE)
##D 
##D # Manually supply an f_meta
##D f_meta <- data.frame(
##D "Proteins" = c(paste0("Mock", 1:3), paste0("Infection", c(1:7)), NA,  "Infection9"),
##D "Lipids" = c(paste0("Mock", 1:3), paste0("Infection", c(1:4)), NA, paste0("Infection", c(6:9))),
##D "Metabolites" = c(paste0("Mock", 1:3), paste0("Infection", c(1:9))),
##D "Condition" = c(rep("A", 3), rep("B", 9))
##D )
##D 
##D mymetab_object <- edata_transform(metab_object, "log2")
##D mymetab_object <- normalize_global(mymetab_object, "all", "median", apply_norm = TRUE)
##D 
##D as.multiData(mylipid_object, pro_object, mymetab_object, f_meta = f_meta)
##D # remove samples that are not common across all data.
##D as.multiData(mylipid_object, pro_object, mymetab_object, f_meta = f_meta, sample_intersect = TRUE)
## End(Not run)




cleanEx()
nameEx("as.nmrData")
### * as.nmrData

flush(stderr()); flush(stdout())

### Name: as.nmrData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.nmrData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("nmr_edata_identified")
##D data("nmr_fdata_identified")
##D mynmrData <- as.nmrData(e_data = nmr_edata_identified,
##D                         f_data = nmr_fdata_identified,
##D                         edata_cname = "Metabolite",
##D                         fdata_cname = "SampleID",
##D                         data_type = "identified")
## End(Not run)




cleanEx()
nameEx("as.pepData")
### * as.pepData

flush(stderr()); flush(stdout())

### Name: as.pepData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.pepData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_edata")
##D data("pep_fdata")
##D data("pep_emeta")
##D mypepData <- as.pepData(e_data = pep_edata,
##D                         e_meta = pep_emeta,
##D                         f_data = pep_fdata,
##D                         edata_cname = "Mass_Tag_ID",
##D                         fdata_cname = "SampleID",
##D                         emeta_cname = "Mass_Tag_ID")
## End(Not run)




cleanEx()
nameEx("as.proData")
### * as.proData

flush(stderr()); flush(stdout())

### Name: as.proData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.proData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pro_edata")
##D data("pro_fdata")
##D myproData <- as.proData(e_data = pro_edata,
##D                         f_data = pro_fdata,
##D                         edata_cname = "Reference",
##D                         fdata_cname = "SampleID",
##D                         is_normalized = TRUE)
## End(Not run)




cleanEx()
nameEx("as.seqData")
### * as.seqData

flush(stderr()); flush(stdout())

### Name: as.seqData
### Title: Convert Data to Appropriate pmartR Class
### Aliases: as.seqData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("seq_edata")
##D data("seq_fdata")
##D data("seq_emeta")
##D mypepData <- as.seqData(e_data = seq_edata,
##D                         e_meta = seq_emeta,
##D                         f_data = seq_fdata,
##D                         edata_cname = "Seq_Tag_ID",
##D                         fdata_cname = "SampleID",
##D                         emeta_cname = "Seq_Tag_ID")
## End(Not run)




cleanEx()
nameEx("as.trelliData")
### * as.trelliData

flush(stderr()); flush(stdout())

### Name: as.trelliData
### Title: Generate an object from omicsData and/or statRes objects to pass
###   to trelliscope building functions
### Aliases: as.trelliData

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D 
##D # Generate an example e_meta file for lipid data 
##D lipid_emeta <- data.frame("LipidCommonName" = lipid_edata_pos$LipidCommonName, 
##D   "LipidFamily" = lipid_edata_pos$LipidCommonName %>% as.character() %>% 
##D     strsplit("(", fixed = TRUE) %>% lapply(function(el) {el[1]}) %>% unlist())
##D     
##D # Extend fdata to have two infection groups 
##D lipid_fdata_pos2 <- data.frame("Sample_Name" = lipid_fdata_pos$Sample_Name, 
##D   "Condition" = c(rep("Mock", 3), rep("Infection_A", 4), rep("Infection_B", 4)),
##D   "Weight" = runif(11))
##D 
##D # Build lipid data object
##D lipids <- as.lipidData(e_data = lipid_edata_pos, f_data = lipid_fdata_pos2, e_meta = lipid_emeta,
##D   edata_cname = "LipidCommonName", fdata_cname = "Sample_Name", 
##D   emeta_cname = "LipidCommonName")
##D 
##D # Transform the data
##D omicsData <- edata_transform(omicsData = lipids, data_scale = "log2")
##D 
##D # Group the data by condition
##D omicsData <- group_designation(omicsData = omicsData, main_effects = c("Condition"))
##D 
##D # Apply the IMD ANOVA filter
##D imdanova_Filt <- imdanova_filter(omicsData = omicsData)
##D omicsData <- applyFilt(filter_object = imdanova_Filt, omicsData = omicsData, min_nonmiss_anova=2)
##D 
##D # Normalize my pepData
##D omicsData <- normalize_global(omicsData, "subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE, "backtransform" = TRUE)
##D 
##D # Implement the IMD ANOVA method and compute all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
##D statRes <- imd_anova(omicsData = omicsData, test_method = 'combined')
##D 
##D # Generate the trelliData object 
##D trelliData2 <- as.trelliData(omicsData = omicsData)
##D trelliData3 <- as.trelliData(statRes = statRes)   
##D trelliData4 <- as.trelliData(omicsData = omicsData, statRes = statRes)
## End(Not run)




cleanEx()
nameEx("as.trelliData.edata")
### * as.trelliData.edata

flush(stderr()); flush(stdout())

### Name: as.trelliData.edata
### Title: Generate an object from edata to pass to trelliscope building
###   functions
### Aliases: as.trelliData.edata

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("lipid_edata_pos")
##D 
##D # Simple example 
##D trelliData <- as.trelliData.edata(e_data = lipid_edata_pos,
##D                                   edata_cname = "LipidCommonName",
##D                                   omics_type = "lipidData")
##D                                   
## End(Not run)   
          



cleanEx()
nameEx("bpquant")
### * bpquant

flush(stderr()); flush(stdout())

### Name: bpquant
### Title: Runs BP-Quant
### Aliases: bpquant

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D 
##D mypepData <- group_designation(omicsData = pep_object,
##D                                main_effects = c("Condition"))
##D mypepData = edata_transform(mypepData, "log2")
##D 
##D imdanova_Filt <- imdanova_filter(omicsData = mypepData)
##D mypepData <- applyFilt(filter_object = imdanova_Filt,
##D                        omicsData = mypepData,
##D                        min_nonmiss_anova=2)
##D 
##D imd_anova_res <- imd_anova(omicsData = mypepData,
##D                            test_method = 'comb',
##D                            pval_adjust='bon')
##D 
##D result = bpquant(statRes = imd_anova_res, pepData = mypepData)
##D 
## End(Not run)




cleanEx()
nameEx("combine_lipidData")
### * combine_lipidData

flush(stderr()); flush(stdout())

### Name: combine_lipidData
### Title: Combines two omicsData objects with identical sample
###   information.
### Aliases: combine_lipidData

### ** Examples

library(pmartR)
library(pmartRdata)

obj_1 <- pmartRdata::lipid_object_neg
obj_2 <- pmartRdata::lipid_object_pos

# de-dupe any duplicate edata identifiers
obj_2$e_data[,get_edata_cname(obj_2)] <- paste0("obj_2_", obj_2$e_data[,get_edata_cname(obj_2)])

combine_object <- combine_lipidData(obj_1, obj_2)

# preprocess and group the data and keep filters/grouping structure

obj_1 <- edata_transform(obj_1, "log2")
obj_1 <- normalize_global(obj_1, "all", "median", apply_norm = TRUE)
obj_2 <- edata_transform(obj_2, "log2")
obj_2 <- normalize_global(obj_2, "all", "median", apply_norm = TRUE)

obj_1 <- group_designation(obj_1, "Condition")
obj_2 <- group_designation(obj_2, "Condition")

obj_1 <- applyFilt(molecule_filter(obj_1), obj_1, min_num = 2)
obj_2 <- applyFilt(cv_filter(obj_2),obj_2, cv_thresh = 60)

combine_lipidData(obj_1, obj_2, retain_groups = TRUE, retain_filters = TRUE)




cleanEx()
nameEx("combine_techreps")
### * combine_techreps

flush(stderr()); flush(stdout())

### Name: combine_techreps
### Title: Combine technical replicates of an omicsData object
### Aliases: combine_techreps

### ** Examples

library(pmartR)
library(pmartRdata)

data(techrep_pep_object)

pep_object_averaged = combine_techreps(techrep_pep_object)




cleanEx()
nameEx("complete_mols")
### * complete_mols

flush(stderr()); flush(stdout())

### Name: complete_mols
### Title: Identify biomolecules with no missing values across samples
### Aliases: complete_mols

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_pepData)
##D complete_peps <- complete_mols(e_data = pep_pepData$e_data, edata_id = attr(pep_pepData, "cnames")$edata_cname)
##D 
## End(Not run)




cleanEx()
nameEx("cor_result")
### * cor_result

flush(stderr()); flush(stdout())

### Name: cor_result
### Title: Correlation matrix of biomolecule data
### Aliases: cor_result

### ** Examples

## Not run: 
##D library(pmartR)
##D 
##D data(pep_object)
##D 
##D my_correlation <- cor_result(omicsData = pep_object)
## End(Not run)




cleanEx()
nameEx("custom_filter")
### * custom_filter

flush(stderr()); flush(stdout())

### Name: custom_filter
### Title: Custom Filter
### Aliases: custom_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("metab_object")
##D to_filter <- custom_filter(metab_object, e_data_remove = "fumaric acid",
##D                            f_data_remove = "Infection1")
##D summary(to_filter)
##D to_filter2 <- custom_filter(metab_object, e_data_remove = "fumaric acid")
##D summary(to_filter2)
## End(Not run)




cleanEx()
nameEx("custom_sampnames")
### * custom_sampnames

flush(stderr()); flush(stdout())

### Name: custom_sampnames
### Title: Creates custom sample names for plots
### Aliases: custom_sampnames

### ** Examples

## Not run: 
##D data(pep_object)
##D results = custom_sampnames(pep_object, firstn = 5)
## End(Not run)




cleanEx()
nameEx("cv_filter")
### * cv_filter

flush(stderr()); flush(stdout())

### Name: cv_filter
### Title: Filter Based on Pooled Coefficient of Variation (CV) Values
### Aliases: cv_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D pep_object2 <- group_designation(omicsData = pep_object,
##D                                  main_effects = "Condition")
##D to_filter <- cv_filter(omicsData = pep_object2, use_groups = TRUE)
##D summary(to_filter, cv_threshold = 30)
## End(Not run)




cleanEx()
nameEx("dim_reduction")
### * dim_reduction

flush(stderr()); flush(stdout())

### Name: dim_reduction
### Title: Reduce Dimension of Data for Exploratory Data Analysis
### Aliases: dim_reduction

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D lipid_object <- edata_transform(omicsData = lipid_object,
##D                                 data_scale="log2")
##D lipid_object <- group_designation(omicsData = lipid_object,
##D                                   main_effects = "Condition")
##D pca_lipids <- dim_reduction(omicsData = lipid_object)
##D plot(pca_lipids)
##D summary(pca_lipids)
## End(Not run)




cleanEx()
nameEx("dot-is_edata")
### * dot-is_edata

flush(stderr()); flush(stdout())

### Name: .is_edata
### Title: Test if a file is an edata file
### Aliases: .is_edata

### ** Examples

## Not run: 
##D 
##D .is_edata(pmartRdata::lipid_edata)
##D 
## End(Not run)



cleanEx()
nameEx("edata_replace")
### * edata_replace

flush(stderr()); flush(stdout())

### Name: edata_replace
### Title: Replace Values Equal to x with y
### Aliases: edata_replace

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(metab_object)
##D metab_object2 <- edata_replace(omicsData = metab_object, x=0, y=NA)
## End(Not run)




cleanEx()
nameEx("edata_summary")
### * edata_summary

flush(stderr()); flush(stdout())

### Name: edata_summary
### Title: Creates a list of six Data Frames, one for each summarizing
###   metric
### Aliases: edata_summary

### ** Examples

## Not run: 
##D library(pmartRdata)
##D 
##D data(lipid_object)
##D 
##D lipid_object2 <- edata_transform(omicsData = lipid_object,
##D                                  data_scale = "log2")
##D 
##D lipid_object2 <- group_designation(omicsData = lipid_object,
##D                                    main_effects = "Condition")
##D 
##D edata_summary(omicsData = lipid_object2, by = "sample", groupvar = NULL)
## End(Not run)

formula'
formula'



cleanEx()
nameEx("edata_transform")
### * edata_transform

flush(stderr()); flush(stdout())

### Name: edata_transform
### Title: Apply a Transformation to the Data
### Aliases: edata_transform

### ** Examples

dontrun{
library(pmartRdata)
data(metab_object)
metab_object2 <- edata_transform(omicsData = metab_object, data_scale="log2")
attr(metab_object2, "data_info")$data_scale
}



cleanEx()
nameEx("fit_surv")
### * fit_surv

flush(stderr()); flush(stdout())

### Name: fit_surv
### Title: Basic survival analysis function
### Aliases: fit_surv

### ** Examples

dontrun{
library(MSomicsSTAT)
library(OvarianPepdataBP)

#Basic analysis without covariates
attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
sfit <- fit_surv(tcga_ovarian_pepdata_bp)
plot(sfit)

#Add some covariate information
attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status", covariates = "g__initial_pathologic_diagnosis_method_g1")
sfit <- fit_surv(tcga_ovarian_pepdata_bp)
plot(sfit,col=c(1,2))
}



cleanEx()
nameEx("get_comparisons")
### * get_comparisons

flush(stderr()); flush(stdout())

### Name: get_comparisons
### Title: Return comparisons of statRes object
### Aliases: get_comparisons

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D 
##D my_prodata = group_designation(omicsData = pro_object,
##D                                main_effects = c("Condition"))
##D 
##D imdanova_Filt = imdanova_filter(omicsData = my_prodata)
##D 
##D my_prodata = applyFilt(filter_object = imdanova_Filt,
##D                        omicsData = my_prodata,
##D                        min_nonmiss_anova=2)
##D                        
##D imd_anova_res = imd_anova(omicsData = my_prodata,
##D                           test_method = 'comb',
##D                           pval_adjust='bon')
##D 
##D result = get_comparisons(imd_anova_res)
## End(Not run)




cleanEx()
nameEx("get_data_class")
### * get_data_class

flush(stderr()); flush(stdout())

### Name: get_data_class
### Title: Return data_class of statRes or trellData object
### Aliases: get_data_class

### ** Examples

## Not run: 
##D library(pmartRdata)
##D 
##D my_prodata = group_designation(omicsData = pro_object,
##D                                main_effects = c("Condition"))
##D 
##D imdanova_Filt = imdanova_filter(omicsData = my_prodata)
##D 
##D my_prodata = applyFilt(filter_object = imdanova_Filt,
##D                        omicsData = my_prodata,
##D                        min_nonmiss_anova=2)
##D 
##D imd_anova_res = imd_anova(omicsData = my_prodata,
##D                           test_method = 'comb',
##D                           pval_adjust='bon')
##D 
##D result = get_data_class(imd_anova_res)
## End(Not run)




cleanEx()
nameEx("get_spans_params")
### * get_spans_params

flush(stderr()); flush(stdout())

### Name: get_spans_params
### Title: Gets the parameters for the highest ranked methods from spans.
### Aliases: get_spans_params

### ** Examples


library(pmartR)
library(pmartRdata)

data(pep_object)

# data must be log transformed and grouped
myobject <- edata_transform(pep_object, data_scale = "log2")
myobject <- group_designation(myobject, main_effects = "Condition")

spans_result <- spans_procedure(myobject)

# a list of the parameters for any normalization procedure with the best SPANS score
best_params <- get_spans_params(spans_result)

# extract the arguments from the first list element
subset_fn = best_params[[1]]$subset_fn
norm_fn = best_params[[1]]$norm_fn
params = best_params[[1]]$params

# pass arguments to normalize global
norm_object <- normalize_global(omicsData = myobject, subset_fn = subset_fn, norm_fn = norm_fn, params = params)




cleanEx()
nameEx("group_designation")
### * group_designation

flush(stderr()); flush(stdout())

### Name: group_designation
### Title: Creates Attribute of omicsData object containing Group
###   Membership Based on Specified Main Effects
### Aliases: group_designation

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D lipid_object2 <- group_designation(omicsData = lipid_object,
##D                                    main_effects = "Condition")
##D attr(lipid_object2, "group_DF")
## End(Not run)




cleanEx()
nameEx("gtest_filter")
### * gtest_filter

flush(stderr()); flush(stdout())

### Name: gtest_filter
### Title: Identifies peptides to be filtered out in preparation for
###   IMD-ANOVA.
### Aliases: gtest_filter

### ** Examples

## Not run: 
##D library(pmartR)
##D data(pep_object)
##D 
##D pep_object2 <- group_designation(omicsData = pep_object,
##D                                  main_effects = "Condition")
##D 
##D nonmissing_result <- nonmissing_per_group(omicsData = pep_object2)
##D 
##D to_filter <- gtest_filter(nonmiss_per_group = nonmissing_result,
##D                           min_nonmiss_gtest = 3)
## End(Not run)




cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



cleanEx()
nameEx("imd_anova")
### * imd_anova

flush(stderr()); flush(stdout())

### Name: imd_anova
### Title: Tests for a qualitative and quantitative difference between
###   groups using IMD and ANOVA, respectively
### Aliases: imd_anova

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D #Transform the data
##D mypepData <- edata_transform(omicsData = pep_object, data_scale = "log2")
##D 
##D #Group the data by condition
##D mypepData <- group_designation(omicsData = mypepData, main_effects = c("Condition"))
##D 
##D #Apply the IMD ANOVA filter
##D imdanova_Filt <- imdanova_filter(omicsData = mypepData)
##D mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
##D 
##D #Implement the IMD ANOVA method and compute all pairwise comparisons (i.e. leave the `comparisons` argument NULL)
##D anova_res <- imd_anova(omicsData = mypepData, test_method = 'anova')
##D imd_res <- imd_anova(omicsData = mypepData, test_method = 'gtest')
##D imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='bon')
##D 
##D #Test with really big dataset
##D library(OvarianPepdataBP)
##D tcga_ovarian_pepdata_bp <- as.pepData(e_data = tcga_ovarian_pepdata_bp$e_data, f_data = tcga_ovarian_pepdata_bp$f_data, e_meta = tcga_ovarian_pepdata_bp$e_meta)
##D tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("race"))
##D ovarian_res <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, test_method = 'anova')
##D #Tukey adjustment is super slow right now because "ptukey" is super slow, not sure how to fix that
##D ovarian_res_holm <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'holm', test_method='gtest')
##D #Dunnett adjustment, should give an error because dunnett correction shouldn't be applied for all pairwise comparisons
##D ovarian_res_dunnett <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, pval_adjust = 'dunnett', test_method='combined')
##D 
##D #Test really big dataset, two factors
##D tcga_ovarian_pepdata_bp <- group_designation(omicsData = tcga_ovarian_pepdata_bp, main_effects = c("vital_status","neoplasm_histologic_grade"))
##D ovarian_res_twofac <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, test_method='comb')
##D 
##D #Same but only test main effects (Dead vs Alive, G2 vs G3)
##D comp_df <- data.frame(Control=c("Alive","G2"), Test=c("Dead","G3"))
##D ovarian_res_twofac_main_effects <- imd_anova(omicsData = tcga_ovarian_pepdata_bp, comparisons = comp_df, test_method='comb')
## End(Not run)




cleanEx()
nameEx("imdanova_filter")
### * imdanova_filter

flush(stderr()); flush(stdout())

### Name: imdanova_filter
### Title: IMD-ANOVA filter object
### Aliases: imdanova_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D pep_pepData2 <- group_designation(omicsData = pep_object,
##D                                   main_effects = "Condition")
##D to_filter <- imdanova_filter(omicsData = pep_pepData2)
##D summary(to_filter, min_nonmiss_anova = 2)
## End(Not run)




cleanEx()
nameEx("los")
### * los

flush(stderr()); flush(stdout())

### Name: los
### Title: Identify Features from the Top L Order Statistics
### Aliases: los

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_pepData)
##D pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
##D pep_subset <- los(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, pmartR_groupDF = attr(pep_pepData2, "group_DF"))
## End(Not run)




cleanEx()
nameEx("mad_transform")
### * mad_transform

flush(stderr()); flush(stdout())

### Name: mad_transform
### Title: Median Absolute Deviation Transformation
### Aliases: mad_transform

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
##D lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, pmartR_groupDF = attr(lipid_object2, "group_DF"))
##D lipids_mad <- mad_transform(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
## End(Not run)



cleanEx()
nameEx("median_center")
### * median_center

flush(stderr()); flush(stdout())

### Name: median_center
### Title: Median Center Transformation
### Aliases: median_center

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D lipid_object2 <- group_designation(omicsData = lipid_object, main_effects = "Condition")
##D lipid_subset <- los(e_data = lipid_object2$e_data, edata_id = attr(lipid_object2, "cnames")$edata_cname, pmartR_groupDF = attr(lipid_object2, "group_DF"))
##D lipids_median <- median_center(e_data = lipid_object2$e_data, edata_id = attr(lipid_object, "cnames")$edata_cname, feature_subset = lipid_subset, backtransform = TRUE)
## End(Not run)




cleanEx()
nameEx("missingval_result")
### * missingval_result

flush(stderr()); flush(stdout())

### Name: missingval_result
### Title: Creates an object of class naRes
### Aliases: missingval_result

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("lipid_object")
##D data("metab_object")
##D 
##D result = missingval_result(lipid_object)
##D result2 = missingval_result(metab_object)
##D 
## End(Not run)




cleanEx()
nameEx("molecule_filter")
### * molecule_filter

flush(stderr()); flush(stdout())

### Name: molecule_filter
### Title: Molecule filter object
### Aliases: molecule_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D to_filter <- molecule_filter(omicsData = pep_object)
##D summary(to_filter, min_num = 2)
## End(Not run)




cleanEx()
nameEx("normalize_global")
### * normalize_global

flush(stderr()); flush(stdout())

### Name: normalize_global
### Title: Calculate Normalization Parameters
### Aliases: normalize_global

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D lipid_object <- edata_transform(omicsData = lipid_object,
##D                                 data_scale="log2")
##D lipid_object <- group_designation(omicsData = lipid_object,
##D                                   main_effects = "Condition")
##D norm_object <- normalize_global(omicsData = lipid_object,
##D                                 subset_fn = "all",
##D                                 norm_fn = "median")
##D norm_data <- normalize_global(omicsData = lipid_object,
##D                               subset_fn = "all",
##D                               norm_fn = "median",
##D                               apply_norm = TRUE,
##D                               backtransform = TRUE)
## End(Not run)




cleanEx()
nameEx("normalize_isobaric")
### * normalize_isobaric

flush(stderr()); flush(stdout())

### Name: normalize_isobaric
### Title: Normalize an object of class isobaricpepData
### Aliases: normalize_isobaric

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(isobaric_object)
##D 
##D isobaric_object = edata_transform(isobaric_object, "log2")
##D isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set",
##D                                    apply_norm = TRUE,
##D                                    channel_cname = "iTRAQ.Channel",
##D                                    refpool_channel = "116")
##D 
##D # alternate specification: #
##D data(isobaric_object)
##D 
##D isobaric_object = edata_transform(isobaric_object, "log2")
##D isobaric_norm = normalize_isobaric(isobaric_object, exp_cname = "Set",
##D                                    apply_norm = TRUE,
##D                                    refpool_cname = "Reference",
##D                                    refpool_notation = "Yes")
##D 
## End(Not run)




cleanEx()
nameEx("normalize_loess")
### * normalize_loess

flush(stderr()); flush(stdout())

### Name: normalize_loess
### Title: Loess Normalization
### Aliases: normalize_loess

### ** Examples

dontrun{
library(pmartR)
library(pmartRdata)
pep_object = edata_transform(pep_object, "log2")
result = normalize_loess(pep_object)

}




cleanEx()
nameEx("normalize_nmr")
### * normalize_nmr

flush(stderr()); flush(stdout())

### Name: normalize_nmr
### Title: Normalize an object of class nmrData
### Aliases: normalize_nmr

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(nmr_object_identified)
##D 
##D nmr_object = edata_transform(nmr_object_identified, "log2")
##D nmr_norm = normalize_nmr(nmr_object, apply_norm = TRUE,
##D                          metabolite_name = "unkm1.53")
##D 
##D # alternate specification: #
##D data(nmr_object_identified)
##D 
##D nmr_object = edata_transform(nmr_object, "log2")
##D nmr_norm = normalize_nmr(nmr_object, apply_norm = TRUE,
##D                          sample_property_cname = "Concentration")
##D 
## End(Not run)




cleanEx()
nameEx("normalize_quantile")
### * normalize_quantile

flush(stderr()); flush(stdout())

### Name: normalize_quantile
### Title: Quantile Normalization
### Aliases: normalize_quantile

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(lipid_object)
##D norm_data <- normalize_quantile(omicsData = lipid_object)
## End(Not run)




cleanEx()
nameEx("plot-RNAFilt")
### * plot-RNAFilt

flush(stderr()); flush(stdout())

### Name: plot.RNAFilt
### Title: plot.RNAFilt
### Aliases: plot.RNAFilt

### ** Examples

## Not run: 
##D data(seq_object)
##D seqfilt <- total_count_filter(pep_object)
##D plot(seqfilt, min_count = 5)
## End(Not run)




cleanEx()
nameEx("plot-cvFilt")
### * plot-cvFilt

flush(stderr()); flush(stdout())

### Name: plot.cvFilt
### Title: plot.cvFilt
### Aliases: plot.cvFilt

### ** Examples

## Not run: 
##D library(pmartRdata)
##D
##D data(pep_object)
##D 
##D pep_object <- group_designation(omicsData = pep_object,
##D                                 main_effects = "Condition")
##D 
##D cvfilt <- cv_filter(pep_object)
##D 
##D plot(cvfilt, cv_threshold = 20)
##D plot(cvfilt, cv_threshold = 10, log_scale = FALSE)
## End(Not run)




cleanEx()
nameEx("plot-dataRes")
### * plot-dataRes

flush(stderr()); flush(stdout())

### Name: plot.dataRes
### Title: Plots an object of class dataRes
### Aliases: plot.dataRes

### ** Examples

## Not run: 
##D library (pmartRdata)
##D 
##D data(lipid_object)
##D 
##D lipid_object = edata_transform(lipid_object, "log2")
##D 
##D result = edata_summary(omicsData = lipid_object,
##D                        by = "molecule",
##D                        groupvar = "Condition")
##D 
##D plot(result)
## End(Not run)




cleanEx()
nameEx("plot-isobaricnormRes")
### * plot-isobaricnormRes

flush(stderr()); flush(stdout())

### Name: plot.isobaricnormRes
### Title: Plots an object of class isobaricnormRes
### Aliases: plot.isobaricnormRes

### ** Examples

## Not run: 
##D library(pmartRdata)
##D 
##D data(isobaric_object)
##D 
##D isobaric_object = edata_transform(isobaric_object, "log2")
##D 
##D result = normalize_isobaric(isobaric_object, exp_cname = "Set",
##D                             apply_norm = FALSE,
##D                             channel_cname = "iTRAQ.Channel",
##D                             refpool_channel = "116")
##D 
##D plot(result)
## End(Not run)




cleanEx()
nameEx("plot-moleculeFilt")
### * plot-moleculeFilt

flush(stderr()); flush(stdout())

### Name: plot.moleculeFilt
### Title: plot.moleculeFilt
### Aliases: plot.moleculeFilt

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_object)
##D molfilt <- molecule_filter(pep_object)
##D plot(molfilt, min_num = 5)
##D plot(molfilt, min_num = 3, cumulative = FALSE)
## End(Not run)




cleanEx()
nameEx("plot-naRes")
### * plot-naRes

flush(stderr()); flush(stdout())

### Name: plot.naRes
### Title: Plots an object of class naRes
### Aliases: plot.naRes

### ** Examples

## Not run: 
##D library(pmartRdata)
##D 
##D data("lipid_object")
##D 
##D result<- missingval_result(lipid_object)
##D 
##D plot(result, plog_type = "bar", x_lab_angle = 50)
## End(Not run)




cleanEx()
nameEx("plot-nmrnormRes")
### * plot-nmrnormRes

flush(stderr()); flush(stdout())

### Name: plot.nmrnormRes
### Title: Plots an object of class nmrnormRes
### Aliases: plot.nmrnormRes

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(nmr_object_identified)
##D 
##D nmr_object = edata_transform(nmr_object_identified, "log2")
##D nmr_norm = normalize_nmr(nmr_object,
##D                          apply_norm = FALSE,
##D                          metabolite_name = "unkm1.53")
##D plot(nmr_norm)
##D 
##D # alternate specification: #
##D data(nmr_object_identified)
##D 
##D nmr_object = edata_transform(nmr_object, "log2")
##D nmr_norm = normalize_nmr(nmr_object,
##D                          apply_norm = FALSE,
##D                          sample_property_cname = "Concentration")
##D plot(nmr_norm)
## End(Not run)




cleanEx()
nameEx("plot-proteomicsFilt")
### * plot-proteomicsFilt

flush(stderr()); flush(stdout())

### Name: plot.proteomicsFilt
### Title: plot.proteomicsFilt
### Aliases: plot.proteomicsFilt

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_object)
##D profilt <- proteomics_filter(pep_object)
##D plot(profilt, min_num_peps = 5)
## End(Not run)




cleanEx()
nameEx("plot-totalCountFilt")
### * plot-totalCountFilt

flush(stderr()); flush(stdout())

### Name: plot.totalCountFilt
### Title: plot.totalCountFilt
### Aliases: plot.totalCountFilt

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(seq_object)
##D seqfilt <- total_count_filter(pep_object)
##D plot(seqfilt, min_count = 5)
## End(Not run)




cleanEx()
nameEx("plot.statRes")
### * plot.statRes

flush(stderr()); flush(stdout())

### Name: plot.statRes
### Title: Plotting function for 'statRes' objects
### Aliases: plot.statRes

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D #Transform the data
##D 
##D #Group the data by condition
##D myproData <- group_designation(omicsData = pro_object,
##D                                main_effects = c("Condition"))
##D 
##D #Apply the IMD ANOVA filter
##D imdanova_Filt <- imdanova_filter(omicsData = myproData)
##D myproData <- applyFilt(filter_obj = imdanova_Filt,
##D                        omicsData = myproData,
##D                        min_nonmiss_anova=2)
##D 
##D #Implement the IMD ANOVA method and compuate all pairwise comparisons
##D #(i.e. leave the `comparisons` argument NULL)
##D anova_res <- imd_anova(omicsData = myproData, test_method = 'anova')
##D plot(anova_res)
##D plot(anova_res, plot_type = "volcano")
##D 
##D imd_res <- imd_anova(omicsData = myproData, test_method = 'gtest')
##D plot(imd_res)
##D plot(imd_res, plot_type = "gheatmap")
##D # using arguments of internal functions:
##D plot(imd_res,
##D      plot_type = "gheatmap",
##D      color_low = "red",
##D      color_high = "green")
##D 
##D imd_anova_res <- imd_anova(omicsData = myproData,
##D                            test_method = 'comb',
##D                            pval_adjust='bon')
##D plot(imd_anova_res, bw_theme = TRUE)
##D plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE)
##D 
## End(Not run)




cleanEx()
nameEx("plot_km")
### * plot_km

flush(stderr()); flush(stdout())

### Name: plot_km
### Title: Basic survival analysis plot
### Aliases: plot_km

### ** Examples

library(MSomicsSTAT)
library(OvarianPepdataBP)
attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
plot_km(omicsData = tcga_ovarian_pepdata_bp)

#Add covariates to "survDF" attribute
attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status", covariates = "age_at_initial_pathologic_diagnosis")
plot_km(omicsData = tcga_ovarian_pepdata_bp)



cleanEx()
nameEx("ppp")
### * ppp

flush(stderr()); flush(stdout())

### Name: ppp
### Title: Identify Proportion of the Peptides Present (PPP) Biomolecules
### Aliases: ppp

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_pepData)
##D pep_subset <- ppp(e_data = pep_pepData$e_data, edata_id = attr(pep_pepData, "cnames")$edata_cname)
##D 
##D pep_subset <- ppp(e_data = pep_pepData$e_data, edata_id = attr(pep_pepData, "cnames")$edata_cname, percent = 0.6)
## End(Not run)




cleanEx()
nameEx("ppp_rip")
### * ppp_rip

flush(stderr()); flush(stdout())

### Name: ppp_rip
### Title: Identify Proportion of Peptides Present (PPP) and Rank Invariant
###   Peptides (RIP)
### Aliases: ppp_rip

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_pepData)
##D pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
##D pep_subset <- ppp_rip(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, samp_id = attr(pep_pepData2, "cnames")$fdata_cname, pmartR_groupDF = attr(pep_pepData2, "group_DF"))
## End(Not run)




cleanEx()
nameEx("protein_quant")
### * protein_quant

flush(stderr()); flush(stdout())

### Name: protein_quant
### Title: protein_quant wrapper function
### Aliases: protein_quant

### ** Examples

## Not run: 
##D library(pmartR)
##D library(pmartRdata)
##D 
##D mypepData <- group_designation(omicsData = pep_object, main_effects = c("Condition"))
##D mypepData = edata_transform(mypepData, "log2")
##D 
##D imdanova_Filt <- imdanova_filter(omicsData = mypepData)
##D mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
##D 
##D imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='bon')
##D 
##D isoformRes = bpquant(statRes = imd_anova_res, pepData = mypepData)
##D 
##D #case where isoformRes is NULL:
##D results<- protein_quant(pepData = mypepData, method = 'rollup', combine_fn = 'median', isoformRes = NULL)
##D 
##D #case where isoformRes is provided:
##D results2 = protein_quant(pepData = mypepData, method = 'rollup', combine_fn = 'mean', isoformRes = isoformRes)
## End(Not run)




cleanEx()
nameEx("proteomics_filter")
### * proteomics_filter

flush(stderr()); flush(stdout())

### Name: proteomics_filter
### Title: Proteomics filter object
### Aliases: proteomics_filter

### ** Examples

## Not run: 
##D library(pmartR)
##D data("pep_object")
##D my_filter <- proteomics_filter(omicsData = pep_object)
##D summary(my_filter, min_num_peps = 3)
## End(Not run)




cleanEx()
nameEx("qrollup")
### * qrollup

flush(stderr()); flush(stdout())

### Name: qrollup
### Title: Applies qrollup function
### Aliases: qrollup

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_object)
##D result = qrollup(pepData = pep_object, qrollup_thresh = 2)
## End(Not run)




cleanEx()
nameEx("report_dataRes")
### * report_dataRes

flush(stderr()); flush(stdout())

### Name: report_dataRes
### Title: Creates a data frame displaying multiple metrics
### Aliases: report_dataRes

### ** Examples

dontrun{
library(pmartRdata)
data(lipid_object)
lipid_object2 <- edata_transform(omicsData = lipid_object, data_scale = "log2")

dataRes_sample = edata_summary(omicsData = lipid_object2, groupvar = NULL, by = "sample")
report_dataRes(dataRes_sample)
}





cleanEx()
nameEx("rip")
### * rip

flush(stderr()); flush(stdout())

### Name: rip
### Title: Identify Rank-Invariant Peptides
### Aliases: rip

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_pepData)
##D pep_pepData2 <- group_designation(omicsData = pep_pepData, main_effects = "Condition")
##D pep_subset <- rip(e_data = pep_pepData2$e_data, edata_id = attr(pep_pepData2, "cnames")$edata_cname, pmartR_groupDF = attr(pep_pepData2, "group_DF"))
## End(Not run)




cleanEx()
nameEx("rmd_conversion")
### * rmd_conversion

flush(stderr()); flush(stdout())

### Name: rmd_conversion
### Title: Conversion between log2(RMD) and p-value
### Aliases: rmd_conversion

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(metab_object)
##D metab_object2 <- edata_transform(omicsData = metab_object,
##D                                  data_scale = "log2")
##D metab_object3 <- group_designation(omicsData = metab_object2,
##D                                    main_effects = "Condition")
##D rmd_results <- rmd_filter(omicsData = metab_object3,
##D                           metrics=c("MAD", "Skewness", "Correlation"))
##D rmd_conversion(log2rmd = rmd_results$Log2.md, df=3)
##D 
##D rmd_conversion(pval = .0001, df = 3)
##D rmd_conversion(log2rmd = 4.5, df = 3)
## End(Not run)




cleanEx()
nameEx("rmd_filter")
### * rmd_filter

flush(stderr()); flush(stdout())

### Name: rmd_filter
### Title: RMD Runs
### Aliases: rmd_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(metab_object)
##D metab_object2 <- edata_transform(omicsData = metab_object, data_scale = "log2")
##D metab_object3 <- group_designation(omicsData = metab_object2, main_effects = "Condition")
##D rmd_results <- rmd_filter(omicsData = metab_object3, metrics=c("MAD", "Skewness", "Correlation"))
##D rmd_results <- rmd_filter(omicsData = metab_object2)
##D 
##D data(pep_pepData)
##D pep_pepData2 <- edata_transform(omicsData = pep_object, data_scale = "log2")
##D pep_pepData3 <- group_designation(omicsData = pep_pepData2, main_effects = "Condition")
##D rmd_results <- rmd_filter(omicsData = pep_pepData3)
## End(Not run)




cleanEx()
nameEx("rrollup")
### * rrollup

flush(stderr()); flush(stdout())

### Name: rrollup
### Title: Applies rrollup function
### Aliases: rrollup

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_object)
##D result = rrollup(pepData = pep_object)
## End(Not run)




cleanEx()
nameEx("spans_procedure")
### * spans_procedure

flush(stderr()); flush(stdout())

### Name: spans_procedure
### Title: Calculate SPANS Score for a Number of Normalization Methods
### Aliases: spans_procedure

### ** Examples

library(pmartR)
library(pmartRdata)

pep_object
pep_object <- edata_transform(pep_object, data_scale = "log2")
pep_object <- group_designation(pep_object, main_effects = "Condition")

## default parameters
# spans_res <- spans_procedure(pep_object)

## specify only certain subset and normalization functions
# spans_res <- spans_procedure(pep_object, norm_fn = c("median", "zscore"), subset_fn = c("all", "los", "ppp"))

## specify parameters for supplied subset functions, notice ppp_rip takes a vector of two numeric arguments.
# spans_res <- spans_procedure(pep_object, subset_fn = c("all", "los", "ppp"), params = list(los = list(0.25, 0.5), ppp = list(0.15, 0.25)))
# spans_res <- spans_procedure(pep_object, subset_fn = c("all", "rip", "ppp_rip"), params = list(rip = list(0.3, 0.4), ppp_rip = list(c(0.15, 0.5), c(0.25, 0.5))))




cleanEx()
nameEx("summary-isobaricnormRes")
### * summary-isobaricnormRes

flush(stderr()); flush(stdout())

### Name: summary-pmartR
### Title: Summarizes an object of class isobaricnormRes
### Aliases: summary-pmartR summary.isobaricnormRes

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(isobaric_object)
##D 
##D isobaric_object = edata_transform(isobaric_object, "log2")
##D result = normalize_isobaric(isobaric_object, exp_cname = "Set", apply_norm = FALSE, channel_cname = "iTRAQ.Channel", refpool_channel = "116")
##D 
##D summary(result)
## End(Not run)




cleanEx()
nameEx("summary-nmrnormRes")
### * summary-nmrnormRes

flush(stderr()); flush(stdout())

### Name: summary-pmartR
### Title: Summarizes an object of class nmrnormRes
### Aliases: summary-pmartR summary.nmrnormRes

### ** Examples

dontrun{
library(pmartRdata)
data(nmr_object_identified)

nmr_object = edata_transform(nmr_object_identified, "log2")
nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, metabolite_name = "unkm1.53")
summary(nmr_norm)

# alternate specification: #
data(nmr_object_identified)

nmr_object = edata_transform(nmr_object, "log2")
nmr_norm = normalize_nmr(nmr_object, apply_norm = FALSE, sample_property_cname = "Concentration")
summary(nmr_norm)

}




cleanEx()
nameEx("summary-omicsData")
### * summary-omicsData

flush(stderr()); flush(stdout())

### Name: summary-omicsData
### Title: Produce a basic summary of a pmartR omicsData S3 Object
### Aliases: summary-omicsData summary.pepData summary.proData
###   summary.lipidData summary.metabData summary.nmrData summary.seqData

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("pep_object")
##D summary(pep_object)
## End(Not run)




cleanEx()
nameEx("summary-pmartR-results")
### * summary-pmartR-results

flush(stderr()); flush(stdout())

### Name: summary-pmartR-results
### Title: Produce summaries of various results objects from pmartR
###   functions.
### Aliases: summary-pmartR-results summary.normRes summary.SPANSRes
###   summary.dimRes summary.corRes

### ** Examples

dontrun{
library(pmartRdata)
data("pep_object")

pep_object <- group_designation(pep_object, main_effects = 'Condition')

norm_result <- normalize_global(pep_object, norm_fn='median', subset_fn='all')
summary(norm_result)

spans_results <- spans_procedure(pep_object)
summary(spans_results)
}




cleanEx()
nameEx("summary-trelliData")
### * summary-trelliData

flush(stderr()); flush(stdout())

### Name: summary-trelliData
### Title: Summarizes potential plotting options for a trelliData object
### Aliases: summary-trelliData summary.trelliData

### ** Examples

## Not run: 
##D 
##D library(dplyr)
##D 
##D # Use an edata example. Build with as.trelliData.edata.
##D summary(trelliData)
##D summary(trelliData %>% trelli_panel_by("LipidCommonName"))
##D summary(trelliData %>% trelli_panel_by("Sample"))
##D 
##D # Use an omicsData example. Build with as.trelliData.
##D summary(trelliData2)
##D 
##D # Use a statRes example. Build with as.trelliData. 
##D summary(trelliData3)
##D 
## End(Not run)




cleanEx()
nameEx("summary.RNAFilt")
### * summary.RNAFilt

flush(stderr()); flush(stdout())

### Name: summary.RNAFilt
### Title: Produce a basic summary of a RNAFilt object
### Aliases: summary.RNAFilt

### ** Examples

dontrun{
library(pmartRdata)
data("pep_object")
myfilter <- molecule_filter(pep_object)
summary(myfilter)
summary(myfilter, min_num = 2)
}




cleanEx()
nameEx("summary.moleculeFilt")
### * summary.moleculeFilt

flush(stderr()); flush(stdout())

### Name: summary.moleculeFilt
### Title: Produce a basic summary of a molecule_filter object
### Aliases: summary.moleculeFilt

### ** Examples

dontrun{
library(pmartRdata)
data("pep_object")
myfilter <- molecule_filter(pep_object)
summary(myfilter)
summary(myfilter, min_num = 2)
}




cleanEx()
nameEx("summary_km")
### * summary_km

flush(stderr()); flush(stdout())

### Name: summary_km
### Title: Basic survival analysis summary
### Aliases: summary_km

### ** Examples

dontrun{
library(OvarianPepdataBP)
attr(tcga_ovarian_pepdata_bp,"survDF") <- list(t_death = "survival_time",ind_death = "vital_status")
#No percent is provided so the entire object is returned
summary_km(tcga_ovarian_pepdata_bp)

#Percent is provided so corresponding time point is returned
summary_km(tcga_ovarian_pepdata_bp, .4)
}



cleanEx()
nameEx("total_count_filter")
### * total_count_filter

flush(stderr()); flush(stdout())

### Name: total_count_filter
### Title: Total Count filter object
### Aliases: total_count_filter

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data("seq_object")
##D to_filter <- total_count_filter(omicsData = seq_object)
##D summary(to_filter, min_num = 2)
## End(Not run)




cleanEx()
nameEx("trelli_abundance_boxplot")
### * trelli_abundance_boxplot

flush(stderr()); flush(stdout())

### Name: trelli_abundance_boxplot
### Title: Boxplot trelliscope building function for abundance data
### Aliases: trelli_abundance_boxplot

### ** Examples

## Not run: 
##D 
##D ## Build the abundance boxplot with an edata file. Generate trelliData in as.trelliData.edata
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName") %>% 
##D    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
##D trelli_panel_by(trelliData = trelliData, panel = "Sample") %>% trelli_abundance_boxplot()
##D 
##D ## Build the abundance boxplot with an omicsData object. Generate trelliData in as.trelliData
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName") %>% 
##D    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidFamily") %>% trelli_abundance_boxplot()
##D     
##D ## Build the abundance boxplot with an omicsData and statRes object. Generate trelliData in as.trelliData.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidCommonName") %>%
##D    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10)
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% trelli_abundance_boxplot()
##D    
##D ## Other options include modifying the ggplot  
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName") %>% 
##D    trelli_abundance_boxplot(test_mode = TRUE, test_example = 1:10, 
##D      ggplot_params = c("ylab('')", "ylim(c(2,20))"))
##D 
##D ## Or making the plot interactive 
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% trelli_abundance_boxplot(interactive = TRUE)
##D 
## End(Not run)




cleanEx()
nameEx("trelli_abundance_heatmap")
### * trelli_abundance_heatmap

flush(stderr()); flush(stdout())

### Name: trelli_abundance_heatmap
### Title: Heatmap trelliscope building function for abundance data
### Aliases: trelli_abundance_heatmap

### ** Examples

## Not run: 
##D 
##D ## Build the abundance heatmap with an omicsData object with emeta variables. Generate trelliData in as.trelliData.
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidFamily") %>% 
##D    trelli_abundance_heatmap(test_mode = TRUE, test_example = 1:3)
##D    
##D ## Users can modify the plotting function with ggplot parameters and interactivity, 
##D ## and can also select certain cognostics.     
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
##D    trelli_abundance_heatmap(test_mode = TRUE, test_example = 1:5, 
##D      ggplot_params = c("ylab('')", "xlab('')"), interactive = TRUE, cognostics = c("mean", "median"))  
##D    
## End(Not run)





cleanEx()
nameEx("trelli_abundance_histogram")
### * trelli_abundance_histogram

flush(stderr()); flush(stdout())

### Name: trelli_abundance_histogram
### Title: Histogram trelliscope building function for abundance data
### Aliases: trelli_abundance_histogram

### ** Examples

## Not run: 
##D 
##D ## Build the abundance histogram with an edata file. Generate trelliData in as.trelliData.edata
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName") %>% 
##D    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
##D 
##D ## Build the abundance histogram with an omicsData object. Generate trelliData in as.trelliData
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName") %>% 
##D    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
##D     
##D ## Build the abundance histogram with an omicsData and statRes object. Generate trelliData in as.trelliData.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidCommonName") %>%
##D    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10)
##D    
##D ## Users can modify the plotting function with ggplot parameters and interactivity, 
##D ## and can also select certain cognostics.     
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName") %>% 
##D    trelli_abundance_histogram(test_mode = TRUE, test_example = 1:10, 
##D      ggplot_params = c("ylab('')", "xlab('Abundance')"), interactive = TRUE,
##D      cognostics = c("mean", "median"))  
##D    
## End(Not run)





cleanEx()
nameEx("trelli_foldchange_bar")
### * trelli_foldchange_bar

flush(stderr()); flush(stdout())

### Name: trelli_foldchange_bar
### Title: Bar chart trelliscope building function for fold_change
### Aliases: trelli_foldchange_bar

### ** Examples

## Not run: 
##D 
##D ## Build fold_change bar plot with statRes data grouped by edata_colname.
##D trelli_panel_by(trelliData = trelliData3, panel = "LipidCommonName") %>% 
##D   trelli_foldchange_bar(test_mode = TRUE, test_example = 1:10, p_value_test = TRUE)
##D   
##D ## Or make the plot interactive  
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidCommonName") %>% 
##D   trelli_foldchange_bar(test_mode = TRUE, test_example = 1:10, p_value_test = TRUE, interactive = TRUE) 
##D    
## End(Not run)




cleanEx()
nameEx("trelli_foldchange_boxplot")
### * trelli_foldchange_boxplot

flush(stderr()); flush(stdout())

### Name: trelli_foldchange_boxplot
### Title: Boxplot trelliscope building function for fold_change
### Aliases: trelli_foldchange_boxplot

### ** Examples

## Not run: 
##D  
##D 
##D ## Build fold_change box plot with statRes data grouped by edata_colname.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
##D   trelli_foldchange_boxplot(p_value_test = TRUE)
##D 
## End(Not run)




cleanEx()
nameEx("trelli_foldchange_heatmap")
### * trelli_foldchange_heatmap

flush(stderr()); flush(stdout())

### Name: trelli_foldchange_heatmap
### Title: Heatmap trelliscope building function for fold_change
### Aliases: trelli_foldchange_heatmap

### ** Examples

## Not run: 
##D  
##D 
##D ## Build fold_change bar plot with statRes data grouped by edata_colname.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
##D   trelli_foldchange_heatmap()
##D 
## End(Not run)




cleanEx()
nameEx("trelli_foldchange_volcano")
### * trelli_foldchange_volcano

flush(stderr()); flush(stdout())

### Name: trelli_foldchange_volcano
### Title: Volcano trelliscope building function for fold_change
### Aliases: trelli_foldchange_volcano

### ** Examples

## Not run: 
##D  
##D 
##D ## Build fold_change bar plot with statRes data grouped by edata_colname.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily") %>% 
##D   trelli_foldchange_volcano(comparison = "Mock_vs_Infection_A")
##D 
## End(Not run)




cleanEx()
nameEx("trelli_missingness_bar")
### * trelli_missingness_bar

flush(stderr()); flush(stdout())

### Name: trelli_missingness_bar
### Title: Bar chart trelliscope building function for missing data
### Aliases: trelli_missingness_bar

### ** Examples

## Not run: 
##D 
##D ## Build the missingness bar plot with an edata file. Generate trelliData in as.trelliData.edata
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName") %>% 
##D   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
##D trelli_panel_by(trelliData = trelliData, panel = "Sample") %>% trelli_missingness_bar()
##D 
##D ## Build the missingness bar plot with an omicsData object. Generate trelliData in as.trelliData
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName") %>% 
##D   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
##D 
##D ## Build the missingness bar plot with a statRes object. Generate trelliData in as.trelliData
##D trelli_panel_by(trelliData = trelliData3, panel = "LipidCommonName") %>%
##D   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10)
##D 
##D ## Build the missingness bar plot with an omicsData and statRes object. Generate trelliData in as.trelliData.
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidCommonName") %>%
##D   trelli_missingness_bar(test_mode = TRUE, test_example = 1:10) 
##D 
##D ## Or making the plot interactive 
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName") %>% 
##D    trelli_missingness_bar(test_mode = TRUE, test_example = 1:5)
##D    
##D ## Or visualize only count data 
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName") %>% 
##D    trelli_missingness_bar(test_mode = TRUE, test_example = 1:5, cognostics = "n", proportion = FALSE)
##D    
## End(Not run)




cleanEx()
nameEx("trelli_panel_by")
### * trelli_panel_by

flush(stderr()); flush(stdout())

### Name: trelli_panel_by
### Title: Set the "panel_by" variable for a trelliData object
### Aliases: trelli_panel_by

### ** Examples

## Not run: 
##D 
##D library(pmartRdata)
##D library(pmartR)
##D 
##D ## "panel_by" with an edata file. Generate with example code in as.trelliData.edata
##D trelli_panel_by(trelliData = trelliData, panel = "LipidCommonName")
##D trelli_panel_by(trelliData = trelliData, panel = "Sample")
##D 
##D ## "panel_by" with trelliData containing omicsData. Generate with example code in as.trelliData
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidCommonName")
##D trelli_panel_by(trelliData = trelliData2, panel = "LipidFamily")
##D 
##D ## "panel_by" with trelliData containing statRes. Generate with example code in as.trelliData
##D trelli_panel_by(trelliData = trelliData3, panel = "LipidCommonName")
##D 
##D ## "panel_by" with trelliData containing both omicsData and statRes. Generate with example code in as.trelliData
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidFamily")
##D trelli_panel_by(trelliData = trelliData4, panel = "LipidCommonName")
##D trelli_panel_by(trelliData = trelliData4, panel = "Sample_Name")
##D 
## End(Not run)  




cleanEx()
nameEx("write_stat_results")
### * write_stat_results

flush(stderr()); flush(stdout())

### Name: write_stat_results
### Title: Creates a list of three sheets, Normalized Data, DA Test
###   Results, and Only DA biomolecules - for OMICS project
### Aliases: write_stat_results

### ** Examples

dontrun{
library(pmartR)
library(pmartRdata)

myproData <- group_designation(omicsData = pro_object, main_effects = c("Condition"))

imdanova_Filt <- imdanova_filter(omicsData = myproData)
myproData <- applyFilt(filter_object = imdanova_Filt, omicsData = myproData, min_nonmiss_anova=2)

imd_anova_res <- imd_anova(omicsData = myproData, test_method = 'comb', pval_adjust='bon')

result = write_stat_results(omicsData = myproData, statResData = imd_anova_res, refCondition = "Mock")
}




cleanEx()
nameEx("zrollup")
### * zrollup

flush(stderr()); flush(stdout())

### Name: zrollup
### Title: Applies zrollup function
### Aliases: zrollup

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_object)
##D result = zrollup(pepData = pep_object)
## End(Not run)




cleanEx()
nameEx("zscore_transform")
### * zscore_transform

flush(stderr()); flush(stdout())

### Name: zscore_transform
### Title: Z-Score Transformation
### Aliases: zscore_transform

### ** Examples

## Not run: 
##D library(pmartRdata)
##D data(pep_edata)
##D peps_subset <- pep_edata$Mass_Tag_ID
##D peps_zscore <- zscore_transform(e_data = pep_edata, edata_id = "Mass_Tag_ID", feature_subset = peps_subset)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
