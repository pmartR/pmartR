context("class: trelliData & trelli_panel_by")

test_that("as.trelliData and trelli_panel_by returns correct data frames and attributes",{ 
  
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.Rdata',
                   package = 'pmartR'))
  
  # Generate omicsData and statRes: peptide expression data---------------------
  
  # Add more classes to the fdata file 
  fdata$Condition <- c(rep("InfectionA", 5), rep("InfectionB", 4), rep("Mock", 3))
  
  # Create omics data object
  pepOmics <- as.pepData(
    e_data = edata, 
    f_data = fdata,
    e_meta = emeta,
    edata_cname = "Mass_Tag_ID", 
    fdata_cname = "SampleID",
    emeta_cname = "Protein"
  )
  
  # Create statRes object 
  pepOmics <- edata_transform(omicsData = pepOmics, data_scale = "log2")
  pepOmics <- group_designation(omicsData = pepOmics, main_effects = c("Condition"))
  pepOmics <- applyFilt(filter_object = imdanova_filter(omicsData = pepOmics), omicsData = pepOmics, min_nonmiss_anova = 2, remove_singleton_groups = FALSE)
  pepOmics <- normalize_global(pepOmics, "subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE, "backtransform" = TRUE)
  pepStat <- imd_anova(omicsData = pepOmics, test_method = "combined")
  
  # Test: peptide expression data-----------------------------------------------
  
  # There are 3 possible trelliData objects that can be built here: one with 
  # only omicsData, one with only statRes, and one with both. 
  pepTrelli_omics <- as.trelliData(omicsData = pepOmics)
  pepTrelli_stat <- as.trelliData(statRes = pepStat)
  pepTrelli_both <- as.trelliData(omicsData = pepOmics, statRes = pepStat)
  
  ## pepOmic only checks 
  
  # 
  
  

})