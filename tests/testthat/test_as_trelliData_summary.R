context("summary: trelliData edata and trelliData")

test_that("Check that trelliData objects return the correct data frames",{ 
  
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
  
  # Make each of the four objects 
  suppressWarnings({pepSummary1 <- as.trelliData.edata(e_data = edata, edata_cname = "Mass_Tag_ID", omics_type = "pepData") %>% summary()})
  suppressWarnings({pepSummary2 <- as.trelliData(omicsData = pepOmics) %>% summary()})
  suppressWarnings({pepSummary3 <- as.trelliData(statRes = pepStat) %>% summary()})
  suppressWarnings({pepSummary4 <- as.trelliData(omicsData = pepOmics, statRes = pepStat) %>% summary()})
  
  # Test: summary of the trelliData objects-------------------------------------
  
  # Test that each have the correct number of rows
  expect_equal(5, nrow(pepSummary1))
  expect_equal(8, nrow(pepSummary2))
  expect_equal(2, nrow(pepSummary3))
  expect_equal(12, nrow(pepSummary4))
  
  
})