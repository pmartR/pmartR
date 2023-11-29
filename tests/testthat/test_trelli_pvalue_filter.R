context("function: trelli_pvalue_filter")

test_that("trelli_pvalue_filter returns correct data frames and attributes", {
  
  ##################
  ## MS/NMR TESTS ##
  ##################
  
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'
  ))
  
  # Process data and make statRes object
  pep_omicsData <- as.pepData(e_data = edata, 
                              f_data = fdata, 
                              e_meta = emeta,
                              edata_cname = "Mass_Tag_ID", 
                              fdata_cname = "SampleID",
                              emeta_cname = "Protein")
  pep_omicsData <- edata_transform(omicsData = pep_omicsData, data_scale = "log2")
  pep_omicsData <- group_designation(omicsData = pep_omicsData, main_effects = "Condition")
  imdanova_Filt <- imdanova_filter(omicsData = pep_omicsData)
  pep_omicsData <- applyFilt(filter_object = imdanova_Filt, omicsData = pep_omicsData, min_nonmiss_anova = 2)
  pep_omicsData <- normalize_global(pep_omicsData, "subset_fn" = "all", "norm_fn" = "median", "apply_norm" = TRUE, "backtransform" = TRUE)
  pep_statRes <- imd_anova(omicsData = pep_omicsData, test_method = 'combined')
  
  # Test: trelli_pvalue_filter function-----------------------------------------
  
  # The trelliData function paramater must be a trelliData object
  expect_error(trelli_pvalue_filter(trelliData = edata), "trelliData must be of the class trelliData.")
  
  # The trelliData object should have stats results in it
  expect_error(trelli_pvalue_filter(as.trelliData(omicsData = pep_omicsData)), "trelliData must contain a statRes object.")
  
  # Save statRes and statRes/omics object
  pep_trelli3 <- as.trelliData(statRes = pep_statRes)
  pep_trelli4 <- as.trelliData(omicsData = pep_omicsData, statRes = pep_statRes)
  
  # p_value_test must be gtest or anova
  expect_error(trelli_pvalue_filter(pep_trelli3, p_value_test = "bologna"), "p_value_test must be anova, or gtest.")
  
  # 
  
  
})