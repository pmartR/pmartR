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
  
  # Test: trelli_pvalue_filter function checks----------------------------------
  
  # The trelliData function paramater must be a trelliData object
  expect_error(trelli_pvalue_filter(trelliData = edata), "trelliData must be of the class trelliData.")
  
  # The trelliData object should have stats results in it
  expect_error(trelli_pvalue_filter(as.trelliData(omicsData = pep_omicsData)), "trelliData must contain a statRes object.")
  
  # Save statRes and statRes/omics object
  pep_trelli3 <- as.trelliData(statRes = pep_statRes)
  pep_trelli4 <- as.trelliData(omicsData = pep_omicsData, statRes = pep_statRes)
  
  # p_value_test must be gtest or anova
  expect_error(trelli_pvalue_filter(pep_trelli3, p_value_test = "bologna"), "p_value_test must be anova, or gtest.")
  
  # p_value_threshold must be a number 
  expect_error(trelli_pvalue_filter(pep_trelli3, p_value_thresh = "joke"), "p_value_threshold must be a number.")
  
  # comparison should be a string
  expect_error(trelli_pvalue_filter(pep_trelli3, comparison = 1), "comparison must be a string of length 1.")
  
  # and comparison should be a single string
  expect_error(trelli_pvalue_filter(pep_trelli3, comparison = c("Infection_vs_Mock", "Infection_vs_Mock")), "comparison must be a string of length 1.")
  
  # comparison should be in the comparison list
  expect_error(trelli_pvalue_filter(pep_trelli4, comparison = "Mock_vs_Infection"), "Mock_vs_Infection is not an acceptable comparison")
  
  # Test: trelli_pvalue_filter functionality------------------------------------
  
  # Filter by anova with a pvalue of 0.05
  anova_test <- trelli_pvalue_filter(pep_trelli3, comparison = "Infection_vs_Mock", p_value_test = "anova")
  
  # There should be no p-values over 0.05 
  expect_true(any(anova_test$trelliData.stat$p_value_anova > 0.05) == FALSE)
  
  # Filter by gtest with a pvalue of 0.10
  gtest_test <- trelli_pvalue_filter(pep_trelli4, p_value_test = "gtest", p_value_thresh = 0.10)
  
  # There should be no p-values over 0.10
  expect_true(any(gtest_test$trelliData.stat$p_value_gtest > 0.10) == FALSE)
  
  # omicsData and statRes should have the same biomolecules 
  expect_true(all(gtest_test$trelliData.omics$Mass_Tag_ID %in% gtest_test$trelliData.stat$Mass_Tag_ID))
  
  ###################
  ## RNA-SEQ TESTS ##
  ###################
  
  # Load: seqData expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'
  ))
  
  seqData_omics <- as.seqData(e_data = edata, f_data = fdata, edata_cname = "ID_REF", fdata_cname = "Samples")
  seqData_omics <- group_designation(seqData_omics, main_effects = "Tissue")
  seqData_omics <- applyFilt(filter_object = total_count_filter(omicsData = seqData_omics), omicsData = seqData_omics, min_count = 15)
  seqData_stat <- diffexp_seq(omicsData = seqData_omics, method = "voom")
  seqTest3 <- as.trelliData(statRes = seqData_stat)
  seqTest4 <- as.trelliData(omicsData = seqData_omics, statRes = seqData_stat)
  
  # Test: trelli_pvalue_filter functionality------------------------------------
  
  # There should be a message when seqData is used letting users know that anova/gtest
  # is not an option for seqData
  expect_message(trelli_pvalue_filter(seqTest3, p_value_thresh = 0.05), "p_value_test is ignored with seqData.")
  
  # Filter by with a pvalue of 0.05
  anova_test2 <- trelli_pvalue_filter(seqTest3, comparison = "cervix_vs_uterus")
  
  # There should be no p-values over 0.05 
  expect_true(any(anova_test2$trelliData.stat$p_value > 0.05) == FALSE)
  
  # Filter with a pvalue of 0.10
  gtest_test2 <- trelli_pvalue_filter(seqTest4, p_value_thresh = 0.10)
  
  # There should be no p-values over 0.10
  expect_true(any(gtest_test2$trelliData.stat$p_value > 0.10) == FALSE)
  
  # omicsData and statRes should have the same biomolecules 
  expect_true(all(gtest_test2$trelliData.omics$ID_REF %in% gtest_test2$trelliData.stat$ID_REF))
  
})