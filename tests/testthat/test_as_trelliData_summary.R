context("summary: trelliData edata and trelliData")

test_that("trelliData object summaries return the correct data frames", {
  testthat::skip_on_cran()
  ##################
  ## MS/NMR TESTS ##
  ##################
  
  # Load: peptide expression data-----------------------------------------------

  load(system.file('testdata',
    'little_pdata.RData',
    package = 'pmartR'
  ))

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

  # Make pepTrelli objects
  pepTrelli1 <- as.trelliData.edata(e_data = edata, edata_cname = "Mass_Tag_ID", omics_type = "pepData")
  pepTrelli2 <- as.trelliData(omicsData = pepOmics)
  pepTrelli3 <- as.trelliData(statRes = pepStat)
  pepTrelli4 <- as.trelliData(omicsData = pepOmics, statRes = pepStat)

  # Make each of the four objects
  suppressWarnings({
    pepSummary1 <- pepTrelli1 %>% summary()
  })
  suppressWarnings({
    pepSummary2 <- pepTrelli2 %>% summary()
  })
  suppressWarnings({
    pepSummary3 <- pepTrelli3 %>% summary()
  })
  suppressWarnings({
    pepSummary4 <- pepTrelli4 %>% summary()
  })

  # Panel by each options
  suppressWarnings({
    pepSumEdata <- pepTrelli4 %>%
      trelli_panel_by("Mass_Tag_ID") %>%
      summary()
  })
  suppressWarnings({
    pepSumFdata <- pepTrelli4 %>%
      trelli_panel_by("SampleID") %>%
      summary()
  })
  suppressWarnings({
    pepSumEmeta <- pepTrelli4 %>%
      trelli_panel_by("Ref_ID") %>%
      summary()
  })

  # Test: summary of the trelliData objects-------------------------------------

  # Test that each have the correct number of rows
  expect_equal(5, nrow(pepSummary1))
  expect_equal(7, nrow(pepSummary2))
  expect_equal(2, nrow(pepSummary3))
  expect_equal(11, nrow(pepSummary4))
  expect_equal(4, nrow(pepSumEdata))
  expect_equal(2, nrow(pepSumFdata))
  expect_equal(5, nrow(pepSumEmeta))
  
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
  
  
  # Test: seqData expression data-----------------------------------------------
  
  seqSummary <- as.trelliData(omicsData = seqData_omics, statRes = seqData_stat) %>% summary()
  
  # No abundance plots should be suggested in summary
  expect_true(any(grepl("abundance", seqSummary$Plot)) == FALSE)
  
})
