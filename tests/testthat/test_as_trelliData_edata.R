context("class: trelliData edata")

test_that("An edata object passed to trelliData edata returns correct data frames and attributes",{
  
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.Rdata',
                   package = 'pmartR'))
  
  # Test: peptide expression data-----------------------------------------------
  
  pep_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Mass_Tag_ID",
    omics_type = "pepData"
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1800 rows.
  edata_cname <- "Mass_Tag_ID"
  
  expect_equal(
    colnames(pep_trelli_edata$trelliData.omics),
    c(edata_cname, "Sample", "Abundance")
  )
  
  expect_equal(1800, nrow(pep_trelli_edata$trelliData.omics))
  
  # omicsData must be of the class pepData
  expect_equal(class(pep_trelli_edata$omicsData), "pepData")
  
  # omicsData must be originally abundance and transformed to log2, which is the default. 
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$data_scale, "log2")
  
  # The normalization should be global, with an "all" subset function and median 
  # normalization function.
  

  # Also, assert fdata of exactly the number of columns - 1, and absolutely no emeta. 
  expect_equal(nrow(pep_trelli_edata$omicsData$f_data), 12)
  expect_null(pep_trelli_edata$omicsData$e_meta)
  
  # trelliData.stats and statRes must both be NULL since no stats were provided.
  expect_null(pep_trelli_edata$trelliData.stat)
  expect_null(pep_trelli_edata$statRes)
  
  # Check attributes. Fdata column should be listed as the "invented" fdata column.
  # Check potential panel_by_options. Nothing should be in the panel_by_omics or stat,
  # since panel_by should be FALSE. The class will be both trelliData and trelli.edata.
  
  expect_equal(attr(pep_trelli_edata, "fdata_col"), "Sample")
  expect_equal(attr(pep_trelli_edata, "panel_by_options"), c(edata_cname, "Sample"))
  expect_equal(attr(pep_trelli_edata, "panel_by_omics"), NA)
  expect_equal(attr(pep_trelli_edata, "panel_by_stat"), NA)
  expect_false(attr(pep_trelli_edata, "panel_by"))
  
  
  
})