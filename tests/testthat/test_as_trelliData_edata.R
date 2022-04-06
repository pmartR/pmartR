context("class: trelliData edata")

test_that("An edata object passed to trelliData edata returns correct data frames and attributes",{
  
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.Rdata',
                   package = 'pmartR'))
  
  # Test: peptide expression data-----------------------------------------------
  
  # Here, we will test a standard case of passing an expression data function, 
  # with absolutely no changes to the default log transformations and normalizations. 
  pep_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Mass_Tag_ID",
    omics_type = "pepData"
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1800 rows.
  expect_equal(
    colnames(pep_trelli_edata$trelliData.omics),
    c("Mass_Tag_ID", "Sample", "Abundance")
  )
  
  expect_equal(nrow(pep_trelli_edata$trelliData.omics), 1800)
  
  # omicsData must be of the class pepData
  expect_equal(class(pep_trelli_edata$omicsData), "pepData")
  
  # omicsData must be originally abundance and transformed to log2, which is the default. 
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$data_scale, "log2")
  
  # The normalization should be global, with an "all" subset function and median 
  # normalization function.
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$norm_info$norm_type, "global")
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$norm_info$subset_fn, "all")
  expect_equal(attributes(pep_trelli_edata$omicsData)$data_info$norm_info$norm_fn, "median")
  
  # The rest are generic checks that are the same for every data type. Define a function to re-use. 
  run_generic_tests <- function(trelliData, edata_cname) {
    
    # Also, assert fdata of exactly the number of columns - 1, and absolutely no emeta. 
    expect_equal(nrow(trelliData$omicsData$f_data), ncol(trelliData$omicsData$e_data) - 1)
    expect_null(trelliData$omicsData$e_meta)
    
    # trelliData.stats and statRes must both be NULL since no stats were provided.
    expect_null(trelliData$trelliData.stat)
    expect_null(trelliData$statRes)
    
    # Check attributes. Fdata column should be listed as the "invented" fdata column "Sample".
    # Check potential panel_by_options which sould be edata_cname and Sample.
    # Nothing should be in the panel_by_omics or stat, since panel_by should be FALSE. 
    # The class will be both trelliData and trelli.edata.
    expect_equal(attr(trelliData, "fdata_col"), "Sample")
    expect_equal(attr(trelliData, "panel_by_options"), c(edata_cname, "Sample"))
    expect_equal(attr(trelliData, "panel_by_omics"), NA)
    expect_equal(attr(trelliData, "panel_by_stat"), NA)
    expect_false(attr(trelliData, "panel_by"))
    
  }
  
  # Run generic tests that are the same for every e_data object 
  run_generic_tests(pep_trelli_edata, "Mass_Tag_ID")
  
  # Load: protein expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_prdata.Rdata',
                   package = 'pmartR'))
  
  # Test: protein expression data-----------------------------------------------
  
  # Here, we will test a case of passing already normalized data to the as.trelliData.edata function.
  
  pro_trelli_edata <- as.trelliData.edata(
   e_data = edata, 
   edata_cname = "Reference",
   omics_type = "proData",
   data_scale_original = "log2",
   data_scale = "log2",
   is_normalized = TRUE
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1650 rows.
  expect_equal(
    colnames(pro_trelli_edata$trelliData.omics),
    c("Reference", "Sample", "Abundance")
  )
  
  expect_equal(nrow(pro_trelli_edata$trelliData.omics), 1650)
  
  # omicsData must be of the class proData
  expect_equal(class(pro_trelli_edata$omicsData), "proData")
  
  # omicsData must be originally log2 and not transformed again. 
  expect_equal(attributes(pro_trelli_edata$omicsData)$data_info$data_scale_orig, "log2")
  expect_equal(attributes(pro_trelli_edata$omicsData)$data_info$data_scale, "log2")
  
  # Normalization should not have been conducted again. 
  expect_true(attributes(pro_trelli_edata$omicsData)$data_info$norm_info$is_normalized)
  expect_null(attributes(pro_trelli_edata$omicsData)$data_info$norm_info$norm_type)
  expect_null(attributes(pro_trelli_edata$omicsData)$data_info$norm_info$subset_fn)
  expect_null(attributes(pro_trelli_edata$omicsData)$data_info$norm_info$norm_fn)
  
  # Now, run the generic checks
  run_generic_tests(pro_trelli_edata, "Reference")
  
  # Load: lipid expression data-------------------------------------------------
  
  load(system.file('testdata',
                   'lipidData.RData',
                   package = 'pmartR'))
  
  # Here, we will test with log10 transformation and loess normalization with the 
  # affy method and a span of 0.2. 
  lip_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "LipidCommonName",
    omics_type = "lipidData",
    data_scale_original = "abundance",
    data_scale = "log10",
    normalization_fun = "loess",
    normalization_params = list(method = "affy", span = 0.2)
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1606 rows.
  expect_equal(
    colnames(lip_trelli_edata$trelliData.omics),
    c("LipidCommonName", "Sample", "Abundance")
  )
  
  expect_equal(nrow(lip_trelli_edata$trelliData.omics), 1606)
  
  
  
})