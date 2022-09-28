context("class: trelliData edata")

test_that("as.trelliData.edata returns correct data frames and attributes",{
  
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
  
  # Now, test again with a different normalization method 
  pep_trelli_edata2 <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Mass_Tag_ID",
    omics_type = "pepData",
    data_scale_original = "abundance",
    data_scale = "log2",
    normalization_fun = "global",
    normalization_params = list(subset_fn = "ppp", norm_fn = "mean", apply_norm = TRUE,
                                params = list(ppp = 0.75), backtransform = TRUE)
  )
  
  # Test normalization parameters where normalization function should be mean, 
  # subset_function should be ppp, and the parameters list should be 0.75.
  expect_equal(attributes(pep_trelli_edata2$omicsData)$data_info$norm_info$norm_type, "global")
  expect_equal(attributes(pep_trelli_edata2$omicsData)$data_info$norm_info$norm_fn, "mean")
  expect_equal(attributes(pep_trelli_edata2$omicsData)$data_info$norm_info$subset_fn, "ppp")
  expect_equal(attributes(pep_trelli_edata2$omicsData)$data_info$norm_info$subset_params$ppp, 0.75)
  
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
  
  # Test: lipid expression data-------------------------------------------------
  
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
  
  # omicsData must be of the class lipdData
  expect_equal(class(lip_trelli_edata$omicsData), "lipidData")
  
  # omicsData must be originally abundance and transformed to log2, which is the default. 
  expect_equal(attributes(lip_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(lip_trelli_edata$omicsData)$data_info$data_scale, "log10")
  
  # The normalization should be loess, with an "affy" method and span of 0.2. 
  # The only info tracked is that normalization has happened 
  expect_true(attributes(lip_trelli_edata$omicsData)$data_info$norm_info$is_normalized)
  
  # Now, run the generic checks
  run_generic_tests(lip_trelli_edata, "LipidCommonName")
  
  # Load: metabolite expression data--------------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
  # Test: metabolite expression data--------------------------------------------
  
  # Here, we will test with ln transformation and quantile normalization. First, 
  # let's fill in the gaps in edata. 
  edata[is.na(edata)] <- rnorm(length(edata[is.na(edata)]), 20000, 1000)
  
  # Then, let's do the quantile normalization
  metab_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Metabolite",
    omics_type = "metabData",
    data_scale_original = "abundance",
    data_scale = "log",
    normalization_fun = "quantile",
    normalization_params = NULL
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 960 rows.
  expect_equal(
    colnames(metab_trelli_edata$trelliData.omics),
    c("Metabolite", "Sample", "Abundance")
  )
  
  expect_equal(nrow(metab_trelli_edata$trelliData.omics), 960)
  
  # omicsData must be of the class metabData
  expect_equal(class(metab_trelli_edata$omicsData), "metabData")
  
  # omicsData must be originally abundance and transformed to log. 
  expect_equal(attributes(metab_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(metab_trelli_edata$omicsData)$data_info$data_scale, "log")
  
  # The normalization should be quantile
  expect_equal(attributes(metab_trelli_edata$omicsData)$data_info$norm_info$norm_type, "quantile")

  # Now, run the generic checks
  run_generic_tests(metab_trelli_edata, "Metabolite")
  
  # Load: nmr expression data---------------------------------------------------
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Test: nmr expression data---------------------------------------------------
  
  # Let skip transformation. Normalization should be skipped as well because it's NMR data. 
  nmr_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Metabolite",
    omics_type = "nmrData",
    data_scale_original = "abundance",
    data_scale = "abundance"
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1558 rows.
  expect_equal(
    colnames(nmr_trelli_edata$trelliData.omics),
    c("Metabolite", "Sample", "Abundance")
  )
  
  expect_equal(nrow(nmr_trelli_edata$trelliData.omics), 1558)
  
  # omicsData must be of the class nmrData
  expect_equal(class(nmr_trelli_edata$omicsData), "nmrData")
  
  # omicsData must be originally abundance and not transformed.
  expect_equal(attributes(nmr_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(nmr_trelli_edata$omicsData)$data_info$data_scale, "abundance")
  
  # No normalization should have occurred. 
  expect_false(attributes(nmr_trelli_edata$omicsData)$data_info$norm_info$is_normalized)
  
  # Now, run the generic checks
  run_generic_tests(nmr_trelli_edata, "Metabolite")
  
  # Load: isobaric expression data----------------------------------------------
  
  load(system.file('testdata',
                   'little_isodata.RData',
                   package = 'pmartR'))
  
  # Test: isobaric expression data----------------------------------------------
  
  # Let's log transform. Normalization should be skipped as well because it's isobaric data. 
  iso_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Peptide",
    omics_type = "isobaricpepData",
    data_scale_original = "abundance",
    data_scale = "log2"
  )
  
  # Check the 4 resulting data.frames. First, trelliData.omics should exist and contain 
  # 3 colummns: edata_cname, Sample, and Abundance. It should also have 1800 rows.
  expect_equal(
    colnames(iso_trelli_edata$trelliData.omics),
    c("Peptide", "Sample", "Abundance")
  )
  
  expect_equal(nrow(iso_trelli_edata$trelliData.omics), 1800)
  
  # omicsData must be of the class nmrData
  expect_equal(class(iso_trelli_edata$omicsData), c("isobaricpepData", "pepData"))
  
  # omicsData must be originally abundance and not transformed.
  expect_equal(attributes(iso_trelli_edata$omicsData)$data_info$data_scale_orig, "abundance")
  expect_equal(attributes(iso_trelli_edata$omicsData)$data_info$data_scale, "log2")
  
  # No normalization should have occurred. 
  expect_false(attributes(iso_trelli_edata$omicsData)$data_info$norm_info$is_normalized)
  
  # Now, run the generic checks
  run_generic_tests(iso_trelli_edata, "Peptide")
  
  # Test: Input checking for is_edata-------------------------------------------
  
  # Here, I will catch specific errors with the is_edata function that didn't fit
  # as a test anywhere else in this script. 
  
  load(system.file('testdata',
                   'little_pdata.Rdata',
                   package = 'pmartR'))
  
  # Catch no edata (passed a null)
  expect_message(.is_edata(NULL), "edata must be a data.frame.")
  
  # Catch a wrong data type type 
  expect_message(.is_edata(as.matrix(edata)),  "edata must be a data.frame.")
  
  # Test: Input checking for as.trelliData.edata--------------------------------
  
  # Here, I will run specific tests that should all fail to ensure the parameter
  # validation portion of the code is running correctly. 
  
  # Anything besides edata should fail 
  expect_error(
    suppressMessages({
      as.trelliData.edata(
        e_data = fdata, 
        edata_cname = "Mass_Tag_ID",
        omics_type = "pepData"
      )
    })
  )
  
  expect_error(
    suppressMessages({
      as.trelliData.edata(
        e_data = emeta, 
        edata_cname = "Mass_Tag_ID",
        omics_type = "pepData"
      )
    })
  )
  
  # A fake omics type should trigger an error 
  expect_error(
    as.trelliData.edata(
      e_data = edata, 
      edata_cname = "Mass_Tag_ID",
      omics_type = "volitomeData"
    ),
    "volitomeData is not an acceptable omics_type."
  )
  
  # Now trigger unacceptable data_scales
  expect_error(
    as.trelliData.edata(
      e_data = edata, 
      edata_cname = "Mass_Tag_ID",
      omics_type = "pepData",
      data_scale_original = "log3"
    ),
    "log3 is not an acceptable data scale."
  )
  expect_error(
    as.trelliData.edata(
      e_data = edata, 
      edata_cname = "Mass_Tag_ID",
      omics_type = "pepData",
      data_scale = "log3"
    ),
    "log3 is not an acceptable data scale."
  )
  
  # Now trigger a wrong normalization choice
  expect_error(
    as.trelliData.edata(
      e_data = edata, 
      edata_cname = "Mass_Tag_ID",
      omics_type = "pepData",
      normalization_fun = "reference_pool"
    ),
    "reference_pool is not an acceptable normalization function type."
  )
  
  # Include a normalization example where the user doesn't apply the normalization
  expect_error(
    as.trelliData.edata(
      e_data = edata, 
      edata_cname = "Mass_Tag_ID",
      omics_type = "pepData",
      normalization_params = list(subset_fn = "all", norm_fn = "median", 
                                  apply_norm = FALSE, backtransform = TRUE)
    ),
    "apply_norm must be TRUE to apply normalization parameters."
  )
  
})