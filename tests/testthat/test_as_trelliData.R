context("class: trelliData & trelli_panel_by")

test_that("as.trelliData and trelli_panel_by returns correct data frames and attributes", {
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

  # Test: peptide expression data-----------------------------------------------

  # There are 3 possible trelliData objects that can be built here: one with
  # only omicsData, one with only statRes, and one with both.
  pepTrelli_omics <- as.trelliData(omicsData = pepOmics)
  pepTrelli_stat <- as.trelliData(statRes = pepStat)
  pepTrelli_both <- as.trelliData(omicsData = pepOmics, statRes = pepStat)

  ## pepOmic only checks

  # trelliData should have just the omicsData
  expect_equal(pepTrelli_omics$omicsData, pepOmics)
  expect_null(pepTrelli_omics$trelliData.stat)
  expect_null(pepTrelli_omics$statRes)

  # trelliData.omics should have edata_cname, fdata_cname, abundance, the main effects, and all emeta categories
  edata_cname <- get_edata_cname(pepOmics)
  fdata_cname <- get_fdata_cname(pepOmics)
  emeta_cols <- colnames(emeta)[colnames(emeta) != edata_cname]
  expect_equal(colnames(pepTrelli_omics$trelliData.omics), c(edata_cname, fdata_cname, "Abundance", "Group", emeta_cols))

  # Check trelliData attributes. The only fdata column should always be
  # fdata_cname. Emeta columns should include everything but the edata_cname.
  # Panel by options will include emeta_cname, and the aforementioned columns.
  # All panel_by properties should be FALSE.
  expect_equal(attr(pepTrelli_omics, "fdata_col"), fdata_cname)
  expect_equal(attr(pepTrelli_omics, "emeta_col"), emeta_cols)
  expect_equal(attr(pepTrelli_omics, "panel_by_options"), c(edata_cname, fdata_cname, emeta_cols))
  expect_false(attr(pepTrelli_omics, "panel_by"))
  expect_equal(attr(pepTrelli_omics, "panel_by_omics"), NA)
  expect_equal(attr(pepTrelli_omics, "panel_by_stat"), NA)

  # Test some trelli_panel_by options
  pepTrelli_omics_edatacname <- trelli_panel_by(pepTrelli_omics, edata_cname)
  pepTrelli_omics_fdatacname <- trelli_panel_by(pepTrelli_omics, fdata_cname)
  pepTrelli_omics_emetacol <- trelli_panel_by(pepTrelli_omics, emeta_cols[1])

  # After paneling by an acceptable choice, the number of plots should equal the
  # number of unique entries in that category.
  expect_equal(
    pepTrelli_omics_edatacname$trelliData.omics %>% nrow(),
    pepTrelli_omics$trelliData.omics[[edata_cname]] %>% unique() %>% length()
  )

  expect_equal(
    pepTrelli_omics_fdatacname$trelliData.omics %>% nrow(),
    pepTrelli_omics$trelliData.omics[[fdata_cname]] %>% unique() %>% length()
  )

  expect_equal(
    pepTrelli_omics_emetacol$trelliData.omics %>% nrow(),
    pepTrelli_omics$trelliData.omics[[emeta_cols[1]]] %>% unique() %>% length()
  )

  # Check that the correct attributes were changed after paneling. Only "panel_by_omics"
  # should contain the panel_by variable. The rest (besides panel_by flag) should be the same.
  expect_true(attr(pepTrelli_omics_edatacname, "panel_by"))
  expect_true(attr(pepTrelli_omics_fdatacname, "panel_by"))
  expect_true(attr(pepTrelli_omics_emetacol, "panel_by"))
  expect_equal(attr(pepTrelli_omics_edatacname, "panel_by_omics"), edata_cname)
  expect_equal(attr(pepTrelli_omics_fdatacname, "panel_by_omics"), fdata_cname)
  expect_equal(attr(pepTrelli_omics_emetacol, "panel_by_omics"), emeta_cols[1])
  expect_equal(attr(pepTrelli_omics_edatacname, "panel_by_stat"), NA)
  expect_equal(attr(pepTrelli_omics_fdatacname, "panel_by_stat"), NA)
  expect_equal(attr(pepTrelli_omics_emetacol, "panel_by_stat"), NA)
  expect_equal(
    attributes(pepTrelli_omics)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_omics_edatacname)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )
  expect_equal(
    attributes(pepTrelli_omics_fdatacname)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_omics_emetacol)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )

  ## pepStat only checks

  # trelliData should have just the statRes
  expect_equal(pepTrelli_stat$statRes, pepStat)
  expect_null(pepTrelli_stat$trelliData.omics)
  expect_null(pepTrelli_stat$omicsData)

  # Without an emeta file, the options with just statRes are limited. The only
  # columns should be edata_cname, Comparison, p_value columns, and fold_change.
  expect_equal(
    colnames(pepTrelli_stat$trelliData.stat),
    c(edata_cname, "Comparison", "p_value_anova", "p_value_gtest", "fold_change")
  )

  # Check trelliData attributes. The only panel_by option should be edata_cname.
  # "panel_by" should be FALSE and no panelling variables listed.
  expect_equal(attr(pepTrelli_stat, "panel_by_options"), edata_cname)
  expect_false(attr(pepTrelli_stat, "panel_by"))
  expect_equal(attr(pepTrelli_stat, "panel_by_omics"), NA)
  expect_equal(attr(pepTrelli_stat, "panel_by_stat"), NA)

  # Test the only panel by option
  pepTrelli_stat_edatacname <- trelli_panel_by(pepTrelli_stat, edata_cname)

  # After paneling by an acceptable choice, the number of plots should equal the
  # number of unique entries in that category.
  expect_equal(
    pepTrelli_stat_edatacname$trelliData.stat %>% nrow(),
    pepTrelli_stat$trelliData.stat[[edata_cname]] %>% unique() %>% length()
  )

  # Check that the correct attributes were changed after paneling. Only "panel_by_stat"
  # should contain the panel_by variable. The rest (besides panel_by flag) should be the same.
  expect_true(attr(pepTrelli_stat_edatacname, "panel_by"))
  expect_equal(attr(pepTrelli_stat_edatacname, "panel_by_omics"), NA)
  expect_equal(attr(pepTrelli_stat_edatacname, "panel_by_stat"), edata_cname)
  expect_equal(
    attributes(pepTrelli_stat)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_stat_edatacname)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )

  ## Both omicsData and statRes checks

  # trelliData should have both omicsData and statRes
  expect_equal(pepTrelli_both$omicsData, pepOmics)
  expect_equal(pepTrelli_both$statRes, pepStat)

  # Now, check the trelliData.omics object. Should be the same as when we just supplied omicsData.
  expect_equal(colnames(pepTrelli_both$trelliData.omics), c(edata_cname, fdata_cname, "Abundance", "Group", emeta_cols))

  # The trelliData.stat object should hav ethe same categories as when we just supplied
  # statRes, with the emeta cols.
  expect_equal(
    colnames(pepTrelli_both$trelliData.stat),
    c(edata_cname, "Comparison", "p_value_anova", "p_value_gtest", "fold_change", emeta_cols)
  )

  # Check trelliData attributes. The only fdata column should always be
  # fdata_cname. Emeta columns should include everything but the edata_cname.
  # Panel by options will include emeta_cname, and the aforementioned columns.
  # All panel_by properties should be FALSE.
  expect_equal(attr(pepTrelli_both, "fdata_col"), fdata_cname)
  expect_equal(attr(pepTrelli_both, "emeta_col"), emeta_cols)
  expect_equal(attr(pepTrelli_both, "panel_by_options"), c(edata_cname, fdata_cname, emeta_cols))
  expect_false(attr(pepTrelli_both, "panel_by"))
  expect_equal(attr(pepTrelli_both, "panel_by_omics"), NA)
  expect_equal(attr(pepTrelli_both, "panel_by_stat"), NA)

  # Test the only panel by edata_cname, fdata_cname, and an emeta_coumns
  pepTrelli_both_edatacname <- trelli_panel_by(pepTrelli_both, edata_cname)
  pepTrelli_both_fdatacname <- trelli_panel_by(pepTrelli_both, fdata_cname)
  pepTrelli_both_emetacol <- trelli_panel_by(pepTrelli_both, emeta_cols[2])

  # After paneling by an acceptable choice, the number of plots should equal the
  # number of unique entries in that category across omicsData and statRes.

  # Number of biomolecules should be the same in both omicsData and statRes
  expect_equal(
    pepTrelli_both_edatacname$trelliData.omics %>% nrow(),
    pepTrelli_both$trelliData.omics[[edata_cname]] %>% unique() %>% length()
  )
  expect_equal(
    pepTrelli_both_edatacname$trelliData.stat %>% nrow(),
    pepTrelli_both$trelliData.stat[[edata_cname]] %>% unique() %>% length()
  )

  # Only omicsData can be paneled_by sample
  expect_equal(
    pepTrelli_both_fdatacname$trelliData.omics %>% nrow(),
    pepTrelli_both$trelliData.omics[[fdata_cname]] %>% unique() %>% length()
  )

  # Number of emeta categories should be the same in both omicsData and statRes
  expect_equal(
    pepTrelli_both_emetacol$trelliData.omics %>% nrow(),
    pepTrelli_both$trelliData.omics[[emeta_cols[2]]] %>% unique() %>% length()
  )
  expect_equal(
    pepTrelli_both_emetacol$trelliData.stat %>% nrow(),
    pepTrelli_both$trelliData.stat[[emeta_cols[2]]] %>% unique() %>% length()
  )

  # Check that the correct attributes were changed after paneling. Only "panel_by_omics"
  # should contain the panel_by variable. The rest (besides panel_by flag) should be the same.
  expect_true(attr(pepTrelli_both_edatacname, "panel_by"))
  expect_equal(attr(pepTrelli_both_edatacname, "panel_by_omics"), edata_cname)
  expect_equal(attr(pepTrelli_both_edatacname, "panel_by_stat"), edata_cname)
  expect_equal(
    attributes(pepTrelli_both)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_both_edatacname)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )

  expect_true(attr(pepTrelli_both_fdatacname, "panel_by"))
  expect_equal(attr(pepTrelli_both_fdatacname, "panel_by_omics"), fdata_cname)
  expect_equal(attr(pepTrelli_both_fdatacname, "panel_by_stat"), NA)
  expect_equal(
    attributes(pepTrelli_both)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_both_fdatacname)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )

  expect_true(attr(pepTrelli_both_emetacol, "panel_by"))
  expect_equal(attr(pepTrelli_both_emetacol, "panel_by_omics"), emeta_cols[2])
  expect_equal(attr(pepTrelli_both_emetacol, "panel_by_stat"), emeta_cols[2])
  expect_equal(
    attributes(pepTrelli_both)[c("fdata_col", "emeta_col", "panel_by_options", "class")],
    attributes(pepTrelli_both_emetacol)[c("fdata_col", "emeta_col", "panel_by_options", "class")]
  )

  # Test: Input checking for as.trelliData--------------------------------------

  # Here, I will run specific tests that should all fail to ensure the parameter
  # validation portion of the code is running correctly

  # No data should return an error
  expect_error(as.trelliData(), "At least 1 omicsData or 1 statRes object must be provided.")

  # Test a wrong omics class
  wrongOmic <- pepOmics
  class(wrongOmic) <- "volitalomeData"
  expect_error(
    as.trelliData(omicsData = wrongOmic),
    "volitalomeData is not a supported omicsData class."
  )

  # Test a non-log transformed omicsData
  pepNotTransformed <- as.pepData(
    e_data = edata,
    f_data = fdata,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID"
  )
  expect_error(
    as.trelliData(omicsData = pepNotTransformed),
    "omicsData must be log transformed."
  )

  # Test a non-normalized omicsData
  pepNotNormalized <- edata_transform(pepNotTransformed, "log2")
  expect_error(
    as.trelliData(omicsData = pepNotNormalized),
    "omicsData must be normalized."
  )

  # Create a group designation with singletons to get singleton warning
  fdata_single <- fdata
  fdata_single[12, 2] <- "Mock2"
  pepSingle <- as.pepData(
    e_data = edata,
    f_data = fdata_single,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID"
  )
  pepSingle <- edata_transform(pepSingle, "log2")
  pepSingle <- group_designation(pepSingle, "Condition")
  pepSingle <- normalize_global(
    omicsData = pepSingle,
    subset_fn = "all",
    norm_fn = "median",
    apply_norm = TRUE
  )
  expect_warning(
    as.trelliData(omicsData = pepSingle),
    "singleton groups found in group_designation: Mock2"
  )

  # Pass a wrong object to statRes
  expect_error(
    as.trelliData(statRes = pepOmics),
    "statRes must be an object of the statRes class. See ?imd_anova.",
    fixed = T
  )

  # Mismatch an omics and statRes dataset
  load(system.file('testdata',
    'nmrData.RData',
    package = 'pmartR'
  ))
  nmrData <- as.nmrData(
    e_data = edata,
    f_data = fdata,
    edata_cname = "Metabolite",
    fdata_cname = "SampleID"
  )

  expect_error(
    as.trelliData(omicsData = nmrData, statRes = pepStat),
    "omicsData and statRes are from different datasets."
  )

  # Test: Input checking for trelli_panel_by------------------------------------

  # Here, I will run specific tests that should all fail to ensure the parameter
  # validation portion of the code is running correctly

  # Only trelliData objects can be paneled
  expect_error(
    trelli_panel_by(trelliData = pepOmics, panel = "Mass_Tag_ID"),
    "trelliData must be of the class trelliData from as.trelliData or as.trelliData.edata."
  )

  # trelliData objects can only be paneled by variables in the list
  expect_error(
    trelli_panel_by(trelliData = pepTrelli_both, panel = "Condition"),
    "panel is not an acceptable option. The following can be selected: Mass_Tag_ID, SampleID, Protein, Ref_ID, Peptide_Sequence"
  )

  # Try to panel an already paneled object
  pepTrelliPanel <- trelli_panel_by(pepTrelli_both, panel = "SampleID")
  expect_error(
    trelli_panel_by(pepTrelliPanel, "Mass_Tag_ID"),
    "trelliData has already been paneled by SampleID"
  )

  pepTrelliStatPanel <- trelli_panel_by(pepTrelli_stat, panel = "Mass_Tag_ID")
  expect_error(
    trelli_panel_by(pepTrelliStatPanel, "Mass_Tag_ID"),
    "trelliData has already been paneled by Mass_Tag_ID"
  )

  # Finally, test the warning with the small groups
  tinyNMR <- as.nmrData(e_data = edata[1, 1:2], f_data = fdata[1, ], edata_cname = "Metabolite", fdata_cname = "SampleID")

  expect_warning(
    as.trelliData(tinyNMR) %>% trelli_panel_by("Metabolite"),
    "Grouping by Metabolite results in panels with less than 3 data points in trelliData.omics."
  )
})

test_that("as.trelliData and trelli_panel_by returns correct data frames and attributes for RNA-seq data", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("limma")
  ###################
  ## RNA SEQ TESTS ##
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
  
  # Expect message
  expect_message(as.trelliData(omicsData = seqData_omics))
  
  # Save object
  seqTest2 <- as.trelliData(omicsData = seqData_omics)
  
  # Expect message
  expect_message(as.trelliData(statRes = seqData_stat))
  
  # Save object
  seqTest3 <- as.trelliData(statRes = seqData_stat)
  
  # Expect message
  expect_message(as.trelliData(omicsData = seqData_omics, statRes = seqData_stat))
  
  # Save object
  seqTest4 <- as.trelliData(omicsData = seqData_omics, statRes = seqData_stat)
  
  # Check all object types
  expect_true(inherits(seqTest2, "trelliData.seqData"))
  expect_true(inherits(seqTest3, "trelliData.seqData"))
  expect_true(inherits(seqTest4, "trelliData.seqData"))
  
})
