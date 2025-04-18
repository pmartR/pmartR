context("trelliscope plotting functions")

test_that("trelliPlots check the correct inputs", {
  testthat::skip_on_cran()
  # Quick trelliscope function check for future changes
  #file_coverage(c("~/Git_Repos/pmartR/R/as.trelliData.R", 
  #                "~/Git_Repos/pmartR/R/summary_trelliData.R",
  #                "~/Git_Repos/pmartR/R/trelliPlots.R",
  #                "~/Git_Repos/pmartR/R/trelliPlots_seqData.R",
  #),
  #c("~/Git_Repos/pmartR/tests/testthat/test_as_trelliData_edata.R", 
  #  "~/Git_Repos/pmartR/tests/testthat/test_as_trelliData.R",
  #  "~/Git_Repos/pmartR/tests/testthat/test_trelli_pvalue_filter.R",
  #  "~/Git_Repos/pmartR/tests/testthat/test_as_trelliData_summary.R",
  #  "~/Git_Repos/pmartR/tests/testthat/test_trelliPlots.R",
  #  "~/Git_Repos/pmartR/tests/testthat/test_trelliPlots_seqData.R")
  #) %>% report()
  
  # Load: lipid expression data-------------------------------------------------

  load(system.file('testdata',
    'metaboliteData.RData',
    package = 'pmartR'
  ))

  # Generate: trelliData objects for testing------------------------------------

  # Add another condition to fdata
  fdata$Condition <- c(rep("Mock", 3), rep("InfectionA", 4), rep("InfectionB", 5))

  # Create the omics object
  metabData <- as.metabData(
    e_data = edata,
    edata_cname = "Metabolite",
    f_data = fdata,
    fdata_cname = "SampleID",
    e_meta = emeta,
    emeta_cname = "MClass"
  )

  # Log transform and set main effects
  metabData <- edata_transform(metabData, "log2")
  metabData <- group_designation(metabData, "Condition")

  # Filter for imd_anova and normalize
  suppressMessages(metabData <- applyFilt(imdanova_filter(metabData), metabData, min_nonmiss_anova = 2))
  metabData <- normalize_global(metabData, "all", "median", apply_norm = TRUE, backtransform = TRUE)

  # Now run statistics
  metabStat <- imd_anova(metabData, test_method = "combined")

  metabStat_gtest <- imd_anova(metabData, test_method = "gtest")

  # Create the four trelliData object testers
  mtrelliData1 <- as.trelliData.edata(e_data = edata, edata_cname = "Metabolite", omics_type = "metabData")
  mtrelliData2 <- as.trelliData(omicsData = metabData)
  mtrelliData3 <- as.trelliData(statRes = metabStat)
  mtrelliData4 <- as.trelliData(omicsData = metabData, statRes = metabStat)

  mtrelliData5 <- as.trelliData(omicsData = metabData, statRes = metabStat_gtest)
  
  ds_test <- data.frame(check.names = FALSE, 
    p_value_anova = 0.7,
    fold_change = NaN
  )
  expect_null(pmartR:::determine_significance(ds_test, 0.05, is_seq = FALSE))
  
  # Load RNA-seq data for testing-----------------------------------------------
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'
  ))
  
  emeta <- data.frame(check.names = FALSE, "ID_REF" = edata$ID_REF, "Classification" = paste("Class", rep(1:200, 6)))
  seqData_omics <- as.seqData(e_data = edata, f_data = fdata, edata_cname = "ID_REF", fdata_cname = "Samples", e_meta = emeta, emeta_cname = "Classification")
  seqData_omics <- group_designation(seqData_omics, main_effects = "Tissue")
  seqData_omics <- applyFilt(filter_object = total_count_filter(omicsData = seqData_omics), omicsData = seqData_omics, min_count = 15)
  seqData_stat <- diffexp_seq(omicsData = seqData_omics, method = "voom")
  seqTrelli1 <- as.trelliData.edata(e_data = edata, edata_cname = "ID_REF", omics_type = "seqData")
  seqTrelli2 <- as.trelliData(omicsData = seqData_omics)
  seqTrelli3 <- as.trelliData(statRes = seqData_stat)
  seqTrelli4 <- as.trelliData(omicsData = seqData_omics, statRes = seqData_stat)
  
  # Test: trelli plotting functions---------------------------------------------

  ## trelli_abundance_boxplot---------------------------------------------------

  # The object must be a trelliData object
  expect_error(
    trelli_abundance_boxplot(metabData),
    "trelliData must be of the class trelliData."
  )

  # The object must be paneled with trelli_panel_by
  expect_error(
    trelli_abundance_boxplot(mtrelliData1),
    "trelliData must be paneled with trelli_panel_by."
  )

  # The object must have omicsData
  expect_error(
    mtrelliData3 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(),
    "trelliData must have omicsData for this plotting function."
  )

  # Only acceptable cognostics are allowed
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(cognostics = "test"),
    "Unacceptable cognostic option included. Acceptable options are:  count, mean abundance, median abundance, cv abundance, anova p-value, fold change"
  )

  # ggplot params should be strings
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(ggplot_params = 1),
    "ggplot_params must be a string, vector of strings, or NULL."
  )

  # Interactive parameter should be a true or false
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(interactive = "Yes"),
    "interactive must be a TRUE or FALSE."
  )

  # Include points should be a true or false
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(include_points = "Yes"),
    "include_points must be a TRUE or FALSE."
  )

  # Test mode should be a true or false
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(test_mode = "test"),
    "test_mode must be a TRUE or FALSE"
  )

  # Test example should be a non-zero integer
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(test_mode = TRUE, test_example = 0),
    "test_example should be a non-zero integer."
  )

  # Test example also can't be more than the number of panels
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(test_mode = TRUE, test_example = -2000.3),
    "test_example must be in the range of possibilities, of 1 to 80"
  )

  # Single plot must be true or false
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(single_plot = "test"),
    "single_plot must be a TRUE or FALSE."
  )

  # Expect a message that median and cv abundance have been removed
  suppressWarnings({
    expect_message(
      mtrelliData2 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(test_mode = TRUE, test_example = 1, cognostics = c("median abundance", "cv abundance")),
      "'median abundance' and 'cv abundance' are not permitted when groups have been specified."
    )
  })
  
  # Expect a message that anova p-value and fold change have been removed 
  suppressWarnings({
    expect_message(
      mtrelliData4 %>% trelli_panel_by("MClass") %>% trelli_abundance_boxplot(test_mode = TRUE, test_example = 1, cognostics = c("anova p-value", "fold change")),
      "Please panel by Metabolite to get 'anova p-value' and 'fold change' as cognostics in the trelliscope display."
    )    
  })
  
  # Blank biomolecules should be removed, and if no data, an error produced
  blankExample <- mtrelliData1 %>% trelli_panel_by("Metabolite")
  blankExample$trelliData.omics <- blankExample$trelliData.omics[1, ]
  blankExample$trelliData.omics$Metabolite <- ""

  expect_error(
    suppressMessages(blankExample %>% trelli_abundance_boxplot()),
    "No data to build trelliscope with."
  )

  # Test mode for windows should pass
  expect_equal(
    .getDownloadsFolder(.test_mode = TRUE),
    file.path(dirname("~"), "Downloads")
  )

  # Expect a single plot object to be made
  abun_boxplot <- mtrelliData1 %>%
    trelli_panel_by("Metabolite") %>%
    trelli_abundance_boxplot(single_plot = TRUE)
  expect_true(inherits(abun_boxplot, "ggplot"))

  # Generate a tests folder
  testFolder <- file.path(tempdir(), "/Trelli_Tests")

  # If the folder exists, remove it
  if (file.exists(testFolder)) {
    unlink(testFolder, recursive = TRUE, force = TRUE)
  }

  # Build a trelliscope that tests multiple functions
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>%
    trelli_abundance_boxplot(
      path = file.path(testFolder, "BoxAbundanceTest1"),
      test_mode = TRUE,
      test_example = 2,
      ggplot_params = "xlab('')",
      interactive = TRUE
    ))
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest1")))

  # Build a trelliscope that tries to call stats cognostics without stats data and
  # has one plot
  singlePlot <- mtrelliData1 %>% trelli_panel_by("Metabolite") 
  singlePlot$trelliData.omics <- singlePlot$trelliData.omics[1, ]
  suppressWarnings(singlePlot %>% 
   trelli_abundance_boxplot(path = file.path(testFolder, "BoxAbundanceTest2"), 
                            cognostics = "anova p-value")
  )
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest2")))

  # Build a trelliscope that tries to call stats cognostics without cognostics data
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>%
    trelli_abundance_boxplot(
      path = file.path(testFolder, "BoxAbundanceTest3"),
      test_mode = TRUE,
      test_example = 2,
      cognostics = NULL
    ))
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest3")))
  
  # Build a trelliscope that panels by MClass
  suppressWarnings(mtrelliData5 %>% trelli_panel_by("MClass") %>%
                     trelli_abundance_boxplot(
                       path = file.path(testFolder, "BoxAbundanceTest4"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = "count"
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest4")))
  
  # Build a trelliscope that returns anova p-value and fold change only
  suppressWarnings(mtrelliData5 %>% trelli_panel_by("Metabolite") %>%
                     trelli_abundance_boxplot(
                       path = file.path(testFolder, "BoxAbundanceTest5"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = c("anova p-value", "fold change")
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest5")))
  
  # Build a trelliscope that returns anova p-value, fold change, and count
  suppressWarnings(mtrelliData5 %>% trelli_panel_by("Metabolite") %>%
                     trelli_abundance_boxplot(
                       path = file.path(testFolder, "BoxAbundanceTest6"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = c("count", "anova p-value", "fold change")
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest6")))
  
  # Seq data is not allowed for this dataset
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_abundance_boxplot(),
    "seqData is not permitted for this plotting function. Use trelli_rnaseq_boxplot instead."
  )

  ## trelli_abundance_histogram-------------------------------------------------

  # trelliData must be grouped by edata_cname
  expect_error(
    mtrelliData2 %>% trelli_panel_by("SampleID") %>% trelli_abundance_histogram(),
    "trelliData must be grouped by edata_cname."
  )

  # Expect a single plot object to be made
  abun_histplot <- mtrelliData1 %>%
    trelli_panel_by("Metabolite") %>%
    trelli_abundance_histogram(single_plot = TRUE)
  expect_true(inherits(abun_histplot, "ggplot"))
  

  # Test that trelliscope builds when passed cognostic doesn't exist  
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>% 
     trelli_abundance_histogram(path = file.path(testFolder, "HistAbundanceTest1"),
                                test_mode = TRUE, 
                               test_example = 2,
                               ggplot_params = "xlab('')",
                               cognostics = c("sample count", "cv abundance"),
                               interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest1")))

  # Test that trelliscope builds even with missing cognostics and a single omic
  suppressWarnings(singlePlot %>%
    trelli_abundance_histogram(
      path = file.path(testFolder, "HistAbundanceTest2"),
      cognostics = NULL
    ))
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest2")))

  # Test that trelliscope builds with all cognostics
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>%
    trelli_abundance_histogram(
      path = file.path(testFolder, "HistAbundanceTest3"),
      test_mode = TRUE,
      test_example = 3
    ))
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest3")))
  
  # Seq data is not allowed for this dataset
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_abundance_histogram(),
    "seqData is not permitted for this plotting function. Use trelli_rnaseq_histogram instead."
  )
  
  ## trelli_abundance_heatmap---------------------------------------------------

  # Test that the data has been grouped by an emeta column
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_heatmap(),
    "trelliData must be paneled_by an e_meta column."
  )

  # Return a single interactive plot
  abun_hmplot <- mtrelliData4 %>%
    trelli_panel_by("MClass") %>%
    trelli_abundance_heatmap(single_plot = TRUE, ggplot_params = "xlab('')", interactive = TRUE)
  expect_true(inherits(abun_hmplot, "plotly"))

  # Test that trelliscope builds with all cognostics
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>%
    trelli_abundance_heatmap(
      path = file.path(testFolder, "HmAbundanceTest1"),
      test_mode = TRUE,
      test_example = 2
    ))
  expect_true(file.exists(file.path(testFolder, "HmAbundanceTest1")))

  # Test that the trelliscope still builds without group data and only one plot
  nogroup <- mtrelliData4 %>% trelli_panel_by("MClass")
  attributes(nogroup$omicsData)$group_DF <- NULL
  nogroup$trelliData.omics <- nogroup$trelliData.omics[1, ]
  suppressWarnings(nogroup %>%
    trelli_abundance_heatmap(path = file.path(testFolder, "HmAbundanceTest2")))
  expect_true(file.exists(file.path(testFolder, "HmAbundanceTest2")))
  
  # Seq data is not allowed for this dataset
  expect_error(
    seqTrelli4 %>% trelli_panel_by("Classification") %>% trelli_abundance_heatmap(),
    "seqData is not permitted for this plotting function. Use trelli_rnaseq_heatmap instead."
  )
  
  ## trelli_missingness_bar-----------------------------------------------------

  # Trigger proportion check warning with something outrageous
  expect_error(
    mtrelliData3 %>% trelli_panel_by("Metabolite") %>% trelli_missingness_bar(proportion = "cantelope"),
    "proportion must be a TRUE or FALSE."
  )
  
  # Expect message if no p-values are provided 
  expect_message(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_missingness_bar(cognostics = "g-test p-value", test_mode = T, test_example = 1)
  )
  
  # Expect message if paneled by class
  expect_message(
    mtrelliData5 %>% trelli_panel_by("MClass") %>% trelli_missingness_bar(cognostics = "g-test p-value", test_mode = T, test_example = 1),
    "'g-test p-value' can only be included if the data has been paneled by the biomolecule column 'edata_cname'"
  )

  # Test trelliscope builds with all cognostics
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>% 
    trelli_missingness_bar(path = file.path(testFolder, "MissingTest1"),
                           ggplot_params = "ylab('')",
                           test_mode = TRUE, 
                           test_example = 2,
                           cognostics = "observed proportion",
                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest1")))
  
  # Test trelliscope with just a statRes object 
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_missingness_bar(path = file.path(testFolder, "MissingTest2"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = c("g-test p-value", "total count"))
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest2")))
  
  # Test trelliscope with paneling by a class
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
                     trelli_missingness_bar(path = file.path(testFolder, "MissingTest3"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = "observed count")
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest3")))
  
  # Test trelliscope with just a statRes object 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_missingness_bar(path = file.path(testFolder, "MissingTest4"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = c("g-test p-value"))
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest4")))

  # Build a single plot
  miss_plot <- singlePlot %>% trelli_missingness_bar(single_plot = TRUE)
  expect_true(inherits(miss_plot, "ggplot"))
  
  # Seq data is not allowed for this dataset
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_missingness_bar(),
    "seqData is not permitted for this plotting function. Use trelli_rnaseq_nonzero_bar instead."
  )

  ## trelli_foldchange_bar------------------------------------------------------

  # statRes data is required for this function
  expect_error(
    mtrelliData2 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(),
    "trelliData must have statRes for this plotting function."
  )

  # Requested test example must not be out of the range of possibilities
  expect_error(
    mtrelliData3 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(test_mode = TRUE, test_example = -2000.3),
    "test_example must be in the range of possibilities, of 1 to 65"
  )

  # trelliData must be paneled_by edata_cname
  expect_error(
    mtrelliData4 %>% trelli_panel_by("SampleID") %>% trelli_foldchange_bar(),
    "trelliData must be grouped by edata_cname."
  )

  # Test that the p-value threshold is a numeric
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(p_value_test = "anova", p_value_thresh = "test"),
    "p_value_thresh must be a numeric."
  )

  # Test that the p-value threshold is a number between 0 and 1
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(p_value_test = "anova", p_value_thresh = -0.05),
    "p_value_thresh must be between 0 and 1."
  )
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(p_value_test = "anova", p_value_thresh = 5),
    "p_value_thresh must be between 0 and 1."
  )
  
  # Generate a trelliscope without proportions and a single plot
  singleStatPlot <- mtrelliData4 %>% trelli_panel_by("Metabolite")
  singleStatPlot$trelliData.omics <- singleStatPlot$trelliData.omics[1, ]
  singleStatPlot$trelliData.stat <- singleStatPlot$trelliData.stat[singleStatPlot$trelliData.stat$Metabolite == singleStatPlot$trelliData.omics$Metabolite, ]
  suppressWarnings(singleStatPlot %>%
     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest1"),
                           ggplot_params = "ylab('')",
                           cognostics = "p-value",
                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest1")))
  
  # Generate a trelliscope with proportions and a gtest 
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest2"),
                                           ggplot_params = "ylab('')",
                                           test_mode = TRUE, 
                                           test_example = 2,
                                           cognostics = "fold change",
                                           p_value_thresh = 0,
                                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest2")))
  
  # Generate a trelliscope with proportions 
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest3"),
                                           test_mode = TRUE, 
                                           test_example = 2)
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest3")))

  # Generate a trelliscope with proportions and no p-values
  fc_barplot <- mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(test_mode = TRUE, 
                                           test_example = 2,
                                           p_value_test = "gtest",
                                           p_value_thresh = 0,
                                           single_plot = TRUE)
  expect_true(inherits(fc_barplot, "ggplot"))
  
  # Generate a trelliscope with seq data 
  suppressWarnings(seqTrelli3 %>% trelli_panel_by("ID_REF") %>% 
                     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest4"),
                                           test_mode = TRUE, 
                                           test_example = 2,
                                           cognostics = "p-value")
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest4")))

  ## trelli_foldchange_boxplot--------------------------------------------------

  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_boxplot(),
    "trelliData must be paneled_by an e_meta column."
  )

  # Include points should be a true or false
  expect_error(
    mtrelliData4 %>% trelli_panel_by("MClass") %>% trelli_foldchange_boxplot(include_points = "FALSE"),
    "include_points must be a TRUE or FALSE."
  )
  
  # Run a standard trelli foldchange boxplot with change ggplot parameters 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
   trelli_foldchange_boxplot(path = file.path(testFolder, "boxFoldChangeTest1"),
                             ggplot_params = "xlab('')",
                             test_mode = TRUE, 
                             test_example = 2,
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest1")))

  # Create a second trelliscope from the single data with no include_points and p_value in the cognostics
  singleEmetaPlot <- mtrelliData4 %>% trelli_panel_by("MClass")
  singleEmetaPlot$trelliData.omics <- singleEmetaPlot$trelliData.omics[singleEmetaPlot$trelliData.omics$MClass == "MClass5", ]
  singleEmetaPlot$trelliData.stat <- singleEmetaPlot$trelliData.stat[singleEmetaPlot$trelliData.stat$MClass == "MClass5", ]

  suppressWarnings(singleEmetaPlot %>%
                     trelli_foldchange_boxplot(path = file.path(testFolder, "boxFoldChangeTest2"),
                                               p_value_thresh = 0,
                                               include_points = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest2")))
  
  # Add RNA-seq example
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_foldchange_boxplot(path = file.path(testFolder, "boxFoldChangeTest3"),
                                               test_mode = T,
                                               test_example = 1)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest3")))
  
  # Add RNA-seq example wit no p-value threshold
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_foldchange_boxplot(path = file.path(testFolder, "boxFoldChangeTest4"),
                                               test_mode = T,
                                               test_example = 1,
                                               p_value_thresh = 0)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest4")))
  

  # Generate a single plot
  fc_boxplot <- singleEmetaPlot %>% trelli_foldchange_boxplot(single_plot = TRUE)
  expect_true(inherits(fc_boxplot, "ggplot"))

  ## trelli_foldchange_volcano--------------------------------------------------
  
  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("MClass") %>% trelli_foldchange_volcano(comparison = "Tigers"),
    "is not an acceptable comparison"
  )

  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_volcano(comparison = "Mock_vs_InfectionA"),
    "trelliData must be paneled_by an e_meta column."
  )

  # Generate a volcano trelliscope with a modified ggplot parameter
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
   trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest1"),
                             comparison = "Mock_vs_InfectionA",
                             ggplot_params = "xlab('')",
                             test_mode = TRUE, 
                             test_example = 2,
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest1")))

  # Generate a single plot volcano trelliscope with no p_value
  suppressWarnings(singleEmetaPlot %>%
                     trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest2"), comparison = "Mock_vs_InfectionA",
                                               p_value_thresh = 0)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest2")))
  
  # Add RNA-seq example
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest3"),
                                               test_mode = T,
                                               test_example = 1)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest3")))
  
  # Add RNA-seq example
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest4"),
                                               test_mode = T,
                                               test_example = 1,
                                               p_value_thresh = 0)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest4")))

  # Test the creation of a single plot
  fc_volcano <- singleEmetaPlot %>% trelli_foldchange_volcano(single_plot = TRUE, p_value_test = NULL, comparison = "Mock_vs_InfectionA")
  expect_true(inherits(fc_volcano, "ggplot"))
  
  # Expect error if trying to return a single plot without selecting a comparison
  expect_error(
    mtrelliData4 %>% trelli_panel_by("MClass") %>% trelli_foldchange_volcano(single_plot = T),
    "single_plot will only work if 1 comparison has been selected."
  )
  
  # Expect error if comparison is not a string
  expect_error(
    mtrelliData4 %>% trelli_panel_by("MClass") %>% trelli_foldchange_volcano(single_plot = T, comparison = 2),
    "comparison must be a string."
  )

  ## trelli_foldchange_heatmap--------------------------------------------------

  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_heatmap(),
    "trelliData must be paneled_by an e_meta column."
  )

  # Generate the fold change heatmap trelliscope with a modified ggplot parameter
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>%
    trelli_foldchange_heatmap(
      path = file.path(testFolder, "hmFoldChangeTest1"),
      ggplot_params = "xlab('')",
      test_mode = TRUE,
      test_example = 2,
      interactive = TRUE
  ))
  expect_true(file.exists(file.path(testFolder, "hmFoldChangeTest1")))

  # Generate the fold change heatmap trelliscope with smaller dataset
  suppressWarnings(singleEmetaPlot %>%
    trelli_foldchange_heatmap(path = file.path(testFolder, "hmFoldChangeTest2"), ))
  expect_true(file.exists(file.path(testFolder, "hmFoldChangeTest2")))
  
  # Add RNA-seq example
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_foldchange_heatmap(path = file.path(testFolder, "hmFoldChangeTest3"),
                                               test_mode = T,
                                               p_value_thresh = 0,
                                               test_example = 1)
  )
  expect_true(file.exists(file.path(testFolder, "hmFoldChangeTest3")))

  # Generate a single plot
  fc_heatmap <- singleEmetaPlot %>% trelli_foldchange_heatmap(single_plot = TRUE)
  expect_true(inherits(fc_heatmap, "ggplot"))
  
})  
