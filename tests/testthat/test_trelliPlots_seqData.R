context("trelliscope plotting functions for seqData")

test_that("trelliPlots check the correct inputs for seqData", {
  testthat::skip_on_cran()
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'
  ))
  
  # Test: peptide expression data-----------------------------------------------
  
  # Here, we will test a standard case of passing an expression data function,
  # with absolutely no changes to the default log transformations and normalizations.
  pep_trelli_edata <- as.trelliData.edata(
    e_data = edata,
    edata_cname = "Mass_Tag_ID",
    omics_type = "pepData"
  )
  
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
  
  # Test: trelli rnaseq functions-----------------------------------------------
  
  ## trelli_rnaseq_boxplot------------------------------------------------------
  
  # seqData required
  expect_error(
    trelli_rnaseq_boxplot(pep_trelli_edata),
    "seqData is required for this plotting function. Use trelli_abundance_boxplot instead."
  )
  
  # The object must be a trelliData object
  expect_error(
    trelli_rnaseq_boxplot(seqData_omics),
    "seqData is required for this plotting function. Use trelli_abundance_boxplot instead."
  )
  
  # The object must be paneled with trelli_panel_by
  expect_error(
    trelli_rnaseq_boxplot(seqTrelli1),
    "trelliData must be paneled with trelli_panel_by."
  )
  
  # The object must have omicsData
  expect_error(
    seqTrelli3 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(),
    "trelliData must have omicsData for this plotting function."
  )
  
  # Only acceptable cognostics are allowed
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(cognostics = "test"),
    "Unacceptable cognostic option included. Acceptable options are:  count, mean lcpm, median lcpm, cv lcpm, p-value, fold change"
  )
  
  # ggplot params should be strings
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(ggplot_params = 1),
    "ggplot_params must be a string, vector of strings, or NULL."
  )
  
  # Interactive parameter should be a true or false
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(interactive = "Yes"),
    "interactive must be a TRUE or FALSE."
  )
  
  # Include points should be a true or false
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(include_points = "Yes"),
    "include_points must be a TRUE or FALSE."
  )
  
  # Test mode should be a true or false
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(test_mode = "test"),
    "test_mode must be a TRUE or FALSE"
  )
  
  # Test example should be a non-zero integer
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(test_mode = TRUE, test_example = 0),
    "test_example should be a non-zero integer."
  )
  
  # Test example also can't be more than the number of panels
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(test_mode = TRUE, test_example = -2000.3),
    "test_example must be in the range of possibilities, of 1 to 1200"
  )
  
  # Single plot must be true or false
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_boxplot(single_plot = "test"),
    "single_plot must be a TRUE or FALSE."
  )
  
  # Expect a message that median and cv abundance have been removed
  suppressWarnings({
    expect_message(
      seqTrelli2 %>% trelli_panel_by("ID_REF") %>% 
        trelli_rnaseq_boxplot(
          test_mode = TRUE, test_example = 1, 
          cognostics = c("median lcpm", "cv lcpm"), path = tempdir()
        ),
      "'median lcpm' and 'cv lcpm' are not permitted when groups have been specified."
    )
  })
  
  # Expect a message that p-value and fold change have been removed 
  suppressWarnings({
    expect_message(
      seqTrelli4 %>% trelli_panel_by("Classification") %>% 
        trelli_rnaseq_boxplot(
          test_mode = TRUE, test_example = 1, 
          cognostics = c("p-value", "fold change"), path = tempdir()
        ),
      "Please panel by ID_REF to get 'p-value' and 'fold change' as cognostics in the trelliscope display."
    )    
  })
  
  # Blank biomolecules should be removed, and if no data, an error produced
  blankExample <- seqTrelli1 %>% trelli_panel_by("ID_REF")
  blankExample$trelliData.omics <- blankExample$trelliData.omics[1, ]
  blankExample$trelliData.omics$ID_REF <- ""
  
  expect_error(
    suppressMessages(blankExample %>% trelli_rnaseq_boxplot(path = tempdir())),
    "No data to build trelliscope with."
  )
  
  # Test mode for windows should pass
  expect_equal(
    .getDownloadsFolder(.test_mode = TRUE),
    file.path(dirname("~"), "Downloads")
  )
  
  # Expect a single plot object to be made
  rna_boxplot <- seqTrelli1 %>%
    trelli_panel_by("ID_REF") %>%
    trelli_rnaseq_boxplot(single_plot = TRUE, path = tempdir())
  expect_true(inherits(rna_boxplot, "ggplot"))
  
  # Generate a tests folder
  testFolder <- file.path(tempdir(), "/Trelli_Tests")
  
  # If the folder exists, remove it
  if (file.exists(testFolder)) {
    unlink(testFolder, recursive = TRUE, force = TRUE)
  }
  
  # Build a trelliscope that tests multiple functions
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("ID_REF") %>%
                     trelli_rnaseq_boxplot(
                       path = file.path(testFolder, "BoxSeqTest1"),
                       test_mode = TRUE,
                       test_example = 2,
                       ggplot_params = "xlab('')",
                       interactive = TRUE
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest1")))
  
  # Build a trelliscope that tries to call stats cognostics without stats data and
  # has one plot
  singlePlot <- seqTrelli1 %>% trelli_panel_by("ID_REF") 
  singlePlot$trelliData.omics <- singlePlot$trelliData.omics[1, ]
  suppressWarnings(singlePlot %>% 
                     trelli_rnaseq_boxplot(path = file.path(testFolder, "BoxSeqTest2"), 
                                              cognostics = "p-value")
  )
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest2")))
  
  # Build a trelliscope that tries to call stats cognostics without cognostics data
  suppressWarnings(seqTrelli1 %>% trelli_panel_by("ID_REF") %>%
                     trelli_rnaseq_boxplot(
                       path = file.path(testFolder, "BoxSeqTest3"),
                       test_mode = TRUE,
                       test_example = 2,
                       cognostics = NULL
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest3")))
  
  # Build a trelliscope that panels by Classification
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_rnaseq_boxplot(
                       path = file.path(testFolder, "BoxSeqTest4"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = "count"
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest4")))
  
  # Build a trelliscope that returns  p-value and fold change only
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("ID_REF") %>%
                     trelli_rnaseq_boxplot(
                       path = file.path(testFolder, "BoxSeqTest5"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = c("p-value", "fold change")
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest5")))
  
  # Build a trelliscope that returns p-value, fold change, and count
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("ID_REF") %>%
                     trelli_rnaseq_boxplot(
                       path = file.path(testFolder, "BoxSeqTest6"),
                       test_mode = TRUE,
                       test_example = 1,
                       cognostics = c("count", "p-value", "fold change")
                     ))
  expect_true(file.exists(file.path(testFolder, "BoxSeqTest6")))
  
  ## trelli_rnaseq_histogram-------------------------------------------------
  
  # trelliData must be grouped by edata_cname
  expect_error(
    seqTrelli2 %>% trelli_panel_by("Samples") %>% trelli_rnaseq_histogram(path = tempdir()),
    "trelliData must be grouped by edata_cname."
  )
  
  # Expect a single plot object to be made
  abun_histplot <- seqTrelli1 %>%
    trelli_panel_by("ID_REF") %>%
    trelli_rnaseq_histogram(single_plot = TRUE, path = tempdir())
  expect_true(inherits(abun_histplot, "ggplot"))
  
  
  # Test that trelliscope builds when passed cognostic doesn't exist  
  suppressWarnings(seqTrelli1 %>% trelli_panel_by("ID_REF") %>% 
                     trelli_rnaseq_histogram(path = file.path(testFolder, "HistSeqTest1"),
                                                test_mode = TRUE, 
                                                test_example = 2,
                                                ggplot_params = "xlab('')",
                                                cognostics = c("sample count", "cv lcpm"),
                                                interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "HistSeqTest1")))
  
  # Test that trelliscope builds even with missing cognostics and a single omic
  suppressWarnings(singlePlot %>%
                     trelli_rnaseq_histogram(
                       path = file.path(testFolder, "HistSeqTest2"),
                       cognostics = NULL
                     ))
  expect_true(file.exists(file.path(testFolder, "HistSeqTest2")))
  
  # Test that trelliscope builds with all cognostics
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("ID_REF") %>%
                     trelli_rnaseq_histogram(
                       path = file.path(testFolder, "HistSeqTest3"),
                       test_mode = TRUE,
                       test_example = 3
                     ))
  expect_true(file.exists(file.path(testFolder, "HistSeqTest3")))
  
  ## trelli_rnaseq_heatmap---------------------------------------------------
  
  # Test that the data has been grouped by an emeta column
  expect_error(
    seqTrelli1 %>% trelli_panel_by("ID_REF") %>% trelli_rnaseq_heatmap(path = tempdir()),
    "trelliData must be paneled_by an e_meta column."
  )
  
  # Return a single interactive plot
  abun_hmplot <- seqTrelli4 %>%
    trelli_panel_by("Classification") %>%
    trelli_rnaseq_heatmap(
      single_plot = TRUE, ggplot_params = "xlab('')", 
      interactive = TRUE, path = tempdir()
    )
  expect_true(inherits(abun_hmplot, "plotly"))
  
  # Test that trelliscope builds with all cognostics
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>%
                     trelli_rnaseq_heatmap(
                       path = file.path(testFolder, "HmSeqTest1"),
                       test_mode = TRUE,
                       test_example = 2
                     ))
  expect_true(file.exists(file.path(testFolder, "HmSeqTest1")))
  
  # Test that the trelliscope still builds without group data and only one plot
  nogroup <- seqTrelli4 %>% trelli_panel_by("Classification")
  attributes(nogroup$omicsData)$group_DF <- NULL
  nogroup$trelliData.omics <- nogroup$trelliData.omics[1, ]
  suppressWarnings(nogroup %>%
                     trelli_rnaseq_heatmap(path = file.path(testFolder, "HmSeqTest2")),)
  expect_true(file.exists(file.path(testFolder, "HmSeqTest2")))
  
  ## trelli_rnaseq_nonzero_bar--------------------------------------------------
  
  # Trigger proportion check warning with something outrageous
  expect_error(
    seqTrelli3 %>% trelli_panel_by("ID_REF") %>% 
    trelli_rnaseq_nonzero_bar(proportion = "cantelope", path = tempdir()),
    "proportion must be a TRUE or FALSE."
  )
  
  # Test trelliscope builds with all cognostics
  suppressWarnings(seqTrelli1 %>% trelli_panel_by("ID_REF") %>% 
                     trelli_rnaseq_nonzero_bar(path = file.path(testFolder, "NonzeroTest1"),
                                            ggplot_params = "ylab('')",
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            cognostics = "non-zero proportion",
                                            interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "NonzeroTest1")))
  
  # Test trelliscope with just a statRes object 
  suppressWarnings(seqTrelli3 %>% trelli_panel_by("ID_REF") %>% 
                     trelli_rnaseq_nonzero_bar(path = file.path(testFolder, "NonzeroTest2"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = "total count")
  )
  expect_true(file.exists(file.path(testFolder, "NonzeroTest2")))
  
  # Test trelliscope with paneling by a class
  suppressWarnings(seqTrelli4 %>% trelli_panel_by("Classification") %>% 
                     trelli_rnaseq_nonzero_bar(path = file.path(testFolder, "NonzeroTest3"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = "non-zero count")
  )
  expect_true(file.exists(file.path(testFolder, "NonzeroTest3")))
  
  # Test trelliscope with just a statRes object 
  suppressWarnings(seqTrelli3 %>% trelli_panel_by("ID_REF") %>% 
                     trelli_rnaseq_nonzero_bar(path = file.path(testFolder, "NonzeroTest4"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = NULL)
  )
  expect_true(file.exists(file.path(testFolder, "NonzeroTest4")))
  
  # Build a single plot
  miss_plot <- singlePlot %>% trelli_rnaseq_nonzero_bar(single_plot = TRUE, path = tempdir())
  expect_true(inherits(miss_plot, "ggplot"))

  
})
