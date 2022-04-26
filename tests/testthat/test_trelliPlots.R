context("trelliscope plotting functions")

test_that("trelliPlots check the correct inputs", {
  
  # Load: lipid expression data-------------------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
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
  
  # Create the four trelliData object testers 
  mtrelliData1 <- as.trelliData.edata(e_data = edata, edata_cname = "Metabolite", omics_type = "metabData")
  mtrelliData2 <- as.trelliData(omicsData = metabData)
  mtrelliData3 <- as.trelliData(statRes = metabStat)
  mtrelliData4 <- as.trelliData(omicsData = metabData, statRes = metabStat)
  
  # Test: trelli plotting functions---------------------------------------------
  
  ## trelli_abundance_boxplot
  
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
    "Unacceptable cognostic option included. Acceptable options are:  n, mean, median, sd, skew, p_value, fold_change"
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
  
  # Blank biomolecules should be removed, and if no data, an error produced 
  blankExample <- mtrelliData1 %>% trelli_panel_by("Metabolite")
  blankExample$trelliData.omics <- blankExample$trelliData.omics[1,]
  blankExample$trelliData.omics$Metabolite <- ""
  
  expect_error(
    suppressMessages(blankExample %>% trelli_abundance_boxplot()),
    "No data to build trelliscope with."
  )
  
  # Test mode for windows should pass 
  expect_equal(
    .getDownloadsFolder(.test_mode = TRUE),
    dirname("~")
  )
  
  # Expect a single plot object to be made 
  abun_boxplot <- mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(single_plot = TRUE)
  expect_true(inherits(abun_boxplot, "ggplot"))
  
  # Generate a tests folder
  testFolder <- file.path(.getDownloadsFolder(), "/Trelli_Tests")
  
  # If the folder exists, remove it
  if (file.exists(testFolder)) {unlink(testFolder, recursive = TRUE, force = TRUE)}
  
  # Build a trelliscope that tests multiple functions 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>% 
    trelli_abundance_boxplot(path = file.path(testFolder, "BoxAbundanceTest1"), 
                             test_mode = TRUE, 
                             test_example = 2,
                             ggplot_params = "xlab('')",
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest1")))
  
  # Build a trelliscope that tries to call stats cognostics without stats data and
  # has one plot
  singlePlot <- mtrelliData1 %>% trelli_panel_by("Metabolite") 
  singlePlot$trelliData.omics <- singlePlot$trelliData.omics[1,]
  suppressWarnings(singlePlot %>% 
   trelli_abundance_boxplot(path = file.path(testFolder, "BoxAbundanceTest2"), 
                            cognostics = "p_value")
  )
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest2")))
  
  # Build a trelliscope that tries to call stats cognostics without cognostics data 
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_abundance_boxplot(path = file.path(testFolder, "BoxAbundanceTest3"), 
                                              test_mode = TRUE, 
                                              test_example = 2,
                                              cognostics = NULL)
  )
  expect_true(file.exists(file.path(testFolder, "BoxAbundanceTest3")))
  
  ## trelli_abundance_histogram
  
  # trelliData must be grouped by edata_cname
  expect_error(
    mtrelliData2 %>% trelli_panel_by("SampleID") %>% trelli_abundance_histogram(), 
    "trelliData must be grouped by edata_cname."
  )
  
  # Expect a single plot object to be made 
  abun_histplot <- mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_histogram(single_plot = TRUE)
  expect_true(inherits(abun_histplot, "ggplot"))
  
  
  # Test that trelliscope builds when passed cognostic doesn't exist  
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>% 
    trelli_abundance_histogram(path = file.path(testFolder, "HistAbundanceTest1"),
                               test_mode = TRUE, 
                               test_example = 2,
                               ggplot_params = "xlab('')",
                               cognostics = c("n", "p_value"),
                               interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest1")))
  
  # Test that trelliscope builds even with missing cognostics and a single omic
  suppressWarnings(singlePlot %>% 
                     trelli_abundance_histogram(path = file.path(testFolder, "HistAbundanceTest2"),
                                                cognostics = NULL)
  )
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest2")))
  
  # Test that trelliscope builds with all cognostics 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_abundance_histogram(path = file.path(testFolder, "HistAbundanceTest3"),
                                                test_mode = TRUE, 
                                                test_example = 3)
  )
  expect_true(file.exists(file.path(testFolder, "HistAbundanceTest3")))
  
  ## trelli_abundance_heatmap
  
  # Test that the data has been grouped by an emeta column 
  expect_error(
    mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_heatmap(),
    "trelliData must be paneled_by an e_meta column."
  )
  
  # Return a single interactive plot 
  abun_hmplot <- mtrelliData4 %>% trelli_panel_by("MClass") %>% 
    trelli_abundance_heatmap(single_plot = TRUE, ggplot_params = "xlab('')", interactive = TRUE)
  expect_true(inherits(abun_hmplot, "plotly"))
  
  # Test that trelliscope builds with all cognostics 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
     trelli_abundance_heatmap(path = file.path(testFolder, "HmAbundanceTest1"),
                              test_mode = TRUE, 
                              test_example = 2)
  )
  expect_true(file.exists(file.path(testFolder, "HmAbundanceTest1")))
  
  # Test that the trelliscope still builds without group data and only one plot
  nogroup <- mtrelliData4 %>% trelli_panel_by("MClass")
  attributes(nogroup$omicsData)$group_DF <- NULL
  nogroup$trelliData.omics <- nogroup$trelliData.omics[1,]
  suppressWarnings(nogroup %>% 
     trelli_abundance_heatmap(path = file.path(testFolder, "HmAbundanceTest2"))
  )
  expect_true(file.exists(file.path(testFolder, "HmAbundanceTest2")))
  
  ## trelli_missingness_bar
  
  # Trigger proportion check warning with something outrageous 
  expect_error(
    mtrelliData3 %>% trelli_panel_by("Metabolite") %>% trelli_missingness_bar(proportion = "cantelope"),
    "proportion must be a TRUE or FALSE."
  )
  
  # Test trelliscope builds with all cognostics
  suppressWarnings(mtrelliData1 %>% trelli_panel_by("Metabolite") %>% 
    trelli_missingness_bar(path = file.path(testFolder, "MissingTest1"),
                           ggplot_params = "ylab('')",
                           test_mode = TRUE, 
                           test_example = 2,
                           cognostics = "proportion",
                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest1")))
  
  # Test trelliscope with just a statRes object 
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_missingness_bar(path = file.path(testFolder, "MissingTest2"),
                                            test_mode = TRUE, 
                                            test_example = 2,
                                            proportion = FALSE,
                                            cognostics = "n")
  )
  expect_true(file.exists(file.path(testFolder, "MissingTest2")))
  
  # Build a single plot
  miss_plot <- singlePlot %>% trelli_missingness_bar(single_plot = TRUE)
  expect_true(inherits(miss_plot, "ggplot"))
  
  ## trelli_foldchange_bar
  
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
  
  # Test that the p-value test is an acceptable entry of anova, g-test, both, or null
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(p_value_test = c("g-test", "anova")),
    "p_value_test must be anova, g-test, combined, or NULL."
  )
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_bar(p_value_test = c("test")),
    "p_value_test must be anova, g-test, combined, or NULL."
  )
  
  # Generate a trelliscope without proportions and a single plot
  singleStatPlot <- mtrelliData4 %>% trelli_panel_by("Metabolite")
  singleStatPlot$trelliData.omics <- singleStatPlot$trelliData.omics[1,]
  singleStatPlot$trelliData.stat <- singleStatPlot$trelliData.stat[singleStatPlot$trelliData.stat$Metabolite == singleStatPlot$trelliData.omics$Metabolite,]
  suppressWarnings(singleStatPlot %>%
     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest1"),
                           ggplot_params = "ylab('')",
                           cognostics = "p_value",
                           p_value_test = "anova",
                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest1")))
  
  # Generate a trelliscope with proportions and a g-test 
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest2"),
                                           ggplot_params = "ylab('')",
                                           test_mode = TRUE, 
                                           test_example = 2,
                                           cognostics = "p_value",
                                           p_value_test = "g-test",
                                           interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest2")))
  
  # Generate a trelliscope with proportions and combined p-values
  suppressWarnings(mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(path = file.path(testFolder, "barFoldChangeTest3"),
                                           test_mode = TRUE, 
                                           test_example = 2,
                                           p_value_test = "combined")
  )
  expect_true(file.exists(file.path(testFolder, "barFoldChangeTest3")))
  
  # Generate a trelliscope with proportions and no p-values
  fc_barplot <- mtrelliData3 %>% trelli_panel_by("Metabolite") %>% 
                     trelli_foldchange_bar(test_mode = TRUE, 
                                           test_example = 2,
                                           p_value_test = "anova",
                                           p_value_thresh = 0,
                                           single_plot = TRUE)
  expect_true(inherits(fc_barplot, "ggplot"))
  
  ## trelli_foldchange_boxplot
  
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
                             p_value_test = "anova",
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest1")))
  
  # Create a second trelliscope from the single data with no include_points and p_value in the cognostics
  singleEmetaPlot <- mtrelliData4 %>% trelli_panel_by("MClass")
  singleEmetaPlot$trelliData.omics <- singleEmetaPlot$trelliData.omics[singleEmetaPlot$trelliData.omics$MClass == "MClass5",]
  singleEmetaPlot$trelliData.stat <- singleEmetaPlot$trelliData.stat[singleEmetaPlot$trelliData.stat$MClass == "MClass5",]
  
  suppressWarnings(singleEmetaPlot %>%
                     trelli_foldchange_boxplot(path = file.path(testFolder, "boxFoldChangeTest2"),
                                               p_value_test = NULL,
                                               include_points = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "boxFoldChangeTest2")))
  
  # Generate a single plot
  fc_boxplot <- singleEmetaPlot %>% trelli_foldchange_boxplot(single_plot = T)
  expect_true(inherits(fc_boxplot, "ggplot"))
  
  ## trelli_foldchange_volcano
  
  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_volcano(),
    "trelliData must be paneled_by an e_meta column."
  )
  
  # Generate a volcano trelliscope with a modified ggplot parameter
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
   trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest1"),
                             ggplot_params = "xlab('')",
                             test_mode = TRUE, 
                             test_example = 2,
                             p_value_test = "anova",
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest1")))
  
  # Generate a single plot volcano trelliscope with no p_value
  suppressWarnings(singleEmetaPlot %>%
                     trelli_foldchange_volcano(path = file.path(testFolder, "volFoldChangeTest2"),
                                               p_value_test = NULL)
  )
  expect_true(file.exists(file.path(testFolder, "volFoldChangeTest2")))
  
  # Test the creation of a single plot
  fc_volcano <- singleEmetaPlot %>% trelli_foldchange_volcano(single_plot = T)
  expect_true(inherits(fc_boxplot, "ggplot"))
  
  ## trelli_foldchange_heatmap
  
  # Data must be grouped by an e_meta column
  expect_error(
    mtrelliData4 %>% trelli_panel_by("Metabolite") %>% trelli_foldchange_heatmap(),
    "trelliData must be paneled_by an e_meta column."
  )
  
  # Generate the fold change heatmap trelliscope with a modified ggplot parameter 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("MClass") %>% 
   trelli_foldchange_heatmap(path = file.path(testFolder, "hmFoldChangeTest1"),
                             ggplot_params = "xlab('')",
                             test_mode = TRUE, 
                             test_example = 2,
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "hmFoldChangeTest1")))
  
  # Generate the fold change heatmap trelliscope with smaller dataset
  suppressWarnings(singleEmetaPlot %>%
                     trelli_foldchange_heatmap(path = file.path(testFolder, "hmFoldChangeTest2"),
                  )
  )
  expect_true(file.exists(file.path(testFolder, "hmFoldChangeTest2")))
  
  # Generate a single plot
  fc_heatmap <- singleEmetaPlot %>% trelli_foldchange_heatmap(single_plot = T)
  expect_true(inherits(fc_heatmap, "ggplot"))
  
})  