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
  
  ## Let's try the trelli abundance boxplot
  
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
  
  # Expect a single plot object to be made 
  plot <- mtrelliData1 %>% trelli_panel_by("Metabolite") %>% trelli_abundance_boxplot(single_plot = TRUE)
  expect_true(inherits(plot, "ggplot"))
  
  # Generate a tests folder
  testFolder <- file.path(getDownloadsFolder(), "/Trelli_Tests")
  
  # If the folder exists, remove it
  if (file.exists(testFolder)) {file.remove(testFolder)}
  
  # Make two example trelliscopes, 
  suppressWarnings(mtrelliData4 %>% trelli_panel_by("Metabolite") %>% 
    trelli_abundance_boxplot(path = file.path(testFolder, "AbundanceTest"), 
                             test_mode = TRUE, 
                             test_example = 2,
                             ggplot_params = "xlab('')",
                             interactive = TRUE)
  )
  expect_true(file.exists(file.path(testFolder, "AbundanceTest")))
  
  
  
  
  
})  