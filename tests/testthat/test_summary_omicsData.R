context('summary: omicsData')

test_that('the omicsData summaries are all square', {
  # Create pepData objects to test with and without group info -----------------

  load(system.file('testdata',
    'little_pdata.RData',
    package = 'pmartR'
  ))

  pdata <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Mass_Tag_ID',
    fdata_cname = 'SampleID',
    emeta_cname = 'Protein'
  )

  pdata_gdf <- group_designation(
    omicsData = pdata,
    main_effects = 'Condition'
  )

  # Create pepData objects to test multiple main effects and covariates --------

  # Copy edata so the names of the samples can be changed.
  edata2_2 <- edata

  # Change some of the Infection samples to Mock samples.
  names(edata2_2) <- c(
    "Mass_Tag_ID",
    paste0("Infection", 1:6),
    paste0("Mock", 1:6)
  )

  # Create additional f_data objects with different main effects and covariates.
  fdata2_2 <- fdata

  # Update the sample names in f_data.
  fdata2_2$SampleID <- c(
    paste0("Infection", 1:6),
    paste0("Mock", 1:6)
  )

  # Update the first main effect to account for changing some infection samples to
  # mock samples.
  fdata2_2$Condition <- c(
    rep("Infection", 6),
    rep("Mock", 6)
  )

  fdata2_2$Level <- c(
    "high", "low", "high", "low", "high", "low", "high",
    "high", "low", "low", "low", "high"
  )
  fdata2_2$Gender <- c(
    "M", "F", "M", "F", "M", "F",
    "F", "M", "F", "F", "M", "F"
  )
  fdata2_2$Age <- round(runif(12, min = 19, max = 89), 2)

  pdata_2_2_4 <- as.pepData(
    e_data = edata2_2,
    f_data = fdata2_2,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID"
  )
  pdata_2_2_4 <- group_designation(
    omicsData = pdata_2_2_4,
    main_effects = c("Condition", "Level"),
    covariates = c("Gender", "Age")
  )

  # Create pepData objects to test paired data ---------------------------------

  load(system.file('testdata',
    'little_pairdata.RData',
    package = 'pmartR'
  ))

  pairdata <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Mass_Tag_ID',
    fdata_cname = 'Name',
    emeta_cname = 'Protein'
  )

  pairdata <- group_designation(pairdata,
    pair_id = "PairID",
    pair_group = "Time",
    pair_denom = "0"
  )

  # Carpe unitae testum --------------------------------------------------------

  expect_identical(
    summary(pdata),
    data.frame(
      temp = c("pepData", "12", "150", "83", "341", "0.189"),
      row.names = c(
        "Class", "Unique SampleIDs (f_data)",
        "Unique Mass_Tag_IDs (e_data)", "Unique Proteins (e_meta)",
        "Missing Observations", "Proportion Missing"
      )
    ) %>%
      `names<-`(NULL)
  )

  expect_identical(
    summary(pdata_gdf),
    data.frame(
      temp = c("pepData", "12", "150", "83", "341", "0.189", "9", "3"),
      row.names = c(
        "Class", "Unique SampleIDs (f_data)",
        "Unique Mass_Tag_IDs (e_data)", "Unique Proteins (e_meta)",
        "Missing Observations", "Proportion Missing",
        "Samples per group: Infection", "Samples per group: Mock"
      )
    ) %>%
      `names<-`(NULL)
  )

  expect_identical(
    summary(pdata_2_2_4),
    data.frame(
      temp = c(
        "pepData", "12", "150", "NA", "341", "0.189",
        "3", "3", "3", "3", "Gender, Age"
      ),
      row.names = c(
        "Class", "Unique SampleIDs (f_data)",
        "Unique Mass_Tag_IDs (e_data)", "Rows (e_meta)",
        "Missing Observations", "Proportion Missing",
        "Samples per group: Infection_high",
        "Samples per group: Infection_low",
        "Samples per group: Mock_high",
        "Samples per group: Mock_low",
        "Covariates"
      )
    ) %>%
      `names<-`(NULL)
  )

  expect_identical(
    summary(pairdata),
    data.frame(
      temp = c("pepData", "30", "150", "16", "1184", "0.263", "30", "15"),
      row.names = c(
        "Class", "Unique Names (f_data)",
        "Unique Mass_Tag_IDs (e_data)", "Unique Proteins (e_meta)",
        "Missing Observations", "Proportion Missing",
        "Total Samples", "Pairs"
      )
    ) %>%
      `names<-`(NULL)
  )
})
