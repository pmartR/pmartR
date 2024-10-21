context('summary: filters')

test_that('the filter object summaries are all square', {
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

  pdata <- edata_transform(pdata, "log")

  pdata_gdf <- group_designation(
    omicsData = pdata,
    main_effects = 'Condition'
  )

  # Create pepData objects to test multiple main effects and covariates --------

  # Copy edata so the names of the samples can be changed.
  edata_2 <- edata

  # Change some of the Infection samples to Mock samples.
  names(edata_2) <- c(
    "Mass_Tag_ID",
    paste0("Infection", 1:6),
    paste0("Mock", 1:6)
  )

  # Create additional f_data objects with different main effects and covariates.
  fdata_2 <- fdata

  # Update the sample names in f_data.
  fdata_2$SampleID <- c(
    paste0("Infection", 1:6),
    paste0("Mock", 1:6)
  )

  # Update the first main effect to account for changing some infection samples to
  # mock samples.
  fdata_2$Condition <- c(
    rep("Infection", 6),
    rep("Mock", 6)
  )

  fdata_2$Level <- c(
    "high", "low", "high", "low", "high", "low", "high",
    "high", "low", "low", "low", "high"
  )
  fdata_2$Gender <- c(
    "M", "F", "M", "F", "M", "F",
    "F", "M", "F", "F", "M", "F"
  )
  fdata_2$Age <- round(runif(12, min = 19, max = 89), 2)

  pdata_2 <- as.pepData(
    e_data = edata_2,
    f_data = fdata_2,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID"
  )
  pdata_2 <- edata_transform(pdata_2, "log")
  pdata_2 <- group_designation(
    omicsData = pdata_2,
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
  pairdata <- edata_transform(pairdata, "log")
  pairdata <- group_designation(pairdata,
    pair_id = "PairID",
    pair_group = "Time",
    pair_denom = "0"
  )

  # Carpe unitae testum --------------------------------------------------------

  # Molecule filter ---------------

  expect_equal(
    summary(molecule_filter(pdata), min_num = 2),
    structure(
      list(
        pep_observation_counts = structure(
          list(
            num_observations = 1:12,
            frequency_counts = c(9, 13, 20, 24, 26, 27, 32, 37, 43, 49, 61, 150)
          ),
          row.names = c(
            "(0,1]", "(1,2]", "(2,3]", "(3,4]",
            "(4,5]", "(5,6]", "(6,7]", "(7,8]",
            "(8,9]", "(9,10]", "(10,11]", "(11,12]"
          ),
          class = "data.frame"
        ),
        min_num = 2,
        num_not_filtered = 141,
        num_filtered = 9
      ),
      names = c(
        "pep_observation_counts", "min_num",
        "num_not_filtered", "num_filtered"
      ),
      class = c("moleculeFilterSummary", "list"),
      use_batch = FALSE,
      use_groups = FALSE
    )
  )

  # CV filter ---------------

  expect_equal(
    summary(cv_filter(pdata_gdf), cv_threshold = 72),
    structure(
      list(
        CVs = summary(na.omit(cv_filter(pdata_gdf)$CV)),
        tot_NAs = 10,
        tot_non_NAs = 140,
        filtered_biomolecules = 1
      ),
      class = c("cvFilterSummary", "list")
    )
  )

  # RMD filter ---------------

  expect_equal(
    summary(rmd_filter(pdata_2), pvalue_threshold = 0.01),
    structure(
      list(
        pvalue = summary(rmd_filter(pdata_2)$pvalue),
        metrics = c(
          "MAD", "Kurtosis", "Skewness",
          "Corr", "Proportion_Missing"
        ),
        filtered_samples = "Mock2"
      ),
      names = c("pvalue", "metrics", "filtered_samples"),
      class = c("rmdFilterSummary", "list")
    )
  )

  expect_equal(
    summary(rmd_filter(pairdata), pvalue_threshold = 0.01),
    structure(
      list(
        pvalue = summary(rmd_filter(pairdata)$pvalue),
        metrics = c(
          "MAD", "Kurtosis", "Skewness",
          "Corr", "Proportion_Missing"
        ),
        filtered_samples = c(
          "Mock_0hr_3", "Mock_18hr_3",
          "AM_0hr_2", "AM_18hr_2"
        )
      ),
      names = c("pvalue", "metrics", "filtered_samples"),
      class = c("rmdFilterSummary", "list")
    )
  )

  # Proteomics filter ---------------

  expect_equal(
    summary(proteomics_filter(pdata_gdf), min_num_peps = 2),
    structure(
      list(
        num_per_pep = summary(proteomics_filter(pdata_gdf)$counts_by_pep$n),
        num_per_pro = summary(proteomics_filter(pdata_gdf)$counts_by_pro$n),
        num_pep_filtered = 0,
        num_pro_filtered = 56,
        num_pro_notfiltered = 27,
        num_pep_notfiltered = 150
      ),
      class = c("proteomicsFilterSummary", "list")
    )
  )

  # IMD-ANOVA filter ---------------

  expect_equal(
    summary(imdanova_filter(pdata_2), min_nonmiss_anova = 2),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 25,
        num_not_filtered = 125
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(imdanova_filter(pdata_2), min_nonmiss_gtest = 3),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 29,
        num_not_filtered = 121
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(imdanova_filter(pdata_2),
      min_nonmiss_anova = 2,
      min_nonmiss_gtest = 3
    ),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 23,
        num_not_filtered = 127
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(imdanova_filter(pairdata), min_nonmiss_anova = 2),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 19,
        num_not_filtered = 131
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(imdanova_filter(pairdata), min_nonmiss_gtest = 3),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 0,
        num_not_filtered = 150
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(imdanova_filter(pairdata),
      min_nonmiss_anova = 2,
      min_nonmiss_gtest = 3
    ),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 19,
        num_not_filtered = 131
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(
      imdanova_filter(pdata_2),
      min_nonmiss_anova = 2,
      comparisons = data.frame(check.names = FALSE, 
        Test = c("Infection_high", "Infection_low"),
        Control = rep("Mock_low", 2)
      )
    ),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 33,
        num_not_filtered = 117
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(
      imdanova_filter(pdata_2),
      min_nonmiss_anova = 2,
      comparisons = data.frame(check.names = FALSE, 
        Test = c("Infection_high", "Mock_high"),
        Control = c("Mock_high", "Infection_low")
      )
    ),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 27,
        num_not_filtered = 123
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  expect_equal(
    summary(
      imdanova_filter(pdata_2),
      min_nonmiss_gtest = 3,
      comparisons = data.frame(check.names = FALSE, 
        Test = c("Infection_high", "Infection_low"),
        Control = rep("Mock_low", 2)
      )
    ),
    structure(
      list(
        pep_observation_counts = 150,
        num_filtered = 29,
        num_not_filtered = 121
      ),
      names = c("pep_observation_counts", "num_filtered", "num_not_filtered"),
      class = c("imdanovaFilterSummary", "list")
    )
  )

  # Custom filter ---------------

  ## Test with remove arguments

  # remove e_data
  expect_equal(
    summary(custom_filter(pdata, e_data_remove = "1406")),
    structure(
      data.frame(check.names = FALSE, 
        Filtered = c(0, 1, 0),
        Remaining = c(12, 149, 83),
        Total = c(12, 150, 83)
      ),
      names = c("Filtered", "Remaining", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # remove f_data
  expect_equal(
    summary(custom_filter(pdata, f_data_remove = "Infection4")),
    structure(
      data.frame(check.names = FALSE, 
        Filtered = c(1, 0, 0),
        Remaining = c(11, 150, 83),
        Total = c(12, 150, 83)
      ),
      names = c("Filtered", "Remaining", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # remove e_meta
  expect_equal(
    summary(custom_filter(pdata, e_meta_remove = "ALBU_HUMAN")),
    structure(
      data.frame(check.names = FALSE, 
        Filtered = c(0, 5, 1),
        Remaining = c(12, 145, 82),
        Total = c(12, 150, 83)
      ),
      names = c("Filtered", "Remaining", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # combination
  expect_equal(
    summary(
      custom_filter(
        pdata,
        e_data_remove = "6948840",
        f_data_remove = "Mock2",
        e_meta_remove = "COR1C_HUMAN"
      )
    ),
    structure(
      data.frame(check.names = FALSE, 
        Filtered = c(1, 1, 1),
        Remaining = c(11, 149, 82),
        Total = c(12, 150, 83)
      ),
      names = c("Filtered", "Remaining", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  ## Test with keep arguments

  # e_meta keep
  expect_equal(
    summary(
      custom_filter(
        pdata,
        e_meta_keep = c("ALBU_HUMAN", "RL40_HUMAN", "ALBU_HUMAN")
      )
    ),
    structure(
      data.frame(check.names = FALSE, 
        Discarded = c(0, 141, 80),
        Kept = c(12, 9, 3),
        Total = c(12, 150, 83)
      ),
      names = c("Discarded", "Kept", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # f_data keep
  expect_equal(
    summary(
      custom_filter(
        pdata,
        f_data_keep = c("Infection1", "Infection2", "Mock1", "Mock2", "Mock3")
      )
    ),
    structure(
      data.frame(check.names = FALSE, 
        Discarded = c(7, 0, 0),
        Kept = c(5, 150, 83),
        Total = c(12, 150, 83)
      ),
      names = c("Discarded", "Kept", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # e_data keep
  expect_equal(
    summary(
      custom_filter(
        pdata,
        e_data_keep = c("1024", "4198254", "6637724", "216191"),
        f_data_keep = c("Infection1", "Infection2", "Mock1", "Mock2", "Mock3"),
        e_meta_keep = c("ALBU_HUMAN", "RL40_HUMAN", "ALBU_HUMAN")
      )
    ),
    structure(
      data.frame(check.names = FALSE, 
        Discarded = c(7, 146, 80),
        Kept = c(5, 4, 3),
        Total = c(12, 150, 83)
      ),
      names = c("Discarded", "Kept", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )

  # combination
  expect_equal(
    summary(
      custom_filter(
        pdata,
        e_data_keep = c("1024", "4198254", "6637724", "216191")
      )
    ),
    structure(
      data.frame(check.names = FALSE, 
        Discarded = c(0, 146, 79),
        Kept = c(12, 4, 4),
        Total = c(12, 150, 83)
      ),
      names = c("Discarded", "Kept", "Total"),
      class = c("customFilterSummary", "data.frame"),
      row.names = c(
        "SampleIDs (f_data)", "Mass_Tag_IDs (e_data)",
        "Proteins (e_meta)"
      )
    )
  )
})
