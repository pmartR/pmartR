context('normalize: spans')

test_that('SPANS correctly spans the data', {
  # Load data and create pepData objects ---------------------------------------

  load(system.file('testdata',
    'little_pdata.RData',
    package = 'pmartR'
  ))

  # Create a pepData object with the reduced data set.
  pdata <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID",
    emeta_cname = "Protein"
  )

  # Natural logate the data.
  pdata <- edata_transform(
    omicsData = pdata,
    data_scale = "log"
  )

  # Forge a group_DF attribute for pdata.
  pdata <- group_designation(
    omicsData = pdata,
    main_effects = 'Condition'
  )

  # Create a vector of the samples for the input to spans_make_distribution.
  group <- c(
    rep("Infection", 9),
    rep("Mock", 3)
  )
  group_eggs <- c(
    "Infection", "Mock", "Mock", "Infection", "Infection",
    "Infection", "Infection", "Infection", "Infection", "Mock",
    "Infection", "Infection"
  )

  # Copy the pepData object to create a scrambled version of it.
  pdata_eggs <- pdata

  set.seed(338)

  # Scramble the order of the columns.
  eggs <- sample(1:12, 12)

  # Reorder the sample columns in e_data (keep the ID column as the first
  # column) and reorder the column in f_data. This information will be used when
  # correctly assign group membership when running group_designation().
  pdata_eggs$e_data <- pdata$e_data[, c(1, eggs + 1)]
  pdata_eggs$f_data <- pdata$f_data[eggs, ]

  # Rerun the group designation function on the scrambled pepData object.
  pdata_eggs <- group_designation(
    omicsData = pdata_eggs,
    main_effects = "Condition"
  )

  # Test spans_make_distribution separately ------------------------------------

  sig_inds <- rep(FALSE, 150)
  true_sig_inds <- c(
    2, 6, 8, 9, 10, 14, 15, 24, 26, 32, 40, 41, 44, 47, 48, 66,
    67, 75, 78, 81, 85, 87, 93, 102, 103, 104, 105, 109, 112,
    113, 115, 117, 129, 131, 147
  )
  sig_inds[true_sig_inds] <- TRUE

  nonsig_inds <- rep(FALSE, 150)
  true_nonsig_inds <- c(
    18, 33, 34, 52, 56, 59, 80, 83, 88, 89, 91, 92, 98, 119,
    122, 123, 127, 136, 142, 143, 146
  )
  nonsig_inds[true_nonsig_inds] <- TRUE

  select_n <- c(
    31, 248, 1363, 875, 1242, 1341, 344, 67, 1321, 681, 1194, 1446,
    107, 316, 1316, 1125, 956, 937, 295, 891, 999, 132, 1252, 280
  )

  # Make distribution - original data ---------------

  # Create a list to hold the output for the spans_make_distribution function.
  test_md <- vector(
    mode = "list",
    length = length(select_n)
  )

  for (e in 1:length(select_n)) {
    # Create a background distribution (whatever that is).
    test_md[[e]] <- pmartR:::spans_make_distribution(
      omicsData = pdata,
      group_vector = group,
      norm_fn = c(
        "median",
        "mean",
        "zscore",
        "mad"
      ),
      sig_inds = sig_inds,
      nonsig_inds = nonsig_inds,
      select_n = select_n[[e]]
    )
  }

  # Calculate the means of the background distribution to use for testing
  # purposes. The test uses a threshold on the mean instead of directly
  # comparing vector values.
  mean1 <- mean(c(sapply(test_md, `[[`, 1)))
  mean2 <- mean(c(sapply(test_md, `[[`, 2)))

  # Make sure the background distribution is correct.
  expect_true(mean1 > 1.33 && mean1 < 1.94)
  expect_true(mean2 > -0.74 && mean2 < -0.25)

  # Make distribution - scrambled data ---------------

  # Create a list to hold the output for the spans_make_distribution function.
  test_eggs <- vector(
    mode = "list",
    length = length(select_n)
  )

  for (e in 1:length(select_n)) {
    # Create a background distribution (whatever that is).
    test_eggs[[e]] <- pmartR:::spans_make_distribution(
      omicsData = pdata_eggs,
      group_vector = group_eggs,
      norm_fn = c(
        "median",
        "mean",
        "zscore",
        "mad"
      ),
      sig_inds = sig_inds,
      nonsig_inds = nonsig_inds,
      select_n = select_n[[e]]
    )
  }

  # Calculate the means of the background distribution to use for testing
  # purposes. The test uses a threshold on the mean instead of directly
  # comparing vector values.
  mean1_eggs <- mean(c(sapply(test_eggs, `[[`, 1)))
  mean2_eggs <- mean(c(sapply(test_eggs, `[[`, 2)))

  # Make sure the background distribution is correct.
  expect_true(mean1_eggs > 1.33 && mean1_eggs < 1.94)
  expect_true(mean2_eggs > -0.74 && mean2_eggs < -0.25)

  # Create SPANS data frame standards ------------------------------------------

  # Original data ---------------

  # Forge a data frame standard for the SPANSRes object.
  spandard <- data.frame(
    subset_method = c(
      "rip", "ppp_rip", "rip", "ppp_rip", "rip", "ppp_rip",
      "rip", "all", "ppp_rip", "all", "all", "all", "los",
      "ppp", "ppp", "ppp", "ppp", "los", "los", "los"
    ),
    normaliztion_method = c(
      "median", "median", "mean", "mean", "zscore",
      "zscore", "mad", "mean", "mad", "zscore", "mad",
      "median", "mean", "mean", "zscore", "median",
      "mad", "median", "zscore", "mad"
    ),
    SPANS_score = c(
      0.87, 0.87, 0.87, 0.87, 0.8, 0.8, 0.8, 0.78, 0.78, 0.7, 0.7,
      0.66, 0.55, 0.55, 0.55, 0.53, 0.53, NA, NA, NA
    ),
    parameters = c(
      "0.2", "0.5;0.2", "0.2", "0.5;0.2", "0.2", "0.5;0.2",
      "0.2", "", "0.5;0.2", "", "", "", "0.05", "0.5", "0.5",
      "0.5", "0.5", "0.05", "0.05", "0.05"
    ),
    mols_used_in_norm = c(
      22, 27, 22, 27, 22, 27, 22, 150, 27, 150, 150, 150,
      14, 124, 124, 124, 124, 14, 14, 14
    ),
    passed_selection = c(
      rep(TRUE, 17),
      rep(FALSE, 3)
    )
  )

  # Set up the attributes for the spandard :)
  attributes(spandard) <- list(
    names = c(
      "subset_method", "normalization_method", "SPANS_score",
      "parameters", "mols_used_in_norm", "passed_selection"
    ),
    row.names = 1:20,
    class = c("SPANSRes", "data.frame"),
    method_selection_pvals = data.frame(
      subset_method = spandard$subset_method,
      normalization_method = spandard$normaliztion_method,
      parameters = spandard$parameters,
      location_p_value = c(
        0.926340489251137, 0.926340489251137,
        0.405380556458942, 0.309177044632912,
        0.405380556458942, 0.309177044632912,
        0.926340489251137, 0.116031614870982,
        0.926340489251137, 0.116031614870982,
        0.0522036353413146, 0.0522036353413146,
        0.0522036353413146, 0.0522036353413146,
        0.0522036353413146, 0.0522036353413146,
        0.0522036353413146, 0.020819004906484,
        0.0522036353413146, 0.020819004906484
      ),
      scale_p_value = c(
        NA, NA, NA, NA, 0.405380556458942, 0.517534719643234,
        0.0789944358337225, NA, 0.16551785869747,
        0.405380556458942, 0.781511294998709, NA, NA, NA,
        0.781511294998709, NA, 0.0522036353413146, NA,
        0.0125549185969666, 0.781511294998709
      ),
      F_log_HSmPV = c(
        rep(1, 17),
        rep(NA, 3)
      ),
      F_log_NSmPV = c(
        0.74, 0.74, 0.74, 0.74, 0.6, 0.6, 0.6, 0.56, 0.56, 0.4,
        0.4, 0.32, 0.1, 0.1, 0.1, 0.06, 0.06, NA, NA, NA
      )
    ),
    group_vector = c(rep("Infection", 9), rep("Mock", 3)),
    significant_thresh = 1e-4,
    nonsignificant_thresh = 0.5,
    n_not_significant = 21,
    n_significant = 35,
    location_threshold = 0.05,
    scale_thresh = 0.05
  )

  # Scrambled data ---------------

  # Forge a data frame standard for the scrambled SPANSRes object.
  spandard_eggs <- data.frame(
    subset_method = spandard$subset_method,
    normaliztion_method = spandard$normalization_method,
    SPANS_score = c(
      0.83, 0.83, 0.83, 0.83, 0.71, 0.71, 0.71, 0.66, 0.66, 0.61,
      0.61, 0.58, 0.53, 0.53, 0.53, 0.52, 0.52, NA, NA, NA
    ),
    parameters = spandard$parameters,
    mols_used_in_norm = spandard$mols_used_in_norm,
    passed_selection = spandard$passed_selection
  )

  # Set up the attributes for the spandard :)
  attributes(spandard_eggs) <- list(
    names = c(
      "subset_method", "normalization_method", "SPANS_score",
      "parameters", "mols_used_in_norm", "passed_selection"
    ),
    row.names = 1:20,
    class = c("SPANSRes", "data.frame"),
    method_selection_pvals = data.frame(
      subset_method = spandard_eggs$subset_method,
      normalization_method = spandard_eggs$normaliztion_method,
      parameters = spandard_eggs$parameters,
      location_p_value = attr(
        spandard,
        "method_selection_pvals"
      )$location_p_value,
      scale_p_value = attr(
        spandard,
        "method_selection_pvals"
      )$scale_p_value,
      F_log_HSmPV = attr(
        spandard,
        "method_selection_pvals"
      )$F_log_HSmPV,
      F_log_NSmPV = c(
        0.66, 0.66, 0.66, 0.66, 0.42, 0.42, 0.42, 0.32, 0.32,
        0.22, 0.22, 0.16, 0.06, 0.06, 0.06, 0.04, 0.04, NA, NA,
        NA
      )
    ),
    group_vector = group_eggs,
    significant_thresh = 1e-4,
    nonsignificant_thresh = 0.5,
    n_not_significant = 21,
    n_significant = 35,
    location_threshold = 0.05,
    scale_thresh = 0.05
  )

  # SPANS the data -------------------------------------------------------------

  # Original data ---------------

  # Set a seed because random functions are used within spans_procedure.
  set.seed(6)

  # Run the spans function with the default settings.
  spansiard <- spans_procedure(pdata,
    params = list(
      "los" = list(0.05),
      "ppp" = list(0.5),
      "rip" = list(0.2),
      "ppp_rip" = list(c(0.5, 0.2))
    ),
    n_iter = 50,
    parallel = FALSE,
    verbose = FALSE
  )

  # Start sleuthing!
  expect_s3_class(spansiard, "SPANSRes")
  expect_equal(spansiard, spandard)

  # Scrambled data ---------------

  # Set a seed because random functions are used within spans_procedure.
  set.seed(6)

  # Run SPANS on the scrambled data. The output will be different from the
  # data that is not scrambled because non-NA elements of e_data are randomly
  # selected to calculate the background distribution. Even if the same seed is
  # set the actual values that are selected to calculate the background
  # distribution will be different because the order of the columns is
  # different.
  spans_eggs <- spans_procedure(pdata_eggs,
    params = list(
      "los" = list(0.05),
      "ppp" = list(0.5),
      "rip" = list(0.2),
      "ppp_rip" = list(c(0.5, 0.2))
    ),
    n_iter = 50,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check the output of the scrambled data against the unscrambled data.
  expect_equal(spans_eggs, spandard_eggs)
})
