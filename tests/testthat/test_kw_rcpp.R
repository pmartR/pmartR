context('Kruskal-Wallis test')

test_that('kw_rcpp calculates the correct p-values by group', {
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

  # Forge a group_DF attribute for pdata.
  pdata <- group_designation(
    omicsData = pdata,
    main_effects = 'Condition'
  )

  # Remove the ID column from e_data
  pe_data <- pdata$e_data[, -1]

  # Kruskal-Wallis tests -------------------------------------------------------

  # Create a vector with group information.
  groups <- c(rep("Infection", 9), rep("Mock", 3))

  # Perform KW test in R.
  kw_r <- as.numeric(
    apply(
      pe_data,
      1,
      function(x) tryCatch(
        kruskal.test(
          x = x[which(!is.na(x))],
          g = groups[which(!is.na(x))]
        )$p.value,
        error = function(err) NA
      )
    )
  )

  # Perform KW test in C++.
  kw_cpp <- pmartR:::kw_rcpp(pe_data %>%
    as.matrix(), groups)

  # Compare R and C++ on the field of battle.
  expect_equal(kw_r, kw_cpp)

  set.seed(338)

  # Scramble the order of the columns/groups.
  eggs <- sample(1:12, 12)

  # Scramble the order of the groups.
  kw_eggs <- pmartR:::kw_rcpp(pe_data[, eggs] %>%
    as.matrix(), groups[eggs])

  # Ensure the p-values are the same even if the data/groups are scrambled.
  expect_equal(kw_r, kw_eggs)
})
