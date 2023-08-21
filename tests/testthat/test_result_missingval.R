context('results: missing value')

test_that('missingval_result correctly counts missing values', {
  # Load data and create omicsData objects ---------------------------------------

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

  # Create standards -----------------------------------------------------------

  # Create rowwise counts.
  count_row <- rowSums(is.na(pdata$e_data[, -1]))

  # Create columnwise counts.
  count_col <- colSums(is.na(pdata$e_data[, -1]))

  # Create standard for naRes object.
  standard <- list(
    "na.by.sample" = data.frame("SampleID" = names(pdata$e_data)[-1],
                                "num_NA" = as.numeric(count_col),
                                "num_non_NA" = nrow(pdata$e_data) - as.numeric(count_col),
                                "Condition" = pdata$f_data[, 2]),
    "na.by.molecule" = data.frame("Mass_Tag_ID" = pdata$e_data[, 1],
                                  "num_NA" = as.numeric(count_row),
                                  "num_non_NA" = nrow(pdata$f_data) - as.numeric(count_row)
                                  )
  )

  # Add class and attribute crap to the standard.
  class(standard) <- "naRes"
  attr(standard, "cnames") <- list(
    "edata_cname" = get_edata_cname(pdata),
    "fdata_cname" = get_fdata_cname(pdata)
  )

  # Missing values: pepData object ---------------------------------------------

  # Create an naRes object.
  missingval <- missingval_result(pdata)

  # Sleuth around the naRes object.
  expect_identical(standard, missingval)
})
