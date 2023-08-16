context('edata structure')

test_that('errors are thrown for bad data structures/values', {
  # Load/create wayward data frames --------------------------------------------

  # Read csv created by Excel that contains values divided by zero and infinity.
  sour_grapes <- read.csv(
    system.file('testdata',
      'sour_grapes.csv',
      package = 'pmartR'
    ),
    fileEncoding = 'UTF-8-BOM'
  )

  # Read csv created by Excel that contains values divided by zero, infinity,
  # and NaN.
  sour_grapes_2 <- read.csv(
    system.file('testdata',
      'sour_grapes_2.csv',
      package = 'pmartR'
    ),
    fileEncoding = 'UTF-8-BOM'
  )

  set.seed(8346)

  # Create a data frame with random numbers.
  cantor <- data.frame(
    id = paste0('e', 1:10),
    col1 = rnorm(10),
    col2 = rcauchy(10),
    col3 = rexp(10),
    col4 = runif(10),
    col5 = rbeta(10, 3, 5)
  )

  # Change some values to +-Inf
  cantor[8, 3] <- Inf

  # Forge an f_data data frame for the naughty edata data frames.
  fdata <- data.frame(
    sampleID = paste0('col', 1:5),
    condition = c(
      rep('case', 3),
      rep('control', 2)
    )
  )

  # Test data structure with as.xxx functions ----------------------------------

  # Construct a pepData object with sour_grapes.
  expect_error(
    as.pepData(
      e_data = sour_grapes,
      f_data = fdata,
      edata_cname = 'id',
      fdata_cname = 'sampleID'
    ),
    'Columns 2 and 4 of e_data contain non-numeric values.'
  )

  # Generate a lipidData object with sour_grapes_2.
  expect_error(
    as.lipidData(
      e_data = sour_grapes_2,
      f_data = fdata,
      edata_cname = 'id',
      fdata_cname = 'sampleID'
    ),
    'Column 4 of e_data contains non-numeric values.'
  )



  # forge a metabdata object with cantor.
  expect_error(
    as.metabData(
      e_data = cantor,
      f_data = fdata,
      edata_cname = 'id',
      fdata_cname = 'sampleID'
    ),
    'Column 3 of e_data contains infinite values.'
  )

  cantor[4, 6] <- -Inf

  # forge a metabdata object with cantor.
  expect_error(
    as.metabData(
      e_data = cantor,
      f_data = fdata,
      edata_cname = 'id',
      fdata_cname = 'sampleID'
    ),
    'Columns 3 and 6 of e_data contain infinite values.'
  )
})
