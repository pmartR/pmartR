context('count non-missing values in e_data')

test_that('nonmissing_per_grp correctly counts non-missing values by group', {
  # Load data and prepare omicsData objects ------------------------------------

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

  # Run the group_designation function on pdata.
  pdata <- group_designation(
    omicsData = pdata,
    main_effects = "Condition"
  )

  # Extract the group information from pdata.
  groupDF <- attr(pdata, "group_DF")

  # Scramble the order of the samples.
  set.seed(4)
  scrambled <- sample(names(edata[, -1]))

  # Scramble the order of the columns in e_data and rows in f_data.
  edata_s <- edata[, c(1, match(scrambled, names(edata)[-1]) + 1)]
  fdata_s <- fdata[match(scrambled, names(edata)[-1]), ]

  # Create a pepData object with the scrambled data frames.
  pdata_s <- as.pepData(
    e_data = edata_s,
    f_data = fdata_s,
    edata_cname = "Mass_Tag_ID",
    fdata_cname = "SampleID"
  )

  # Run the group_designation function on pdata_s.
  pdata_s <- group_designation(
    omicsData = pdata_s,
    main_effects = "Condition"
  )

  # Fish out the group information from the scrambled pepData object.
  groupDF_s <- attr(pdata_s, "group_DF")

  # Test nonmissing_per_grp with ordered columns -------------------------------

  # Assemble the inputs to nonmissing_per_grp.
  group_dat <- as.character(groupDF$Group)
  e_data <- pdata$e_data[, -1]

  # Run nonmissing_per_grp with the columns in the original order.
  nonmissing <- pmartR:::nonmissing_per_grp(
    as.matrix(e_data),
    group_dat
  )

  # Check the counts for each column in nonmissing are correct.
  expect_equal(
    nonmissing[, 1],
    as.vector(rowSums(!is.na(pdata$e_data[, 2:10])))
  )
  expect_equal(
    nonmissing[, 2],
    as.vector(rowSums(!is.na(pdata$e_data[, 11:13])))
  )

  # Test nonmissing_per_grp with unordered columns -----------------------------

  # The following tests are to make sure how the group information and the data
  # matrix are ordered works properly. These objects are ordered within the
  # nonmissing_per_group R function before calling the nonmissin_per_grp c++
  # function.

  # Reorder the samples.
  group_dat_s <- as.character(groupDF_s$Group[order(groupDF_s$Group)])

  # First remove the Mass_Tag_ID column
  e_data_s <- pdata_s$e_data[, -1]
  e_data_s <- e_data_s[, order(groupDF_s$Group)]

  # Run nonmissing_per_grp with the column order scrambled and then reordered.
  nonmissing_s <- pmartR:::nonmissing_per_grp(
    as.matrix(e_data_s),
    group_dat_s
  )

  # Check the counts for each column in nonmissing_s match those from before.
  expect_identical(
    nonmissing,
    nonmissing_s
  )
})
