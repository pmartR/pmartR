context('filter by molecule')

test_that('molecule_filter and applyFilt produce the correct output',{

  # Load data and prepare omicsData objects ------------------------------------

  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))

  # Create a pepData object with the reduced data set.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")

  # Create a function to count the number of non-missing values to be used in
  # connection with apply.
  count <- function (x) {

    present <- sum(!is.na(x))

    return (present)

  }

  # Test molecule_filter -------------------------------------------------------

  # Try creating a moleculeFilt object with inappropriate input objects.
  expect_error(molecule_filter(omicsData = edata),
               paste("omicsData must be of class 'pepData', 'proData',",
                     "'metabData', 'lipidData', or 'nmrData'",
                     sep = ' '))

  # Run molecule filter on the reduced data frame.
  filter <- molecule_filter(omicsData = pdata)

  # Review the class for the filter object.
  expect_s3_class(filter,
                  c('moleculeFilt', 'data.frame'))

  # Check the dimensions of filter.
  expect_equal(dim(filter),
               c(150, 2))

  # Count the number of molecules outside of the molecule_filter function
  s_count <- apply(pdata$e_data[, -1], 1, count)

  # Remove the names from the counts vector.
  names(s_count) <- NULL

  # Ensure the row sums are correct.
  expect_equal(filter$Num_Observations,
               s_count)

  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(filter, 'num_samps'),
                   12)

  # Test applyFilt.moleculeFilt ------------------------------------------------

  # Apply the filter to the reduced peptide data set.
  filtered <- applyFilt(filter_object = filter,
                        omicsData = pdata,
                        min_num = 2)

  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered))

  # Find the peptide IDs that will be filtered (without molecule_filter
  # functions).
  s_pepes <- pdata$e_data[which(s_count < 2), 1]

  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered, 'filters')[[1]]$filtered,
               s_pepes)
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered$e_data)),
         prop_missing = (sum(is.na(filtered$e_data)) /
                           prod(dim(filtered$e_data[, -1]))),
         num_samps = ncol(filtered$e_data[, -1]),
         data_types = NULL)
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered$e_meta$Protein)))
  )

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(141, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(141, 4))

  # Test applying the same filter a second time --------------------------------

  # Run molecule filter on the already filtered data set.
  filter_2 <- molecule_filter(omicsData = filtered)

  # Review the class for the filter object.
  expect_s3_class(filter_2,
                  c('moleculeFilt', 'data.frame'))

  # Check the dimensions of filter.
  expect_equal(dim(filter_2),
               c(141, 2))

  # Count the number of molecules outside of the molecule_filter function
  s_count_2 <- apply(filtered$e_data[, -1], 1, count)

  # Remove the names from the count standard.
  names(s_count_2) <- NULL

  # Find the peptide IDs that will be filtered (without molecule_filter
  # functions).
  s_pepes_2 <- filtered$e_data[which(s_count_2 < 3), 1]

  # Ensure the row sums are correct.
  expect_equal(filter_2$Num_Observations,
               s_count_2)

  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(filter, 'num_samps'),
                   12)

  # Test outcome of applying the same filter more than once.
  expect_warning(filtered_2 <- applyFilt(filter_object = filter_2,
                                         omicsData = filtered,
                                         min_num = 3),
                 'A molecule filter has already been applied to this data set.')

  # Check that the filters attribute is a list with two elements.
  expect_true(length(attr(filtered_2, 'filters')) == 2)

  # Examine the first element in the filters attribute.
  expect_equal(attr(filtered_2, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_2, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered_2, 'filters')[[1]]$filtered,
               s_pepes)
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))

  # Examine the second element in the filters attribute.
  expect_equal(attr(filtered_2, 'filters')[[2]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_2, 'filters')[[2]]$threshold,
                   3)
  expect_equal(attr(filtered_2, 'filters')[[2]]$filtered,
               s_pepes_2)
  expect_true(is.na(attr(filtered_2, 'filters')[[2]]$method))

  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_2$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered_2$e_data)),
         prop_missing = (sum(is.na(filtered_2$e_data)) /
                           prod(dim(filtered_2$e_data[, -1]))),
         num_samps = ncol(filtered_2$e_data[, -1]),
         data_types = NULL)
  )

  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_2, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_2$e_meta$Protein)))
  )

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_2$e_data),
               c(137, 13))
  expect_equal(dim(filtered_2$f_data),
               c(12, 2))
  expect_equal(dim(filtered_2$e_meta),
               c(137, 4))

  # Test scenario when nothing is filtered -------------------------------------

  # Apply the filter with a value for min_num that will not filter any rows.
  expect_message(noFilta <- applyFilt(filter_object = filter,
                                      omicsData = pdata,
                                      min_num = 1),
                 paste("No biomolecules were filtered with the value specified",
                       "for the min_num argument.",
                       sep = " "))

  # The output of applyFilt should be the same as the omicsData object used as
  # the input because the filter was not applied. Therefore, the filters
  # attribute should remain how it was before running applyFilt.
  expect_identical(noFilta, pdata)

})
