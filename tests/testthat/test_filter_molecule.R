context('filter by molecule')

test_that('molecule_filter and applyFilt produce the correct output',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'filter_data_mol.RData',
                   package = 'pmartR'))
  
  # Test molecule_filter -------------------------------------------------------
  
  # Create a pepData object with the reduced data set.
  pdata_mol <- as.pepData(e_data = edata_mol,
                          f_data = fdata_mol,
                          e_meta = emeta_mol,
                          edata_cname = 'Mass_Tag_ID',
                          fdata_cname = 'SampleID',
                          emeta_cname = 'Protein')
  
  # Try creating a moleculeFilt object with inappropriate input objects.
  expect_error(molecule_filter(omicsData = edata_mol),
               paste("omicsData must be of class 'pepData', 'proData',",
                     "'metabData', 'lipidData', or 'nmrData'",
                     sep = ' '))
  
  # Run molecule filter on the reduced data frame.
  filter_mol <- molecule_filter(omicsData = pdata_mol)
  
  # Review the class for the filter_mol object.
  expect_s3_class(filter_mol,
                  c('moleculeFilt', 'data.frame'))
  
  # Check the dimensions of filter_mol.
  expect_equal(dim(filter_mol),
               c(150, 2))
  
  # Ensure the row sums are correct.
  expect_identical(filter_mol$Num_Observations,
                   c(1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2,
                     2, 1, 1, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 2, 2, 2, 1,
                     1, 1, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 2, 12, 11, 12, 11, 12,
                     12, 12, 12, 12, 11, 12, 12, 12, 4, 10, 12, 3, 12, 11, 9,
                     12, 12, 8, 6, 3, 10, 12, 5, 12, 12, 12, 12, 3, 12, 12, 3,
                     9, 12, 7, 4, 12, 12, 12, 8, 9, 12, 12, 12, 8, 4, 4, 12, 10,
                     11, 12, 12, 5, 3, 3, 12, 12, 12, 3, 12, 12, 10, 12, 12, 11,
                     12, 12, 12, 12, 12, 11, 12, 9, 12, 12, 12, 11, 12, 12, 12,
                     11, 7, 12, 7, 12, 12, 12, 12, 11, 11, 11, 12, 12, 8, 12,
                     12, 12, 12, 10))
  
  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(filter_mol, 'num_samps'),
                   12)
  
  # Test applyFilt.moleculeFilt ------------------------------------------------
  
  # Apply the filter to the reduced peptide data set.
  filtered <- applyFilt(filter_object = filter_mol,
                        omicsData = pdata_mol,
                        min_num = 2)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_mol, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata_mol, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata_mol),
                   class(filtered))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered, 'filters')[[1]]$filtered,
               c(1024, 15714, 16636, 976139, 6769231, 6769844, 6832528,
                 6901575, 6934326, 6948923, 6949061, 6949902, 6954178,
                 6955444, 6959184, 6959316, 6964584, 6966306, 6967194,
                 6967301, 6967315, 6967481))
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(128, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(128, 4))
  
})
