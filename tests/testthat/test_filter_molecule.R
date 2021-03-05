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
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered, 'data_info')$num_edata,
               128)
  expect_equal(attr(filtered, 'data_info')$num_miss_obs,
               430)
  expect_equal(round(attr(filtered, 'data_info')$prop_missing, 4),
               0.2799)
  expect_equal(attr(filtered, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered, 'meta_info')$meta_data)
  expect_equal(attr(filtered, 'meta_info')$num_emeta,
               75)
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(128, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(128, 4))
  
  # Test applying the same filter a second time --------------------------------
  
  # Run molecule filter on the already filtered data set.
  filter_mol_2 <- molecule_filter(omicsData = filtered)
  
  # Review the class for the filter_mol object.
  expect_s3_class(filter_mol_2,
                  c('moleculeFilt', 'data.frame'))
  
  # Check the dimensions of filter_mol.
  expect_equal(dim(filter_mol_2),
               c(128, 2))
  
  # Ensure the row sums are correct.
  expect_identical(filter_mol_2$Num_Observations,
                   c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                     2, 2, 2, 2, 2, 12, 11, 12, 11, 12, 12, 12, 12, 12, 11, 12,
                     12, 12, 4, 10, 12, 3, 12, 11, 9, 12, 12, 8, 6, 3, 10, 12,
                     5, 12, 12, 12, 12, 3, 12, 12, 3, 9, 12, 7, 4, 12, 12, 12,
                     8, 9, 12, 12, 12, 8, 4, 4, 12, 10, 11, 12, 12, 5, 3, 3, 12,
                     12, 12, 3, 12, 12, 10, 12, 12, 11, 12, 12, 12, 12, 12, 11,
                     12, 9, 12, 12, 12, 11, 12, 12, 12, 11, 7, 12, 7, 12, 12,
                     12, 12, 11, 11, 11, 12, 12, 8, 12, 12, 12, 12, 10))
  
  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(filter_mol, 'num_samps'),
                   12)
  
  # Test outcome of applying the same filter more than once.
  expect_warning(filtered_2 <- applyFilt(filter_object = filter_mol_2,
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
               c(1024, 15714, 16636, 976139, 6769231, 6769844, 6832528,
                 6901575, 6934326, 6948923, 6949061, 6949902, 6954178,
                 6955444, 6959184, 6959316, 6964584, 6966306, 6967194,
                 6967301, 6967315, 6967481))
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))
  
  # Examine the second element in the filters attribute.
  expect_equal(attr(filtered_2, 'filters')[[2]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_2, 'filters')[[2]]$threshold,
                   3)
  expect_equal(attr(filtered_2, 'filters')[[2]]$filtered,
               c(1687, 11083, 6809644, 6831118, 6948918, 6949195, 6949275,
                 6955031, 6955082, 6955096, 6955110, 6955122, 6959166, 6962332,
                 6962373, 6962578, 6964253, 6967197, 6967234, 6967241, 6967247,
                 6967304, 6967435, 6967456, 6967706))
  expect_true(is.na(attr(filtered_2, 'filters')[[2]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered_2, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered_2, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_2, 'data_info')$num_edata,
               103)
  expect_equal(attr(filtered_2, 'data_info')$num_miss_obs,
               180)
  expect_equal(round(attr(filtered_2, 'data_info')$prop_missing, 4),
               0.1456)
  expect_equal(attr(filtered_2, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered_2, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered_2, 'meta_info')$meta_data)
  expect_equal(attr(filtered_2, 'meta_info')$num_emeta,
               59)
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_2$e_data),
               c(103, 13))
  expect_equal(dim(filtered_2$f_data),
               c(12, 2))
  expect_equal(dim(filtered_2$e_meta),
               c(103, 4))
  
})
