context('filter by peptide and protein counts')

test_that('proteomics_filter and applyFilt produce the correct output',{
  
  # Load the peptide data frames -----------------------------------------------
  
  load(system.file('testdata',
                   'filter_data_pro.RData',
                   package = 'pmartR'))
  
  # Test proteomics_filter -----------------------------------------------------
  
  # Create a pepData object with the reduced data set.
  pdata_pro <- as.pepData(e_data = edata_pro,
                          f_data = fdata_pro,
                          e_meta = emeta_pro,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID",
                          emeta_cname = "Protein")
  
  # Try creating a proteomicsFilt object with unholy input objects.
  expect_error(proteomics_filter(omicsData = fdata_pro),
               paste("omicsData must be of class 'pepData'",
                     sep = " "))
  
  # Try creating a proteomicsFilt object without the e_meta data frame.
  expect_error(proteomics_filter(as.pepData(e_data = edata_pro,
                                            f_data = fdata_pro,
                                            edata_cname = "Mass_Tag_ID",
                                            fdata_cname = "SampleID")),
               paste("e_meta must be non-NULL",
                     sep = " "))
  
  # Run proteomics_filter with holy input objects.
  filter_pro <- proteomics_filter(omicsData = pdata_pro)
  
  # Review the class for the filter_pro object.
  expect_s3_class(filter_pro,
                  c('proteomicsFilt', 'data.frame'))
  
  # Check the length of the filter_pro list.
  expect_equal(length(filter_pro), 2)
  
  # Ensure the counts for each peptide are correct.
  expect_identical(filter_pro[[1]],
                   pfStandard[[1]])
  
  # Verify the counts for each protein are correct.
  expect_identical(filter_pro[[2]],
                   pfStandard[[2]])
  
  # Test applyFilt.proteomicsFilt ----------------------------------------------
  
  # Apply the proteomics filter only using the min_num_peps argument.
  filtered <- applyFilt(filter_object = filter_pro,
                        omicsData = pdata_pro,
                        min_num_peps = 2)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_pro, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata_pro, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata_pro),
                   class(filtered))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'proteomicsFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                   data.frame(min_num_peps = 2,
                              degen_peps = as.character(FALSE)))
  expect_equal(attr(filtered, 'filters')[[1]]$filtered$e_meta_remove,
               attr(afStandard, "filters")[[1]]$filtered$emeta_filt)
  expect_equal(attr(filtered, 'filters')[[1]]$filtered$e_data_remove,
               attr(afStandard, "filters")[[1]]$filtered$edata_filt)
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered, 'data_info')$num_edata,
               380)
  expect_equal(attr(filtered, 'data_info')$num_miss_obs,
               635)
  expect_equal(round(attr(filtered, 'data_info')$prop_missing, 4),
               0.1393)
  expect_equal(attr(filtered, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered, 'meta_info')$meta_data)
  expect_equal(attr(filtered, 'meta_info')$num_emeta,
               79)
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(380, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(380, 4))
  
})
