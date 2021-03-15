context('filter by coefficient of variation')

test_that('cv_filter and applyFilt produce the correct output',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'filter_data_cv.RData',
                   package = 'pmartR'))
  
  # Test filter_cv without groups ----------------------------------------------
  
  # Create a pepData object with the reduced data set.
  pdata_cv <- as.pepData(e_data = edata_cv,
                         f_data = fdata_cv,
                         e_meta = emeta_cv,
                         edata_cname = 'Mass_Tag_ID',
                         fdata_cname = 'SampleID',
                         emeta_cname = 'Protein')
  
  # Try creating a cvFilt object with an untoward input object.
  expect_error(cv_filter(omicsData = fdata_cv),
               paste("omicsData must be of class 'pepData', 'proData',",
                     "'metabData', 'lipidData', or 'nmrData'",
                     sep = ' '))
  
  # Run cv filter on the reduced data frame.
  filter_cv <- cv_filter(omicsData = pdata_cv)
  
  # Review the class for the cv_filter object.
  expect_s3_class(filter_cv,
                  c('cvFilt', 'data.frame'))
  
  # Check the dimensions of cv_filter.
  expect_equal(dim(filter_cv),
               c(150, 2))
  
  # Ensure the row sums are correct.
  expect_identical(round(filter_cv$CV_pooled, 3),
                   c(53.413, 44.024, 48.222, 51.462, 44.7, 55.745, 
                     80.581, 64.386, 44.421, 53.135, 28.859, 46.544, 
                     34.115, 48.714, 22.616, 105.602, 57.335, 12.66, 
                     7.951, 46.164, 72.192, 56.355, 54.41, 64.163, 
                     44.241, 75.603, 6.917, 59.268, 69.866, 33.917, 
                     35.984, 40.864, 53.374, 39.86, 8.21, 57.073, 
                     55.883, 23.939, 42.402, 49.515, 106.378, 60.911, 
                     52.559, 67.488, 31.147, 60.425, 44.164, 27.938, 
                     34.481, 24.016, 91.029, 77.308, 126.354, 33.327, 
                     69.152, 28.834, 48.358, 58.589, 70.466, 41.951, 
                     55.988, 39.315, 48.917, 77.449, 40.776, 61.521, 
                     74.333, 24.448, 57.545, 57.277, 33.157, 56.268, 
                     48.047, 50.174, 42.885, 43.122, 43.826, 72.966, 
                     34.606, 40.572, 62.338, 28.218, 63.366, 37.894, 
                     35.297, 28.522, 34.866, 26.911, 46.264, 63.554, 
                     25.862, 39.274, 54.486, 50.965, 66.826, 45.046, 
                     41.899, 41.899, 42.382, 38.602, 34.281, 46.25, 
                     50.955, 36.156, 31.934, 47.73, 43.722, 59.005, 
                     56.884, 50.58, 36.924, 73.401, 27.896, 47.313, 
                     43.419, 38.397, 27.599, 35.151, 31.254, 171.535, 
                     197.314, 171.27, 171.27, 186.315, 152.481, 225.596, 
                     157.757, 209.967, 213.983, 203.991, 209.497, 154.497, 
                     199.676, 162.603, 150.238, 151.622, 166.633, 162.054, 
                     174.03, 154.332, 171.218, 154.864, 164.243, 157.345, 
                     159.167, 160.373, 151.62, 178.301, 160.972, 153.229))
  
  # Inspect the attributes of the cvFilt object.
  expect_identical(attr(filter_cv, 'max_x_val'),
                   200)
  expect_equal(attr(filter_cv, 'tot_nas'),
               0)
  expect_false(attr(filter_cv, 'pooled'))
  
  # Test filter_cv with groups -------------------------------------------------
  
  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata_cv,
                                 main_effects = 'Condition')
  
  # Run cv filter on the reduced data frame with groups added.
  filter_gdf <- cv_filter(omicsData = pdata_gdf)
  
  # Review the class for the filter_cv object.
  expect_s3_class(filter_gdf,
                  c('cvFilt', 'data.frame'))
  
  # Check the dimensions of filter_cv.
  expect_equal(dim(filter_gdf),
               c(150, 2))
  
  # Ensure the row sums are correct.
  expect_identical(round(filter_gdf$CV_pooled, 3),
                   c(36.012, 32.768, 31.93, 47.891, 20.448, 33.369, 42.747, 
                     36.528, 19.768, 26.027, 28.859, 29.821, 20.875, 19.629, 
                     22.616, 44.259, 36.152, 12.66, 7.951, 29.045, 52.619, 
                     38.169, 26.964, 24.18, 34.24, 57.684, 6.917, 22.252, 
                     20.755, 30.455, 32.856, 28.478, 54.257, 28.876, 8.21, 
                     31.248, 33.82, 23.939, 29.86, 23.593, 31.444, 15.693, 
                     28.477, 23.84, 27.968, 11.561, 44.164, 27.526, 27.125, 
                     23.237, 54.696, 27.129, NaN, 33.327, 69.152, 21.461, 
                     49.465, 14.271, 49.843, 16.571, 31.854, 39.315, 48.917, 
                     45.667, 34.881, 22.933, 74.333, 24.193, 20.719, 36.68, 
                     29.829, 23.409, 28.334, 51.737, 42.851, 27.862, 41.238, 
                     20.086, 31.641, 42.73, 62.338, 28.521, 49.657, 25.383, 
                     16.76, 23.977, 30.713, 22.888, 42.954, 38.729, 20.915, 
                     39.274, 18.885, 26.121, 23.757, 25.38, 13.901, 13.901, 
                     22.142, 28.448, 30.159, 46.25, 24.448, 14.45, 27.964, 
                     24.852, 33.451, 25.986, 27.423, 41.608, 32.924, 41.563, 
                     26.448, 40.992, 31.513, 31.325, 26.42, 33.803, 28.45, 
                     71.728, 66.075, 63.178, 63.178, 94.102, 61.05, 149.144, 
                     64.717, 92.595, 89.103, 118.727, 156.935, 48.823, 92.817, 
                     77.15, 62.025, 90.059, 79.209, 59.814, 60.821, 45.93, 
                     67.457, 71.802, 100.432, 62.53, 63.508, 97.046, 59.292, 
                     54.676, 46.204, 66.356))
  
  # Inspect the attributes of the cvFilt object.
  expect_equal(round(attr(filter_gdf, 'max_x_val'), 3),
               134.544)
  expect_equal(attr(filter_gdf, 'tot_nas'),
               1)
  expect_true(attr(filter_gdf, 'pooled'))
  
  # Test applyFilt without groups ----------------------------------------------
  
  # Apply the filter without groups to the reduced peptide data set.
  filtered <- applyFilt(filter_object = filter_cv,
                        omicsData = pdata_cv,
                        cv_threshold = 150)

  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_cv, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata_cv, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata_cv),
                   class(filtered))

  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'cvFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                   150)
  expect_equal(attr(filtered, 'filters')[[1]]$filtered,
               c(6949002, 6955078, 6959151, 6959152, 6959154, 6967199, 6967232,
                 6967250, 6967253, 6967261, 6967276, 6967283, 6967345, 6967513,
                 6976251, 6976332, 6976356, 6976383, 6976619, 7039620, 7040132,
                 7045499, 7105633, 7217077, 7322754, 7439172, 7440084, 7451135,
                 7537762, 7537769, 7590799))
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(attr(filtered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered, 'data_info')$num_edata,
               119)
  expect_equal(attr(filtered, 'data_info')$num_miss_obs,
               228)
  expect_equal(round(attr(filtered, 'data_info')$prop_missing, 4),
               0.1597)
  expect_equal(attr(filtered, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered, 'data_info')$data_types)

  # Explore the meta_info attribute.
  expect_true(attr(filtered, 'meta_info')$meta_data)
  expect_equal(attr(filtered, 'meta_info')$num_emeta,
               65)

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(119, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(119, 4))
  
  # Test applyFilt with groups -------------------------------------------------
  
  # Apply the filter without groups to the reduced peptide data set.
  filtered_gdf <- applyFilt(filter_object = filter_gdf,
                            omicsData = pdata_gdf,
                            cv_threshold = 150)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_cv, 'cnames'),
                   attr(filtered_gdf, 'cnames'))
  expect_identical(attr(pdata_cv, 'check.names'),
                   attr(filtered_gdf, 'check.names'))
  expect_identical(class(pdata_cv),
                   class(filtered_gdf))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_gdf, 'filters')[[1]]$type,
               'cvFilt')
  expect_identical(attr(filtered_gdf, 'filters')[[1]]$threshold,
                   150)
  expect_equal(attr(filtered_gdf, 'filters')[[1]]$filtered,
               6967283)
  expect_true(is.na(attr(filtered_gdf, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered_gdf, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered_gdf, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_gdf, 'data_info')$num_edata,
               149)
  expect_equal(attr(filtered_gdf, 'data_info')$num_miss_obs,
               295)
  expect_equal(round(attr(filtered_gdf, 'data_info')$prop_missing, 4),
               0.165)
  expect_equal(attr(filtered_gdf, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered_gdf, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered_gdf, 'meta_info')$meta_data)
  expect_equal(attr(filtered_gdf, 'meta_info')$num_emeta,
               85)
  
  # Inspect the filtered_gdf e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_gdf$e_data),
               c(149, 13))
  expect_equal(dim(filtered_gdf$f_data),
               c(12, 2))
  expect_equal(dim(filtered_gdf$e_meta),
               c(149, 4))
  
})
