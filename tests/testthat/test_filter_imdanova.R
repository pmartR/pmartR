context('filter by IMD-ANOVA')

test_that('imdanova_filter and applyFilt produce the correct output',{
  
  # Load data and standards and prepare omicsData objects ----------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  load(system.file('testdata',
                   'imdanova_standards.RData',
                   package = 'pmartR'))
  
  # Create a pepData object with the reduced data set.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = "Mass_Tag_ID",
                      fdata_cname = "SampleID",
                      emeta_cname = "Protein")
  
  # The group_designation function will be run on pdata further down the script.
  
  # Create a pepData object without e_meta.
  no_emeta <- as.pepData(e_data = edata,
                         f_data = fdata,
                         edata_cname = "Mass_Tag_ID",
                         fdata_cname = "SampleID")
  
  # Run the group_designation function on pepData without e_meta.
  no_emeta <- group_designation(omicsData = no_emeta,
                                main_effects = "Condition")
  
  # Forge a pepData object with only one sample from Mock. This will create a
  # singleton group when group designating the data. sg: singleton group.
  pdata_sg <- as.pepData(e_data = edata[, 1:11],
                         f_data = fdata[1:10, ],
                         e_meta = emeta,
                         edata_cname = "Mass_Tag_ID",
                         fdata_cname = "SampleID",
                         emeta_cname = "Protein")
  
  # Run the group_designation function on pdata_sg.
  pdata_sg <- group_designation(omicsData = pdata_sg,
                                main_effects = "Condition")
  
  # Test imdanova_filter no singleton groups -----------------------------------
  
  # Try creating a imdanovaFilt object with unholy input objects.
  expect_error(imdanova_filter(omicsData = emeta),
               paste("",
                     sep = " "))
  
  # Try creating a imdanovaFilt object without the group_DF attribute.
  expect_error(imdanova_filter(pdata),
               paste("omicsData must contain attribute information for",
                     "'group_DF'. See documentation for group_designation",
                     "function for more information.",
                     sep = " "))
  
  # Run the group_designation function on pdata.
  pdata <- group_designation(omicsData = pdata,
                             main_effects = "Condition")
  
  # Run imdanova_filter with virtuous input objects.
  filter <- imdanova_filter(omicsData = pdata)
  
  # Review the class for the filter object.
  expect_s3_class(filter,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter object.
  expect_identical(attr(filter, "group_sizes"),
                   attr(ifStandard, "group_sizes"))
  expect_identical(attr(filter, "nonsingleton_groups"),
                   attr(ifStandard, "nonsingleton_groups"))
  
  # Check the dimensions of the filter data frame.
  expect_equal(dim(filter),
               c(150, 3))
  
  # Ensure the counts for the infection group are correct.
  expect_identical(filter[, 2],
                   ifStandard[, 2])
  
  # Verify the counts for the mock group are correct.
  expect_identical(filter[, 3],
                   ifStandard[, 3])
  
  # Run imdanova_filter without e_meta. ne: no e_meta
  filter_ne <- imdanova_filter(omicsData = no_emeta)
  
  # Review the class for the filter_ne object.
  expect_s3_class(filter_ne,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter_ne object.
  expect_identical(attr(filter_ne, "group_sizes"),
                   attr(ifStandard, "group_sizes"))
  expect_identical(attr(filter_ne, "nonsingleton_groups"),
                   attr(ifStandard, "nonsingleton_groups"))
  
  # Check the dimensions of the filter_ne data frame.
  expect_equal(dim(filter_ne),
               c(150, 3))
  
  # Ensure the counts for the infection group are correct.
  expect_identical(filter_ne[, 2],
                   ifStandard[, 2])
  
  # Verify the counts for the mock group are correct.
  expect_identical(filter_ne[, 3],
                   ifStandard[, 3])
  
  # Test imdanova_filter with singletons ---------------------------------------
  
  # Run imdanova_filter with virtuous input objects.
  filter_sg <- imdanova_filter(omicsData = pdata_sg)
  
  # Review the class for the filter_sg object.
  expect_s3_class(filter_sg,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter_sg object.
  expect_identical(attr(filter_sg, "group_sizes"),
                   attr(ifStandard_sg, "group_sizes"))
  expect_identical(attr(filter_sg, "nonsingleton_groups"),
                   attr(ifStandard_sg, "nonsingleton_groups"))
  
  # Check the dimensions of the filter data frame.
  expect_equal(dim(filter_sg),
               c(150, 2))
  
  # Ensure the counts for the infection group are correct.
  expect_identical(filter_sg[, 2],
                   ifStandard_sg[, 2])
  
  # Test applyFilt anova -------------------------------------------------------
  
  # Filter the reduced pepData object using anova.
  aFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata,
                         min_nonmiss_anova = 3, 
                         min_nonmiss_gtest = NULL, 
                         remove_singleton_groups = FALSE)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(aFiltered, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(aFiltered, 'check.names'))
  expect_identical(class(pdata),
                   class(aFiltered))
  
  # Investigate the filters attribute.
  expect_equal(attr(aFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(aFiltered, 'filters')[[1]]$threshold,
                   attr(afStandard, 'filters')[[1]]$threshold)
  expect_identical(attr(aFiltered, 'filters')[[1]]$filtered,
                   attr(afStandard, 'filters')[[1]]$filtered)
  expect_equal(attr(aFiltered, 'filters')[[1]]$method,
               "anova")
  
  # Examinate the data_info attribute.
  expect_equal(attr(aFiltered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(aFiltered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(aFiltered, 'data_info')$num_edata,
               110)
  expect_equal(attr(aFiltered, 'data_info')$num_miss_obs,
               50)
  expect_equal(round(attr(aFiltered, 'data_info')$prop_missing, 4),
               0.0379)
  expect_equal(attr(aFiltered, 'data_info')$num_samps,
               12)
  expect_null(attr(aFiltered, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(aFiltered, 'meta_info')$meta_data)
  expect_equal(attr(aFiltered, 'meta_info')$num_emeta,
               61)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(aFiltered, "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter)
  expect_identical(attr(aFiltered, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter, "group_sizes"))
  expect_identical(attr(aFiltered, "imdanova")$test_with_anova,
                   attr(afStandard, "imdanova")$test_with_anova)
  expect_null(attr(aFiltered, "imdanova")$test_with_gtest)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_identical(dim(aFiltered$e_data),
                   dim(afStandard$e_data))
  expect_identical(dim(aFiltered$f_data),
                   dim(afStandard$f_data))
  expect_identical(dim(aFiltered$e_meta),
                   dim(afStandard$e_meta))
  
  # Test applyFilt gtest -------------------------------------------------------
  
  # Filter the reduced pepData object using gtest.
  gFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata,
                         min_nonmiss_anova = NULL, 
                         min_nonmiss_gtest = 3, 
                         remove_singleton_groups = FALSE)
  
  # Investigate the filters attribute.
  expect_equal(attr(gFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(gFiltered, 'filters')[[1]]$threshold,
                   attr(gfStandard, 'filters')[[1]]$threshold)
  expect_identical(attr(gFiltered, 'filters')[[1]]$filtered,
                   attr(gfStandard, 'filters')[[1]]$filtered)
  expect_equal(attr(gFiltered, 'filters')[[1]]$method,
               "gtest")
  
  # Examinate the data_info attribute.
  expect_equal(attr(gFiltered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(gFiltered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(gFiltered, 'data_info')$num_edata,
               136)
  expect_equal(attr(gFiltered, 'data_info')$num_miss_obs,
               194)
  expect_equal(round(attr(gFiltered, 'data_info')$prop_missing, 4),
               0.1189)
  expect_equal(attr(gFiltered, 'data_info')$num_samps,
               12)
  expect_null(attr(gFiltered, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(gFiltered, 'meta_info')$meta_data)
  expect_equal(attr(gFiltered, 'meta_info')$num_emeta,
               75)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(gFiltered, "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter)
  expect_identical(attr(gFiltered, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter, "group_sizes"))
  expect_null(attr(gFiltered, "imdanova")$test_with_anova)
  expect_identical(attr(gFiltered, "imdanova")$test_with_gtest,
                   attr(gfStandard, "imdanova")$test_with_gtest)
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_identical(dim(gFiltered$e_data),
                   dim(gfStandard$e_data))
  expect_identical(dim(gFiltered$f_data),
                   dim(gfStandard$f_data))
  expect_identical(dim(gFiltered$e_meta),
                   dim(gfStandard$e_meta))
  
  # Test applyFilt combined ----------------------------------------------------
  
  # Filter the reduced pepData object using both anova and gtest.
  bFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata,
                         min_nonmiss_anova = 3, 
                         min_nonmiss_gtest = 3, 
                         remove_singleton_groups = FALSE)
  
  # Investigate the filters attribute.
  expect_equal(attr(bFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(bFiltered, 'filters')[[1]]$threshold,
                   attr(bfStandard, 'filters')[[1]]$threshold)
  expect_identical(attr(bFiltered, 'filters')[[1]]$filtered,
                   attr(bfStandard, 'filters')[[1]]$filtered)
  expect_equal(attr(bFiltered, 'filters')[[1]]$method,
               "combined")
  
  # Examinate the data_info attribute.
  expect_equal(attr(bFiltered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(bFiltered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(bFiltered, 'data_info')$num_edata,
               136)
  expect_equal(attr(bFiltered, 'data_info')$num_miss_obs,
               194)
  expect_equal(round(attr(bFiltered, 'data_info')$prop_missing, 4),
               0.1189)
  expect_equal(attr(bFiltered, 'data_info')$num_samps,
               12)
  expect_null(attr(bFiltered, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(bFiltered, 'meta_info')$meta_data)
  expect_equal(attr(bFiltered, 'meta_info')$num_emeta,
               75)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(bFiltered, "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter)
  expect_identical(attr(bFiltered, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter, "group_sizes"))
  expect_identical(attr(bFiltered, "imdanova")$test_with_anova,
                   attr(bfStandard, "imdanova")$test_with_anova)
  expect_identical(attr(bFiltered, "imdanova")$test_with_gtest,
                   attr(bfStandard, "imdanova")$test_with_gtest)
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_identical(dim(bFiltered$e_data),
                   dim(bfStandard$e_data))
  expect_identical(dim(bFiltered$f_data),
                   dim(bfStandard$f_data))
  expect_identical(dim(bFiltered$e_meta),
                   dim(bfStandard$e_meta))
  
  # Test applyFilt no e_meta ---------------------------------------------------
  
  # Filter the reduced pepData object using anova without an e_meta data frame.
  filtered_ne <- applyFilt(filter_object = filter_ne,
                           omicsData = no_emeta,
                           min_nonmiss_anova = 3, 
                           min_nonmiss_gtest = NULL, 
                           remove_singleton_groups = FALSE)
  
  # Investigate the filters attribute.
  expect_equal(attr(filtered_ne, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(filtered_ne, 'filters')[[1]]$threshold,
                   attr(afStandard, 'filters')[[1]]$threshold)
  expect_identical(attr(filtered_ne, 'filters')[[1]]$filtered,
                   attr(afStandard, 'filters')[[1]]$filtered)
  expect_equal(attr(filtered_ne, 'filters')[[1]]$method,
               "anova")
  
  # Examinate the data_info attribute.
  expect_equal(attr(filtered_ne, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered_ne, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_ne, 'data_info')$num_edata,
               110)
  expect_equal(attr(filtered_ne, 'data_info')$num_miss_obs,
               50)
  expect_equal(round(attr(filtered_ne, 'data_info')$prop_missing, 4),
               0.0379)
  expect_equal(attr(filtered_ne, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered_ne, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_false(attr(filtered_ne, 'meta_info')$meta_data)
  expect_null(attr(filtered_ne, 'meta_info')$num_emeta)
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_identical(dim(aFiltered$e_data),
                   dim(afStandard$e_data))
  expect_identical(dim(aFiltered$f_data),
                   dim(afStandard$f_data))
  expect_null(filtered_ne$e_meta)
  
  # Test applyFilt with singletons ---------------------------------------------
  
  # Filter the reduced pepData object with removing singleton groups.
  expect_error(suppressMessages(applyFilt(filter_object = filter_sg,
                                          omicsData = pdata_sg,
                                          min_nonmiss_anova = 3, 
                                          min_nonmiss_gtest = NULL, 
                                          remove_singleton_groups = TRUE)),
               paste("An IMD-ANOVA filter cannot be used because there is only",
                     "one non-singleton group.",
                     sep = " "))

})
