context('filter by IMD-ANOVA')

test_that('imdanova_filter and applyFilt produce the correct output',{
  
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
  
  # Group designate the pdata object.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = "Condition")
  
  # Create a pepData object without e_meta.
  no_emeta <- as.pepData(e_data = edata,
                         f_data = fdata,
                         edata_cname = "Mass_Tag_ID",
                         fdata_cname = "SampleID")
  
  # Run the group_designation function on pepData without e_meta.
  no_emeta <- group_designation(omicsData = no_emeta,
                                main_effects = "Condition")
  
  # Forge a pepData object with only one sample from Mock. This will create a
  # singleton group when group designating the data.
  # sg_1: singleton group with one non-singleton group.
  pdata_sg_1 <- as.pepData(e_data = edata[, 1:11],
                           f_data = fdata[1:10, ],
                           e_meta = emeta,
                           edata_cname = "Mass_Tag_ID",
                           fdata_cname = "SampleID",
                           emeta_cname = "Protein")
  
  # Run the group_designation function on pdata_sg_1.
  pdata_sg_1 <- group_designation(omicsData = pdata_sg_1,
                                  main_effects = "Condition")
  
  # Copy the pepData object created above. This object will have three groups
  # (two Infection and one Mock). The Mock group will be reduced to a singleton
  # group later. sg: singleton group.
  pdata_sg <- pdata
  
  # Break the Infection group into two groups.
  pdata_sg$f_data$Condition <- c(rep("Infection1", 4),
                                 rep("Infection2", 5),
                                 rep("Mock", 3))
  
  # Remove two of the Mock samples.
  pdata_sg$e_data <- pdata_sg$e_data[, -c(12, 13)]
  pdata_sg$f_data <- pdata_sg$f_data[-c(11, 12), ]
  
  # Run the group_designation function on the singleton group pepData object.
  pdata_sg <- group_designation(omicsData = pdata_sg,
                                main_effects = "Condition")
  
  # Create count standards -----------------------------------------------------
  
  # Create a function to count the number of non-missing values to be used in
  # connection with apply.
  count <- function (x) {
    
    present <- sum(!is.na(x))
    
    return (present)
    
  }
  
  # Fish out the indices for the Infection and Mock groups (no singletons).
  idx_i <- which(attr(pdata_gdf, "group_DF")$Group == "Infection")
  idx_m <- which(attr(pdata_gdf, "group_DF")$Group == "Mock")
  
  # Count the non-missing values for the Infection group and remove the names
  # from the vector.
  count_i <-  apply(pdata_gdf$e_data[, idx_i + 1], 1, count)
  names(count_i) <- NULL
  
  # Count the non-missing values for the Mock group and remove the names from
  # the vector.
  count_m <-  apply(pdata_gdf$e_data[, idx_m + 1], 1, count)
  names(count_m) <- NULL
  
  # Fish out the indices for the Infection and Mock groups (with singletons).
  idx_i1 <- which(attr(pdata_sg, "group_DF")$Group == "Infection1")
  idx_i2 <- which(attr(pdata_sg, "group_DF")$Group == "Infection2")
  idx_m2 <- which(attr(pdata_sg, "group_DF")$Group == "Mock")
  
  # Count the non-missing values for the Infection groups and remove the names
  # from the vector (with singletons).
  count_i1 <-  apply(pdata_sg$e_data[, idx_i1 + 1], 1, count)
  count_i2 <-  apply(pdata_sg$e_data[, idx_i2 + 1], 1, count)
  names(count_i1) <- NULL
  names(count_i2) <- NULL
  
  # Count the non-missing values for the Mock group and remove the names from
  # the vector (with singletons).
  count_m2 <-  apply(pdata_sg$e_data[, idx_m2 + 1, drop = FALSE], 1, count)
  names(count_m2) <- NULL
  
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
  
  # Run imdanova_filter with virtuous input objects.
  filter <- imdanova_filter(omicsData = pdata_gdf)
  
  # Review the class for the filter object.
  expect_s3_class(filter,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter object.
  expect_equal(attr(filter, "group_sizes"),
               data.frame(Group = c("Infection", "Mock"),
                          n_group = c(9, 3)))
  expect_identical(attr(filter, "nonsingleton_groups"),
                   c("Infection", "Mock"))
  
  # Check the dimensions of the filter data frame.
  expect_equal(dim(filter),
               c(150, 3))
  
  # Ensure the counts for the infection group are correct.
  expect_identical(filter[, 2],
                   count_i)
  
  # Verify the counts for the mock group are correct.
  expect_identical(filter[, 3],
                   count_m)
  
  # Run imdanova_filter without e_meta. ne: no e_meta
  filter_ne <- imdanova_filter(omicsData = no_emeta)
  
  # Review the class for the filter_ne object.
  expect_s3_class(filter_ne,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter_ne object.
  expect_equal(attr(filter_ne, "group_sizes"),
               data.frame(Group = c("Infection", "Mock"),
                          n_group = c(9, 3)))
  expect_identical(attr(filter_ne, "nonsingleton_groups"),
                   c("Infection", "Mock"))
  
  # Check the dimensions of the filter_ne data frame.
  expect_equal(dim(filter_ne),
               c(150, 3))
  
  # Ensure the counts for the infection group are correct.
  expect_identical(filter_ne[, 2],
                   count_i)
  
  # Verify the counts for the mock group are correct.
  expect_identical(filter_ne[, 3],
                   count_m)
  
  # Test imdanova_filter with singletons ---------------------------------------
  
  # Run imdanova_filter with virtuous input objects.
  filter_sg_1 <- imdanova_filter(omicsData = pdata_sg_1)
  
  # Review the class for the filter_sg object.
  expect_s3_class(filter_sg_1,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter_sg object.
  expect_equal(attr(filter_sg_1, "group_sizes"),
               data.frame(Group = c("Infection", "Mock"),
                          n_group = c(9, 1)))
  expect_identical(attr(filter_sg_1, "nonsingleton_groups"),
                   "Infection")
  
  # Check the dimensions of the filter data frame.
  expect_equal(dim(filter_sg_1),
               c(150, 2))
  
  # Ensure the counts for the Infection group are correct.
  expect_identical(filter_sg_1[, 2],
                   count_i)

  # Run imdanova_filter with virtuous input objects and one singleton group.
  filter_sg <- imdanova_filter(omicsData = pdata_sg)
  
  # Review the class for the filter_sg object.
  expect_s3_class(filter_sg,
                  c('imdanovaFilt', 'data.frame'))
  
  # Inspect the attributes of the filter_sg object.
  expect_equal(attr(filter_sg, "group_sizes"),
               data.frame(Group = c("Infection1", "Infection2", "Mock"),
                          n_group = c(4, 5, 1)))
  expect_identical(attr(filter_sg, "nonsingleton_groups"),
                   c("Infection1", "Infection2"))
  
  # Check the dimensions of the filter data frame.
  expect_equal(dim(filter_sg),
               c(150, 3))
  
  # Ensure the counts for the Infection groups are correct.
  expect_identical(filter_sg[, 2],
                   count_i1)
  expect_identical(filter_sg[, 3],
                   count_i2)
  
  # Create filter standards ----------------------------------------------------
  
  # No singleton groups ---------------
  
  size <- attributes(filter)$group_sizes
  
  # Create the anova, gtest, and combined standards.
  s_anova <- anova_filter(nonmiss_per_group = list(nonmiss_totals = filter,
                                                   group_sizes = size),
                          min_nonmiss_anova = 3,
                          cname_id =  get_edata_cname(pdata_gdf))
  
  s_gtest <- gtest_filter(nonmiss_per_group = list(nonmiss_totals = filter,
                                                   group_sizes = size),
                          omicsData = pdata_gdf,
                          min_nonmiss_gtest = 3)
  
  s_both <- intersect(s_anova, s_gtest)
  
  # With singleton groups ---------------
  
  size_sg <- attributes(filter_sg)$group_sizes
  
  # Create the anova, gtest, and combined standards.
  s_anova_sg <- anova_filter(
    nonmiss_per_group = list(nonmiss_totals = filter_sg,
                             group_sizes = size_sg),
    min_nonmiss_anova = 3,
    cname_id =  get_edata_cname(pdata_sg)
  )
  
  s_gtest_sg <- gtest_filter(
    nonmiss_per_group = list(nonmiss_totals = filter_sg,
                             group_sizes = size_sg),
    omicsData = pdata_sg,
    min_nonmiss_gtest = 3
  )
  
  s_both_sg <- intersect(s_anova_sg, s_gtest_sg)
  
  # Test applyFilt anova no singletons -----------------------------------------
  
  # Filter the reduced pepData object using anova.
  aFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata_gdf,
                         min_nonmiss_anova = 3, 
                         min_nonmiss_gtest = NULL, 
                         remove_singleton_groups = FALSE)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_gdf, 'cnames'),
                   attr(aFiltered, 'cnames'))
  expect_identical(attr(pdata_gdf, 'check.names'),
                   attr(aFiltered, 'check.names'))
  expect_identical(class(pdata_gdf),
                   class(aFiltered))
  
  # Investigate the filters attribute.
  expect_equal(attr(aFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(aFiltered, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = NA))
  expect_identical(attr(aFiltered, 'filters')[[1]]$filtered,
                   s_anova)
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
                   setdiff(x = pdata_gdf$e_data[, 1],
                           y = s_anova))
  expect_null(attr(aFiltered, "imdanova")$test_with_gtest)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(aFiltered$e_data),
               c(110, 13))
  expect_equal(dim(aFiltered$f_data),
               c(12, 2))
  expect_equal(dim(aFiltered$e_meta),
               c(110, 4))
  
  # Test applyFilt gtest no singletons -----------------------------------------
  
  # Filter the reduced pepData object using gtest.
  gFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata_gdf,
                         min_nonmiss_anova = NULL, 
                         min_nonmiss_gtest = 3, 
                         remove_singleton_groups = FALSE)
  
  # Investigate the filters attribute.
  expect_equal(attr(gFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(gFiltered, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = NA,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(gFiltered, 'filters')[[1]]$filtered,
                   s_gtest)
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
                   setdiff(x = pdata_gdf$e_data[, 1],
                           y = s_gtest))
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(gFiltered$e_data),
               c(136, 13))
  expect_equal(dim(gFiltered$f_data),
               c(12, 2))
  expect_equal(dim(gFiltered$e_meta),
               c(136, 4))
  
  # Test applyFilt combined no singletons --------------------------------------
  
  # Filter the reduced pepData object using both anova and gtest.
  bFiltered <- applyFilt(filter_object = filter,
                         omicsData = pdata_gdf,
                         min_nonmiss_anova = 3, 
                         min_nonmiss_gtest = 3, 
                         remove_singleton_groups = FALSE)
  
  # Investigate the filters attribute.
  expect_equal(attr(bFiltered, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(bFiltered, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(bFiltered, 'filters')[[1]]$filtered,
                   s_both)
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
                   setdiff(x = pdata_gdf$e_data[, 1],
                           y = s_anova))
  expect_identical(attr(bFiltered, "imdanova")$test_with_gtest,
                   setdiff(x = pdata_gdf$e_data[, 1],
                           y = s_gtest))
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(bFiltered$e_data),
               c(136, 13))
  expect_equal(dim(bFiltered$f_data),
               c(12, 2))
  expect_equal(dim(bFiltered$e_meta),
               c(136, 4))
  
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
                   attr(aFiltered, 'filters')[[1]]$threshold)
  expect_identical(attr(filtered_ne, 'filters')[[1]]$filtered,
                   attr(aFiltered, 'filters')[[1]]$filtered)
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
  expect_identical(dim(filtered_ne$e_data),
                   dim(aFiltered$e_data))
  expect_identical(dim(filtered_ne$f_data),
                   dim(aFiltered$f_data))
  expect_null(filtered_ne$e_meta)

  # Test applyFilt anova with singletons ---------------------------------------
  
  # Filter the reduced pepData object with removing singleton groups.
  expect_error(suppressMessages(applyFilt(filter_object = filter_sg_1,
                                          omicsData = pdata_sg_1,
                                          min_nonmiss_anova = 3, 
                                          min_nonmiss_gtest = NULL, 
                                          remove_singleton_groups = TRUE)),
               paste("An IMD-ANOVA filter cannot be used because there is only",
                     "one non-singleton group.",
                     sep = " "))
  
  # Filter the reduced pepData object using anova.
  expect_message(aFiltered_sg <- applyFilt(filter_object = filter_sg,
                                           omicsData = pdata_sg,
                                           min_nonmiss_anova = 3, 
                                           min_nonmiss_gtest = NULL, 
                                           remove_singleton_groups = TRUE),
                 paste("You have specified remove_single_groups = TRUE, so we",
                       "have removed the following sample\\(s\\) that",
                       "correspond to singleton groups prior to proceeding",
                       "with the IMD-ANOVA filter: Mock1",
                       sep = " "))
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, 'cnames'),
                   attr(aFiltered_sg, 'cnames'))
  expect_identical(attr(pdata_sg, 'check.names'),
                   attr(aFiltered_sg, 'check.names'))
  expect_identical(class(pdata_sg),
                   class(aFiltered_sg))
  
  # Investigate the filters attribute.
  expect_equal(attr(aFiltered_sg, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(aFiltered_sg, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = NA))
  expect_identical(attr(aFiltered_sg, 'filters')[[1]]$filtered,
                   list(samples = "Mock1",
                        biomolecules = s_anova_sg))
  expect_equal(attr(aFiltered_sg, 'filters')[[1]]$method,
               "anova")
  
  # Examinate the data_info attribute.
  expect_equal(attr(aFiltered_sg, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(aFiltered_sg, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(aFiltered_sg, 'data_info')$num_edata,
               113)
  expect_equal(attr(aFiltered_sg, 'data_info')$num_miss_obs,
               23)
  expect_equal(round(attr(aFiltered_sg, 'data_info')$prop_missing, 4),
               0.0226)
  expect_equal(attr(aFiltered_sg, 'data_info')$num_samps,
               9)
  expect_null(attr(aFiltered_sg, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(aFiltered_sg, 'meta_info')$meta_data)
  expect_equal(attr(aFiltered_sg, 'meta_info')$num_emeta,
               65)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(aFiltered_sg, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes")[1:2, ])
  expect_identical(attr(aFiltered_sg, "imdanova")$test_with_anova,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_anova_sg))
  expect_null(attr(aFiltered_sg, "imdanova")$test_with_gtest)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(aFiltered_sg$e_data),
               c(113, 10))
  expect_equal(dim(aFiltered_sg$f_data),
               c(9, 2))
  expect_equal(dim(aFiltered_sg$e_meta),
               c(113, 4))
  
  # Test applyFilt gtest with singletons ---------------------------------------
  
  # Filter the reduced pepData object using gtest.
  expect_message(gFiltered_sg <- applyFilt(filter_object = filter_sg,
                                           omicsData = pdata_sg,
                                           min_nonmiss_anova = NULL, 
                                           min_nonmiss_gtest = 3, 
                                           remove_singleton_groups = TRUE),
                 paste("You have specified remove_single_groups = TRUE, so we",
                       "have removed the following sample\\(s\\) that",
                       "correspond to singleton groups prior to proceeding",
                       "with the IMD-ANOVA filter: Mock1",
                       sep = " "))
  
  # Investigate the filters attribute.
  expect_equal(attr(gFiltered_sg, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(gFiltered_sg, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = NA,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(gFiltered_sg, 'filters')[[1]]$filtered,
                   list(samples = "Mock1",
                        biomolecules = s_gtest_sg))
  expect_equal(attr(gFiltered_sg, 'filters')[[1]]$method,
               "gtest")
  
  # Examinate the data_info attribute.
  expect_equal(attr(gFiltered_sg, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(gFiltered_sg, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(gFiltered_sg, 'data_info')$num_edata,
               128)
  expect_equal(attr(gFiltered_sg, 'data_info')$num_miss_obs,
               87)
  expect_equal(round(attr(gFiltered_sg, 'data_info')$prop_missing, 4),
               0.0755)
  expect_equal(attr(gFiltered_sg, 'data_info')$num_samps,
               9)
  expect_null(attr(gFiltered_sg, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(gFiltered_sg, 'meta_info')$meta_data)
  expect_equal(attr(gFiltered_sg, 'meta_info')$num_emeta,
               71)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(gFiltered_sg, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes")[1:2, ])
  expect_null(attr(gFiltered_sg, "imdanova")$test_with_anova)
  expect_identical(attr(gFiltered_sg, "imdanova")$test_with_gtest,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_gtest_sg))
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(gFiltered_sg$e_data),
               c(128, 10))
  expect_equal(dim(gFiltered_sg$f_data),
               c(9, 2))
  expect_equal(dim(gFiltered_sg$e_meta),
               c(128, 4))
  
  # Test applyFilt combined with singletons ------------------------------------
  
  # Filter the reduced pepData object using both anova and gtest.
  expect_message(bFiltered_sg <- applyFilt(filter_object = filter_sg,
                                           omicsData = pdata_sg,
                                           min_nonmiss_anova = 3, 
                                           min_nonmiss_gtest = 3, 
                                           remove_singleton_groups = TRUE),
                 paste("You have specified remove_single_groups = TRUE, so we",
                       "have removed the following sample\\(s\\) that",
                       "correspond to singleton groups prior to proceeding",
                       "with the IMD-ANOVA filter: Mock1",
                       sep = " "))
  
  # Investigate the filters attribute.
  expect_equal(attr(bFiltered_sg, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(bFiltered_sg, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(bFiltered_sg, 'filters')[[1]]$filtered,
                   list(samples = "Mock1",
                        biomolecules = s_both_sg))
  expect_equal(attr(bFiltered_sg, 'filters')[[1]]$method,
               "combined")
  
  # Examinate the data_info attribute.
  expect_equal(attr(bFiltered_sg, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(bFiltered_sg, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(bFiltered_sg, 'data_info')$num_edata,
               128)
  expect_equal(attr(bFiltered_sg, 'data_info')$num_miss_obs,
               87)
  expect_equal(round(attr(bFiltered_sg, 'data_info')$prop_missing, 4),
               0.0755)
  expect_equal(attr(bFiltered_sg, 'data_info')$num_samps,
               9)
  expect_null(attr(bFiltered_sg, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(bFiltered_sg, 'meta_info')$meta_data)
  expect_equal(attr(bFiltered_sg, 'meta_info')$num_emeta,
               71)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(bFiltered_sg, "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes")[1:2, ])
  expect_identical(attr(bFiltered_sg, "imdanova")$test_with_anova,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_anova_sg))
  expect_identical(attr(bFiltered_sg, "imdanova")$test_with_gtest,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_gtest_sg))
  
  # Inspectate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(bFiltered_sg$e_data),
               c(128, 10))
  expect_equal(dim(bFiltered_sg$f_data),
               c(9, 2))
  expect_equal(dim(bFiltered_sg$e_meta),
               c(128, 4))

  # Test applyFilt with ignoring singleton groups ------------------------------
  
  # Filter the reduced pepData object using anova.
  # _f: remove_singleton_groups is set to FALSE.
  expect_message(aFiltered_sg_f <- applyFilt(filter_object = filter_sg,
                                             omicsData = pdata_sg,
                                             min_nonmiss_anova = 3, 
                                             min_nonmiss_gtest = NULL, 
                                             remove_singleton_groups = FALSE),
                 paste("You have specified remove_singleton_groups = FALSE, so",
                       "the sample\\(s\\) corresponding to the singleton",
                       "group\\(s\\) were not utilized in the IMD-ANOVA filter",
                       "and will be retained in the resulting omicsData",
                       "object.",
                       sep = " "))
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, 'cnames'),
                   attr(aFiltered_sg_f, 'cnames'))
  expect_identical(attr(pdata_sg, 'check.names'),
                   attr(aFiltered_sg_f, 'check.names'))
  expect_identical(class(pdata_sg),
                   class(aFiltered_sg_f))
  
  # Investigate the filters attribute.
  expect_equal(attr(aFiltered_sg_f, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(aFiltered_sg_f, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = NA))
  expect_identical(attr(aFiltered_sg_f, 'filters')[[1]]$filtered,
                   s_anova_sg)
  expect_equal(attr(aFiltered_sg_f, 'filters')[[1]]$method,
               "anova")
  
  # Examinate the data_info attribute.
  expect_equal(attr(aFiltered_sg_f, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(aFiltered_sg_f, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(aFiltered_sg_f, 'data_info')$num_edata,
               113)
  expect_equal(attr(aFiltered_sg_f, 'data_info')$num_miss_obs,
               31)
  expect_equal(round(attr(aFiltered_sg_f, 'data_info')$prop_missing, 4),
               0.0274)
  expect_equal(attr(aFiltered_sg_f, 'data_info')$num_samps,
               10)
  expect_null(attr(aFiltered_sg_f, 'data_info')$data_types)
  
  # Explorate the meta_info attribute.
  expect_true(attr(aFiltered_sg_f, 'meta_info')$meta_data)
  expect_equal(attr(aFiltered_sg_f, 'meta_info')$num_emeta,
               65)
  
  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(aFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)
  expect_identical(attr(aFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes"))
  expect_identical(attr(aFiltered_sg_f, "imdanova")$test_with_anova,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_anova_sg))
  expect_null(attr(aFiltered_sg_f, "imdanova")$test_with_gtest)
  
  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(aFiltered_sg_f$e_data),
               c(113, 11))
  expect_equal(dim(aFiltered_sg_f$f_data),
               c(10, 2))
  expect_equal(dim(aFiltered_sg_f$e_meta),
               c(113, 4))
  
  # Test applyFilt with a molecule filter applied first ------------------------
  
  # Run molecule filter on the reduced data frame.
  filter_mol <- molecule_filter(omicsData = pdata_sg)
  
  # Apply the filter to the reduced peptide data set.
  filtered_mol <- applyFilt(filter_object = filter_mol,
                            omicsData = pdata_sg,
                            min_num = 2)
  
  # Group designate the filtered_mol object.
  filtered_mol_gdf <- group_designation(filtered_mol,
                                        main_effects = "Condition")
  
  # Apply the IMD-ANOVA filter after the molecule filter.
  filter_imd <- imdanova_filter(omicsData = filtered_mol_gdf)
  
  # applyFilt the previously filtered data.
  filtered_imd <- suppressMessages(applyFilt(filter_object = filter_imd,
                                             omicsData = filtered_mol_gdf,
                                             min_nonmiss_anova = 3, 
                                             min_nonmiss_gtest = NULL, 
                                             remove_singleton_groups = TRUE))
  
  # Check the number of filters in the filters attribute is correct.
  expect_equal(length(attr(filtered_imd, "filters")),
               2)
  
  # Ensure the molecule filter attribute is correct.
  expect_identical(attr(filtered_mol, "filters")[[1]],
                   attr(filtered_imd, "filters")[[1]])
  
  # Inspect the filters attribute for the IMD-ANOVA filter
  expect_equal(attr(filtered_imd, 'filters')[[2]]$type,
               'imdanovaFilt')
  expect_identical(attr(filtered_imd, 'filters')[[2]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = NA))
  expect_equal(attr(filtered_imd, 'filters')[[1]]$filtered,
               c(1024, 11939, 15714, 16636, 21168, 976139, 6637724, 6654733,
                 6769231, 6769844, 6832528, 6901575, 6909787, 6934326))
  expect_equal(attr(aFiltered_sg_f, 'filters')[[1]]$method,
               "anova")

  # Check the group_sizes attribute. This is done at the end of the tests
  # because the singleton groups are removed from the group_sizes attribute
  # within the applyFilt function. However, in the test script this has to be
  # done manually after the filter_sg object has been called in the applyFilt
  # function for all tests with singleton groups.
  
  # Remove the singleton group from the filter_sg group_sizes attribute.
  attr(filter_sg, "group_sizes") <- attr(filter_sg, "group_sizes")[1:2, ]
  
  expect_identical(attr(aFiltered_sg,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)
  expect_identical(attr(gFiltered_sg,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)
  expect_identical(attr(bFiltered_sg,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)

})
