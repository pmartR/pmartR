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

  # Create a vector of rows that pass the anova criteria.
  p_anova <- which(rowSums(data.frame(infection = count_i,
                                      mock = count_m) >= 3) < 2)

  # Extract the peptide IDs corresponding to the rows below the anova threshold.
  s_anova <- as.character(pdata_gdf$e_data[p_anova, 1])

  # Fashion a vector of rows that pass the gtest criteria.
  p_gtest <- apply(data.frame(infection = count_i,
                              mock = count_m),
                   1,
                   function (x) max(x) >= 3)

  # Extract the peptide IDs corresponding to the rows below the gtest threshold.
  s_gtest <- as.character(pdata_gdf$e_data[!p_gtest, 1])

  # Find the peptide IDs that occur in both the anova and gtest sets.
  s_both <- intersect(s_anova, s_gtest)

  # With singleton groups ---------------

  size_sg <- attributes(filter_sg)$group_sizes

  # Create a vector of rows that pass the anova criteria for data with
  # singleton groups.
  p_anova_sg <- which(rowSums(data.frame(infection1 = count_i1,
                                         infection2 = count_i2,
                                         mock = count_m2) >= 3) < 2)

  # Extract the peptide IDs corresponding to the rows below the anova threshold.
  s_anova_sg <- as.character(pdata_sg$e_data[p_anova_sg, 1])

  # Fashion a vector of rows that pass the gtest criteria for data with
  # singleton groups.
  p_gtest_sg <- apply(data.frame(infection1 = count_i1,
                                 infection2 = count_i2,
                                 mock = count_m2),
                      1,
                      function (x) max(x) >= 3)

  # Extract the peptide IDs corresponding to the rows below the gtest threshold.
  s_gtest_sg <- as.character(pdata_sg$e_data[!p_gtest_sg, 1])

  # Find the peptide IDs that occur in both the anova and gtest sets.
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
  expect_equal(
    attr(aFiltered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(aFiltered$e_data[, 1])),
         num_miss_obs = sum(is.na(aFiltered$e_data)),
         prop_missing = (sum(is.na(aFiltered$e_data)) /
                           prod(dim(aFiltered$e_data[, -1]))),
         num_samps = ncol(aFiltered$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(aFiltered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(aFiltered$e_meta$Protein)))
  )

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
  expect_equal(
    attr(gFiltered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(gFiltered$e_data[, 1])),
         num_miss_obs = sum(is.na(gFiltered$e_data)),
         prop_missing = (sum(is.na(gFiltered$e_data)) /
                           prod(dim(gFiltered$e_data[, -1]))),
         num_samps = ncol(gFiltered$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(gFiltered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(gFiltered$e_meta$Protein)))
  )

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
  expect_equal(
    attr(bFiltered, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(bFiltered$e_data[, 1])),
         num_miss_obs = sum(is.na(bFiltered$e_data)),
         prop_missing = (sum(is.na(bFiltered$e_data)) /
                           prod(dim(bFiltered$e_data[, -1]))),
         num_samps = ncol(bFiltered$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(bFiltered, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(bFiltered$e_meta$Protein)))
  )

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
  expect_equal(
    attr(filtered_ne, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_ne$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered_ne$e_data)),
         prop_missing = (sum(is.na(filtered_ne$e_data)) /
                           prod(dim(filtered_ne$e_data[, -1]))),
         num_samps = ncol(filtered_ne$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(filtered_ne, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )

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
  expect_equal(
    attr(aFiltered_sg, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(aFiltered_sg$e_data[, 1])),
         num_miss_obs = sum(is.na(aFiltered_sg$e_data)),
         prop_missing = (sum(is.na(aFiltered_sg$e_data)) /
                           prod(dim(aFiltered_sg$e_data[, -1]))),
         num_samps = ncol(aFiltered_sg$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(aFiltered_sg, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(aFiltered_sg$e_meta$Protein)))
  )

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
  expect_equal(
    attr(gFiltered_sg, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(gFiltered_sg$e_data[, 1])),
         num_miss_obs = sum(is.na(gFiltered_sg$e_data)),
         prop_missing = (sum(is.na(gFiltered_sg$e_data)) /
                           prod(dim(gFiltered_sg$e_data[, -1]))),
         num_samps = ncol(gFiltered_sg$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(gFiltered_sg, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(gFiltered_sg$e_meta$Protein)))
  )

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
  expect_equal(
    attr(bFiltered_sg, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(bFiltered_sg$e_data[, 1])),
         num_miss_obs = sum(is.na(bFiltered_sg$e_data)),
         prop_missing = (sum(is.na(bFiltered_sg$e_data)) /
                           prod(dim(bFiltered_sg$e_data[, -1]))),
         num_samps = ncol(bFiltered_sg$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(bFiltered_sg, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(bFiltered_sg$e_meta$Protein)))
  )

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

  # Test applyFilt anova ignoring singleton groups -----------------------------

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
  expect_equal(
    attr(aFiltered_sg_f, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(aFiltered_sg_f$e_data[, 1])),
         num_miss_obs = sum(is.na(aFiltered_sg_f$e_data)),
         prop_missing = (sum(is.na(aFiltered_sg_f$e_data)) /
                           prod(dim(aFiltered_sg_f$e_data[, -1]))),
         num_samps = ncol(aFiltered_sg_f$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(aFiltered_sg_f, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(aFiltered_sg_f$e_meta$Protein)))
  )

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

  # Test applyFilt gtest ignoring singleton groups -----------------------------

  # Filter the reduced pepData object using gtest.
  # _f: remove_singleton_groups is set to FALSE.
  expect_message(gFiltered_sg_f <- applyFilt(filter_object = filter_sg,
                                             omicsData = pdata_sg,
                                             min_nonmiss_anova = NULL,
                                             min_nonmiss_gtest = 3,
                                             remove_singleton_groups = FALSE),
                 paste("You have specified remove_singleton_groups = FALSE, so",
                       "the sample\\(s\\) corresponding to the singleton",
                       "group\\(s\\) were not utilized in the IMD-ANOVA filter",
                       "and will be retained in the resulting omicsData",
                       "object.",
                       sep = " "))

  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, 'cnames'),
                   attr(gFiltered_sg_f, 'cnames'))
  expect_identical(attr(pdata_sg, 'check.names'),
                   attr(gFiltered_sg_f, 'check.names'))
  expect_identical(class(pdata_sg),
                   class(gFiltered_sg_f))

  # Investigate the filters attribute.
  expect_equal(attr(gFiltered_sg_f, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(gFiltered_sg_f, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = NA,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(gFiltered_sg_f, 'filters')[[1]]$filtered,
                   s_gtest_sg)
  expect_equal(attr(gFiltered_sg_f, 'filters')[[1]]$method,
               "gtest")

  # Examinate the data_info attribute.
  expect_equal(
    attr(gFiltered_sg_f, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(gFiltered_sg_f$e_data[, 1])),
         num_miss_obs = sum(is.na(gFiltered_sg_f$e_data)),
         prop_missing = (sum(is.na(gFiltered_sg_f$e_data)) /
                           prod(dim(gFiltered_sg_f$e_data[, -1]))),
         num_samps = ncol(gFiltered_sg_f$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(gFiltered_sg_f, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(gFiltered_sg_f$e_meta$Protein)))
  )

  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(gFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)
  expect_identical(attr(gFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes"))
  expect_identical(attr(gFiltered_sg_f, "imdanova")$test_with_anova,
                   NULL)
  expect_identical(attr(gFiltered_sg_f, "imdanova")$test_with_gtest,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_gtest_sg))

  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(gFiltered_sg_f$e_data),
               c(128, 11))
  expect_equal(dim(gFiltered_sg_f$f_data),
               c(10, 2))
  expect_equal(dim(gFiltered_sg_f$e_meta),
               c(128, 4))

  # Test applyFilt combined ignoring singleton groups --------------------------

  # Filter the reduced pepData object using gtest.
  # _f: remove_singleton_groups is set to FALSE.
  expect_message(bFiltered_sg_f <- applyFilt(filter_object = filter_sg,
                                             omicsData = pdata_sg,
                                             min_nonmiss_anova = 3,
                                             min_nonmiss_gtest = 3,
                                             remove_singleton_groups = FALSE),
                 paste("You have specified remove_singleton_groups = FALSE, so",
                       "the sample\\(s\\) corresponding to the singleton",
                       "group\\(s\\) were not utilized in the IMD-ANOVA filter",
                       "and will be retained in the resulting omicsData",
                       "object.",
                       sep = " "))

  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata_sg, 'cnames'),
                   attr(bFiltered_sg_f, 'cnames'))
  expect_identical(attr(pdata_sg, 'check.names'),
                   attr(bFiltered_sg_f, 'check.names'))
  expect_identical(class(pdata_sg),
                   class(bFiltered_sg_f))

  # Investigate the filters attribute.
  expect_equal(attr(bFiltered_sg_f, 'filters')[[1]]$type,
               'imdanovaFilt')
  expect_identical(attr(bFiltered_sg_f, 'filters')[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 3,
                              min_nonmiss_gtest = 3))
  expect_identical(attr(bFiltered_sg_f, 'filters')[[1]]$filtered,
                   s_both_sg)
  expect_equal(attr(bFiltered_sg_f, 'filters')[[1]]$method,
               "combined")

  # Examinate the data_info attribute.
  expect_equal(
    attr(bFiltered_sg_f, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(bFiltered_sg_f$e_data[, 1])),
         num_miss_obs = sum(is.na(bFiltered_sg_f$e_data)),
         prop_missing = (sum(is.na(bFiltered_sg_f$e_data)) /
                           prod(dim(bFiltered_sg_f$e_data[, -1]))),
         num_samps = ncol(bFiltered_sg_f$e_data[, -1]),
         data_types = NULL)
  )

  # Explorate the meta_info attribute.
  expect_equal(
    attr(bFiltered_sg_f, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(bFiltered_sg_f$e_meta$Protein)))
  )

  # Checkerate the imdanova attribute. This attribute is added to the omicsData
  # object only when the imdanova filter functions are run.
  expect_identical(attr(bFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$nonmiss_totals,
                   filter_sg)
  expect_identical(attr(bFiltered_sg_f,
                        "imdanova")$nonmiss_per_group$group_sizes,
                   attr(filter_sg, "group_sizes"))
  expect_identical(attr(bFiltered_sg_f, "imdanova")$test_with_anova,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_anova_sg))
  expect_identical(attr(bFiltered_sg_f, "imdanova")$test_with_gtest,
                   setdiff(x = pdata_sg$e_data[, 1],
                           y = s_gtest_sg))

  # Inspecticate the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(bFiltered_sg_f$e_data),
               c(128, 11))
  expect_equal(dim(bFiltered_sg_f$f_data),
               c(10, 2))
  expect_equal(dim(bFiltered_sg_f$e_meta),
               c(128, 4))

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

  # Test scenario when nothing is filtered -------------------------------------

  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))

  # Produce a pretty little nmrData object.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'nmrClass')
  nmrdata <- group_designation(omicsData = nmrdata,
                               main_effects = "Condition")

  nmrFilta <- imdanova_filter(omicsData = nmrdata)

  # Apply the filter with a value for min_nonmiss_anova that will not filter any
  # rows.
  expect_message(noFilta <- applyFilt(filter_object = nmrFilta,
                                      omicsData = nmrdata,
                                      min_nonmiss_anova = 2,
                                      remove_singleton_groups = FALSE),
                 paste("No biomolecules were filtered with the values",
                       "specified for the min_nonmiss_anova and/or",
                       "min_nonmiss_gtest arguments.",
                       sep = " "))

  # The output of applyFilt should be the same as the omicsData object used as
  # the input because the filter was not applied. Therefore, the filters
  # attribute should remain how it was before running applyFilt.
  expect_identical(noFilta, nmrdata)

  # Create paired objects ------------------------------------------------------

  load(system.file('testdata',
                   'little_pairdata.RData',
                   package = 'pmartR'))

  # Create a pepData object with the original main effect and pairing variable.
  pairdata <- as.pepData(e_data = edata,
                         f_data = fdata,
                         e_meta = emeta,
                         edata_cname = 'Mass_Tag_ID',
                         fdata_cname = 'Name',
                         emeta_cname = 'Protein')
  pairdata <- edata_transform(pairdata,
                              data_scale = "log")
  pairdata <- group_designation(pairdata,
                                main_effects = "Virus",
                                pairs = "PairID")

  pair_filter <- imdanova_filter(omicsData = pairdata)

  anova_2 <- applyFilt(filter_object = pair_filter,
                       omicsData = pairdata,
                       min_nonmiss_anova = 2,
                       remove_singleton_groups = FALSE)

  gtest_2 <- applyFilt(filter_object = pair_filter,
                       omicsData = pairdata,
                       min_nonmiss_gtest = 2,
                       remove_singleton_groups = FALSE)

  combined_2 <- applyFilt(filter_object = pair_filter,
                          omicsData = pairdata,
                          min_nonmiss_anova = 2,
                          min_nonmiss_gtest = 2,
                          remove_singleton_groups = FALSE)

  differential <- data.frame(
    Mass_Tag_ID = pairdata$e_data$Mass_Tag_ID,
    Mock_1 = pairdata$e_data$Mock_0hr_1 - pairdata$e_data$Mock_18hr_1,
    Mock_2 = pairdata$e_data$Mock_0hr_2 - pairdata$e_data$Mock_18hr_2,
    Mock_3 = pairdata$e_data$Mock_0hr_3 - pairdata$e_data$Mock_18hr_3,
    Mock_4 = pairdata$e_data$Mock_0hr_4 - pairdata$e_data$Mock_18hr_4,
    Mock_5 = pairdata$e_data$Mock_0hr_5 - pairdata$e_data$Mock_18hr_5,
    FM_1 = pairdata$e_data$FM_0hr_1 - pairdata$e_data$FM_18hr_1,
    FM_2 = pairdata$e_data$FM_0hr_2 - pairdata$e_data$FM_18hr_2,
    FM_3 = pairdata$e_data$FM_0hr_3 - pairdata$e_data$FM_18hr_3,
    FM_4 = pairdata$e_data$FM_0hr_4 - pairdata$e_data$FM_18hr_4,
    FM_5 = pairdata$e_data$FM_0hr_5 - pairdata$e_data$FM_18hr_5,
    AM_1 = pairdata$e_data$AM_0hr_1 - pairdata$e_data$AM_18hr_1,
    AM_2 = pairdata$e_data$AM_0hr_2 - pairdata$e_data$AM_18hr_2,
    AM_3 = pairdata$e_data$AM_0hr_3 - pairdata$e_data$AM_18hr_3,
    AM_4 = pairdata$e_data$AM_0hr_4 - pairdata$e_data$AM_18hr_4,
    AM_5 = pairdata$e_data$AM_0hr_5 - pairdata$e_data$AM_18hr_5,
    row.names = NULL
  )

  count_diff <- data.frame(
    Mass_Tag_ID = differential$Mass_Tag_ID,
    count_mock = rowSums(!is.na(differential[, 2:6])),
    count_fm = rowSums(!is.na(differential[, 7:11])),
    count_am = rowSums(!is.na(differential[, 12:16]))
  )

  # Holy paired IMD-ANOVA filter tests, Batman ---------------------------------

  # Filter object tests ---------------

  expect_equal(pair_filter$AM,
               unname(rowSums(!is.na(pairdata$e_data[, 22:31]))))
  expect_equal(pair_filter$FM,
               unname(rowSums(!is.na(pairdata$e_data[, 12:21]))))
  expect_equal(pair_filter$Mock,
               unname(rowSums(!is.na(pairdata$e_data[, 2:11]))))
  expect_equal(
    attributes(pair_filter),
    list(
     names = c("Mass_Tag_ID", "AM", "FM", "Mock"),
     class = c("imdanovaFilt", "data.frame"),
     row.names = 1:150,
     group_sizes = data.frame(
       Group = c("AM", "FM", "Mock"),
       n_group = rep(10, 3)
     ),
     nonsingleton_groups = c("AM", "FM", "Mock")
    )
  )

  # Apply filter tests: ANOVA ---------------

  # filters attribute
  expect_equal(attr(anova_2, "filters")[[1]]$type,
               "imdanovaFilt")
  expect_identical(attr(anova_2, "filters")[[1]]$threshold,
                   data.frame(min_nonmiss_anova = 2,
                              min_nonmiss_gtest = NA))
  expect_identical(
    attr(anova_2, "filters")[[1]]$filtered,
    as.character(count_diff[which(rowSums(count_diff[, -1] >= 2) < 2), 1])
  )
  expect_equal(attr(anova_2, "filters")[[1]]$method,
               "anova")

  # data_info attribute
  expect_equal(
    attr(anova_2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(anova_2$e_data[, 1])),
         num_miss_obs = sum(is.na(anova_2$e_data)),
         prop_missing = (sum(is.na(anova_2$e_data)) /
                           prod(dim(anova_2$e_data[, -1]))),
         num_samps = ncol(anova_2$e_data[, -1]),
         data_types = NULL)
  )

  # meta_info attribute
  expect_equal(
    attr(anova_2, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(anova_2$e_meta$Protein)))
  )

  # imdanova attribute
  expect_equal(
    attr(anova_2, "imdanova")$nonmiss_per_group$nonmiss_totals$Mock,
    rowSums(!is.na(differential[, 2:6]))
  )
  expect_equal(
    attr(anova_2, "imdanova")$nonmiss_per_group$nonmiss_totals$FM,
    rowSums(!is.na(differential[, 7:11]))
  )
  expect_equal(
    attr(anova_2, "imdanova")$nonmiss_per_group$nonmiss_totals$AM,
    rowSums(!is.na(differential[, 12:16]))
  )
  expect_equal(attr(anova_2, "imdanova")$nonmiss_per_group$group_sizes,
               data.frame(Group = c("AM", "FM", "Mock"),
                          n_group = rep(5, 3)))
  expect_identical(attr(anova_2, "imdanova")$test_with_anova,
                   count_diff[which(rowSums(count_diff[, -1] >= 2) >= 2), 1])
  expect_null(attr(anova_2, "imdanova")$test_with_gtest)

  # Dimensions of e_data, f_data, and e_meta
  expect_equal(dim(anova_2$e_data),
               c(104, 31))
  expect_equal(dim(anova_2$f_data),
               c(30, 5))
  expect_equal(dim(anova_2$e_meta),
               c(119, 6))

  # Apply filter tests: G-test ---------------

  # filters attribute
  expect_equal(attr(gtest_2, "filters")[[1]]$type,
               "imdanovaFilt")
  expect_identical(attr(gtest_2, "filters")[[1]]$threshold,
                   data.frame(min_nonmiss_anova = NA,
                              min_nonmiss_gtest = 2))
  expect_identical(
    attr(gtest_2, "filters")[[1]]$filtered,
    as.character(count_diff[which(apply(count_diff[, -1], 1, max) < 2), 1])
  )
  expect_equal(attr(gtest_2, "filters")[[1]]$method,
               "gtest")

  # data_info attribute
  expect_equal(
    attr(gtest_2, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "log",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(gtest_2$e_data[, 1])),
         num_miss_obs = sum(is.na(gtest_2$e_data)),
         prop_missing = (sum(is.na(gtest_2$e_data)) /
                           prod(dim(gtest_2$e_data[, -1]))),
         num_samps = ncol(gtest_2$e_data[, -1]),
         data_types = NULL)
  )

  # meta_info attribute
  expect_equal(
    attr(gtest_2, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(gtest_2$e_meta$Protein)))
  )

  # imdanova attribute
  expect_equal(
    attr(gtest_2, "imdanova")$nonmiss_per_group$nonmiss_totals$Mock,
    rowSums(!is.na(differential[, 2:6]))
  )
  expect_equal(
    attr(gtest_2, "imdanova")$nonmiss_per_group$nonmiss_totals$FM,
    rowSums(!is.na(differential[, 7:11]))
  )
  expect_equal(
    attr(gtest_2, "imdanova")$nonmiss_per_group$nonmiss_totals$AM,
    rowSums(!is.na(differential[, 12:16]))
  )
  expect_equal(attr(gtest_2, "imdanova")$nonmiss_per_group$group_sizes,
               data.frame(Group = c("AM", "FM", "Mock"),
                          n_group = rep(5, 3)))
  expect_identical(attr(gtest_2, "imdanova")$test_with_gtest,
                   count_diff[which(apply(count_diff[, -1], 1, max) >= 2), 1])
  expect_null(attr(gtest_2, "imdanova")$test_with_anova)

  # Dimensions of e_data, f_data, and e_meta
  expect_equal(dim(gtest_2$e_data),
               c(124, 31))
  expect_equal(dim(gtest_2$f_data),
               c(30, 5))
  expect_equal(dim(gtest_2$e_meta),
               c(144, 6))

})
