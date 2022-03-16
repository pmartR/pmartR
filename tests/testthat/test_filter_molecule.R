context('filter by molecule')

test_that('molecule_filter and applyFilt produce the correct output',{

  # Load data and prepare omicsData objects ------------------------------------

  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # add group and batch to fdata
  fdata$Condition <- c(rep("Infection",6),rep("Mock",6))
  fdata$Batch <- rep(seq(1:2),6)

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
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
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
               c(12, 3))
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
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
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
               c(12, 3))
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
  
  # Show errors if batch or group not specified and we use use_batch and use_groups
  expect_error(molecule_filter(omicsData = pdata,use_batch = TRUE),
               paste("omicsData must have batch_id specified if use_batch = TRUE"))
  expect_error(molecule_filter(omicsData = pdata,use_groups = TRUE),
                paste("omicsData must have groups specified if use_groups = TRUE"))
  
  pdata_gp <- group_designation(pdata,main_effects = "Condition",
                                batch_id = "Batch")

  # Test scenario when batch is used for molecule filter -----------------------
  batch_filter <- molecule_filter(pdata_gp,use_batch = TRUE)
  
  # Review the class for the filter object.
  expect_s3_class(batch_filter,
                  c('moleculeFilt', 'data.frame'))
  
  # Check the dimensions of filter.
  expect_equal(dim(batch_filter),
               c(150, 2))
  
  # ensure the number of observations are correct by calculating by hand as well
  s_count <- pdata_gp$e_data %>%
    tidyr::pivot_longer(cols = -Mass_Tag_ID, names_to = "SampleID",values_to = "value") %>%
    dplyr::left_join(pdata$f_data, by = "SampleID") %>%
    dplyr::group_by(Mass_Tag_ID, Batch) %>%
    dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
    dplyr::group_by(Mass_Tag_ID) %>%
    dplyr::summarise(min_num_obs = min(num_obs),.groups = "keep") %>%
    dplyr::ungroup()
  s_count <- dplyr::pull(s_count[,2])
  
  # Ensure the row sums are correct.
  expect_equal(batch_filter$Num_Observations,
               s_count)
  
  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(batch_filter, 'num_samps'),
                   12)
  
  # Apply the filter to the reduced peptide data set.
  filtered_b <- applyFilt(filter_object = batch_filter,
                        omicsData = pdata_gp,
                        min_num = 2)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_b, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_b, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_b))
  
  # Find the peptide IDs that will be filtered (without molecule_filter
  # functions).
  s_pepes <- pdata_gp$e_data[which(s_count < 2), 1]

  # Examine the filters attribute.
  expect_equal(attr(filtered_b, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_b, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered_b, 'filters')[[1]]$filtered,
               s_pepes)
  expect_true(is.na(attr(filtered_b, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_b, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_b$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered_b$e_data)),
         prop_missing = (sum(is.na(filtered_b$e_data)) /
                           prod(dim(filtered_b$e_data[, -1]))),
         num_samps = ncol(filtered_b$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_b, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_b$e_meta$Protein)))
  )
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_b$e_data),
               c(126, 13))
  expect_equal(dim(filtered_b$f_data),
               c(12, 3))
  expect_equal(dim(filtered_b$e_meta),
               c(126, 4))
  
  # Test scenario when group is used for molecule filter -----------------------
  group_filter <- molecule_filter(pdata_gp,use_groups = TRUE)
  
  # Review the class for the filter object.
  expect_s3_class(group_filter,
                  c('moleculeFilt', 'data.frame'))
  
  # Check the dimensions of filter.
  expect_equal(dim(group_filter),
               c(150, 2))
  
  # ensure the number of observations are correct by calculating by hand as well
  s_count <- pdata_gp$e_data %>%
    tidyr::pivot_longer(cols = -Mass_Tag_ID, names_to = "SampleID",values_to = "value") %>%
    dplyr::left_join(attr(pdata_gp,"group_DF"), by = "SampleID") %>%
    dplyr::group_by(Mass_Tag_ID, Group) %>%
    dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
    dplyr::group_by(Mass_Tag_ID) %>%
    dplyr::summarise(min_num_obs = min(num_obs),.groups = "keep") %>%
    dplyr::ungroup()
  s_count <- dplyr::pull(s_count[,2])
  
  # Ensure the row sums are correct.
  expect_equal(group_filter$Num_Observations,
               s_count)
  
  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(group_filter, 'num_samps'),
                   12)
  
  # Apply the filter to the reduced peptide data set.
  filtered_g <- applyFilt(filter_object = group_filter,
                          omicsData = pdata_gp,
                          min_num = 2)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_g, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_g, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_g))
  
  # Find the peptide IDs that will be filtered (without molecule_filter
  # functions).
  s_pepes <- pdata_gp$e_data[which(s_count < 2), 1]
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_g, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_g, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered_g, 'filters')[[1]]$filtered,
               s_pepes)
  expect_true(is.na(attr(filtered_g, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_g, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_g$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered_g$e_data)),
         prop_missing = (sum(is.na(filtered_g$e_data)) /
                           prod(dim(filtered_g$e_data[, -1]))),
         num_samps = ncol(filtered_g$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_g, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_g$e_meta$Protein)))
  )
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_g$e_data),
               c(127, 13))
  expect_equal(dim(filtered_g$f_data),
               c(12, 3))
  expect_equal(dim(filtered_g$e_meta),
               c(127, 4))
  
  # Test scenario when batch and group used for molecule filter ----------------
  
  bg_filter <- molecule_filter(pdata_gp,use_batch = TRUE, use_groups = TRUE)
  
  # Review the class for the filter object.
  expect_s3_class(bg_filter,
                  c('moleculeFilt', 'data.frame'))
  
  # Check the dimensions of filter.
  expect_equal(dim(bg_filter),
               c(150, 2))
  
  # ensure the number of observations are correct by calculating by hand as well
  s_count <- pdata_gp$e_data %>%
    tidyr::pivot_longer(cols = -Mass_Tag_ID, names_to = "SampleID",values_to = "value") %>%
    dplyr::left_join(attr(pdata_gp,"group_DF"), by = "SampleID") %>%
    dplyr::left_join(pdata_gp$f_data, by = "SampleID") %>%
    dplyr::group_by(Mass_Tag_ID, Group, Batch) %>%
    dplyr::summarise(num_obs = sum(!is.na(value)),.groups = "keep") %>%
    dplyr::group_by(Mass_Tag_ID) %>%
    dplyr::summarise(min_num_obs = min(num_obs),.groups = "keep") %>%
    dplyr::ungroup()
  s_count <- dplyr::pull(s_count[,2])
  
  # Ensure the row sums are correct.
  expect_equal(bg_filter$Num_Observations,
               s_count)
  
  # Inspect the number of samples in the moleculeFilt object.
  expect_identical(attr(bg_filter, 'num_samps'),
                   12)
  
  # Apply the filter to the reduced peptide data set.
  filtered_bg <- applyFilt(filter_object = bg_filter,
                          omicsData = pdata_gp,
                          min_num = 2)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_bg, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_bg, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_bg))
  
  # Find the peptide IDs that will be filtered (without molecule_filter
  # functions).
  s_pepes <- pdata_gp$e_data[which(s_count < 2), 1]
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_bg, 'filters')[[1]]$type,
               'moleculeFilt')
  expect_identical(attr(filtered_bg, 'filters')[[1]]$threshold,
                   2)
  expect_equal(attr(filtered_bg, 'filters')[[1]]$filtered,
               s_pepes)
  expect_true(is.na(attr(filtered_bg, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_bg, "data_info"),
    list(data_scale_orig = "abundance",
         data_scale = "abundance",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_bg$e_data[, 1])),
         num_miss_obs = sum(is.na(filtered_bg$e_data)),
         prop_missing = (sum(is.na(filtered_bg$e_data)) /
                           prod(dim(filtered_bg$e_data[, -1]))),
         num_samps = ncol(filtered_bg$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_bg, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(filtered_bg$e_meta$Protein)))
  )
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_bg$e_data),
               c(109, 13))
  expect_equal(dim(filtered_bg$f_data),
               c(12, 3))
  expect_equal(dim(filtered_bg$e_meta),
               c(109, 4))


})
