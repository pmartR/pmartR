context('filter by coefficient of variation')

test_that('RNA_filter and applyFilt produce the correct output',{

  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  # Create a seqData object with the reduced data set.
  pdata <- as.seqData(e_data = edata,
                      f_data = fdata,
                      edata_cname = 'ID_REF',
                      fdata_cname = 'Samples'
  )
  
  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = c('Treatment', "Tissue"))
  
  # Copy the seqData object created above. This object will have three groups
  # (two Infection and one Mock). The Mock group will be reduced to a singleton
  # group later. sg: singleton group.
  pdata_sg <- pdata
  
  # Remove nine of the samples.
  rm_inds = c(11:19)
  rm_samps = pdata_sg$f_data[rm_inds, get_fdata_cname(pdata_sg)]
  rm_inds_edata = which(colnames(pdata_sg$e_data) %in% rm_samps)
  
  pdata_sg$e_data <- pdata_sg$e_data[, -rm_inds_edata]
  pdata_sg$f_data <- pdata_sg$f_data[-rm_inds, ]
  
  # Run the group_designation function on the singleton group seqData object.
  pdata_sg_gdf <- group_designation(omicsData = pdata_sg,
                                    main_effects = c('Treatment', "Tissue"))
  
  set.seed(38)
  
  # Scramble the order of the order of the samples.
  scrambled <- sample(1:(ncol(edata) - 1), (ncol(edata) - 1))
  
  # Generate a seqData object with the columns of just edata scrambled.
  eggs <- as.seqData(e_data = edata[, c(1, scrambled + 1)],
                     f_data = fdata,
                     edata_cname = 'ID_REF',
                     fdata_cname = 'Samples'
  )
  
  # Add groupies to the scrambled seqData object.
  eggs_gdf <- group_designation(omicsData = eggs,
                                main_effects = c('Treatment', "Tissue"))
  
  # Test RNA_filter without groups ----------------------------------------------
  
  # Try creating a RNAFilt object with an untoward input object.
  expect_error(RNA_filter(omicsData = fdata),
               paste("omicsData must be of class 'seqData'",
                     sep = ' '))
  
  # Calculate the RNA with R functions. The first ID column of edata is removed.
  use_r_lib <- as.numeric(apply(pdata$e_data[ -1], 2, sum, pooled = FALSE))
  use_r_nonzero <- as.numeric(apply(pdata$e_data[ -1] != 0, 2, sum, pooled = FALSE))
  
  # Run RNA_filter on the reduced data frame.
  filter <- RNA_filter(omicsData = pdata)
  
  # Review the class for the filter object.
  expect_s3_class(filter,
                  c('RNAFilt', 'data.frame'))
  
  # Check the dimensions of filter.
  expect_equal(dim(filter),
               c(40, 4))
  
  # Ensure the RNA values are correct.
  expect_equal(round(sort(filter$LibrarySize), 3),
               round(sort(use_r_lib), 3))
  
  expect_equal(round(sort(filter$NonZero), 3),
               round(sort(use_r_nonzero), 3))
  
  expect_equal(round(sort(filter$ProportionNonZero), 3),
               round(sort(use_r_nonzero/nrow(edata)), 3))
  
  # Test RNA_filter with groups -------------------------------------------------
  
  # Run RNA_filter on the reduced data frame with groups added.
  filter_gdf <- RNA_filter(omicsData = pdata_gdf)
  
  # Review the class for the filter_gdf object.
  expect_s3_class(filter_gdf,
                  c('RNAFilt', 'data.frame'))
  
  # Check the dimensions of filter_gdf.
  expect_equal(dim(filter_gdf),
               c(40, 4))
  
  # Ensure the RNA values are correct.
  expect_equal(round(sort(filter_gdf$LibrarySize), 3),
               round(sort(use_r_lib), 3))
  
  expect_equal(round(sort(filter_gdf$NonZero), 3),
               round(sort(use_r_nonzero), 3))
  
  expect_equal(round(sort(filter_gdf$ProportionNonZero), 3),
               round(sort(use_r_nonzero/nrow(edata)), 3))
  
  # Run the RNA filter on the data with scrambled columns.
  filter_eggs <- RNA_filter(omicsData = eggs_gdf)
  
  expect_equal(filter_gdf,  filter_eggs)
  
  # Test RNA_filter with singleton groups ---------------------------------------
  
  # Run the RNA filter on the data with scrambled columns.
  filter_sg_gdf <- RNA_filter(omicsData = pdata_sg_gdf)
  
  # Calculate the RNA with R functions. The first ID column of edata is removed.
  use_r_lib <- as.numeric(apply(pdata_sg_gdf$e_data[ -1], 2, sum, pooled = FALSE))
  use_r_nonzero <- as.numeric(apply(pdata_sg_gdf$e_data[ -1] != 0, 2, sum, pooled = FALSE))
  
  
  # Review the class for the filter_sg_gdf object.
  expect_s3_class(filter_sg_gdf,
                  c('RNAFilt', 'data.frame'))
  
  # Check the dimensions of filter_sg_gdf.
  expect_equal(dim(filter_sg_gdf),
               c(31, 4))
  
  # Ensure the RNA values are correct.
  expect_equal(round(sort(filter_sg_gdf$LibrarySize), 3),
               round(sort(use_r_lib), 3))
  
  expect_equal(round(sort(filter_sg_gdf$NonZero), 3),
               round(sort(use_r_nonzero), 3))
  
  expect_equal(round(sort(filter_sg_gdf$ProportionNonZero), 3),
               round(sort(use_r_nonzero/nrow(edata)), 3))
  
  # Test applyFilt without groups ----------------------------------------------
  
  # Apply the filter without groups to the reduced peptide data set.
  filtered <- applyFilt(filter_object = filter,
                        omicsData = pdata,
                        size_library = 10000,
                        min_nonzero = 885)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'RNAFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold$size_library,
                   10000)
  expect_identical(attr(filtered, 'filters')[[1]]$threshold$min_nonzero,
                   885)
  expect_equal(sort(attr(filtered, 'filters')[[1]]$filtered),
               sort(c("cervix_hCG_10", "cervix_PBS_8",  "uterus_hCG_1",  
                 "uterus_PBS_8",  "uterus_PBS_R3")))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered, "data_info"),
    list(data_scale_orig = "counts",
         data_scale = "counts",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered$e_data[, 1])),
         num_zero_obs = sum(filtered$e_data == 0),
         prop_zeros = (sum(filtered$e_data == 0) /
                           prod(dim(filtered$e_data[, -1]))),
         num_samps = ncol(filtered$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(1200, 36))
  expect_equal(dim(filtered$f_data),
               c(35, 4))
  
  # Test applyFilt with groups -------------------------------------------------
  
  # Apply the filter with groups to the reduced peptide data set.
  filtered_gdf <- applyFilt(filter_object = filter_gdf,
                            omicsData = pdata,
                            size_library = 10000,
                            min_nonzero = 885)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_gdf, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_gdf, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_gdf))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_gdf, 'filters')[[1]]$type,
               'RNAFilt')
  expect_identical(attr(filtered_gdf, 'filters')[[1]]$threshold$size_library,
                   10000)
  expect_identical(attr(filtered_gdf, 'filters')[[1]]$threshold$min_nonzero,
                   885)
  expect_equal(sort(attr(filtered_gdf, 'filters')[[1]]$filtered),
               sort(c("cervix_hCG_10", "cervix_PBS_8",  "uterus_hCG_1",  
                 "uterus_PBS_8",  "uterus_PBS_R3")))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_gdf, "data_info"),
    list(data_scale_orig = "counts",
         data_scale = "counts",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_gdf$e_data[, 1])),
         num_zero_obs = sum(filtered_gdf$e_data == 0),
         prop_zeros = (sum(filtered_gdf$e_data == 0) /
                         prod(dim(filtered_gdf$e_data[, -1]))),
         num_samps = ncol(filtered_gdf$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_gdf, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Inspect the filtered_gdf e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_gdf$e_data),
               c(1200, 36))
  expect_equal(dim(filtered_gdf$f_data),
               c(35, 4))
  
  # Test applyFilt with singleton groups ---------------------------------------
  
  # Apply the filter with groups to the reduced peptide data set.
  filtered_sg_gdf <- applyFilt(filter_object = filter_sg_gdf,
                               omicsData = pdata_sg_gdf,
                               size_library = 10000,
                               min_nonzero = 885)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_sg_gdf, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_sg_gdf, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_sg_gdf))
  
  # Examine the filters attribute.

  expect_equal(attr(filtered_sg_gdf, 'filters')[[1]]$type,
               'RNAFilt')
  expect_identical(attr(filtered_sg_gdf, 'filters')[[1]]$threshold$size_library,
                   10000)
  expect_identical(attr(filtered_sg_gdf, 'filters')[[1]]$threshold$min_nonzero,
                   885)
  expect_equal(sort(attr(filtered_sg_gdf, 'filters')[[1]]$filtered),
               sort(c("cervix_hCG_10", "cervix_PBS_8",  "uterus_hCG_1")))
  
  expect_true(is.na(attr(filtered_sg_gdf, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(
    attr(filtered_sg_gdf, "data_info"),
    list(data_scale_orig = "counts",
         data_scale = "counts",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(filtered_sg_gdf$e_data[, 1])),
         num_zero_obs = sum(filtered_sg_gdf$e_data == 0),
         prop_zeros = (sum(filtered_sg_gdf$e_data == 0) /
                         prod(dim(filtered_sg_gdf$e_data[, -1]))),
         num_samps = ncol(filtered_sg_gdf$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Explore the meta_info attribute.
  expect_equal(
    attr(filtered_sg_gdf, "meta_info"),
    list(meta_data = FALSE,
         num_emeta = NULL)
  )
  
  # Inspect the filtered_sg_gdf e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_sg_gdf$e_data),
               c(1200, 29))
  expect_equal(dim(filtered_sg_gdf$f_data),
               c(28, 4))
  
  # Test scenario when nothing is filtered -------------------------------------

  # Apply the filter with a value for RNA_threshold that will not filter any
  # rows.
  expect_error(noFilta <- applyFilt(filter_object = filter_sg_gdf,
                                      omicsData = pdata_sg_gdf,
                                      size_library = max(na.omit(filter_sg_gdf$LibrarySize)) + 1),
                                      "size_library must be integer of length 1 less than max")
  
  # Expect warning if data has already been filtered ---------------------------
  
  filter1 <- RNA_filter(omicsData = pdata)
  filter2 <- RNA_filter(omicsData = pdata)
  filter3 <- RNA_filter(omicsData = pdata)
  
  filtered1 <- applyFilt(filter_object = filter1,
                        omicsData = pdata,
                        size_library = 10000,
                        min_nonzero = 885)
  warnings <- capture_warnings(
    filtered2 <- applyFilt(filter_object = filter2,
                      omicsData = filtered1,
                      size_library = 10000,
                      min_nonzero = 885)
  )
  
  expect_match(
    warnings,
    "An RNA filter has already been applied to this data set.",
    all = FALSE
  )
  expect_match(
    warnings,
    paste0(
      "Specified samples cervix_PBS_8, cervix_hCG_10, uterus_PBS_8, uterus_PBS",
      "_R3, uterus_hCG_1 were not found in the data\\."
    ),
    all = FALSE
  )
  
  expect_identical(filtered1$e_data, filtered2$e_data)
})
