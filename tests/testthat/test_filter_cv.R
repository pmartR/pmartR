context('filter by coefficient of variation')

test_that('cv_filter and applyFilt produce the correct output',{
  
  # Functions for calculating the CV -------------------------------------------
  
  cv.calc <- function(data, pooled) {
    
    # calculate mean of observations #
    ybar <- mean(data, na.rm = TRUE)
    
    # calculate standard deviation of observations #
    ystd <- sd(data, na.rm = TRUE)
    
    # Check if the pooled CV will be calculated.
    if (pooled) {
      
      # Check if ybar is a number and ystd is NA. If it is then there is only
      # one value for the current group and the CV should be 0.
      if (is.numeric(ybar) && is.na(ystd)) {
        
        sampcv <- 0
        
      } else {
        
        # calculate the sample cv value #
        sampcv <- ystd / ybar
        
        # Multiply the CV value by 100. (It's for fun)
        sampcv <- sampcv * 100
        
      }
      
    } else {
      
      # calculate the sample cv value #
      sampcv <- ystd / ybar
      
      # Multiply the CV value by 100. (It's for fun)
      sampcv <- sampcv * 100
      
    }
    
    # return cv values #
    return(sampcv)
    
  }
  
  # Create a function to count the number of non-missing values to be used in
  # connection with apply.
  count <- function (x) {
    
    present <- sum(!is.na(x))
    
    return (present)
    
  }
  
  # This function only calculates the pooled CV for the peptide data in the
  # little_pdata.RData file. It assumes the peptide ID column is the first
  # column in e_data. The order of the remaining columns can be changed though.
  cv.calc.pooled <- function (x, pooled, rm_single) {
    
    # Go fishin for the group attribute.
    groupDF <- attr(x, "group_DF")
    
    # Check if singleton groups should be removed.
    if (rm_single) {
      
      # Fetch the names of the non-singleton groups.
      group_names <- attributes(groupDF)$nonsingleton_groups
      
      # Determine the number of non-singleton groups.
      n_groups <- length(group_names)
      
    } else {
      
      # Fish out the unique group names.
      group_names <- unique(groupDF$Group)
      
      # Determine the number of groups present.
      n_groups <- length(group_names)
      
    }
    
    # Create a list that will hold the row indices for each sample that belongs
    # to each group in f_data. Note: a one needs to be added when using these
    # indices to correctly subset e_data because the ID column is the first
    # column of the test data.
    indices <- vector(mode = "list",
                      length = n_groups)
    
    # Loop through each group to find the indices of the samples that belong to
    # each group.
    for (e in 1:n_groups) {
      
      indices[[e]] <- which(groupDF$Group == group_names[[e]])
      
    }
    
    # Create a list that will hold the counts of non-missing values for each
    # group. Each element of n_nonmissing will be a vector with a length equal
    # to the number of rows in e_data. The elements of the vector will be the
    # counts of non-missing values for the current group. The list will have a
    # length equal to the number of groups.
    n_nonmissing <- vector(mode = "list",
                           length = n_groups)
    
    # Loop through each group.
    for (e in 1:n_groups) {
      
      # Count the non-missing values for the eth group.
      n_nonmissing[[e]] <- apply(x$e_data[, indices[[e]] + 1], 1, count)
      
    }
    
    # Create a list that will hold the CV for each group. Each element of cv
    # will be a vector with a length equal to the number of rows in e_data. The
    # elements of the vector will be the CV for each biomolecule calculated by
    # group. The length of the list is equal to the number of groups.
    cv <- vector(mode = "list",
                 length = n_groups)
    
    # Loop through each group.
    for (e in 1:n_groups) {
      
      # Calculate the CV for the eth group.
      cv[[e]] <- apply(x$e_data[, indices[[e]] + 1], 1, cv.calc, pooled = TRUE)
      
    }
    
    # Create a list that will hold the indices of the CV values that are equal to
    # zero.
    cv_zero <- vector(mode = "list",
                      length = n_groups)
    
    # Loop through each group and find the indices corresponding to a group with
    # a CV of zero. These will be used to change the counts from 1 to 0. This is
    # done because the CV will only depend on the group(s) that have more than
    # one non-missing value.
    for (e in 1:n_groups) {
      
      # Extricate the indices of cv whose value is zero.
      cv_zero[[e]] <- which(cv[[e]] == 0)
      
      # Change any elements in n_nonmissing[[e]] to zero that correspond to the
      # elements in cv[[e]] that are zero. These elements correspond to groups
      # that only have one non-missing value and should be ignored when
      # calculating the CV.
      n_nonmissing[[e]][cv_zero[[e]]] <- 0
      
    }
    
    # Fashion a list to hold the numerator of the pooled CV for each group
    # within each biomolecule. For example, if there are 150 rows in e_data then
    # each element in numerator (as well as each element in n_nonmissing and cv)
    # will have a length of 150. These 150 elements correspond to the numerator
    # of the pooled CV for each biomolecule in e_data.
    numerator <- vector(mode = "list",
                        length = n_groups)
    
    # Loop through each group and calculate the numerator of the pooled CV.
    for (e in 1:n_groups) {
      
      # Calculate the numerator of the pooled CV for each group.
      numerator[[e]] <- n_nonmissing[[e]] * cv[[e]]
      
    }
    
    # Calculate the pooled cv. The Reduce function performs element-wise
    # calculations across the vectors in numerator and n_nonmissing. For
    # example, if numerator <- list(1:3, 4:6) then Reduce(`+`, numerator) will
    # produce the vector 5, 7, 9.
    cv_pooled <- Reduce(`+`, numerator) / Reduce(`+`, n_nonmissing)
    
    # Return the pooled CV!!!
    return (cv_pooled)
    
  }
  
  # Load the reduced peptide data frames ---------------------------------------
  
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
  
  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = 'Condition')
  
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
  pdata_sg_gdf <- group_designation(omicsData = pdata_sg,
                                    main_effects = "Condition")
  
  # Test cv_filter without groups ----------------------------------------------
  
  # Try creating a cvFilt object with an untoward input object.
  expect_error(cv_filter(omicsData = fdata),
               paste("omicsData must be of class 'pepData', 'proData',",
                     "'metabData', 'lipidData', or 'nmrData'",
                     sep = ' '))
  
  # Calculate the CV with R functions. The first column of edata is removed
  # because it contains the peptide IDs.
  use_r <- apply(pdata$e_data[ -1], 1, cv.calc, pooled = FALSE)
  
  # Remove the names from the use_r vector.
  names(use_r) <- NULL
  
  # Run cv_filter on the reduced data frame.
  filter <- cv_filter(omicsData = pdata)
  
  # Review the class for the filter object.
  expect_s3_class(filter,
                  c('cvFilt', 'data.frame'))
  
  # Check the dimensions of filter.
  expect_equal(dim(filter),
               c(150, 2))
  
  # Ensure the CV values are correct.
  expect_equal(round(filter$CV, 3),
               round(use_r, 3))
  
  # Inspect the attributes of the cvFilt object.
  expect_identical(round(attr(filter, 'max_x_val'), 3),
                   106.068)
  expect_equal(attr(filter, 'tot_nas'),
               9)
  expect_false(attr(filter, 'pooled'))
  
  # Log transmogrify the data.
  pdata_log <- edata_transform(omicsData = pdata,
                               data_scale = "log")
  
  # CV filter the log transmogrified data.
  filter_log <- cv_filter(omicsData = pdata_log)
  
  # Check that the cv_filter object is the same before and after log
  # transmogrification.
  expect_equal(filter, filter_log)
  
  # Test cv_filter with groups -------------------------------------------------
  
  # Calculate the pooled CV with R funcitons.
  use_r_pooled <- cv.calc.pooled(x = pdata_gdf,
                                 pooled = TRUE,
                                 rm_single = FALSE)
  
  # Remove the names from use_r_pooled.
  names(use_r_pooled) <- NULL
  
  # Run cv_filter on the reduced data frame with groups added.
  filter_gdf <- cv_filter(omicsData = pdata_gdf)
  
  # Review the class for the filter_gdf object.
  expect_s3_class(filter_gdf,
                  c('cvFilt', 'data.frame'))
  
  # Check the dimensions of filter_gdf.
  expect_equal(dim(filter_gdf),
               c(150, 2))
  
  # Ensure the CV values are correct.
  expect_equal(round(filter_gdf$CV, 3),
               round(use_r_pooled, 3))
  
  # Inspect the attributes of the cvFilt object.
  expect_equal(round(attr(filter_gdf, 'max_x_val'), 3),
               69.295)
  expect_equal(attr(filter_gdf, 'tot_nas'),
               10)
  expect_true(attr(filter_gdf, 'pooled'))
  
  # Test cv_filter with singleton groups ---------------------------------------
  
  # Calculate the pooled CV with R funcitons for data with singleton groups.
  use_r_sg <- cv.calc.pooled(x = pdata_sg_gdf,
                             pooled = TRUE,
                             rm_single = TRUE)
  
  # Remove the names from use_r_sg.
  names(use_r_sg) <- NULL
  
  # Run cv_filter on the data with singleton groups.
  expect_warning(filter_sg_gdf <- cv_filter(omicsData = pdata_sg_gdf,
                                            use_groups = TRUE),
                 paste("Grouping information is being utilized when",
                       "calculating the CV, and there are group\\(s\\)",
                       "consisting of a single sample. The singleton",
                       "group\\(s\\) will be ignored by this filter.",
                       sep = " "))
  
  # Review the class for the filter_sg_gdf object.
  expect_s3_class(filter_sg_gdf,
                  c('cvFilt', 'data.frame'))
  
  # Check the dimensions of filter_sg_gdf.
  expect_equal(dim(filter_sg_gdf),
               c(150, 2))
  
  # Ensure the CV values are correct.
  expect_equal(round(filter_sg_gdf$CV, 3),
               round(use_r_sg, 3))
  
  # Inspect the attributes of the cvFilt object.
  expect_equal(round(attr(filter_sg_gdf, 'max_x_val'), 3),
               71.332)
  expect_equal(attr(filter_sg_gdf, 'tot_nas'),
               20)
  expect_true(attr(filter_sg_gdf, 'pooled'))
  
  # Test applyFilt without groups ----------------------------------------------
  
  # Apply a filter with the threshold being the max CV value (no peptides should
  # be filtered).
  expect_identical(applyFilt(filter_object = filter,
                             omicsData = pdata,
                             cv_threshold = max(na.omit(filter$CV)))$e_data,
                   pdata$e_data)
  
  # Apply the filter without groups to the reduced peptide data set.
  filtered <- applyFilt(filter_object = filter,
                        omicsData = pdata,
                        cv_threshold = 90)

  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered))

  # Examine the filters attribute.
  expect_equal(attr(filtered, 'filters')[[1]]$type,
               'cvFilt')
  expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                   90)
  expect_equal(attr(filtered, 'filters')[[1]]$filtered,
               c(11055, 6701524, 6781846, 6809644, 6948899))
  expect_true(is.na(attr(filtered, 'filters')[[1]]$method))

  # Investigate the data_info attribute.
  expect_equal(attr(filtered, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered, 'data_info')$num_edata,
               145)
  expect_equal(attr(filtered, 'data_info')$num_miss_obs,
               318)
  expect_equal(round(attr(filtered, 'data_info')$prop_missing, 4),
               0.1828)
  expect_equal(attr(filtered, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered, 'data_info')$data_types)

  # Explore the meta_info attribute.
  expect_true(attr(filtered, 'meta_info')$meta_data)
  expect_equal(attr(filtered, 'meta_info')$num_emeta,
               81)

  # Inspect the filtered e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered$e_data),
               c(145, 13))
  expect_equal(dim(filtered$f_data),
               c(12, 2))
  expect_equal(dim(filtered$e_meta),
               c(145, 4))
  
  # Test applyFilt with groups -------------------------------------------------
  
  # Apply the filter with groups to the reduced peptide data set.
  filtered_gdf <- applyFilt(filter_object = filter_gdf,
                            omicsData = pdata_gdf,
                            cv_threshold = 60)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_gdf, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_gdf, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_gdf))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_gdf, 'filters')[[1]]$type,
               'cvFilt')
  expect_identical(attr(filtered_gdf, 'filters')[[1]]$threshold,
                   60)
  expect_equal(attr(filtered_gdf, 'filters')[[1]]$filtered,
               c(6850636, 6948820, 6948835, 6948904))
  expect_true(is.na(attr(filtered_gdf, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered_gdf, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered_gdf, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_gdf, 'data_info')$num_edata,
               146)
  expect_equal(attr(filtered_gdf, 'data_info')$num_miss_obs,
               317)
  expect_equal(round(attr(filtered_gdf, 'data_info')$prop_missing, 4),
               0.1809)
  expect_equal(attr(filtered_gdf, 'data_info')$num_samps,
               12)
  expect_null(attr(filtered_gdf, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered_gdf, 'meta_info')$meta_data)
  expect_equal(attr(filtered_gdf, 'meta_info')$num_emeta,
               80)
  
  # Inspect the filtered_gdf e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_gdf$e_data),
               c(146, 13))
  expect_equal(dim(filtered_gdf$f_data),
               c(12, 2))
  expect_equal(dim(filtered_gdf$e_meta),
               c(146, 4))
  
  # Test applyFilt with singleton groups ---------------------------------------
  
  # Apply the filter with groups to the reduced peptide data set.
  filtered_sg_gdf <- applyFilt(filter_object = filter_sg_gdf,
                               omicsData = pdata_sg_gdf,
                               cv_threshold = 60)
  
  # Ensure the class and attributes that shouldn't have changed didn't change.
  expect_identical(attr(pdata, 'cnames'),
                   attr(filtered_sg_gdf, 'cnames'))
  expect_identical(attr(pdata, 'check.names'),
                   attr(filtered_sg_gdf, 'check.names'))
  expect_identical(class(pdata),
                   class(filtered_sg_gdf))
  
  # Examine the filters attribute.
  expect_equal(attr(filtered_sg_gdf, 'filters')[[1]]$type,
               'cvFilt')
  expect_identical(attr(filtered_sg_gdf, 'filters')[[1]]$threshold,
                   60)
  expect_equal(attr(filtered_sg_gdf, 'filters')[[1]]$filtered,
               c(6948820, 6948835, 6948904))
  expect_true(is.na(attr(filtered_sg_gdf, 'filters')[[1]]$method))
  
  # Investigate the data_info attribute.
  expect_equal(attr(filtered_sg_gdf, 'data_info')$data_scale,
               'abundance')
  expect_false(attr(filtered_sg_gdf, 'data_info')$norm_info$is_normalized,
               FALSE)
  expect_equal(attr(filtered_sg_gdf, 'data_info')$num_edata,
               147)
  expect_equal(attr(filtered_sg_gdf, 'data_info')$num_miss_obs,
               280)
  expect_equal(round(attr(filtered_sg_gdf, 'data_info')$prop_missing, 4),
               0.1905)
  expect_equal(attr(filtered_sg_gdf, 'data_info')$num_samps,
               10)
  expect_null(attr(filtered_sg_gdf, 'data_info')$data_types)
  
  # Explore the meta_info attribute.
  expect_true(attr(filtered_sg_gdf, 'meta_info')$meta_data)
  expect_equal(attr(filtered_sg_gdf, 'meta_info')$num_emeta,
               81)
  
  # Inspect the filtered_sg_gdf e_data, f_data, and e_meta data frames.
  expect_equal(dim(filtered_sg_gdf$e_data),
               c(147, 11))
  expect_equal(dim(filtered_sg_gdf$f_data),
               c(10, 2))
  expect_equal(dim(filtered_sg_gdf$e_meta),
               c(147, 4))
  
})
