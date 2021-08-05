context('calculate pooled CV')

test_that('pooled_cv_rcpp correctly calculates the CV by group',{
  
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
    
    # Create a list that will hold the indices of the CV values that are equal
    # to zero.
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
  
  # Extract e_data and remove the ID column.
  e_data <- pdata_gdf$e_data[, -1]
  
  # Fish out the group names from the group_DF attribute.
  groups <- attr(pdata_gdf, "group_DF")$Group
  
  # Test pooled_cv_rcpp with ordered rows --------------------------------------
  
  # Calculate the CV with C++ functions.
  cv_cpp <- pooled_cv_rcpp(mtr = as.matrix(e_data),
                           group = groups) * 100
  
  # Calculate the CV with R functions.
  cv_r <- cv.calc.pooled(pdata_gdf,
                         pooled = TRUE,
                         rm_single = FALSE)
  names(cv_r) <- NULL
  
  # Check the number of NaN entries and where they occur.
  expect_equal(which(is.na(cv_cpp)),
               c(1, 25, 27, 35, 53, 54, 59, 61, 68, 73))
  
  # See if R and C++ agree.
  expect_identical(round(cv_cpp, 3),
                   round(cv_r, 3))
  
  # Test pooled_cv_rcpp with unordered rows ------------------------------------
  
  set.seed(38)
  eggs <- sample(1:12, 12)
  
  # Scramble the order of the e_data columns and f_data rows. _s: scrambled
  edata_s <- edata[, c(1, eggs + 1)]
  fdata_s <- fdata[eggs, ]
  
  # Create a pepData object with the reduced data set.
  pdata_s <- as.pepData(e_data = edata_s,
                        f_data = fdata_s,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID")
  
  # Forge a group_DF attribute for pdata.
  pdata_s_gdf <- group_designation(omicsData = pdata_s,
                                   main_effects = 'Condition')
  
  # Extract e_data and remove the ID column.
  e_data_s <- pdata_s_gdf$e_data[, -1]
  
  # Go fishin for the group names from the group_DF attribute.
  groups_s <- attr(pdata_s_gdf, "group_DF")$Group
  
  # Calculate the CV with C++ functions when the order of the columns is
  # scrambled.
  cv_cpp_s <- pooled_cv_rcpp(mtr = as.matrix(e_data_s),
                             group = groups_s) * 100
  
  # Calculate the CV with R functions.
  cv_r_s <- cv.calc.pooled(pdata_s_gdf,
                           pooled = TRUE,
                           rm_single = FALSE)
  names(cv_r_s) <- NULL
  
  # Check the number of NaN entries and where they occur.
  expect_equal(which(is.na(cv_cpp_s)),
               c(1, 25, 27, 35, 53, 54, 59, 61, 68, 73))
  
  # See if R and C++ are playing nicely.
  expect_identical(round(cv_cpp_s, 3), round(cv_r_s, 3))
  
})
