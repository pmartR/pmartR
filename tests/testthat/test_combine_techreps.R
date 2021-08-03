context("combine technical replicates")

test_that("combine_techreps properly aggregates technical replicates", {
  
  # Load the data and create a pepData object ----------------------------------
  
  load(system.file('testdata',
                   'little_techdata.RData',
                   package = 'pmartR'))
  
  # Construct a pepData object.
  tdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'RunID',
                      techrep_cname = 'TECH_REP')
  
  # Create edata standard ------------------------------------------------------
  
  # Create a counter that will be used to fill the standard data frame.
  counter <- 1
  
  # Copy the first 33 rows of edata: one row for the ID column and the remaining
  # 32 rows for the averaged technical replicates.
  standard <- tdata$e_data[, 1:33]
  
  # Loop through each set of technical replicates and take their mean.
  for (e in seq(from = 2, to = 64, by = 2)) {
    
    # Update the counter in order to add the tech rep average to the correct
    # column of the standard.
    counter <- counter + 1
    
    # Average the technical replicates.
    standard[, counter] <- rowMeans(tdata$e_data[, e:(e + 1)],
                                    na.rm = TRUE)
    
  }
  
  # Convert NaNs to NAs to make the standard directly comparable to the output
  # from combine_techreps().
  for (e in 2:33) standard[[e]][is.nan(standard[[e]])] <- NA
  
  # Rename the columns of the standard to the unique names in the TECH_REP
  # column in f_data.
  colnames(standard)[2:33] <- unique(fdata$TECH_REP)
  
  # Create a list to hold the names of the technical replicates. This list will
  # become the standard for the tech_rep_info attribute.
  the_list <- vector(mode = "list",
                     length = 32)
  
  # Start a counter for filling in the_list.
  counter <- 0
  
  # Loop through each tech rep and extract the RunIDs.
  for (e in unique(tdata$f_data$TECH_REP)) {
    
    # Update the counter so the correct technical replicates are placed in the
    # proper order of the list.
    counter <- counter + 1
    
    the_list[[counter]] <- tdata$f_data$RunID[tdata$f_data$TECH_REP == e]
    
  }
  
  # Name the elements of the list with the technical replicate number.
  names(the_list) <- unique(fdata$TECH_REP)
  
  # tech_replicate the heck out of the data ------------------------------------
  
  # Use default settings for combine_techreps (the current data won't allow for
  # any other settings to be used).
  techie <- combine_techreps(tdata)
  
  # Put on my sleuth hat and start the investigation!
  expect_identical(techie$e_data, standard)
  expect_equal(dim(techie$f_data),
               c(32, 3))
  expect_null(attr(techie, "cnames")$techrep_cname)
  expect_identical(get_fdata_cname(techie), "TECH_REP")
  expect_identical(attr(techie, "tech_rep_info")$combine_method, "mean")
  expect_identical(attr(techie, "tech_rep_info")$tech_reps_by_sample,
                   the_list)
  expect_identical(attr(techie, "data_info")$num_edata,
                   nrow(techie$e_data))
  expect_identical(attr(techie, "data_info")$num_miss_obs,
                   sum(is.na(techie$e_data)))
  expect_identical(attr(techie, "data_info")$prop_missing,
                   (sum(is.na(techie$e_data)) /
                      prod(dim(techie$e_data[, -1]))))
  expect_identical(attr(techie, "data_info")$num_edata,
                   length(unique(techie$e_data[, 1])))
  
  # Test for errors when bio_sample_names does not fulfill the 1-to-1 mapping
  # requirement with the tech_rep column.
  expect_error(combine_techreps(tdata, bio_sample_names = "asdf"),
               paste("Specified display name column was not found in",
                     "f_data or was the same as fdata_cname",
                     sep = " "))
  expect_error(combine_techreps(tdata, bio_sample_names = "DILUTION"),
               paste("Specified display name column did not have a",
                     "one-to-one correspondence with the techrep ID column",
                     sep = " "))
  expect_error(combine_techreps(tdata, bio_sample_names = "FACTOR"),
               paste("Specified display name column did not have a",
                     "one-to-one correspondence with the techrep ID column",
                     sep = " "))
  expect_error(combine_techreps(tdata, bio_sample_names = "RunID"),
               paste("Specified display name column was not found in",
                     "f_data or was the same as fdata_cname",
                     sep = " "))
  
})
