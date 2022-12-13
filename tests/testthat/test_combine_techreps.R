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
  
  # Test: seqData ------------------------------------------------------
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  # Construct a seqData object with the edata, fdata, and emeta data frames.
  seqdata <- as.seqData(e_data = edata,
                        f_data = fdata,
                        edata_cname = 'ID_REF',
                        fdata_cname = 'Samples'
  )
  
  f_data_temp <- fdata
  f_data_temp$techrep <- paste(sapply(c("A", "B", "C", "D"), rep, 10),  rep(1:5, 8))
  f_data_temp$bio_samp_names <- paste(sapply(c("crock", "o", "dile", "rock"), rep, 10),  rep(1:5, 8))
  f_data_temp$Label <- NULL
  
  temp_omics <- as.seqData(e_data = edata, edata_cname = 'ID_REF',
                           f_data = f_data_temp, fdata_cname = "Samples",
                           techrep_cname = "techrep", check.names = F)
  
  vector_biosamp <- f_data_temp$bio_samp_names
  
  ## check equal function
  check_equal <- function(og_obj, ct_obj, cmbn_fn, tr_col = NULL){
    
    tr_col <- ifelse(is.null(tr_col), attr(og_obj, "cname")$techrep_cname, tr_col)
    
    reps <- split(og_obj$f_data[[get_fdata_cname(og_obj)]], 
                  og_obj$f_data[[tr_col]])
    purrr::map(1:length(reps), function(n){
      cols <- reps[[n]]
      x <- expect_identical(
        as.numeric(apply(og_obj$e_data[cols], 1, cmbn_fn)),
        ct_obj$e_data[[names(reps)[n]]])
      return()
    })
    return(invisible(NULL))
  }
  
  ## Trigger class errors
  expect_error(combine_techreps(fdata), 
               "omicsData must be of class 'pepData', 'proData'")
  
  expect_error(combine_techreps(temp_omics, bio_sample_names = 5),
               "bio_sample_names must be a character string specifying a column in f_data")
  
  expect_error(combine_techreps(temp_omics, bio_sample_names = c("tiger",vector_biosamp[-1])),
  "character vector of sample names does not have the same number of names as the number of biological samples")
  
  temp_omics_grp_bad <- temp_omics %>% group_designation("Treatment")
  temp_omics_grp_bad$f_data$Treatment[[5]] <- "Tiger"
  
  expect_warning(combine_techreps(temp_omics_grp_bad), 
                                  "The following columns in f_data have multiple values across technical replicates")
  
  temp_comb1 <- combine_techreps(temp_omics)
  expect_identical(unique(temp_omics$f_data$techrep), temp_comb1$f_data$techrep)
  expect_identical(unique(temp_omics$f_data$techrep), colnames(temp_comb1$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb1, sum) ## Check equal sums
  
  temp_comb2 <- combine_techreps(temp_omics, combine_fn = "mean")
  expect_identical(unique(temp_omics$f_data$techrep), temp_comb2$f_data$techrep)
  expect_identical(unique(temp_omics$f_data$techrep), colnames(temp_comb2$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb2, mean)
  
  temp_comb3 <- combine_techreps(temp_omics, bio_sample_names = "bio_samp_names")
  expect_identical(unique(temp_omics$f_data$bio_samp_names), temp_comb3$f_data$bio_samp_names)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), colnames(temp_comb3$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb3, sum, tr_col = "bio_samp_names")
  
  temp_comb4 <- combine_techreps(temp_omics, bio_sample_names = vector_biosamp)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), temp_comb4$f_data$bio_samp_names)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), colnames(temp_comb4$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb4, sum, tr_col = "bio_samp_names")
  
  temp_omics <- group_designation(temp_omics, main_effects = c("Tissue", "Treatment"))
  
  temp_comb5 <- combine_techreps(temp_omics)
  expect_identical(unique(temp_omics$f_data$techrep), temp_comb5$f_data$techrep)
  expect_identical(unique(temp_omics$f_data$techrep), colnames(temp_comb5$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb5, sum) ## Check equal sums
  
  temp_comb6 <- combine_techreps(temp_omics, combine_fn = "mean")
  expect_identical(unique(temp_omics$f_data$techrep), temp_comb6$f_data$techrep)
  expect_identical(unique(temp_omics$f_data$techrep), colnames(temp_comb6$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb6, mean)
  
  temp_comb7 <- combine_techreps(temp_omics, bio_sample_names = "bio_samp_names")
  expect_identical(unique(temp_omics$f_data$bio_samp_names), temp_comb7$f_data$bio_samp_names)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), colnames(temp_comb7$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb7, sum, tr_col = "bio_samp_names")
  
  temp_comb8 <- combine_techreps(temp_omics, bio_sample_names = vector_biosamp)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), temp_comb8$f_data$bio_samp_names)
  expect_identical(unique(temp_omics$f_data$bio_samp_names), colnames(temp_comb8$e_data)[-1])
  check_equal(og_obj = temp_omics, ct_obj = temp_comb8, sum, tr_col = "bio_samp_names")
  
  
  ## No tech reps, differing #
  err_1 <- "This object did not have technical replicates specified."
  err_2 <- "Differing number of technical replicates per sample"
  
  ## Errors expected in techrep designation
  testthat::expect_error(combine_techreps(seqdata), err_1) ## no tech reps
  
  f_data_temp <- fdata
  f_data_temp$techrep <- paste(sapply(c("A", "B", "C", "D"), rep, 10),  rep(1:5, 8))
  f_data_temp$bio_samp_names <- paste(sapply(c("crock", "o", "dile", "rock"), rep, 10),  rep(1:5, 8))
  f_data_temp$Label <- NULL
  
  
  f_data_temp2 <- f_data_temp
  f_data_temp2$techrep[5] <- NA
  
  temp_omics <- as.seqData(e_data = edata, edata_cname = 'ID_REF',
                           f_data = f_data_temp2, fdata_cname = "Samples",
                           techrep_cname = "techrep")
  
  testthat::expect_error(combine_techreps(temp_omics), err_2) ## unequal tech reps
  
  
})
