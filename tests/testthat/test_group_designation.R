context('populate the group_DF attribute')

test_that('the correct group data frame and attributes are created',{

  # Load the reduced peptide data frames ---------------------------------------

  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))

  # Fabricate objects and tests for 1 main effect ------------------------------

  # Construct a pepData object with the edata, fdata, and emeta data frames.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')

  # Forge a group_DF attribute for pdata.
  pdata_gdf <- group_designation(omicsData = pdata,
                                 main_effects = 'Condition')

  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf$e_data),
               c(150, 13))
  expect_equal(dim(pdata_gdf$f_data),
               c(12, 2))
  expect_equal(dim(pdata_gdf$e_meta),
               c(150, 4))

  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf, 'group_DF')),
               data.frame(SampleID = as.character(fdata$SampleID),
                          Group = as.character(fdata$Condition),
                          stringsAsFactors = FALSE))

  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf, 'group_DF'))$main_effects,
               'Condition')
  expect_equal(attributes(attr(pdata_gdf, 'group_DF'))$nonsingleton_groups,
               c('Infection', 'Mock'))
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$pairs)
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$batch_id)
  expect_null(attributes(attr(pdata_gdf, 'group_DF'))$time_course)

  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata, 'cnames'),
                   attr(pdata_gdf, 'cnames'))
  expect_identical(attr(pdata, 'data_info'),
                   attr(pdata_gdf, 'data_info'))
  expect_identical(attr(pdata, 'meta_info'),
                   attr(pdata_gdf, 'meta_info'))
  expect_identical(attr(pdata, 'filters'),
                   attr(pdata_gdf, 'filters'))

  # Generate objects and tests for 2 main effects ------------------------------

  # Add another column to fdata for a second main effect.
  fdata_2 <- data.frame(fdata,
                        Intensity = c('low', 'low', 'high', 'low', 'high',
                                      'high', 'high', 'high', 'low', 'none',
                                      'none', 'none'),
                        stringsAsFactors = FALSE)

  # Produce a pepData object with the new f_data data frame.
  pdata_2 <- as.pepData(e_data = edata,
                        f_data = fdata_2,
                        e_meta = emeta,
                        edata_cname = 'Mass_Tag_ID',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'Protein')

  # Run group_designation with two main effects.
  pdata_gdf_2 <- group_designation(omicsData = pdata_2,
                                   main_effects = c('Condition', 'Intensity'))

  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf_2$e_data),
               c(150, 13))
  expect_equal(dim(pdata_gdf_2$f_data),
               c(12, 3))
  expect_equal(dim(pdata_gdf_2$e_meta),
               c(150, 4))

  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf_2, 'group_DF')),
               data.frame(SampleID = as.character(fdata_2$SampleID),
                          Group = as.character(paste(fdata_2$Condition,
                                                     fdata_2$Intensity,
                                                     sep = '_')),
                          Condition = as.character(fdata_2$Condition),
                          Intensity = as.character(fdata_2$Intensity),
                          stringsAsFactors = FALSE))

  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$main_effects,
               c('Condition', 'Intensity'))
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$nonsingleton_groups,
               c('Infection_high', 'Infection_low', 'Mock_none'))
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$pairs)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$batch_id)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$time_course)

  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata_2, 'cnames'),
                   attr(pdata_gdf_2, 'cnames'))
  expect_identical(attr(pdata_2, 'data_info'),
                   attr(pdata_gdf_2, 'data_info'))
  expect_identical(attr(pdata_2, 'meta_info'),
                   attr(pdata_gdf_2, 'meta_info'))
  expect_identical(attr(pdata_2, 'filters'),
                   attr(pdata_gdf_2, 'filters'))

  # Add another column to fdata for a second main effect with some NA values.
  fdata_2 <- data.frame(fdata,
                        Intensity = c('low', 'low', 'high', 'low', 'high',
                                      'high', 'high', 'high', 'low', NA, NA,
                                      NA),
                        stringsAsFactors = FALSE)

  # Produce a pepData object with the new f_data data frame.
  pdata_2 <- as.pepData(e_data = edata,
                        f_data = fdata_2,
                        e_meta = emeta,
                        edata_cname = 'Mass_Tag_ID',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'Protein')

  # Run group_designation with two main effects.
  expect_warning(pdata_gdf_2 <- group_designation(omicsData = pdata_2,
                                                  main_effects = c('Condition',
                                                                   'Intensity')),
                 paste('The following 3 sample\\(s\\) has/have been removed',
                       'from the dataset due to missing group information:',
                       'Mock1, Mock2, Mock3',
                       sep = ' '))

  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gdf_2$e_data),
               c(150, 10))
  expect_equal(dim(pdata_gdf_2$f_data),
               c(9, 3))
  expect_equal(dim(pdata_gdf_2$e_meta),
               c(150, 4))

  # Examinate the group_DF data frame.
  expect_equal(data.frame(attr(pdata_gdf_2, 'group_DF')),
               data.frame(SampleID = as.character(fdata_2$SampleID[1:9]),
                          Group = as.character(paste(fdata_2$Condition[1:9],
                                                     fdata_2$Intensity[1:9],
                                                     sep = '_')),
                          Condition = as.character(fdata_2$Condition[1:9]),
                          Intensity = as.character(fdata_2$Intensity[1:9]),
                          stringsAsFactors = FALSE))

  # Inspecticate the attributes of the group_DF data frame.
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$main_effects,
               c('Condition', 'Intensity'))
  expect_equal(attributes(attr(pdata_gdf_2, 'group_DF'))$nonsingleton_groups,
               c('Infection_high', 'Infection_low'))
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$covariates)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$pairs)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$batch_id)
  expect_null(attributes(attr(pdata_gdf_2, 'group_DF'))$time_course)

  # Checkerate the data_info attribute. The elements below should not be the
  # same as the original pepData object because the Mock samples were removed.
  expect_equal(get_data_info(pdata_gdf_2)$num_miss_obs,
               267)
  expect_equal(round(get_data_info(pdata_gdf_2)$prop_missing,
                     4),
               0.1978)
  expect_equal(get_data_info(pdata_gdf_2)$num_samps,
               9)

  # The following elements of data_info should be the same as the original
  # pepData object because they were not affected by removing the Mock samples.
  expect_identical(get_data_scale_orig(pdata),
                   get_data_scale_orig(pdata_gdf_2))
  expect_identical(get_data_scale(pdata),
                   get_data_scale(pdata_gdf_2))
  expect_identical(get_data_info(pdata)$norm_info$is_normalized,
                   get_data_info(pdata_gdf_2)$norm_info$is_normalized)
  expect_identical(get_data_info(pdata)$num_edata,
                   get_data_info(pdata_gdf_2)$num_edata)
  expect_identical(get_data_info(pdata)$data_types,
                   get_data_info(pdata_gdf_2)$data_types)

  # Create objects and tests for covariates ------------------------------------

  # Copy edata so the names of the samples can be changed.
  edata_3 <- edata

  # Change some of the Infection samples to Mock samples.
  names(edata_3) <- c("Mass_Tag_ID",
                       paste0("Infection", 1:6),
                       paste0("Mock", 1:6))

  # Create additional f_data objects with different main effects and covariates.
  fdata_3 <- fdata

  # Update the sample names in f_data.
  fdata_3$SampleID <- c(paste0("Infection", 1:6),
                         paste0("Mock", 1:6))

  # Update the first main effect to account for changing some infection samples
  # to mock samples.
  fdata_3$Condition <- c(rep("Infection", 6),
                          rep("Mock", 6))

  # Add a second main effect and two covariates.
  fdata_3$Level <- c("high", "low", "high", "low", "high", "low", "high",
                     "high", "low", "low", "low", "high")
  set.seed(720)
  fdata_3$Gender <- sample(c("F", "M"), 12, replace = TRUE)
  fdata_3$Age <- round(runif(12, min = 19, max = 89), 2)

  # Create a pepData object and run group_designation.
  pdata_3 <- as.pepData(e_data = edata_3,
                        f_data = fdata_3,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID")
  pdata_gdf_3 <- group_designation(omicsData = pdata_3,
                                   main_effects = c("Condition", "Level"),
                                   covariates = c("Gender", "Age"))

  # Investigate the covariates attribute.
  expect_identical(attr(attr(pdata_gdf_3, "group_DF"), "covariates"),
                   pdata_3$f_data[, c(1, 4, 5)])

  # Use a covariate that is a factor.
  fdata_3.2 <- fdata_3

  # Convert from a factor to a character.
  fdata_3.2$Gender <- factor(fdata_3$Gender,
                             levels = c("F", "M"))

  # Create a pepData object and run group_designation.
  pdata_3.2 <- as.pepData(e_data = edata_3,
                          f_data = fdata_3,
                          edata_cname = "Mass_Tag_ID",
                          fdata_cname = "SampleID")

  # Make a feeble attempt at changing the covariate to something other than a
  # character vector. Quickly give up and let fate overtake me as I realize it
  # is impossible.
  pdata_gdf_3.2 <- group_designation(omicsData = pdata_3.2,
                                     main_effects = c("Condition", "Level"),
                                     covariates = c("Gender", "Age"),
                                     cov_type = c("logical", "numeric"))

  # Exploricate the covariates attribute.
  expect_type(pdata_gdf_3.2$f_data$Gender, "character")
  expect_type(pdata_gdf_3.2$f_data$Age, "double")
  expect_type(attr(attr(pdata_gdf_3.2, "group_DF"), "covariates")$Gender,
              "character")
  expect_type(attr(attr(pdata_gdf_3.2, "group_DF"), "covariates")$Age,
              "double")
  expect_null(attributes(attr(pdata_gdf_3.2, 'group_DF'))$batch_id)
  
  # Correctly specify that the Gender covariate should be a character vector.
  pdata_gdf_3.4 <- group_designation(omicsData = pdata_3.2,
                                     main_effects = c("Condition", "Level"),
                                     covariates = c("Gender", "Age"),
                                     cov_type = c("character", "numeric"))

  # Interrogate the covariates attribute.
  expect_type(pdata_gdf_3.4$f_data$Gender, "character")
  expect_type(pdata_gdf_3.4$f_data$Age, "double")
  expect_type(attr(attr(pdata_gdf_3.4, "group_DF"), "covariates")$Gender,
              "character")
  expect_type(attr(attr(pdata_gdf_3.4, "group_DF"), "covariates")$Age,
              "double")
  expect_null(attributes(attr(pdata_gdf_3.4, 'group_DF'))$batch_id)
  

  # Assemble another f_data object with 0/1 representing m/f. This will be a
  # numeric vector that should be changed to a character vector after running
  # group_designation.
  fdata_4 <- fdata_3
  set.seed(85)
  fdata_4$Gender <- sample(0:1, 12, replace = TRUE)
  fdata_4$Age <- fdata_3$Age

  # Create a pepData object and run group_designation with the cov_type input.
  pdata_4 <- as.pepData(e_data = edata_3,
                        f_data = fdata_4,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID")
  pdata_gdf_4 <- group_designation(omicsData = pdata_4,
                                   main_effects = c("Condition", "Level"),
                                   covariates = c("Gender", "Age"),
                                   cov_type = c("numeric", "numeric"))

  # Examinate the covariates attribute.
  expect_identical(attr(attr(pdata_gdf_4, "group_DF"), "covariates"),
                   pdata_4$f_data[, c(1, 4, 5)])

  # Change the Gender covariate to a character vector.
  pdata_gdf_4.3 <- group_designation(omicsData = pdata_4,
                                     main_effects = c("Condition", "Level"),
                                     covariates = c("Gender", "Age"),
                                     cov_type = c("character", "numeric"))

  # Checkipate the type of the covariates before and after group_designation.
  expect_type(pdata_4$f_data$Gender, "integer")
  expect_type(pdata_4$f_data$Age, "double")
  expect_type(pdata_gdf_4.3$f_data$Gender, "character")
  expect_type(pdata_gdf_4.3$f_data$Age, "double")
  expect_type(attr(attr(pdata_gdf_4.3, "group_DF"), "covariates")$Gender,
              "character")
  expect_type(attr(attr(pdata_gdf_4.3, "group_DF"), "covariates")$Age,
              "double")
  expect_null(attributes(attr(pdata_gdf_4.3, 'group_DF'))$batch_id)
  
  # Try to change the covariates to something other than a character vector.
  pdata_gdf_4.5 <- group_designation(omicsData = pdata_4,
                                     main_effects = c("Condition", "Level"),
                                     covariates = c("Gender", "Age"),
                                     cov_type = c("logical", "raw"))

  # Diagnosticate the group_DF attribute. BAM!
  expect_type(pdata_gdf_4.5$f_data$Gender, "character")
  expect_type(pdata_gdf_4.5$f_data$Age, "character")
  expect_type(attr(attr(pdata_gdf_4.5, "group_DF"), "covariates")$Gender,
              "character")
  expect_type(attr(attr(pdata_gdf_4.5, "group_DF"), "covariates")$Age,
              "character")
  expect_null(attributes(attr(pdata_gdf_4.5, 'group_DF'))$batch_id)

  # Make sure pmart correctly brains the user with an error when the covariates
  # and cov_type vectors are not the same length.
  expect_error(group_designation(omicsData = pdata_4,
                                 main_effects = c("Condition", "Level"),
                                 covariates = c("Gender", "Age"),
                                 cov_type = c("character")))
  
  # Create objects and tests for batch_id -------------------------------------
  # Copy edata so the names of the samples can be changed.
  edata_5 <- edata
  
  # Change some of the Infection samples to Mock samples.
  names(edata_5) <- c("Mass_Tag_ID",
                      paste0("Infection", 1:6),
                      paste0("Mock", 1:6))
  
  # Create additional f_data objects with different main effects and covariates.
  fdata_5 <- fdata
  
  # Update the sample names in f_data.
  fdata_5$SampleID <- c(paste0("Infection", 1:6),
                        paste0("Mock", 1:6))
  
  # Update the first main effect to account for changing some infection samples
  # to mock samples.
  fdata_5$Condition <- c(rep("Infection", 6),
                         rep("Mock", 6))
  
  # Add a second main effect and two covariates.
  fdata_5$Level <- c("high", "low", "high", "low", "high", "low", "high",
                     "high", "low", "low", "low", "high")
  set.seed(720)
  fdata_5$Gender <- sample(c("F", "M"), 12, replace = TRUE)
  fdata_5$Age <- round(runif(12, min = 19, max = 89), 2)
  
  # Add a batch id
  fdata_5$BatchID <- rep(seq(1:2),6)
  fdata_5$BatchName <- rep(c("Batch1","Batch2"),6)
  
  # Create a pepData object and run group_designation.
  pdata_5 <- as.pepData(e_data = edata_5,
                        f_data = fdata_5,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID")
  
  # run a very simple model with just one main effect and batch id
  pdata_gdf_5 <- group_designation(omicsData = pdata_5,
                                   main_effects = "Condition",
                                   batch_id = "BatchID")
  # Investigate the batch_id attribute
  expect_identical(attr(attr(pdata_gdf_5, "group_DF"), "batch_id"),
                   pdata_5$f_data[,c(1,6)])
  
  # create a more complex model with two main effects and two covariates and batch_id
  pdata_gdf_5.1 <- group_designation(omicsData = pdata_5,
                                   main_effects = c("Condition", "Level"),
                                   covariates = c("Gender", "Age"),
                                   batch_id = "BatchID")

  # Investigate the covariates attribute.
  expect_identical(attr(attr(pdata_gdf_5.1, "group_DF"), "covariates"),
                   pdata_5$f_data[, c(1, 4, 5)])
  # Investigate the batch_id attribute
  expect_identical(attr(attr(pdata_gdf_5.1, "group_DF"), "batch_id"),
                   pdata_5$f_data[,c(1,6)])
  
  # Use a batch_id that is a factor rather than number
  pdata_gdf_5.2 <- group_designation(omicsData = pdata_5,
                                     main_effects = c("Condition", "Level"),
                                     covariates = c("Gender", "Age"),
                                     batch_id = "BatchName")
  
  # Exploricate the covariates attribute.
  # Investigate the batch_id attribute
  expect_identical(attr(attr(pdata_gdf_5.2, "group_DF"), "batch_id"),
                   pdata_5$f_data[,c(1,7)])
  # should be a character
  expect_type(attr(attr(pdata_gdf_5.2, "group_DF"),"batch_id")$BatchName,"character")

  # Create objects for paired data ---------------------------------------------

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

  # Create additional an additional main effect.
  pairdata_2_me <- pairdata
  pairdata_2_me$f_data$subclass <- c(
    rep(c("one", "one", "two", "two", "one"), 2),
    rep(c("two", "two", "two", "one", "one"), 2),
    rep(c("two", "two", "one", "one", "two"), 2)
  )

  # Generate a covariate.
  set.seed(12)
  pairdata_cov <- pairdata
  pairdata_cov$f_data$age <- c(
    rep(sample(runif(5)), 2),
    rep(sample(runif(5)), 2),
    rep(sample(runif(5)), 2)
  )

  # Create a scenario where one of the main effects does not match between the
  # two pairs.
  pairdata_bad_me <- pairdata_2_me
  pairdata_bad_me$f_data$subclass[[7]] <- "five"

  # Make a situation where the covariate is not the same for a pair.
  pairdata_bad_cov <- pairdata_cov
  pairdata_bad_cov$f_data$age[[29]] <- rcauchy(1)

  # Holy paired group designation tests, Batman --------------------------------

  # One main effect.
  expect_identical(
    attributes(
      attr(group_designation(pairdata,
                             main_effects = "Virus",
                             pairs = "PairID"),
           "group_DF")
    ),
    list(names = c("Name", "Group"),
         class = "data.frame",
         row.names = 1:30,
         main_effects = "Virus",
         pairs = "PairID",
         nonsingleton_groups = c("AM", "FM", "Mock"))
  )

  # Two main effects.
  expect_identical(
    attributes(
      attr(group_designation(pairdata_2_me,
                             main_effects = c("Virus", "subclass"),
                             pairs = "PairID"),
           "group_DF")
    ),
    list(names = c("Name", "Group", "Virus", "subclass"),
         class = "data.frame",
         row.names = 1:30,
         main_effects = c("Virus", "subclass"),
         pairs = "PairID",
         nonsingleton_groups = c("AM_one", "AM_two", "FM_one",
                                 "FM_two", "Mock_one", "Mock_two"))
  )

  # No main effects.
  expect_identical(
    attributes(
      attr(group_designation(pairdata,
                             pairs = "PairID"),
           "group_DF")
    ),
    list(names = c("Name", "Group"),
         class = "data.frame",
         row.names = 1:30,
         main_effects = "no_main_effect",
         pairs = "PairID",
         nonsingleton_groups = "zzzz")
  )

  expect_identical(
    attributes(
      attr(group_designation(pairdata_cov,
                             main_effects = "Virus",
                             covariates = "age",
                             pairs = "PairID"),
           "group_DF")
    ),
    list(names = c("Name", "Group"),
         class = "data.frame",
         row.names = 1:30,
         main_effects = "Virus",
         covariates = data.frame(
           Name = fdata$Name,
           age = pairdata_cov$f_data$age
         ),
         pairs = "PairID",
         nonsingleton_groups = c("AM", "FM", "Mock"))
  )

  expect_error(
      group_designation(pairdata_bad_me,
                             main_effects = c("Virus", "subclass"),
                             pairs = "PairID"),
    paste("The following samples have main effects that differ between",
          "pairs: Mock_0hr_2 and Mock_18hr_2",
          sep = " ")
  )

  expect_error(
    group_designation(pairdata_bad_cov,
                      main_effects = "Virus",
                      covariates = "age",
                      pairs = "PairID"),
    paste("The following samples have covariates that differ between",
          "pairs: AM_0hr_4 and AM_18hr_4",
          sep = " ")
  )

})
 