context('edata summary')

test_that('edata_summary correctly summarizes the data',{
  
  # Load the data and create a pepData object ----------------------------------
  
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Generate a second main effect.
  fdata2 <- data.frame(fdata,
                       Intensity = c('low', 'low', 'high', 'low', 'high',
                                     'high', 'high', 'high', 'low', 'none',
                                     'none', 'none'),
                       stringsAsFactors = FALSE)
  
  # Construct a pepData object.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata2,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')
  
  # Create test standards ------------------------------------------------------
  
  # Convert sample names into a character vector because in pdata they are
  # factors but in the output from edata_summary they are not.
  samp_char <- as.character(pdata$f_data[, 1])
  
  # Make the summary data frame standard for the summary by sample.
  stan_samp <- list(
    mean = data.frame(sample = samp_char,
                      mean = colMeans(pdata$e_data[, -1],
                                      na.rm = TRUE),
                      row.names = 1:12),
    sd = data.frame(sample = samp_char,
                    sd = apply(pdata$e_data[, -1], 2,
                               sd, na.rm = TRUE),
                    row.names = 1:12),
    median = data.frame(sample = samp_char,
                        median = apply(pdata$e_data[, -1], 2,
                                       median, na.rm = TRUE),
                        row.names = 1:12),
    pct_obs = data.frame(sample = samp_char,
                         pct_obs = apply(pdata$e_data[, -1], 2,
                                         function (x) {
                                           sum(!is.na(x)) / length(x)
                                         }),
                         row.names = 1:12),
    min = data.frame(sample = samp_char,
                     min = apply(pdata$e_data[, -1], 2,
                                 min, na.rm = TRUE),
                     row.names = 1:12),
    max = data.frame(sample = samp_char,
                     max = apply(pdata$e_data[, -1], 2,
                                 max, na.rm = TRUE),
                     row.names = 1:12)
  )
  
  # Add necessary attributes to the standard.
  attr(stan_samp, "by") <- "sample"
  attr(stan_samp, "cnames") <- list(edata_cname = "Mass_Tag_ID",
                                    fdata_cname = "SampleID")
  attr(stan_samp, "data_scale") <- "abundance"
  
  # Make the sample standard classy.
  class(stan_samp) <- "dataRes"
  
  # Make the summary data frame standard for the summary by molecule no main
  # effects.
  stan_mole <- list(
    mean = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                      mean = rowMeans(pdata$e_data[, -1],
                                      na.rm = TRUE)),
    sd = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                    sd = apply(pdata$e_data[, -1], 1,
                               sd, na.rm = TRUE)),
    median = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                        median = apply(pdata$e_data[, -1], 1,
                                       median, na.rm = TRUE)),
    pct_obs = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                         pct_obs = apply(pdata$e_data[, -1], 1,
                                         function (x) {
                                           sum(!is.na(x)) / length(x)
                                         })),
    min = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     min = apply(pdata$e_data[, -1], 1,
                                 min, na.rm = TRUE)),
    max = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     max = apply(pdata$e_data[, -1], 1,
                                 max, na.rm = TRUE))
  )
  
  # Add necessary attributes to the standard.
  attr(stan_mole, "by") <- "molecule"
  attr(stan_mole, "cnames") <- list(edata_cname = "Mass_Tag_ID",
                                    fdata_cname = "SampleID")
  attr(stan_mole, "data_scale") <- "abundance"
  
  # Make the molecule standard classy.
  class(stan_mole) <- "dataRes"
  
  # Make the summary data frame standard for the summary by molecule one main
  # effect.
  stan_mole_1 <- list(
    n_per_grp = data.frame(Group = factor(x = c("Infection", "Mock"),
                                          levels = c("Infection", "Mock")),
                           count = as.integer(c(9, 3))),
    mean = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                      Infection = rowMeans(pdata$e_data[, 2:10],
                                           na.rm = TRUE),
                      Mock = rowMeans(pdata$e_data[, 11:13],
                                      na.rm = TRUE),
                      row.names = 1:150),
    sd = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                    Infection = apply(pdata$e_data[, 2:10], 1,
                                      sd, na.rm = TRUE),
                    Mock = apply(pdata$e_data[, 11:13], 1,
                                 sd, na.rm = TRUE),
                    row.names = 1:150),
    median = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                        Infection = apply(pdata$e_data[, 2:10], 1,
                                          median, na.rm = TRUE),
                        Mock = apply(pdata$e_data[, 11:13], 1,
                                     median, na.rm = TRUE),
                        row.names = 1:150),
    pct_obs = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                         Infection = apply(pdata$e_data[, 2:10], 1,
                                           function (x) {
                                             sum(!is.na(x)) / length(x)
                                             }),
                         Mock = apply(pdata$e_data[, 11:13], 1,
                                      function (x) {
                                        sum(!is.na(x)) / length(x)
                                      }),
                         row.names = 1:150),
    min = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     Infection = apply(pdata$e_data[, 2:10], 1,
                                       function(x){if(all(is.na(x))){
                                         min(x)
                                       }else{
                                         min(x, na.rm = T)
                                       }}),
                     Mock = apply(pdata$e_data[, 11:13], 1,
                                  function(x){if(all(is.na(x))){
                                    min(x)
                                  }else{
                                    min(x, na.rm = T)
                                  }}),
                     row.names = 1:150),
    max = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     Infection = apply(pdata$e_data[, 2:10], 1,
                                       function(x){if(all(is.na(x))){
                                         max(x)
                                       }else{
                                         max(x, na.rm = T)
                                       }}),
                     Mock = apply(pdata$e_data[, 11:13], 1,
                                  function(x){if(all(is.na(x))){
                                    max(x)
                                  }else{
                                    max(x, na.rm = T)
                                  }}),
                     row.names = 1:150)
  )
  
  # Replace NaNs with NAs in the mean data frame so the standard can be compared
  # to the actual output from edata_summary.
  stan_mole_1[[2]][is.nan(stan_mole_1[[2]][, 2]), 2] <- NA
  stan_mole_1[[2]][is.nan(stan_mole_1[[2]][, 3]), 3] <- NA
  
  # Add necessary attributes to the standard.
  attr(stan_mole_1, "by") <- "molecule"
  attr(stan_mole_1, "groupvar") <- "Condition"
  attr(stan_mole_1, "cnames") <- list(edata_cname = "Mass_Tag_ID",
                                    fdata_cname = "SampleID")
  attr(stan_mole_1, "data_scale") <- "abundance"
  
  # Make the molecule standard classy.
  class(stan_mole_1) <- "dataRes"
  
  # Make the summary data frame standard for the summary by molecule two main
  # effects.
  stan_mole_2 <- list(
    n_per_grp = data.frame(Group = c("Infection_low",
                                     "Infection_high",
                                     "Mock_none"),
                           count = as.integer(c(4, 5, 3)),
                           row.names = c(2, 1, 3)),
    mean = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                      Infection_high = rowMeans(pdata$e_data[, c(4, 6:9)],
                                                na.rm = TRUE),
                      Infection_low = rowMeans(pdata$e_data[, c(2, 3, 5, 10)],
                                               na.rm = TRUE),
                      Mock_none = rowMeans(pdata$e_data[, 11:13],
                                      na.rm = TRUE),
                      row.names = 1:150),
    sd = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                    Infection_high = apply(pdata$e_data[, c(4, 6:9)], 1,
                                           sd, na.rm = TRUE),
                    Infection_low = apply(pdata$e_data[, c(2, 3, 5, 10)], 1,
                                          sd, na.rm = TRUE),
                    Mock_none = apply(pdata$e_data[, 11:13], 1,
                                 sd, na.rm = TRUE),
                    row.names = 1:150),
    median = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                        Infection_high = apply(pdata$e_data[, c(4, 6:9)], 1,
                                               median, na.rm = TRUE),
                        Infection_low = apply(pdata$e_data[, c(2, 3, 5, 10)], 1,
                                              median, na.rm = TRUE),
                        Mock_none = apply(pdata$e_data[, 11:13], 1,
                                     median, na.rm = TRUE),
                        row.names = 1:150),
    pct_obs = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                         Infection_high = apply(pdata$e_data[, c(4, 6:9)], 1,
                                                function (x) {
                                                  sum(!is.na(x)) / length(x)
                                                }),
                         Infection_low = apply(pdata$e_data[, c(2, 3, 5, 10)],
                                               1,
                                               function (x) {
                                                 sum(!is.na(x)) / length(x)
                                               }),
                         Mock_none = apply(pdata$e_data[, 11:13], 1,
                                           function (x) {
                                             sum(!is.na(x)) / length(x)
                                           }),
                         row.names = 1:150),
    min = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     Infection_high = apply(pdata$e_data[, c(4, 6:9)], 1,
                                            function(x){if(all(is.na(x))){
                                              min(x)
                                            }else{
                                              min(x, na.rm = T)
                                            }}),
                     Infection_low = apply(pdata$e_data[, c(2, 3, 5, 10)], 1,
                                           function(x){if(all(is.na(x))){
                                             min(x)
                                           }else{
                                             min(x, na.rm = T)
                                           }}),
                     Mock_none = apply(pdata$e_data[, 11:13], 1,
                                       function(x){if(all(is.na(x))){
                                         min(x)
                                       }else{
                                         min(x, na.rm = T)
                                       }}),
                     row.names = 1:150),
    max = data.frame(Mass_Tag_ID = pdata$e_data$Mass_Tag_ID,
                     Infection_high = apply(pdata$e_data[, c(4, 6:9)], 1,
                                            function(x){if(all(is.na(x))){
                                              max(x)
                                            }else{
                                              max(x, na.rm = T)
                                            }}),
                     Infection_low = apply(pdata$e_data[, c(2, 3, 5, 10)], 1,
                                           function(x){if(all(is.na(x))){
                                             max(x)
                                           }else{
                                             max(x, na.rm = T)
                                           }}),
                     Mock_none = apply(pdata$e_data[, 11:13], 1,
                                       function(x){if(all(is.na(x))){
                                         max(x)
                                       }else{
                                         max(x, na.rm = T)
                                       }}),
                     row.names = 1:150)
  )
  
  # Convert the row.names attribute of the group data frame to match the actual
  # output of the edata_summary function.
  class(attributes(stan_mole_2[[1]])$row.names) <- "integer"
  
  # Replace NaNs with NAs in the mean data frame so the standard can be compared
  # to the actual output from edata_summary.
  stan_mole_2[[2]][is.nan(stan_mole_2[[2]][, 2]), 2] <- NA
  stan_mole_2[[2]][is.nan(stan_mole_2[[2]][, 3]), 3] <- NA
  stan_mole_2[[2]][is.nan(stan_mole_2[[2]][, 4]), 4] <- NA
  
  # Add necessary attributes to the standard.
  attr(stan_mole_2, "by") <- "molecule"
  attr(stan_mole_2, "groupvar") <- c("Condition", "Intensity")
  attr(stan_mole_2, "cnames") <- list(edata_cname = "Mass_Tag_ID",
                                      fdata_cname = "SampleID")
  attr(stan_mole_2, "data_scale") <- "abundance"
  
  # Make the molecule standard classy.
  class(stan_mole_2) <- "dataRes"
  
  # Holy edata_summary tests, Batman! ------------------------------------------
  
  # Summarize by sample, no groupvar.
  samp <- edata_summary(pdata, by = "sample", groupvar = NULL)
  
  # Summarize by molecule, no groupvar.
  mole <- edata_summary(pdata, by = "molecule", groupvar = NULL)
  
  # Summarize by molecule, one groupvar.
  expect_warning(mole_1 <- edata_summary(pdata,
                                         by = "molecule",
                                         groupvar = "Condition"))
  
  # Summarize by molecule, two groupvars.
  expect_warning(mole_2 <- edata_summary(pdata,
                                         by = "molecule",
                                         groupvar = c("Condition",
                                                      "Intensity")))
  
  # Sleuth around each object.
  expect_identical(samp, stan_samp)
  expect_identical(mole, stan_mole)
  expect_identical(mole_1, stan_mole_1)
  expect_identical(mole_2, stan_mole_2)
  
})
