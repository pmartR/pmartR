context('KW/AoV tests on normalization parameters')

test_that('independece tests output correct p-values based on group info',{
  
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
  
  # Log the heck out of the data.
  pdata <- edata_transform(pdata, "log")
  
  # Run the group_designation function on pdata.
  pdata_g <- group_designation(omicsData = pdata,
                               main_effects = "Condition")
  
  # Extract the group information from pdata.
  groupDF <- attr(pdata, "group_DF")
  
  # Normalize the data globally. Don't apply the normalization. 
  pdata_n <- normalize_global(pdata_g, subset_fn = "all",
                              norm_fn = "zscore",
                              apply_norm = FALSE,
                              backtransform = FALSE)
  
  # Globally normalize the data and apply the normalization
  pdata_an <- normalize_global(pdata_g, subset_fn = "all",
                               norm_fn = "zscore",
                               apply_norm = TRUE,
                               backtransform = FALSE)
  
  # Scramble the order of the samples.
  set.seed(4)
  scrambled <- sample(names(edata[, -1]))
  
  # Scramble the order of the columns in e_data and rows in f_data.
  edata_s <- edata[, c(1, match(scrambled, names(edata)[-1]) + 1)]
  
  # Create a pepData object with the scrambled data frames.
  pdata_s <- as.pepData(e_data = edata_s,
                        f_data = fdata,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID")
  
  # Log the heck out of the scrambled omicsData object.
  pdata_s <- edata_transform(pdata_s, "log")
  
  # Run the group_designation function on pdata_s.
  pdata_s <- group_designation(omicsData = pdata_s,
                               main_effects = "Condition")
  
  # Fish out the group information from the scrambled pepData object.
  groupDF_s <- attr(pdata_s, "group_DF")
  
  # Globally normalize the scrambled data and apply the normalization
  pdata_s <- normalize_global(pdata_s, subset_fn = "all",
                              norm_fn = "median",
                              apply_norm = TRUE,
                              backtransform = FALSE)
  
  # Calculate p-value standards ------------------------------------------------
  
  # Compute standard for omicsData object ---------------
  
  kw_loc_r <- kruskal.test(
    as.numeric(attributes(pdata_an)$data_info$norm_info$params$norm_location) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  kw_sca_r <- kruskal.test(
    as.numeric(attributes(pdata_an)$data_info$norm_info$params$norm_scale) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  standard_kw_1 <- list(p_location = kw_loc_r$p.value,
                        p_scale = kw_sca_r$p.value)
  
  aov_loc_r <- aov(
    as.numeric(attributes(pdata_an)$data_info$norm_info$params$norm_location) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  aov_sca_r <- aov(
    as.numeric(attributes(pdata_an)$data_info$norm_info$params$norm_scale) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  standard_aov_1 <- list(p_location = summary(aov_loc_r)[[1]]$`Pr(>F)`[[1]],
                         p_scale = summary(aov_sca_r)[[1]]$`Pr(>F)`[[1]])
  
  # Compute standard for normRes object ---------------
  
  kw_loc_r <- kruskal.test(
    as.numeric(pdata_n$parameters$normalization$location) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  kw_sca_r <- kruskal.test(
    as.numeric(pdata_n$parameters$normalization$scale) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  standard_kw_2 <- list(p_location = kw_loc_r$p.value,
                        p_scale = kw_sca_r$p.value)
  
  aov_loc_r <- aov(
    as.numeric(pdata_n$parameters$normalization$location) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  aov_sca_r <- aov(
    as.numeric(pdata_n$parameters$normalization$scale) ~
      c(rep("Infection", 9), rep("Mock", 3))
  )
  
  standard_aov_2 <- list(p_location = summary(aov_loc_r)[[1]]$`Pr(>F)`[[1]],
                         p_scale = summary(aov_sca_r)[[1]]$`Pr(>F)`[[1]])
  
  # Compute standard for scrambled omicsData object ---------------
  
  kw_loc_r <- kruskal.test(
    as.numeric(attributes(pdata_s)$data_info$norm_info$params$norm_location) ~
      c("Infection", "Mock", "Infection", "Mock", "Infection", "Infection",
        "Infection", "Infection", "Infection",  "Infection", "Mock",
        "Infection")
  )
  
  standard_kw_3 <- list(p_location = kw_loc_r$p.value,
                        p_scale = NULL)
  
  aov_loc_r <- aov(
    as.numeric(attributes(pdata_s)$data_info$norm_info$params$norm_location) ~
      c("Infection", "Mock", "Infection", "Mock", "Infection", "Infection",
        "Infection", "Infection", "Infection",  "Infection", "Mock",
        "Infection")
  )
  
  standard_aov_3 <- list(p_location = summary(aov_loc_r)[[1]]$`Pr(>F)`[[1]],
                         p_scale = NULL)
  
  # Test the heck out of the normRes tests -------------------------------------
  
  kw_1 <- normRes_tests(pdata_an, test_fn = "kw")
  expect_equal(kw_1, standard_kw_1)
  aov_1 <- normRes_tests(pdata_an, test_fn = "anova")
  expect_equal(aov_1, standard_aov_1)
  
  kw_2 <- normRes_tests(pdata_n, test_fn = "kw")
  expect_equal(kw_2, standard_kw_2)
  aov_2 <- normRes_tests(pdata_n, test_fn = "anova")
  expect_equal(aov_2, standard_aov_2)
  
  kw_3 <- normRes_tests(pdata_s, test_fn = "kw")
  expect_equal(kw_3, standard_kw_3)
  aov_3 <- normRes_tests(pdata_s, test_fn = "anova")
  expect_equal(aov_3, standard_aov_3)
  
  # Test argument checks ----------------------------------------------------
  
  # file_coverage("R/normRes_tests.R", "tests/testthat/test_normRes_tests.R")
  # 
  # expect_error(normRes_tests(pdata),
  #              "Normalization has not been run on this data")
  # expect_error(normRes_tests(normalize_global(pdata,
  #                                             subset_fn = "all",
  #                                             norm_fn = "median")),
  #              "No grouping structure present in object")
  
})
