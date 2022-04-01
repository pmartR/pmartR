context('normalize: nmr')

test_that('normalize_nmr produces the correct output',{
  
  # Load the nmr data frames ---------------------------------------------------
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Create an nmrData object.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = "Metabolite",
                        fdata_cname = "SampleID",
                        emeta_cname = "Bin")
  
  # Natural logate the data.
  lognmrdata <- edata_transform(omicsData = nmrdata,
                                data_scale = "log")
  
  # Take the natural log of the Concentration column in f_data. This will be
  # used for an nmrData object whose original scale is the log scale.
  logfdata <- fdata
  logfdata$Concentration <- log(fdata$Concentration)
  
  # Create an nmrData object whose original scale is the log scale.
  lnmrdata <- as.nmrData(e_data = lognmrdata$e_data,
                         f_data = logfdata,
                         e_meta = emeta,
                         edata_cname = "Metabolite",
                         fdata_cname = "SampleID",
                         emeta_cname = "Bin",
                         data_scale = "log")
  
  # Abundate the logged nmrData object.
  anmrdata <- edata_transform(omicsData = lnmrdata,
                              data_scale = "abundance")
  
  # Calculate the normalization standards --------------------------------------
  
  # Grab all of the row indices with the reference metabolite.
  ref_idx <- which(edata[, "Metabolite"] == "unkm1.53")
  
  # Extract the reference metabolite rows from edata.
  ref_spec1 <- as.numeric(edata[ref_idx, -1])
  
  # Natural logate the reference metabolite.
  ref_spec1_log <- log(ref_spec1)
  
  # Normalize the data in an nmr fashion to the reference metabolite.
  norm_spec1 <- (nmrdata$e_data[-ref_idx, -1] /
                   rep(ref_spec1, each = 37))
  
  # Back transmogrify spec1 abundance scale.
  norm_spec1_bt <- norm_spec1 * median(ref_spec1)
  
  # Normalize the log data in an nmr fashion to the reference metabolite.
  norm_spec1_log <- (lognmrdata$e_data[-ref_idx, -1] -
                       rep(ref_spec1_log, each = 37))
  
  # Back transmute spec1 log scale.
  norm_spec1_log_bt <- norm_spec1_log + median(ref_spec1_log)
  
  # Pluck out the normalizing values from f_data like an annoying chin whisker.
  ref_spec2 <- fdata$Concentration
  
  # Take the natural log of the concentration values.
  ref_spec2_log <- log(ref_spec2)
  
  # Normalize the data in an nmr fashion to the concentration values.
  norm_spec2 <- (nmrdata$e_data[, -1] /
                   rep(ref_spec2, each = 38))
  
  # Back transmogrify spec2 abundance scale.
  norm_spec2_bt <- norm_spec2 * median(ref_spec2)
  
  # Normalize the log data in an nmr fashion to the concentration values.
  norm_spec2_log <- (lognmrdata$e_data[, -1] -
                       rep(ref_spec2_log, each = 38))
  
  # Back transmute spec2 log scale.
  norm_spec2_log_bt <- norm_spec2_log + median(ref_spec2_log)
  
  # Test normalize_nmr: nmrnormRes ---------------------------------------------
  
  # Abundance scale ---------------
  
  # Use the first specification for identifying the reference samples.
  spec1 <- normalize_nmr(nmrdata,
                         apply_norm = FALSE,
                         metabolite_name = "unkm1.53")
  
  # Use the second specification for identifying the reference samples.
  spec2 <- normalize_nmr(nmrdata,
                         apply_norm = FALSE,
                         sample_property_cname = "Concentration")
  
  # Inspect the class of the two specifications.
  expect_s3_class(spec1, c("nmrnormRes", "list"))
  expect_s3_class(spec2, c("nmrnormRes", "list"))
  
  # Put on our sleuth hat and investigate the nmrnormRes attributes.
  expect_identical(attr(spec1, "cnames"),
                   list(edata_cname = "Metabolite",
                        fdata_cname = "SampleID"))
  expect_identical(attr(spec1, "nmr_info"),
                   list(metabolite_name = "unkm1.53",
                        sample_property_cname = NULL))
  expect_identical(attr(spec2, "cnames"),
                   list(edata_cname = "Metabolite",
                        fdata_cname = "SampleID"))
  expect_identical(attr(spec2, "nmr_info"),
                   list(metabolite_name = NULL,
                        sample_property_cname = "Concentration"))
  
  # Sleuth around the list elements.
  expect_identical(spec1$Sample,
                   fdata$SampleID)
  expect_identical(spec1$Metabolite,
                   "unkm1.53")
  expect_identical(spec1$value,
                   ref_spec1)
  expect_identical(spec2$Sample,
                   fdata$SampleID)
  expect_identical(spec2$Property,
                   "Concentration")
  expect_identical(spec2$value,
                   fdata$Concentration)
  
  # Log scale ---------------
  
  # Use the first specification for identifying the reference samples.
  spec1log <- normalize_nmr(lognmrdata,
                            apply_norm = FALSE,
                            metabolite_name = "unkm1.53")
  
  # Use the second specification for identifying the reference samples.
  spec2log <- normalize_nmr(lognmrdata,
                            apply_norm = FALSE,
                            sample_property_cname = "Concentration")
  
  # Examine the class of the two specifications.
  expect_s3_class(spec1log, c("nmrnormRes", "list"))
  expect_s3_class(spec2log, c("nmrnormRes", "list"))
  
  # Scrutinize the nmrnormRes objects after they have been logated. The objects
  # created with the second specification should be identical whether or not a
  # log transformation has bee applied.
  expect_identical(attributes(spec1log), attributes(spec1))
  expect_identical(spec1log$Sample, spec1$Sample)
  expect_identical(spec1log$Metabolite, spec1$Metabolite)
  expect_identical(spec1log$value, ref_spec1_log)
  expect_identical(spec2log, spec2)
  
  # Test normalize_nmr: apply_norm, no backtransmogrification ------------------
  
  # Abundance scale ---------------
  
  # Use the first specification for identifying the reference samples.
  expect_message(
    spec1 <- normalize_nmr(nmrdata,
                           apply_norm = TRUE,
                           metabolite_name = "unkm1.53"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Use the second specification for identifying the reference samples.
  expect_message(
    spec2 <- normalize_nmr(nmrdata,
                           apply_norm = TRUE,
                           sample_property_cname = "Concentration"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Compare the normalization to the standards.
  expect_identical(spec1$e_data[, -1],
                   norm_spec1)
  expect_identical(spec2$e_data[, -1],
                   norm_spec2)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec1, "nmrData")
  expect_s3_class(spec2, "nmrData")
  
  # Call in the big dogs (Holmes and Watson) to investigate the nmrData
  # attributes.
  expect_identical(attr(spec1, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec1, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 37,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec1, "nmr_info"),
    list(metabolite_name = "unkm1.53",
         sample_property_cname = NULL,
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using metabolite_name:",
                                              "unkm1.53",
                                              sep = " "),
                          norm_params = ref_spec1))
  )
  expect_equal(attr(spec1, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 37))
  expect_identical(attr(spec1, "filters"), list())
  expect_identical(attr(spec2, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec2, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2))
  )
  expect_equal(attr(spec2, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2, "filters"), list())
  
  # Log scale ---------------
  
  # Use the first specification for identifying the reference samples.
  expect_message(
    spec1_log <- normalize_nmr(lognmrdata,
                               apply_norm = TRUE,
                               metabolite_name = "unkm1.53"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Use the second specification for identifying the reference samples.
  expect_message(
    spec2_log <- normalize_nmr(lognmrdata,
                               apply_norm = TRUE,
                               sample_property_cname = "Concentration"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Compare the normalization to the standards.
  expect_identical(spec1_log$e_data[, -1],
                   norm_spec1_log)
  expect_equal(spec2_log$e_data[, -1],
               norm_spec2_log)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec1_log, "nmrData")
  expect_s3_class(spec2_log, "nmrData")
  
  # Call in the big dogs (Holmes and Watson) to investigate the nmrData
  # attributes.
  expect_identical(attr(spec1_log, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec1_log, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 37,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec1_log, "nmr_info"),
    list(metabolite_name = "unkm1.53",
         sample_property_cname = NULL,
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using metabolite_name:",
                                              "unkm1.53",
                                              sep = " "),
                          norm_params = ref_spec1_log))
  )
  expect_equal(attr(spec1_log, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 37))
  expect_identical(attr(spec1_log, "filters"), list())
  expect_identical(attr(spec2_log, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_log, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec2_log, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2_log))
  )
  expect_equal(attr(spec2_log, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_log, "filters"), list())
  
  # Test normalize_nmr: apply_norm, with backtransmogrification ----------------
  
  # Abundance scale ---------------
  
  # Use the first specification for identifying the reference samples.
  spec1_bt <- normalize_nmr(nmrdata,
                            apply_norm = TRUE,
                            backtransform = TRUE,
                            metabolite_name = "unkm1.53")
  
  # Use the second specification for identifying the reference samples
  spec2_bt <- normalize_nmr(nmrdata,
                            apply_norm = TRUE,
                            backtransform = TRUE,
                            sample_property_cname = "Concentration")
  
  # Compare the normalization to the standards.
  expect_identical(spec1_bt$e_data[, -1],
                   norm_spec1_bt)
  expect_equal(spec2_bt$e_data[, -1],
               norm_spec2_bt)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec1_bt, "nmrData")
  expect_s3_class(spec2_bt, "nmrData")
  
  # Call in the big dogs (Holmes and Watson) to investigate the nmrData
  # attributes.
  expect_identical(attr(spec1_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec1_bt, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 37,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec1_bt, "nmr_info"),
    list(metabolite_name = "unkm1.53",
         sample_property_cname = NULL,
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using metabolite_name:",
                                              "unkm1.53",
                                              sep = " "),
                          norm_params = ref_spec1))
  )
  expect_equal(attr(spec1_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 37))
  expect_identical(attr(spec1_bt, "filters"), list())
  expect_identical(attr(spec2_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_bt, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec2_bt, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2))
  )
  expect_equal(attr(spec2_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_bt, "filters"), list())
  
  # Log scale ---------------
  
  # Use the first specification for identifying the reference samples.
  spec1_log_bt <- normalize_nmr(lognmrdata,
                                apply_norm = TRUE,
                                backtransform = TRUE,
                                metabolite_name = "unkm1.53")
  
  # Use the second specification for identifying the reference samples
  spec2_log_bt <- normalize_nmr(lognmrdata,
                                apply_norm = TRUE,
                                backtransform = TRUE,
                                sample_property_cname = "Concentration")
  
  # Compare the normalization to the standards.
  expect_identical(spec1_log_bt$e_data[, -1],
                   norm_spec1_log_bt)
  expect_equal(spec2_log_bt$e_data[, -1],
               norm_spec2_log_bt)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec1_log_bt, "nmrData")
  expect_s3_class(spec2_log_bt, "nmrData")
  
  # Call in the big dogs (Holmes and Watson) to investigate the nmrData
  # attributes.
  expect_identical(attr(spec1_log_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec1_log_bt, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 37,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec1_log_bt, "nmr_info"),
    list(metabolite_name = "unkm1.53",
         sample_property_cname = NULL,
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using metabolite_name:",
                                              "unkm1.53",
                                              sep = " "),
                          norm_params = ref_spec1_log))
  )
  expect_equal(attr(spec1_log_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 37))
  expect_identical(attr(spec1_log_bt, "filters"), list())
  expect_identical(attr(spec2_log_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_log_bt, "data_info"),
               list(data_scale_orig = "abundance",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(
    attr(spec2_log_bt, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2_log))
  )
  expect_equal(attr(spec2_log_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_log_bt, "filters"), list())
  
  # Test normalize_nmr: when original data scale is log ------------------------
  
  # No backtransform ---------------
  
  # Use the second specification for identifying the reference samples
  expect_message(
    spec2_l <- normalize_nmr(lnmrdata,
                             apply_norm = TRUE,
                             backtransform = FALSE,
                             sample_property_cname = "Concentration"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Use the second specification for identifying the reference samples
  expect_message(
    spec2_a <- normalize_nmr(anmrdata,
                             apply_norm = TRUE,
                             backtransform = FALSE,
                             sample_property_cname = "Concentration"),
    paste("backtransform is set to FALSE. Examine the",
          "distribution of your data to ensure this is",
          "reasonable.",
          sep = " ")
  )
  
  # Compare the normalization to the standards.
  expect_identical(spec2_l$e_data[, -1],
                   norm_spec2_log)
  expect_equal(spec2_a$e_data[, -1],
               norm_spec2)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec2_l, "nmrData")
  expect_s3_class(spec2_a, "nmrData")
  
  # Sleuth around the nmr_info attribute.
  expect_equal(
    attr(spec2_l, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2_log))
  )
  expect_equal(
    attr(spec2_a, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = FALSE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2))
  )
  
  # Poke around the remaining attributes.
  expect_identical(attr(spec2_l, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_l, "data_info"),
               list(data_scale_orig = "log",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(attr(spec2_l, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_l, "filters"), list())
  
  expect_identical(attr(spec2_a, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_a, "data_info"),
               list(data_scale_orig = "log",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(attr(spec2_a, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_a, "filters"), list())
  
  # With backtransform ---------------
  
  # Use the second specification for identifying the reference samples
  spec2_l_bt <- normalize_nmr(lnmrdata,
                              apply_norm = TRUE,
                              backtransform = TRUE,
                              sample_property_cname = "Concentration")
  
  # Use the second specification for identifying the reference samples
  spec2_a_bt <- normalize_nmr(anmrdata,
                              apply_norm = TRUE,
                              backtransform = TRUE,
                              sample_property_cname = "Concentration")
  
  # Compare the normalization to the standards.
  expect_identical(spec2_l_bt$e_data[, -1],
                   norm_spec2_log_bt)
  expect_equal(spec2_a_bt$e_data[, -1],
               norm_spec2_bt)
  
  # Have a looksie at the class of the output for the two specifications.
  expect_s3_class(spec2_l_bt, "nmrData")
  expect_s3_class(spec2_a_bt, "nmrData")
  
  # Sleuth around the nmr_info attribute.
  expect_equal(
    attr(spec2_l_bt, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2_log))
  )
  expect_equal(
    attr(spec2_a_bt, "nmr_info"),
    list(metabolite_name = NULL,
         sample_property_cname = "Concentration",
         norm_info = list(is_normalized = TRUE,
                          backtransform = TRUE,
                          norm_method = paste("nmrObject was normalized",
                                              "using sample property:",
                                              "Concentration",
                                              sep = " "),
                          norm_params = ref_spec2))
  )
  
  # Poke around the remaining attributes.
  expect_identical(attr(spec2_l_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_l_bt, "data_info"),
               list(data_scale_orig = "log",
                    data_scale = "log",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(attr(spec2_l_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_l_bt, "filters"), list())
  
  expect_identical(attr(spec2_a_bt, "cnames"),
                   list(edata_cname = "Metabolite",
                        emeta_cname = "Bin",
                        fdata_cname = "SampleID",
                        techrep_cname = NULL))
  expect_equal(attr(spec2_a_bt, "data_info"),
               list(data_scale_orig = "log",
                    data_scale = "abundance",
                    norm_info = list(is_normalized = FALSE),
                    num_edata = 38,
                    num_miss_obs = 0,
                    prop_missing = 0,
                    num_samps = 41,
                    data_types = NULL,
                    batch_info = list(is_bc = FALSE)))
  expect_equal(attr(spec2_a_bt, "meta_info"),
               list(meta_data = TRUE,
                    num_emeta = 38))
  expect_identical(attr(spec2_a_bt, "filters"), list())
  
  # Test mutate_fdata ----------------------------------------------------------
  
  # Inspect output when the original data scale is log.
  expect_equal(pmartR:::mutate_fdata(ds = "log2",
                            ds_orig = "log",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log(fdata$Concentration)),
               log2(fdata$Concentration))
  expect_equal(pmartR:::mutate_fdata(ds = "log10",
                            ds_orig = "log",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log(fdata$Concentration)),
               log10(fdata$Concentration))
  
  # Inspect output when the original data scale is log2.
  expect_equal(pmartR:::mutate_fdata(ds = "log",
                            ds_orig = "log2",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log2(fdata$Concentration)),
               log(fdata$Concentration))
  expect_equal(pmartR:::mutate_fdata(ds = "log10",
                            ds_orig = "log2",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log2(fdata$Concentration)),
               log10(fdata$Concentration))
  
  # Inspect output when the original data scale is log10.
  expect_equal(pmartR:::mutate_fdata(ds = "log",
                            ds_orig = "log10",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log10(fdata$Concentration)),
               log(fdata$Concentration))
  expect_equal(pmartR:::mutate_fdata(ds = "log2",
                            ds_orig = "log10",
                            is_log = TRUE,
                            is_log_orig = TRUE,
                            sample_property = log10(fdata$Concentration)),
               log2(fdata$Concentration))
  
})
