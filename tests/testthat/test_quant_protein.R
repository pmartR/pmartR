context('quantitation: protein')

test_that('each rollup method correctly quantifies proteins',{
  
  # Load data and prepare pepData objects --------------------------------------
  
  # Load peptide data.
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'pmartR'))
  
  # Construct a pepData object.
  pdata <- as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')
  
  # Log transmute the peptide data.
  pdata <- edata_transform(pdata, "log")
  
  # Group designate pdata.
  pdata2 <- group_designation(omicsData = pdata,
                              main_effects = "Condition")
  
  # Apply an IMD-ANOVA filter to pdata2 and save as pdata3.
  pdata3 <- applyFilt(filter_object = imdanova_filter(omicsData = pdata2),
                      omicsData = pdata2,
                      min_nonmiss_anova = 2)
  
  # Apply a filter to pdata2.
  pdata2 <- applyFilt(filter_object = molecule_filter(pdata2),
                      omicsData = pdata2)
  
  # Run some statisiticalness on the filtered data.
  inova <- imd_anova(omicsData = pdata3,
                     test_method = 'comb',
                     pval_adjust = 'bon')
  
  # Run bpquant pdata with the IMD-ANOVA output as the statRes object.
  bayes <- bpquant(statRes = inova, pepData = pdata3, parallel = FALSE)
  
  # Create objects that will be used throughout --------------------------------
  
  # Extricate the the group_DF attribute
  groupies <- attr(pdata2, "group_DF")
  
  # Merge e_data and e_meta (no filter) by "Mass_Tag_ID".
  merged <- merge(x = pdata$e_meta[, c("Mass_Tag_ID", "Protein")],
                  y = pdata$e_data,
                  by = "Mass_Tag_ID",
                  all.x = FALSE,
                  all.y = TRUE)
  
  # Merge e_data and e_meta (with filter) by "Mass_Tag_ID".
  merged2 <- merge(x = pdata2$e_meta[, c("Mass_Tag_ID", "Protein")],
                   y = pdata2$e_data,
                   by = "Mass_Tag_ID",
                   all.x = FALSE,
                   all.y = TRUE)
  
  # Extract unique proteins from e_meta (no filter).
  unique_proteins <- unique(as.character(pdata$e_meta$Protein))
  
  # Extract unique proteins from e_meta (with filter).
  unique_proteins2 <- unique(as.character(pdata2$e_meta$Protein))
  
  # Calculate standards - rollup -----------------------------------------------
  
  # No filters or group_DF ---------------
  
  # Standard for rollup - median.
  stan_med <- as.proData(
    e_data = merged %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(dplyr::across(.fns = median, na.rm = TRUE)) %>%
      data.frame(),
    f_data = pdata$f_data,
    e_meta = pdata$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::distinct(Protein, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_med, "data_info")$data_scale_orig <- "abundance"
  
  # Standard for rollup - mean.
  stan_mea <- as.proData(
    e_data = merged %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(dplyr::across(.fns = combine_fn_mean)) %>%
      data.frame(),
    f_data = pdata$f_data,
    e_meta = pdata$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::select(Protein, peps_per_pro, n_peps_used, Peptide_Sequence) %>%
      dplyr::distinct(.) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_mea, "data_info")$data_scale_orig <- "abundance"
  
  # With filters and group_DF ---------------
  
  # Standard for rollup - median.
  stan_med2 <- as.proData(
    e_data = merged2 %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(dplyr::across(.fns = median, na.rm = TRUE)) %>%
      data.frame(),
    f_data = pdata2$f_data,
    e_meta = pdata2$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::distinct(Protein, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_med2, "data_info")$data_scale_orig <- "abundance"
  
  # Update the group_DF attribute.
  attr(stan_med2, "group_DF") <- groupies
  
  # Quantitate without isoformRes - rollup -------------------------------------
  
  # No filters or group_DF ---------------
  
  # Quantitate using median.
  pq_med <- protein_quant(pepData = pdata,
                          method = 'rollup',
                          combine_fn = 'median',
                          isoformRes = NULL)
  
  # Quantitate using mean.
  pq_mea <- protein_quant(pepData = pdata,
                          method = 'rollup',
                          combine_fn = 'mean',
                          isoformRes = NULL,
                          emeta_cols = c("Peptide_Sequence"))
  
  # Ensure the output from protein_quant matches the standards.
  expect_identical(pq_med, stan_med)
  expect_identical(pq_mea, stan_mea)
  
  # With filters and group_DF ---------------
  
  # Quantitate using median.
  pq_med2 <- protein_quant(pepData = pdata2,
                           method = 'rollup',
                           combine_fn = 'median',
                           isoformRes = NULL)
  
  # Ensure the output from protein_quant matches the standards.
  expect_identical(pq_med2, stan_med2)
  
  # Calculate standards - rrollup ----------------------------------------------
  
  ratio <- function (cData, fn) {
    
    # Check the number of rows for cData. If there is only one row skip the
    # ratio calculation.
    if (nrow(cData) != 1) {
      
      # Filter the reference data frame by the current protein.
      ref_na <- rowSums(is.na(cData))
      ref_med <- apply(cData, 1, median, na.rm = TRUE)
      reference <- cData[ref_med == max(ref_med[ref_na == min(ref_na)]), ]
      
      # Calculate the ratio of each protein to its reference.
      ratiod <- rep(as.numeric(reference),
                    each = nrow(cData)) - cData
      
      # Compute the row-wise median. This will be used to scale each sample.
      mediand <- apply(ratiod, 1, median, na.rm = TRUE)
      
      # Scale the input data (cData).
      scaled <- cData + rep(mediand, ncol(cData))
      
      # Take the mean/median of the scaled data.
      protein <- scaled %>%
        dplyr::summarize(dplyr::across(.fns = fn))
      
    } else {
      
      protein <- cData
      
    }
    
    return (protein)
    
  }
  
  # Apply rrollup median.
  temp_med <-  merged %>%
    dplyr::select(-Mass_Tag_ID) %>%
    dplyr::nest_by(Protein) %>%
    dplyr::mutate(ratio = list(ratio(cData = data,
                                     fn = combine_fn_median)))
  
  # Produce an e_data object to create the rrollup standard.
  edata_med <- data.frame(temp_med$Protein,
                          data.table::rbindlist(temp_med$ratio))
  names(edata_med)[1] <- "Protein"
  edata_med <- edata_med[match(unique_proteins, edata_med$Protein), ]
  rownames(edata_med) <- NULL
  
  # Standard for rrollup - median.
  stan_rr_med <- as.proData(
    e_data = edata_med,
    f_data = pdata$f_data,
    e_meta = pdata$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::distinct(Protein, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_rr_med, "data_info")$data_scale_orig <- "abundance"
  
  # Quantitate without isoformRes - rrollup ------------------------------------
  
  # Quantitate using rrollup - median.
  rr_med <- protein_quant(pepData = pdata,
                          method = 'rrollup',
                          combine_fn = 'median',
                          use_parallel = FALSE,
                          isoformRes = NULL)
  
  # Compare the output to the standards.
  expect_identical(stan_rr_med, rr_med)
  
  # Calculate standards - qrollup ----------------------------------------------
  
  qtile <- function (cData, fn, qthold) {

    # Check the number of rows for cData. If there is only one row skip the
    # ratio calculation.
    if (nrow(cData) != 1) {
      
      # Compute rowwise means.
      r_mean <- rowMeans(cData, na.rm = TRUE)
      
      # Compute the threshold.
      threshold <- quantile(r_mean, probs = qthold, na.rm = TRUE)
      
      # Only keep the rows above the threshold.
      protein <- cData[r_mean >= threshold, ] %>%
        # Count the number of pepes (rows).
        dplyr::mutate(n_pepes = nrow(.)) %>%
        # Take the mean/median of the data.
        dplyr::summarize(dplyr::across(.fns = fn))
        
      
    } else {
      
      protein <- cData %>%
        # Add a count for the number of pepes (rows).
        dplyr::mutate(n_pepes = 1)
      
    }
    
    return (protein)
    
  }
  
  # Apply qrollup median.
  temp_med <-  merged %>%
    dplyr::select(-Mass_Tag_ID) %>%
    dplyr::nest_by(Protein) %>%
    dplyr::mutate(qtile = list(qtile(cData = data,
                                     fn = combine_fn_median,
                                     qthold = 0.4)))
  
  # Produce an e_data object to create the qrollup standard.
  edata_med <- data.frame(temp_med$Protein,
                          data.table::rbindlist(temp_med$qtile))
  names(edata_med)[1] <- "Protein"
  edata_med <- edata_med[match(unique_proteins, edata_med$Protein), ]
  rownames(edata_med) <- NULL
  
  # Match the protein counts with proteins in e_meta.
  qr_pepes <- dplyr::inner_join(x = pdata$e_meta,
                                y = edata_med[, c(1, 14)])
  
  # Standard for qrollup - median.
  stan_qr_med <- as.proData(
    e_data = edata_med[, -14],
    f_data = pdata$f_data,
    e_meta = pdata$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::distinct(Protein, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Update e_meta with the qrollup pepe counts.
  stan_qr_med$e_meta$n_peps_used <- dplyr::distinct(qr_pepes,
                                                    Protein, n_pepes)$n_pepes
  
  # Set the original data scale to abundance.
  attr(stan_qr_med, "data_info")$data_scale_orig <- "abundance"
  
  # Quantitate without isoformRes - qrollup ------------------------------------
  
  # Quantitate using qrollup - median.
  qr_med <- protein_quant(pepData = pdata,
                          method = 'qrollup',
                          combine_fn = 'median',
                          qrollup_thresh = 0.4,
                          use_parallel = FALSE,
                          isoformRes = NULL)
  
  # Compare the output to the standards.
  expect_identical(stan_qr_med, qr_med)
  
  # Calculate standards - zrollup ----------------------------------------------
  
  zombie <- function (cData, fn) {
    
    # Determine the number of rows for the current protein.
    num_peps <- nrow(cData)
    
    protein <- matrix(NA, nrow = 1, ncol =  ncol(cData))
    
    # Compute row-wise median and standard deviation.
    mds <- apply(cData, 1, median, na.rm = TRUE)
    sds <- apply(cData, 1, sd, na.rm = TRUE)
    
    # Create a vector with the row-wise medians and standard deviations repeated
    # the same number of times as the number of columns. For example, if the
    # median vector is 1, 2, 3 and there are two columns then the repeated
    # vector will be 1, 2, 3, 1, 2, 3.
    vector_mds <- rep(mds, ncol(cData))
    vector_sds <- rep(sds, ncol(cData))
    
    # Shift each peptide by its median and scale it by its standard deviation.
    scaled_pepes <- apply((cData - vector_mds) / vector_sds, 2, fn)
    
    protein[1, ] <- scaled_pepes
    protein <- data.frame(protein)
    names(protein) <- names(cData)
    
    return (protein)
    
  }
  
  # Run both filters on pdata (without group information). This data will be
  # used when calculating the standards for zrollup because a proteomics and
  # molecule filter are applied within the zrollup function when there are
  # singleton groups or proteins with a single peptide mapping to them.
  # Apply a filter to pdata2.
  pdataz <- applyFilt(filter_object = molecule_filter(pdata),
                      omicsData = pdata)
  pdataz <- applyFilt(proteomics_filter(pdataz),
                      pdataz,
                      min_num_peps = 2)
  
  # Merge e_data and e_meta (with both filters) by "Mass_Tag_ID".
  mergedz <- merge(x = pdataz$e_meta[, c("Mass_Tag_ID", "Protein")],
                   y = pdataz$e_data,
                   by = "Mass_Tag_ID",
                   all.x = FALSE,
                   all.y = TRUE)
  
  # Extract unique proteins from e_meta (with both filters).
  unique_proteinz <- unique(as.character(pdataz$e_meta$Protein))

  # Apply qrollup median.
  temp_med <-  mergedz %>%
    dplyr::select(-Mass_Tag_ID) %>%
    dplyr::nest_by(Protein) %>%
    dplyr::mutate(zombie = list(zombie(cData = data,
                                       fn = combine_fn_median)))
  
  # Produce an e_data object to create the qrollup standard.
  edata_med <- data.frame(temp_med$Protein,
                          data.table::rbindlist(temp_med$zombie))
  names(edata_med)[1] <- "Protein"
  edata_med <- edata_med[match(unique_proteinz, edata_med$Protein), ]
  rownames(edata_med) <- NULL
  
  # Standard for zrollup - median.
  stan_zr_med <- as.proData(
    e_data = edata_med,
    f_data = pdataz$f_data,
    e_meta = pdataz$e_meta %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein) %>%
      dplyr::mutate(peps_per_pro = dplyr::n()) %>%
      dplyr::mutate(n_peps_used = peps_per_pro) %>%
      dplyr::distinct(Protein, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein",
    fdata_cname = "SampleID",
    emeta_cname = "Protein",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_zr_med, "data_info")$data_scale_orig <- "abundance"
  
  # Quantitate without isoformRes - zrollup ------------------------------------
  
  # Quantitate using zrollup - median.
  zr_med <- protein_quant(pepData = pdata,
                          method = 'zrollup',
                          combine_fn = 'median',
                          single_observation = TRUE,
                          single_pep = TRUE,
                          use_parallel = FALSE,
                          isoformRes = NULL)
  
  # Compare the output to the standards.
  expect_identical(stan_zr_med, zr_med)
  
  # Calculate standards with isoformRes - rollup -------------------------------
  
  # Find peptides in both pdata3 and bayes.
  pepes <- which(
    pdata3$e_data[, "Mass_Tag_ID"] %in%
      attr(bayes, "isoformRes_subset")[, "Mass_Tag_ID"]
  )
  
  # Combine e_data and e_meta with output from the bayes object.
  merged_bayes <- merge(x = attr(bayes, "isoformRes_subset"),
                        y = pdata3$e_data[pepes, ],
                        by = "Mass_Tag_ID",
                        all.x = FALSE,
                        all.y = TRUE)
  
  # Count number of peptides per protein in the original e_meta object.
  pepes_per_pro <- emeta %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(peps_per_pro = dplyr::n()) %>%
    dplyr::select(Protein, peps_per_pro)
  
  # Count number of peptides per protein in the original e_meta object and keep
  # the column containing the peptide sequences.
  pepes_per_pro2 <- emeta %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(peps_per_pro = dplyr::n()) %>%
    dplyr::select(Protein, peps_per_pro, Peptide_Sequence)
  
  # Standard for rollup - median.
  stan_med_bayes <- as.proData(
    e_data = merged_bayes %>%
      dplyr::select(-Mass_Tag_ID, -Protein) %>%
      dplyr::group_by(Protein_Isoform) %>%
      dplyr::mutate(dplyr::across(.fns = median, na.rm = TRUE)) %>%
      data.frame(),
    f_data = pdata3$f_data,
    e_meta = attr(bayes, "isoformRes_subset") %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein_Isoform) %>%
      dplyr::mutate(n_peps_used = dplyr::n()) %>%
      dplyr::left_join(pepes_per_pro2, by = "Protein") %>%
      dplyr::relocate(n_peps_used, .after = peps_per_pro) %>%
      dplyr::distinct(Protein, Protein_Isoform, peps_per_pro,
                      n_peps_used, Peptide_Sequence) %>%
      data.frame(),
    edata_cname = "Protein_Isoform",
    fdata_cname = "SampleID",
    emeta_cname = "Protein_Isoform",
    data_scale = "log"
  )
  
  # Set the original data scale to abundance.
  attr(stan_med_bayes, "data_info")$data_scale_orig <- "abundance"
  
  # Quantitate with isoformRes - rollup ----------------------------------------
  
  # Rollup the pepes to proteins with isoformRes and pquant median.
  pq_med_bayes <- protein_quant(pepData = pdata,
                                method = 'rollup',
                                combine_fn = 'median',
                                use_parallel = FALSE,
                                isoformRes = bayes,
                                emeta_cols = c("Peptide_Sequence"))
  
  # Compare the output to the standards.
  expect_identical(stan_med_bayes, pq_med_bayes)
  
  # Calculate standards with isoformRes - qrollup ------------------------------
  
  # Apply qrollup median.
  edata_med <-  merged_bayes %>%
    dplyr::select(-Mass_Tag_ID, -Protein) %>%
    dplyr::nest_by(Protein_Isoform) %>%
    dplyr::mutate(qtile = list(qtile(cData = data,
                                     fn = combine_fn_median,
                                     qthold = 0.4))) %>%
    dplyr::select(Protein_Isoform, qtile) %>%
    tidyr::unnest(cols = c(Protein_Isoform, qtile)) %>%
    data.frame()
  rownames(edata_med) <- NULL
  
  # Match the protein counts with proteins in e_meta.
  qr_pepes <- dplyr::inner_join(x = attr(bayes, "isoformRes_subset"),
                                y = edata_med[, c(1, 14)])
  
  # Standard for qrollup - median.
  stan_qr_med_bayes <- as.proData(
    e_data = edata_med[, -14],
    f_data = pdata3$f_data,
    e_meta = attr(bayes, "isoformRes_subset") %>%
      dplyr::select(-Mass_Tag_ID) %>%
      dplyr::group_by(Protein_Isoform) %>%
      dplyr::mutate(n_peps_used = dplyr::n()) %>%
      dplyr::left_join(pepes_per_pro, by = "Protein") %>%
      dplyr::relocate(n_peps_used, .after = peps_per_pro) %>%
      dplyr::distinct(Protein, Protein_Isoform, peps_per_pro, n_peps_used) %>%
      data.frame(),
    edata_cname = "Protein_Isoform",
    fdata_cname = "SampleID",
    emeta_cname = "Protein_Isoform",
    data_scale = "log"
  )
  
  # Update e_meta with the qrollup pepe counts.
  stan_qr_med_bayes$e_meta$n_peps_used <- dplyr::distinct(
    qr_pepes, Protein, n_pepes
  )$n_pepes
  
  # Set the original data scale to abundance.
  attr(stan_qr_med_bayes, "data_info")$data_scale_orig <- "abundance"
  
  # Quantitate with isoformRes - qrollup ---------------------------------------
  
  # Quantitate using qrollup - median.
  qr_med_bayes <- protein_quant(pepData = pdata,
                                method = 'qrollup',
                                combine_fn = 'median',
                                qrollup_thresh = 0.4,
                                use_parallel = FALSE,
                                isoformRes = bayes)
  
  # Reorder the rows in the standard because the nest_by function creates a
  # different order than the qrollup function.
  stan_qr_med_bayes$e_data <- stan_qr_med_bayes$e_data[match(
    qr_med_bayes$e_data$Protein_Isoform,
    stan_qr_med_bayes$e_data$Protein_Isoform
  ), ]
  rownames(stan_qr_med_bayes$e_data) <- NULL
  
  # Compare the output to the standards.
  expect_identical(stan_qr_med_bayes, qr_med_bayes)
  
})
