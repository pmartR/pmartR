context('filter by coefficient of variation')

test_that('total_count_filter and applyFilt produce the correct output',{
 
  # Load the reduced seqdata data frames ---------------------------------------
  
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
  pdata_sg$e_data <- pdata_sg$e_data[, -c(11:19)]
  pdata_sg$f_data <- pdata_sg$f_data[-c(11:19), ]
  
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
  
  obj_list <- list(pdata, pdata_gdf, pdata_sg, pdata_sg_gdf, eggs_gdf)
  
  # One time tests
  # Try creating a total_countFilt object with an untoward input object.
  expect_error(total_count_filter(omicsData = fdata),
               paste("omicsData must be of class 'seqData'",
                     sep = ' '))
  
  # Repeat for all object types
  purrr::map(1:length(obj_list), function(n){
    tester <- obj_list[[n]]
    
    #### Filter tests ####
    
    # Run total_count_filter on the reduced data frame.
    filter <- total_count_filter(omicsData = tester)
    
    # Review the class for the filter object.
    expect_s3_class(filter,
                    c('total_countFilt', 'data.frame'))
    
    # Check the dimensions of filter.
    expect_equal(dim(filter),
                 c(1200, 2))
    
    # Calculate the total_count with R functions. The first column of edata is removed
    # because it contains the peptide IDs.
    use_r <- as.data.frame(apply(tester$e_data[ -1], 1, sum))
    use_r <- use_r[[1]]
    
    temp_data <- tester$e_data[ -1]
    samp_sum <- apply(temp_data, 
                      2, 
                      sum,
                      na.rm = TRUE) + 1
    
    # divide adjusted (ensure non-zero) counts by library size
    div_sum <- sweep((temp_data + .5), 2, samp_sum, `/`)
    
    # Apply per million multiplier and log2
    temp_data <- log2(div_sum * 10^6)
    
    temp_data <- cbind(tester$e_data[1], temp_data)
    
    use_r_lcpm <- reshape2::melt(temp_data, id.var = "ID_REF", 
                                 value.name = "lcpm", variable.name = "Samples")
    use_r_lcpm$ID_REF <- as.character(use_r_lcpm$ID_REF)
    use_r_lcpm$Samples <- as.character(use_r_lcpm$Samples)
    
    use_r_lcpm <- dplyr::arrange(use_r_lcpm, 
                                    !!rlang::sym(get_edata_cname(tester)),
                                    !!rlang::sym(get_fdata_cname(tester)), 
                                    lcpm)
    
    # Ensure the total_count values are correct.
    expect_equal(sort(filter$Total_Counts),
                 sort(use_r))
    
    # Inspect the attributes of the total_countFilt object.
    expect_identical(attr(filter, "e_data_lcpm"),
                     use_r_lcpm)
    
    #### Apply filt tests ####
    
    # Apply the filter without groups to the reduced peptide data set.
    filtered <- applyFilt(filter_object = filter,
                          omicsData = tester,
                          min_count = 2)
    
    # Ensure the class and attributes that shouldn't have changed didn't change.
    expect_identical(attr(tester, 'cnames'),
                     attr(filtered, 'cnames'))
    expect_identical(attr(tester, 'check.names'),
                     attr(filtered, 'check.names'))
    expect_identical(class(tester),
                     class(filtered))
    
    # Examine the filters attribute.
    expect_equal(attr(filtered, 'filters')[[1]]$type,
                 'totalCountFilt')
    expect_identical(attr(filtered, 'filters')[[1]]$threshold,
                     2)
    
    expect_equal(sort(as.character(attr(filtered, 'filters')[[1]]$filtered)),
                 sort(as.character(filter$ID_REF[filter$Total_Counts < 2]))) ############### change after autofilter corrections
    expect_equal(attr(filtered, 'filters')[[1]]$method, NA)
    
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
    
    ## less samps for sg groups
    n_samps <- ifelse(n %in% c(3, 4), 31, 40)
    n_bios <- ifelse(n %in% c(3, 4), 1061, 1070) ## extra is filtered
    
    # Inspect the filtered e_data, f_data, and e_meta data frames.
    expect_equal(dim(filtered$e_data),
                 c(n_bios, n_samps + 1))
    
    expect_equal(dim(filtered$f_data),
                 c(n_samps, 4))
    
    expect_error(applyFilt(filter_object = filter,
                           omicsData = tester,
                           min_count = 1),
                   "min_count must be an integer greater than or equal to 2"
                   )
    
    
  })
  
  # Object Specific tests
  
  # Run the RNA filter on the data with scrambled columns.
  filter_gdf <- total_count_filter(omicsData = pdata_gdf)
  filter_eggs <- total_count_filter(omicsData = eggs_gdf)
  
  # Sleuth around the scrambled data
  expect_equal(filter_gdf,  filter_eggs)

  # Expect warning if data has already been filtered ---------------------------
  
  filter1 <- total_count_filter(omicsData = pdata)
  filter2 <- total_count_filter(omicsData = pdata)
  filter3 <- total_count_filter(omicsData = pdata)
  filtered1 <- applyFilt(filter1,
                         pdata,
                         min_count = 2)
  warnings <- capture_warnings(
    filtered2 <- applyFilt(filter2,
                           filtered1,
                           min_count = 2)
  )
  expect_match(
    warnings,
    "A total count filter has already been applied to this data set\\.",
    all = FALSE
  )
  expect_match(
    warnings,
    paste0(
      "Specified biomolecules Name_100040016, Name_100041811, Name_10031",
      "6710, Name_100316712, Name_100504323, Name_100526527, Name_100529075, ",
      "Name_114600, Name_114640, Name_116911, Name_11922, Name_12951, Name_12",
      "958, Name_140489, Name_14237, Name_14380, Name_15228, Name_16637, Name",
      "_16814, Name_171206, Name_171245, Name_171253, Name_171262, Name_18328",
      ", Name_207151, Name_210155, Name_210876, Name_215854, Name_216783, Nam",
      "e_21907, Name_22040, Name_22109, Name_228731, Name_229550, Name_230025",
      ", Name_235256, Name_23793, Name_238568, Name_245492, Name_245509, Name",
      "_257872, Name_258086, Name_258166, Name_258181, Name_258201, Name_2582",
      "25, Name_258249, Name_258276, Name_258330, Name_258355, Name_258361, N",
      "ame_258431, Name_258444, Name_258461, Name_258488, Name_258497, Name_2",
      "58536, Name_258570, Name_258621, Name_258624, Name_258628, Name_258652",
      ", Name_258684, Name_258716, Name_258747, Name_258752, Name_258767, Nam",
      "e_258818, Name_258853, Name_258865, Name_258893, Name_258899, Name_258",
      "904, Name_258917, Name_258944, Name_259005, Name_259008, Name_259019, ",
      "Name_259024, Name_259045, Name_259072, Name_27140, Name_27423, Name_31",
      "9579, Name_329404, Name_381816, Name_382244, Name_384696, Name_385328,",
      " Name_387170, Name_387230, Name_387344, Name_387349, Name_387354, Name",
      "_404288, Name_433449, Name_435207, Name_435951, Name_50774, Name_50929",
      ", Name_546648, Name_546912, Name_60505, Name_624765, Name_627537, Name",
      "_637898, Name_667240, Name_667262, Name_669149, Name_670464, Name_6740",
      "2, Name_70993, Name_71057, Name_73353, Name_73360, Name_73470, Name_73",
      "5267, Name_735283, Name_735285, Name_73887, Name_73972, Name_74441, Na",
      "me_74529, Name_74987, Name_75470, Name_75718, Name_75861, Name_76405, ",
      "Name_78471, Name_93968 were not found in e_data\\."
    ),
    all = FALSE
  )
  
  expect_identical(filtered1$e_data, filtered2$e_data)
  
  warnings <- capture_warnings(
    filtered3 <- applyFilt(filter3,
                           filtered1,
                           min_count = 3)
  )
  expect_match(
    warnings,
    "A total count filter has already been applied to this data set\\.",
    all = FALSE
  )
  expect_match(
    warnings,
    paste0(
      "Specified biomolecules Name_100040016, Name_100041811, Name_10031",
      "6710, Name_100316712, Name_100504323, Name_100526527, Name_100529075, ",
      "Name_114600, Name_114640, Name_116911, Name_11922, Name_12951, Name_12",
      "958, Name_140489, Name_14237, Name_14380, Name_15228, Name_16637, Name",
      "_16814, Name_171206, Name_171245, Name_171253, Name_171262, Name_18328",
      ", Name_207151, Name_210155, Name_210876, Name_215854, Name_216783, Nam",
      "e_21907, Name_22040, Name_22109, Name_228731, Name_229550, Name_230025",
      ", Name_235256, Name_23793, Name_238568, Name_245492, Name_245509, Name",
      "_257872, Name_258086, Name_258166, Name_258181, Name_258201, Name_2582",
      "25, Name_258249, Name_258276, Name_258330, Name_258355, Name_258361, N",
      "ame_258431, Name_258444, Name_258461, Name_258488, Name_258497, Name_2",
      "58536, Name_258570, Name_258621, Name_258624, Name_258628, Name_258652",
      ", Name_258684, Name_258716, Name_258747, Name_258752, Name_258767, Nam",
      "e_258818, Name_258853, Name_258865, Name_258893, Name_258899, Name_258",
      "904, Name_258917, Name_258944, Name_259005, Name_259008, Name_259019, ",
      "Name_259024, Name_259045, Name_259072, Name_27140, Name_27423, Name_31",
      "9579, Name_329404, Name_381816, Name_382244, Name_384696, Name_385328,",
      " Name_387170, Name_387230, Name_387344, Name_387349, Name_387354, Name",
      "_404288, Name_433449, Name_435207, Name_435951, Name_50774, Name_50929",
      ", Name_546648, Name_546912, Name_60505, Name_624765, Name_627537, Name",
      "_637898, Name_667240, Name_667262, Name_669149, Name_670464, Name_6740",
      "2, Name_70993, Name_71057, Name_73353, Name_73360, Name_73470, Name_73",
      "5267, Name_735283, Name_735285, Name_73887, Name_73972, Name_74441, Na",
      "me_74529, Name_74987, Name_75470, Name_75718, Name_75861, Name_76405, ",
      "Name_78471, Name_93968 were not found in e_data\\."
    ),
    all = FALSE
  )
  
  expect_false(identical(filtered1, filtered3))
})
