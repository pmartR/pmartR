testthat::context("test pmartR::edata_transform")

##### Load and log-transform data #####

# Load data #
mypepData <- pmartRdata::pep_object
myproData <- pmartRdata::pro_object
myisoData <- pmartRdata::isobaric_object
mytechData <- pmartRdata::techrep_pep_object
mylipData <- pmartRdata::lipid_object
mymetData <- pmartRdata::metab_object

# Log transform untransformed data #
myisoData   <- pmartR::edata_transform(myisoData, "log2")
mylipData      <- pmartR::edata_transform(mylipData, "log2")
mymetData      <- pmartR::edata_transform(mymetData, "log2")
mypepData        <- pmartR::edata_transform(mypepData, "log2")
mytechData        <- pmartR::edata_transform(mytechData, "log2")

##### Check carried over attributes ######

testthat::test_that("Initial transformation carries all attribute similarities except data scale",{
  
  testthat::expect_match(
    all.equal(attributes(pmartRdata::pep_object), attributes(mypepData)),
    "data_scale")
  testthat::expect_equal(
    length(all.equal(attributes(pmartRdata::pep_object), attributes(mypepData))),
    1)
  testthat::expect_match(
    all.equal(attributes(pmartRdata::lipid_object), attributes(mylipData)),
    "data_scale")
  testthat::expect_equal(
    length(all.equal(attributes(pmartRdata::lipid_object), attributes(mylipData))),
    1)
  testthat::expect_match(
    all.equal(attributes(pmartRdata::isobaric_object), attributes(myisoData)),
    "data_scale")
  testthat::expect_equal(
    length(all.equal(attributes(pmartRdata::isobaric_object), attributes(myisoData))),
    1)
  testthat::expect_match(
    all.equal(attributes(pmartRdata::metab_object), attributes(mymetData)),
    "data_scale")
  testthat::expect_equal(
    length(all.equal(attributes(pmartRdata::metab_object), attributes(mymetData))),
    1)
  testthat::expect_match(
    all.equal(attributes(pmartRdata::techrep_pep_object), attributes(mytechData)),
    "data_scale")
  testthat::expect_equal(
    length(all.equal(attributes(pmartRdata::techrep_pep_object), attributes(mytechData))),
    1)

  # No Transformation for protein, NA
})

##### Check error throwing ######

testthat::test_that("edata_transform error throwing",{
  
  testthat::expect_error(pmartR::edata_transform("gfhjlkdsfhsd", "log2"), 
                         "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData")
  testthat::expect_error(pmartR::edata_transform(mypepData, "hgfghdhjlhld"), 
                         "is not a valid option for 'data_scale'. See details of as.pepData for specifics.")
  testthat::expect_error(pmartR::edata_transform(myproData, "log2"), 
                         "Data is already on log2 scale.")
  
  x <- pmartR::edata_transform(myproData, "log10")
  testthat::expect_error(pmartR::edata_transform(x, "log10"), 
                         "Data is already on")
  x <- pmartR::edata_transform(myproData, "log")
  testthat::expect_error(pmartR::edata_transform(x, "log"),
                         "Data is already on")
  x <- pmartR::edata_transform(myproData, "abundance")
  testthat::expect_error(pmartR::edata_transform(x, "abundance"), 
                         "Data is already on")

})
  

##### Fill possible object attributes ######

# isobaric norm #
myisoData <- pmartR::normalize_isobaric(myisoData, apply_norm = T, channel_cname = "iTRAQ.Channel", refpool_channel = "116")
myisoData <- pmartR::normalize_isobaric(myisoData, apply_norm = T, refpool_cname = "Reference", refpool_notation = "Yes")

# Set appropriate group designations #
myisoData    <- pmartR::group_designation(myisoData, "Group")
mylipData       <- pmartR::group_designation(mylipData, "Condition")
mymetData       <- pmartR::group_designation(mymetData, "Condition")
mypepData         <- pmartR::group_designation(mypepData, "Condition")
myproData         <- pmartR::group_designation(myproData, "Condition")
mytechData <- pmartR::group_designation(mytechData,
                                                c("FACTOR", "DILUTION"))

# Set filters #
myisoData <- pmartR::applyFilt(pmartR::imdanova_filter(myisoData),
                                      myisoData,
                                      min_nonmiss_anova = 2,
                                      min_nonmiss_gtest = 3)
mylipData <- pmartR::applyFilt(pmartR::imdanova_filter(mylipData),
                                  mylipData,
                                  min_nonmiss_anova = 2,
                                  min_nonmiss_gtest = 3)
mymetData <- pmartR::applyFilt(pmartR::imdanova_filter(mymetData),
                                  mymetData,
                                  min_nonmiss_anova = 2,
                                  min_nonmiss_gtest = 3)
mypepData <- pmartR::applyFilt(pmartR::imdanova_filter(mypepData),
                                mypepData,
                                min_nonmiss_anova = 2,
                                min_nonmiss_gtest = 3)
myproData <- pmartR::applyFilt(pmartR::imdanova_filter(myproData),
                                myproData,
                                min_nonmiss_anova = 2,
                                min_nonmiss_gtest = 3)
mytechData <- pmartR::applyFilt(pmartR::imdanova_filter(mytechData),
                                mytechData,
                                min_nonmiss_anova = 2,
                                min_nonmiss_gtest = 3)

# Normalize #

myisoData <- pmartR::normalize_global(omicsData = myisoData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)
mylipData <- pmartR::normalize_global(omicsData = mylipData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)
mymetData <- pmartR::normalize_global(omicsData = mymetData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)
mypepData <- pmartR::normalize_global(omicsData = mypepData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)
myproData <- pmartR::normalize_global(omicsData = myproData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)
mytechData <- pmartR::normalize_global(omicsData = mytechData, subset_fn = "all", 
                              norm_fn = "median", apply_norm = TRUE, 
                              backtransform = TRUE)

##### Test processed data #####

datalist <- list(mypepData, myproData, myisoData, 
                  mytechData, mylipData, mymetData)
  ##### Dataframe population #####
testthat::test_that("Edata dimensions, f_data, and e_meta are identical after transformation",{
  
  purrr::map(datalist, function(omicsData){
    
    transData <- omicsData # For comparing transformed data dimensions
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      transData <- pmartR::edata_transform(transData, "abundance")
    }
    
    # Check e_data dimensions #
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(pmartR::edata_transform(transData, "log")$e_data))
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(pmartR::edata_transform(transData, "log2")$e_data))
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(pmartR::edata_transform(transData, "log10")$e_data))
    
    
    # Check f_data df #
    testthat::expect_identical(
      omicsData$f_data,
      pmartR::edata_transform(transData, "log")$f_data)
    testthat::expect_identical(
      omicsData$f_data,
      pmartR::edata_transform(transData, "log2")$f_data)
    testthat::expect_identical(
      omicsData$f_data,
      pmartR::edata_transform(transData, "log10")$f_data)
    
    # Check e_meta df #
    if (!is.null(omicsData$e_meta)){
      testthat::expect_identical(
        omicsData$e_meta,
        pmartR::edata_transform(transData, "log")$e_meta)
      testthat::expect_identical(
        omicsData$e_meta,
        pmartR::edata_transform(transData, "log2")$e_meta)
      testthat::expect_identical(
        omicsData$e_meta,
        pmartR::edata_transform(transData, "log10")$e_meta)
    }
    
    
    ## Check abundance setting ##
    transData <- pmartR::edata_transform(transData, "log10")
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(pmartR::edata_transform(transData, "abundance")$e_data))
    testthat::expect_identical(
      omicsData$f_data,
      pmartR::edata_transform(transData, "abundance")$f_data)
    if (!is.null(omicsData$e_meta)){
      testthat::expect_identical(
        omicsData$e_meta,
        pmartR::edata_transform(transData, "abundance")$e_meta)
    }
  })
})


testthat::test_that("Edata transformed values are as expected",{
  
  purrr::map(datalist, function(omicsData){
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      omicsData <- pmartR::edata_transform(omicsData, "abundance")
    }
    # Check values #
    testthat::expect_identical(
      log10(omicsData$e_data[2:ncol(omicsData$e_data)]),
      pmartR::edata_transform(omicsData, "log10")$e_data[2:ncol(omicsData$e_data)])
    testthat::expect_identical(
      log(omicsData$e_data[2:ncol(omicsData$e_data)]),
      pmartR::edata_transform(omicsData, "log")$e_data[2:ncol(omicsData$e_data)])
    testthat::expect_identical(
      log2(omicsData$e_data[2:ncol(omicsData$e_data)]),
      pmartR::edata_transform(omicsData, "log2")$e_data[2:ncol(omicsData$e_data)])

    omicsData <- pmartR::edata_transform(omicsData, "log10")
    testthat::expect_identical(
      10^(omicsData$e_data[2:ncol(omicsData$e_data)]),
      pmartR::edata_transform(omicsData, "abundance")$e_data[2:ncol(omicsData$e_data)])
  })
})

  ##### Attribute population #####
testthat::test_that("All object attributes are maintained across transformations",{
  
  purrr::map(datalist, function(omicsData){
    
    transData <- omicsData # For comparing transformed data attributes
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      transData <- pmartR::edata_transform(transData, "abundance")
    }
    
    ## Check unchanging attributes; changed attributes proccesed below ##
      # data_scale expected to change -> grep removes data_info #
      # e_meta populates in names in edata_transform even if NULL -> grep removes names #
    
    testthat::expect_identical(
      grep("data_scale|e_data", 
           attributes(omicsData), value = TRUE, invert = TRUE),
      grep("data_scale|e_data", 
           attributes(pmartR::edata_transform(transData, "log")), value = TRUE, invert = TRUE))
    
    testthat::expect_identical(
      grep("data_scale|e_data",
           attributes(omicsData), value = TRUE, invert = TRUE),
      grep("data_scale|e_data",
           attributes(pmartR::edata_transform(transData, "log2")), value = TRUE, invert = TRUE))
    
    testthat::expect_identical(
      grep("data_scale|e_data",
           attributes(omicsData), value = TRUE, invert = TRUE),
      grep("data_scale|e_data",
           attributes(pmartR::edata_transform(transData, "log10")), value = TRUE, invert = TRUE))
    
    ## Test data_info attribute ##
      ## num_emeta registers even if emeta is NULL, data_scale expected to change;
      ## Unlisting removes "NULL" num_emeta for comparison, [-1] removes data_scale
    testthat::expect_identical(
      unlist(attr(omicsData, "data_info"))[-1],
      unlist(attr(pmartR::edata_transform(transData, "log"), "data_info"))[-1])
    testthat::expect_identical(
      unlist(attr(omicsData, "data_info"))[-1],
      unlist(attr(pmartR::edata_transform(transData, "log2"), "data_info"))[-1])
    testthat::expect_identical(
      unlist(attr(omicsData, "data_info"))[-1],
      unlist(attr(pmartR::edata_transform(transData, "log10"), "data_info"))[-1])

    
    ## Test names
    if(is.null(omicsData$e_meta)){
      # emeta is null, but names is populated. -3 removes emeta from transdata attribute
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log"), "names"))[-3])
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log2"), "names"))[-3])
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log10"), "names"))[-3])
    } else {
      # emeta is populated, correct names for both
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log"), "names")))
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log2"), "names")))
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "log10"), "names")))
    }
    
    # Check abundance settings for all attributes #
    transData <- pmartR::edata_transform(transData, "log10")
    
    testthat::expect_identical(
      grep("data_scale|e_data",
           attributes(omicsData),
           value = TRUE, invert = TRUE),
      grep("data_scale|e_data",
           attributes(pmartR::edata_transform(transData, "abundance")),
           value = TRUE, invert = TRUE))
    
    testthat::expect_identical(
      unlist(attr(omicsData, "data_info"))[-1],
      unlist(attr(pmartR::edata_transform(transData, "abundance"), "data_info"))[-1])
    
    if(is.null(omicsData$e_meta)){
      # emeta is null, but names is populated. -3 removes emeta from transdata attribute
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "abundance"), "names"))[-3])
    } else {
      # emeta is populated, correct names for both
      testthat::expect_identical(
        unlist(attr(omicsData, "names")),
        unlist(attr(pmartR::edata_transform(transData, "abundance"), "names")))
    }

    
    # Double Transform attribute check #
    
    testthat::expect_true(
      all.equal(
        attributes(pmartR::edata_transform(omicsData, "log")),
        attributes(pmartR::edata_transform(pmartR::edata_transform(omicsData, "log10"), "log"))))
    testthat::expect_true(
      all.equal(
        attributes(pmartR::edata_transform(omicsData, "log10")),
        attributes(pmartR::edata_transform(pmartR::edata_transform(omicsData, "log"), "log10"))))
    testthat::expect_true(
      all.equal(
        attributes(pmartR::edata_transform(omicsData, "abundance")),
        attributes(pmartR::edata_transform(pmartR::edata_transform(omicsData, "log"), "abundance"))))
    testthat::expect_true(
      all.equal(
        attributes(pmartR::edata_transform(transData, "log2")),
        attributes(pmartR::edata_transform(pmartR::edata_transform(transData, "log"), "log2"))))
  })
})

