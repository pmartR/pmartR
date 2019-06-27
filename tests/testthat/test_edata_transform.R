context("test edata_transform")

# Load data #
mypepData <- pmartRdata::pep_object
myproData <- pmartRdata::pro_object
myisoData <- pmartRdata::isobaric_object
mytechData <- pmartRdata::techrep_pep_object
mylipData <- pmartRdata::lipid_object
mymetData <- pmartRdata::metab_object

# Set list of objects #
datalist <- list(mypepData, myproData, myisoData, 
                 mytechData, mylipData, mymetData)



# processed data #
datalist2 <- 


# Remove checknames problems #
datalist <- purrr::map(datalist, function(omicsData){
  attr(omicsData, "check.names") <- FALSE
  return(omicsData)
})

testthat::test_that("Edata dimensions, f_data, and e_meta are identical after transformation",{
  
  purrr::map(datalist, function(omicsData){
    
    transData <- omicsData # For comparing transformed data dimensions
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      transData <- edata_transform(transData, "abundance")
    }
    
    # Check e_data dimensions #
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(edata_transform(transData, "log")$e_data))
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(edata_transform(transData, "log2")$e_data))
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(edata_transform(transData, "log10")$e_data))
    
    
    # Check f_data df #
    testthat::expect_identical(
      omicsData$f_data,
      edata_transform(transData, "log")$f_data)
    testthat::expect_identical(
      omicsData$f_data,
      edata_transform(transData, "log2")$f_data)
    testthat::expect_identical(
      omicsData$f_data,
      edata_transform(transData, "log10")$f_data)
    
    # Check e_meta df #
    if (!is.null(omicsData$e_meta)){
      testthat::expect_identical(
        omicsData$e_meta,
        edata_transform(transData, "log")$e_meta)
      testthat::expect_identical(
        omicsData$e_meta,
        edata_transform(transData, "log2")$e_meta)
      testthat::expect_identical(
        omicsData$e_meta,
        edata_transform(transData, "log10")$e_meta)
    }
    
    
    ## Check abundance setting ##
    transData <- edata_transform(transData, "log10")
    testthat::expect_identical(
      dim(omicsData$e_data),
      dim(edata_transform(transData, "abundance")$e_data))
    testthat::expect_identical(
      omicsData$f_data,
      edata_transform(transData, "abundance")$f_data)
    if (!is.null(omicsData$e_meta)){
      testthat::expect_identical(
        omicsData$e_meta,
        edata_transform(transData, "abundance")$e_meta)
    }
  })
})


testthat::test_that("Edata transformed values are as expected",{
  
  purrr::map(datalist, function(omicsData){
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      omicsData <- edata_transform(omicsData, "abundance")
    }
    # Check values #
    testthat::expect_identical(
      log10(omicsData$e_data[2:ncol(omicsData$e_data)]),
      edata_transform(omicsData, "log10")$e_data[2:ncol(omicsData$e_data)])
    testthat::expect_identical(
      log(omicsData$e_data[2:ncol(omicsData$e_data)]),
      edata_transform(omicsData, "log")$e_data[2:ncol(omicsData$e_data)])
    testthat::expect_identical(
      log2(omicsData$e_data[2:ncol(omicsData$e_data)]),
      edata_transform(omicsData, "log2")$e_data[2:ncol(omicsData$e_data)])

    omicsData <- edata_transform(omicsData, "log10")
    testthat::expect_identical(
      10^(omicsData$e_data[2:ncol(omicsData$e_data)]),
      edata_transform(omicsData, "abundance")$e_data[2:ncol(omicsData$e_data)])
  })
})


testthat::test_that("Edata attributes are maintained across transformations",{
  
  purrr::map(datalist, function(omicsData){
    
    transData <- omicsData # For comparing transformed data attributes
    
    # Set to no scaling #
    if(pmartR:::get_data_scale(omicsData) != "abundance"){
      transData <- edata_transform(transData, "abundance")
    }
    
    # Check attributes #
    testthat::expect_identical(
      grep("data_scale|cname", 
           attributes(omicsData), 
           value = TRUE, invert = TRUE),
      grep("data_scale|cname", 
           attributes(edata_transform(transData, "log")), 
           value = TRUE, invert = TRUE))
    testthat::expect_identical(
      grep("data_scale|cname", 
           attributes(omicsData), 
           value = TRUE, invert = TRUE),
      grep("data_scale|cname", 
           attributes(edata_transform(transData, "log2")), 
           value = TRUE, invert = TRUE))
    testthat::expect_identical(
      grep("data_scale|cname", 
           attributes(omicsData), 
           value = TRUE, invert = TRUE),
      grep("data_scale|cname", 
           attributes(edata_transform(transData, "log10")), 
           value = TRUE, invert = TRUE))
    
    #### Data_norm and median center of data_info is not passed on from protein object
    # testthat::expect_identical(
    #   unlist(attr(omicsData, "data_info"))[-1],
    #   unlist(attr(edata_transform(transData, "log"), "data_info"))[-1])
    # testthat::expect_identical(
    #   unlist(attr(omicsData, "data_info"))[-1],
    #   unlist(attr(edata_transform(transData, "log2"), "data_info"))[-1])
    # testthat::expect_identical(
    #   unlist(attr(omicsData, "data_info"))[-1],
    #   unlist(attr(edata_transform(transData, "log10"), "data_info"))[-1])
    
    ##### Currently seperate to account for tech_rep cname abscence in objects
    testthat::expect_identical(
      unlist(attr(omicsData, "cnames")),
      unlist(attr(edata_transform(transData, "log"), "cnames")))
    testthat::expect_identical(
      unlist(attr(omicsData, "cnames")),
      unlist(attr(edata_transform(transData, "log2"), "cnames")))
    testthat::expect_identical(
      unlist(attr(omicsData, "cnames")),
      unlist(attr(edata_transform(transData, "log10"), "cnames")))
    
    # Check abundance settings #
    transData <- edata_transform(transData, "log10")
    testthat::expect_identical(
      grep("data_scale|cname", 
           attributes(omicsData), 
           value = TRUE, invert = TRUE),
      grep("data_scale|cname", 
           attributes(edata_transform(transData, "abundance")), 
           value = TRUE, invert = TRUE))
    testthat::expect_identical(
      unlist(attr(omicsData, "cnames")),
      unlist(attr(edata_transform(transData, "abundance"), "cnames")))
    
    #### Data_norm and median center of data_info is not passed on from protein object
    testthat::expect_identical(
      unlist(attr(omicsData, "data_info"))[-1],
      unlist(attr(edata_transform(transData, "abundance"), "data_info"))[-1])
    
  })
})
  
