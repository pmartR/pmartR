context('class: seqData')


### TO do: add emeta

test_that('diffexp_seq returns the correct data frame and attributes',{
  
  # Load the reduced peptide data frames ---------------------------------------
  
  load(system.file('testdata',
                   'little_seqdata.RData',
                   package = 'pmartR'))
  
  # Run as.seqData with agreeable data frames ----------------------------------
  
  # Construct a seqData object with the edata, fdata, and emeta data frames.
  seqdata <- as.seqData(e_data = edata,
                        f_data = fdata,
                        edata_cname = 'ID_REF',
                        fdata_cname = 'Samples'
  )
  
  ## add new stuff to dfs
  ## Always factor levels
  seqdata$f_data$Pair <- paste0("Pair_", rep(1:20, 2))
  seqdata$f_data$Pair_group <- paste0("Potato_", c(rep(1, 20), rep(2, 20)))
  
  ## Alt with numeric main effect
  # seqdata$f_data$group2 <- rnorm(40)
  
  seqdata$f_data$covar <- rep(rnorm(20), 2) ## numeric, redundant with pair id
  seqdata$f_data$covar2 <- rep(c(1,1, 1,1, rnorm(16)), 2) ## numeric, not redundant with pair id
  seqdata$f_data$covar3 <- rep(c(1,1, 1,1, 5:20), 2) ## factor, not redundant with pair id
  seqdata$f_data$covar4 <- as.numeric(as.factor(seqdata$f_data$Treatment))  ## factor, redundant with group
  
  expect_error(diffexp_seq(seqdata), 
               "Running group_designation is required")
  
  ## Factor main effect
  seqdata_grp <- group_designation(seqdata, main_effects = c("Tissue", "Treatment"))
  
  ## numeric main effects ## -> this is automatically read as a factor level
  # seqdata_grp2 <- suppressWarnings(group_designation(seqdata, main_effects = c("group2")))
  
  ## Factor main effects and covar
  seqdata_grp3 <- group_designation(seqdata, main_effects = c("Tissue", "Treatment"),
                                    covariates = "covar")
  
  seqdata_grp4 <- group_designation(seqdata, main_effects = c("Treatment"),
                                    covariates = "covar4")
  
  ## Treatment and paired groupp ## works
  seqdata_grp5 <- group_designation(seqdata, main_effects = c("Treatment"),
                                    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1")
  
  ## Test rejection of redundant -- not considered redundant since pair_group is changed for processing
  seqdata_grp6 <- group_designation(seqdata, main_effects = c("Treatment"), covariates = "covar",
                                    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1")
  
  ## Non-redundant ## works
  seqdata_grp7 <- group_designation(seqdata, main_effects = c("Treatment"), covariates = "covar2",
                                    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1")
  
  ## pairs only
  seqdata_grp8 <- group_designation(seqdata,
                                    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1")
  
  
  list_of_omics <- list(seqdata_grp, #seqdata_grp2, 
                        seqdata_grp3, seqdata_grp4,
                        seqdata_grp5, seqdata_grp6,
                        seqdata_grp7, seqdata_grp8
                        )
  
  #### Argument errors ####
  
  expect_error(diffexp_seq(seqdata_grp, method = "fhk"),
               "method must a single character string of length 1 in 'edgeR', 'DESeq2', or 'voom'")
  
  expect_error(diffexp_seq(seqdata_grp$f_data, method = "edgeR"), 
               "omicsData must be of class 'seqData'")
  
  expect_error(diffexp_seq(seqdata_grp, p_adjust = "tiger"), 
               "p_adjust must a single character string of")
  
  expect_error(diffexp_seq(seqdata_grp, p_cutoff = "tiger"), 
               "p_cutoff must numeric of length 1")
  
  # Custom comparisons ---------------
  comparisons <- data.frame(
    Test = "uterus_hCG",
    Control = "uterus_PBS"
  )
  
  bad_comparisons <- data.frame(
    Test = "uterus_hC",
    Control = "uterus_PBS"
  )

  err_comps <- "Invalid comparisons given"
  res <- purrr::map(c("edgeR", "DESeq2", "voom"), function(method){
    expect_error(diffexp_seq(seqdata_grp, method = method, comparisons = bad_comparisons), err_comps)
  })
  
  ## Ok versions ##
  run_all <- purrr::map(c("edgeR", "DESeq2", "voom"), function(method){
    
    ## basic
    basic_list <- purrr::map(1:length(list_of_omics), function(n){
      
      omics <- list_of_omics[[n]]
      
      ## warnings tested in get_group_formula as well as minor warnings:
      # character -> factor, rowname/colname coef magic
      suppressWarnings(diffexp_seq(omics, method = method))
    })
    
    ## pvalue adjust
    p_adjust_none <-  list(diffexp_seq(seqdata_grp, method = method, p_adjust = "none"))

    p_adjust_holm <- list(diffexp_seq(seqdata_grp, method = method, p_adjust = "holm"))
    
    ## custom list
    comps <- list(diffexp_seq(seqdata_grp, method = method, comparisons = comparisons))
    
    
    return(c(basic_list, p_adjust_none, p_adjust_holm, comps))
  })
  
  names(run_all) <- c("edgeR", "DESeq2", "voom")
  

  ## check IDs on all dfs for consistency with edata
  res <- purrr::map(run_all, function(method_res){
    purrr::map(method_res, function(df){
      expect_true(all(df$ID_REF == seqdata_grp$e_data$ID_REF))
    })
  })
  
  
  ## Check dimensions
  
  
  ## Check with previously saved values
  
  ## check using ... in calls
  fitType <- "parametric"
  fittype_change <- diffexp_seq(seqdata_grp, method = "DESeq2", fitType = fitType)
  
  logratioTrim <- 0.2
  robust_change <- diffexp_seq(seqdata_grp, logratioTrim = 0.2)
  
  
})