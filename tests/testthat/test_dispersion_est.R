context('class: seqData')

test_that('dispersion_est returns the correct data frame and attributes',{
  
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
  
  expect_error(dispersion_est(seqdata,  method = "edgeR"), 
               "OmicsData requires group_designation for statistical analysis")
  
  seqdata_grp <- group_designation(seqdata, main_effects = c("Tissue", "Treatment"))
  
  expect_error(dispersion_est(seqdata_grp, method = "fhk"),
               "method must a single character string of length 1 in 'edgeR', 'DESeq2', or 'voom'")
  
  expect_error(dispersion_est(seqdata_grp$f_data, method = "edgeR"), 
               "omicsData must be of class 'seqData'")
  
  edgeR_res <- dispersion_est(seqdata_grp, method = "edgeR")
  deseq2_res <- suppressWarnings(dispersion_est(seqdata_grp, method = "DESeq2"))
  voom_res <- dispersion_est(seqdata_grp, method = "voom")
  
  res <- purrr::map(list(edgeR_res, deseq2_res, voom_res), function(gg){
    expect_true(inherits(gg, "ggplot"))
  })
  
  ## Custom theme warning/error
  edgeR_res <- expect_error(expect_warning(
    dispersion_est(seqdata_grp, method = "edgeR", bw_theme = T, custom_theme = 2),
    "Setting both bw_theme to TRUE and specifying a custom"),
    "custom_theme must be a valid 'theme' object as used in ggplot")
  
  ## plotly result
  edgeR_res_py <- dispersion_est(seqdata_grp, method = "edgeR", interactive = T)
  deseq2_res_py <- suppressWarnings(dispersion_est(seqdata_grp, method = "DESeq2", interactive = T))
  voom_res_py <- dispersion_est(seqdata_grp, method = "voom", interactive = T)
  
  res <- purrr::map(list(edgeR_res_py, deseq2_res_py, voom_res_py), function(gg){
    expect_true(inherits(gg, "plotly"))
  })
  
})
