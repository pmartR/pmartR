context('RMD conversion')

test_that('conversion between log2(RMD) and p-values is seamless',{
  
  # Load data and prepare omicsData objects ------------------------------------
  
  load(system.file('testdata',
                   'metaboliteData.RData',
                   package = 'pmartR'))
  
  mdata <- as.metabData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'MClass')
  
  mdata <- edata_transform(omicsData = mdata,
                           data_scale = "log")
  
  mdata <- group_designation(omicsData = mdata,
                             main_effects = "Condition")
  
  # Holy RMD conversion tests, Batman! -----------------------------------------
  
  rmd_filta <- rmd_filter(omicsData = mdata,
                          metrics = c("MAD", "Skewness", "Correlation"))
  
  expect_identical(
    rmd_conversion(log2rmd = rmd_filta$Log2.md, df = 3),
    1 - pchisq(2^(rmd_filta$Log2.md), df = 3)
  )
  
  expect_identical(
    rmd_conversion(pval = .0001, df = 5),
    log(qchisq(1 - 0.0001, df = 5), base = 2)
  )
  
  expect_identical(
    rmd_conversion(log2rmd = 4.5, df = 3),
    1 - pchisq(2^(4.5), df = 3)
  )
  
})
