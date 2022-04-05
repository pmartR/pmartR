context("class: trelliData edata")

test_that("An edata object passed to trelliData edata returns correct data frames and attributes",{
  
  # Load: peptide expression data-----------------------------------------------
  
  load(system.file('testdata',
                   'little_pdata.Rdata',
                   package = 'pmartR'))
  
  # Test: peptide expression data-----------------------------------------------
  
  pep_trelli_edata <- as.trelliData.edata(
    e_data = edata, 
    edata_cname = "Mass_Tag_ID",
    omics_type = "pepData"
  )
  
  
  
})