context('results: correlation')

test_that('cor_result correctly calcualtes the correlation',{
  
  # Load data and create omicsData objects ---------------------------------------
  
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
  
  # Natural logate the data.
  pdata <- edata_transform(omicsData = pdata,
                           data_scale = "log")
  
  load(system.file('testdata',
                   'little_isodata.RData',
                   package = 'pmartR'))
  
  # Construct an isobaricpepData object.
  isodata <- as.isobaricpepData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Peptide',
                                fdata_cname = 'Sample',
                                emeta_cname = 'Protein')
  
  # Logify the isobaric data.
  isodata <- edata_transform(isodata, "log")
  
  load(system.file('testdata',
                   'nmrData.RData',
                   package = 'pmartR'))
  
  # Produce a nmrData object.
  nmrdata <- as.nmrData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'Metabolite',
                        fdata_cname = 'SampleID',
                        emeta_cname = 'nmrClass')
  
  # Logitate the data nmr data.
  nmrdata <- edata_transform(nmrdata, "log")
  
  # Correlation: pepData object ------------------------------------------------
  
  # Correlation standard.
  standard <- cor(pdata$e_data[, -1],
                      use = "pairwise.complete.obs")
  
  # Add the correct class to the standard.
  class(standard) <- c("corRes", "matrix", "array")
  
  # Add the attributes to the standard.
  attr(standard, "sample_names") <- names(pdata$e_data[, -1])
  attr(standard, "is_normalized") <- FALSE
  attr(standard, "cor_method") <- "pearson"
  
  # Use cor_result to calculate the correlation.
  iocor <- cor_result(pdata)
  
  # Compare the standard to the cor_result output.
  expect_identical(standard, iocor)
  
  # Correlation: isobaricpepData object ----------------------------------------
  
  # Correlation standard.
  standard <- cor(isodata$e_data[, -1],
                  use = "pairwise.complete.obs")
  
  # Add the correct class to the standard.
  class(standard) <- c("corRes", "matrix", "array")
  
  # Add the attributes to the standard.
  attr(standard, "sample_names") <- names(isodata$e_data[, -1])
  attr(standard, "isobaric_norm") <- FALSE
  attr(standard, "is_normalized") <- FALSE
  attr(standard, "cor_method") <- "pearson"
  
  # Use cor_result to calculate the correlation.
  iocor <- cor_result(isodata)
  
  # Compare the standard to the cor_result output.
  expect_identical(standard, iocor)
  
  # Correlation: nmrData object ------------------------------------------------
  
  # Correlation standard.
  standard <- cor(nmrdata$e_data[, -1],
                  use = "pairwise.complete.obs")
  
  # Add the correct class to the standard.
  class(standard) <- c("corRes", "matrix", "array")
  
  # Add the attributes to the standard.
  attr(standard, "sample_names") <- names(nmrdata$e_data[, -1])
  attr(standard, "nmr_norm") <- FALSE
  attr(standard, "is_normalized") <- FALSE
  attr(standard, "cor_method") <- "pearson"
  
  # Use cor_result to calculate the correlation.
  iocor <- cor_result(nmrdata)
  
  # Compare the standard to the cor_result output.
  expect_identical(standard, iocor)
  
})
