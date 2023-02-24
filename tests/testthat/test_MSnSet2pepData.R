context('MSnSet conversion')

test_that('MSnSet objects are correctly converted to pepData objects',{
  
  # Load data from the MSnbase package -----------------------------------------
  
  the_set <- MSnbase::msnset
  
  # Convert an MSnSet object to a pepData object.
  msn_pepe <- MSnSet2pepData(the_set,
                             data_scale = "log2",
                             edata_cname = "UniqueID",
                             fdata_cname = "SampleID",
                             emeta_cname = "UniqueID")
  
  # Holy MSnSet unit tests, Batman! --------------------------------------------
  
  # Ensure the returned data frames are the correct dimension.
  expect_equal(dim(msn_pepe$e_data),
               c(55, 5))
  expect_equal(dim(msn_pepe$f_data),
               c(4, 3))
  expect_equal(dim(msn_pepe$e_meta),
               c(55, 16))
  
  # Confirm the correct attributes are present in the pepData object.
  expect_equal(names(attributes(msn_pepe)),
               c("names", "cnames", "data_info", "meta_info",
                 "filters", "class"))
  
  # Scrutinize the column names attribute.
  expect_equal(attr(msn_pepe, "cnames"),
               list(edata_cname = "UniqueID",
                    emeta_cname = "UniqueID",
                    fdata_cname = "SampleID",
                    techrep_cname = NULL))
  
  # Investigate the elements of the data_info attribute.
  expect_equal(
    attr(msn_pepe, "data_info"),
    list(data_scale_orig = "log2",
         data_scale = "log2",
         norm_info = list(is_normalized = FALSE),
         num_edata = length(unique(msn_pepe$e_data[, 1])),
         num_miss_obs = sum(is.na(msn_pepe$e_data)),
         prop_missing = (sum(is.na(msn_pepe$e_data)) /
                           prod(dim(msn_pepe$e_data[, -1]))),
         num_samps = ncol(msn_pepe$e_data[, -1]),
         data_types = NULL,
         batch_info = list(is_bc = FALSE))
  )
  
  # Inspect the elements of the meta_info attribute.
  expect_equal(
    attr(msn_pepe, "meta_info"),
    list(meta_data = TRUE,
         num_emeta = length(unique(msn_pepe$e_meta$UniqueID)))
  )
  
  # Take a looksie at the filters attribute.
  expect_identical(attr(msn_pepe, "filters"), list())
  
  # Ensure the omicsData object is classy.
  expect_s3_class(msn_pepe, "pepData")
  
})
