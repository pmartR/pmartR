context('principal component analysis')

test_that('PCA produces the correct output',{
  
  # Load data and prepare omicsData objects ------------------------------------
  
  load(system.file('testdata',
                   'lipidData.RData',
                   package = 'pmartR'))
  
  ldata <- as.lipidData(e_data = edata,
                        f_data = fdata,
                        e_meta = emeta,
                        edata_cname = 'LipidCommonName',
                        fdata_cname = 'Sample_Name',
                        emeta_cname = 'LipidClass')
  
  # Log the heck out of the data.
  ldata <- edata_transform(ldata, "log")
  
  # Run the group_designation function on ldata.
  ldata_g <- group_designation(omicsData = ldata,
                               main_effects = "Condition")
  
  # Calculate PCA standards ----------------------------------------------------
  
  # No group_DF attribute ---------------
  
  # Find rows with fewer than two non-NA values. These will be removed before
  # performing the PCA.
  naughty <- which(rowSums(!is.na(ldata$e_data[, -1])) < 2)
  
  # Set a seed before running the pca function because it has random components.
  # This seed will need to match the seed before the dim_reduction function.
  set.seed(338)
  pals <- pcaMethods::pca(
    object = as.matrix(t(ldata$e_data[-naughty, -1])),
    method = "ppca",
    scale = "vector",
    nPcs = 3
  )
  
  standard <- list(SampleID = names(ldata$e_data)[-1],
                   PC1 = as.numeric(pals@scores[, 1]),
                   PC2 = as.numeric(pals@scores[, 2]),
                   PC3 = as.numeric(pals@scores[, 3]))
  
  class(standard) <- "dimRes"
  attr(standard, "group_DF") <- attr(ldata, "group_DF")
  attr(standard, "R2") <- pals@R2
  attr(standard, "row.names") <- names(ldata$e_data[, -1])
  
  # With group_DF attribute ---------------
  
  # Set a seed before running the pca function because it has random components.
  # This seed will need to match the seed before the dim_reduction function.
  set.seed(40)
  pals_g <- pcaMethods::pca(
    object = as.matrix(t(ldata_g$e_data[-naughty, -1])),
    method = "ppca",
    scale = "vector",
    nPcs = 2
  )
  
  standard_g <- list(SampleID = names(ldata_g$e_data)[-1],
                     PC1 = as.numeric(pals_g@scores[, 1]),
                     PC2 = as.numeric(pals_g@scores[, 2]))
  
  class(standard_g) <- "dimRes"
  attr(standard_g, "group_DF") <- attr(ldata_g, "group_DF")
  attr(standard_g, "R2") <- pals_g@R2
  attr(standard_g, "row.names") <- names(ldata_g$e_data[, -1])
  
  # Test the heck out of the dimRes objects ------------------------------------
  
  set.seed(338)
  expect_warning(
    dimmer <- dim_reduction(ldata, k = 3),
    paste("group_designation has not been run on this data and may limit",
          "plotting options",
          sep = " ")
  )
  expect_identical(dimmer, standard)
  
  set.seed(40)
  dimmer_g <- dim_reduction(ldata_g, k = 2)
  expect_identical(dimmer_g, standard_g)
  
})
