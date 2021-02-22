# Generate data to test the as.proData function. -------------------------------

# Load necessary libraries.
library (pmartRdata)

# Load the peptide data objects.
data("pro_edata")
data("pro_fdata")

# Take a small subset of the rows in pro_edata to use for testing.
edata <- pro_edata[1:150, ]

# Keep a copy of the original pro_fdata data frame.
fdata <- pro_fdata

set.seed(338)

# Construct an emeta data frame from edata.
emeta <- data.frame(Reference = edata[, 1],
                    # Make up 38 different protein classes.
                    PClass = sample(x = paste0('PClass', 1:38),
                                    size = 150,
                                    replace = TRUE))

save(edata,
     fdata,
     emeta,
     file = '~/pmartR/inst/testdata/little_prdata.RData')
