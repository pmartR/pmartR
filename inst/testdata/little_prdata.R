# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Generate data to test the as.proData function. -------------------------------

# Load necessary libraries.
library (pmartRdata)

# Load the peptide data objects.
data("pro_edata")
data("pro_fdata")

# Take a small subset of the rows in pro_edata to use for testing.
edata <- pro_edata[1:150, ]

# Remove factors from data columns to reduce the size of the .RData file.
edata$Reference <- as.character(edata$Reference)

# Keep a copy of the original pro_fdata data frame.
fdata <- pro_fdata

set.seed(338)

# Construct an emeta data frame from edata.
emeta <- data.frame(check.names = FALSE, Reference = edata[, 1],
                    # Make up 38 different protein classes.
                    PClass = sample(x = paste0('PClass', 1:38),
                                    size = 150,
                                    replace = TRUE))

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/little_prdata.RData')
