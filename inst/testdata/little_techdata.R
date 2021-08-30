# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Generate data for tests that include tech reps. ------------------------------

# Load necessary libraries.
library (pmartRdata)

# Load the tech rep data objects.
data("techrep_edata")
data("techrep_fdata")

# Take a small subset of the rows in techrep_edata to use for testing.
edata <- techrep_edata[1:150, ]

# Keep a copy of the original techrep_fdata data frame.
fdata <- techrep_fdata

# save(edata,
#      fdata,
#      file = '/pmartR/inst/testdata/little_techdata.RData')
