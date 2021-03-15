# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Forge the data to test the as.lipidData function. ----------------------------

# I am saving these data sets (even though they are exact copies) in case the
# data sets in the pmartRdata package change in the future. If they do, the
# tests will not need to be updated to reflect changes in the new data.

# Load necessary libraries.
library (pmartRdata)

# Load the lipid data objects.
data("lipid_edata")
data("lipid_fdata")

# These data sets are small, allowing us to use the entire data set for testing
# purposes.
edata <- lipid_edata
fdata <- lipid_fdata

# Construct an emeta data frame from edata.
emeta <- data.frame(LipidCommonName = edata[, 1],
                    LipidClass = sub("\\(.*", "", edata[, 1]))
# The sub function extracts all characters before the first (.

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/lipidData.RData')
