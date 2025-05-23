# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Construct the data to test the as.nmrData function. --------------------------

# I am saving these data sets (even though they are exact copies) in case the
# data sets in the pmartRdata package change in the future. If they do, the
# tests will not need to be updated to reflect changes in the new data.

# Load necessary libraries.
library (pmartRdata)

# Load the nmr data objects.
data("nmr_edata_identified")
data("nmr_fdata_identified")
data("nmr_emeta_identified")

# These data sets are small, allowing us to use the entire data set for testing
# purposes.
edata <- nmr_edata_identified
fdata <- nmr_fdata_identified
emeta <- nmr_emeta_identified

# Fashion an nmr "type" column in the emeta data frame.
emeta <- data.frame(check.names = FALSE, emeta,
                    nmrClass = sub("_.*", "", emeta[, 2]))
# The sub function extracts all characters before the first _.

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/nmrData.RData')
