# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Fabricate the data to test the as.metabData function. ------------------------

# I am saving these data sets (even though they are exact copies) in case the
# data sets in the pmartRdata package change in the future. If they do, the
# tests will not need to be updated to reflect changes in the new data.

# Load necessary libraries.
library (pmartRdata)

# Load the metabolite data objects.
data("metab_edata")
data("metab_fdata")

# These data sets are small, allowing us to use the entire data set for testing
# purposes.
edata <- metab_edata
fdata <- metab_fdata

set.seed(22)

# Construct an emeta data frame from edata.
emeta <- data.frame(check.names = FALSE, Metabolite = edata[, 1],
                    # Make up 17 different metabolite classes.
                    MClass = sample(x = paste0('MClass', 1:17),
                                    size = 80,
                                    replace = TRUE))

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/metaboliteData.RData')
