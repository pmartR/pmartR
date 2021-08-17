# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created from the paired data .RData file.

# The following data set no longer exists. It is being reduced to keep the size
# of the package small.
# load("~/pmartR/inst/paired_data_for_test.RData")

# Create the edata object and reduce it.
edata <- pep_data$e_data
edata <- edata[1:150, ]

# Create the emeta object and reduce it.
emeta <- pep_data$e_meta

# Find all rows in emeta that match the Mass_Tag_IDs from edata. These IDs will
# be used to subset the full emeta object.
id_idx <- which(emeta$Mass_Tag_ID %in% edata$Mass_Tag_ID)

emeta <- emeta[id_idx, ]

# Keep a copy of the original f_data object.
fdata <- pep_data$f_data

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/little_pairdata.RData')
