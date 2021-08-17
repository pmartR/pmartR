# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Generate data to test the as.pepData function. -------------------------------

# Load necessary libraries.
library (pmartRdata)

# Load the peptide data objects.
data("pep_edata")
data("pep_fdata")
data("pep_emeta")

# Take a small subset of the rows in pep_edata and pep_emeta to use for testing.
edata <- pep_edata[1:150, ]

# Find all rows in emeta that match the Mass_Tag_IDs from edata. These IDs will
# be used to subset the full emeta object.
id_idx <- which(pep_emeta$Mass_Tag_ID %in% edata$Mass_Tag_ID)

emeta <- pep_emeta[id_idx, ]

# Remove factor levels from emeta variables. This greatly reduces the size of
# the .RData file.
emeta$Protein <- as.character(emeta$Protein)
emeta$Peptide_Sequence <- as.character(emeta$Peptide_Sequence)

# Keep a copy of the original pep_fdata data frame.
fdata <- pep_fdata

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/little_pdata.RData')
