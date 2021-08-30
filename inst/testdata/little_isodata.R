# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Generate data to test the as.isobaricpepData function. -----------------------

# Load necessary libraries.
library (pmartRdata)

# Load the peptide data objects.
data("isobaric_edata")
data("isobaric_fdata")
data("isobaric_emeta")

# Take a small subset of the rows in isobaric_edata and isobaric_emeta to use
# for testing.
edata <- isobaric_edata[1:150, ]
emeta <- isobaric_emeta[1:150, ]

# Amend the first column of emeta to be the same as the first column of edata.
# This ensure that all the peptides in edata are present in emeta. This changes
# the data from the original but it doesn't matter for testing purposes.
emeta[, 1] <- edata[, 1]

# Remove factors from data columns to reduce the size of the .RData file.
edata$Peptide <- as.character(edata$Peptide)
emeta$Peptide <- as.character(emeta$Peptide)
emeta$Protein <- as.character(emeta$Protein)

# Keep a copy of the original pep_fdata data frame.
fdata <- isobaric_fdata

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/little_isodata.RData')
