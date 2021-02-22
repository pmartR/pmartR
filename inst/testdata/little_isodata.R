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
emeta[, 1] <- edata[, 1]

# Keep a copy of the original pep_fdata data frame.
fdata <- isobaric_fdata

save(edata,
     fdata,
     emeta,
     file = '~/pmartR/inst/testdata/little_isodata.RData')