# Generate data to test the as.pepData function. -------------------------------

# Load necessary libraries.
library (pmartRdata)

# Load the peptide data objects.
data("pep_edata")
data("pep_fdata")
data("pep_emeta")

# Take a small subset of the rows in pep_edata and pep_emeta to use for testing.
edata <- pep_edata[1:150, ]
emeta <- pep_emeta[1:150, ]

# Keep a copy of the original pep_fdata data frame.
fdata <- pep_fdata

save(edata,
     fdata,
     emeta,
     file = '~/pmartR/inst/testdata/little_pdata.RData')
