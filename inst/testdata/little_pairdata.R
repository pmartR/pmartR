# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created from the paired data .RData file.

# The following data set may be removed from pmartR later on.
# load("~/pmartR/inst/paired_data_for_test.RData")

library (tidyverse)

# Nab a subset of the data -----------------------------------------------------

# Create the reduced edata object.
edata <- pep_data$e_data[1001:1150, ]

# Find all rows in emeta that match the Mass_Tag_IDs from edata. These IDs will
# be used to subset the full emeta object.
idx <- which(pep_data$e_meta$Mass_Tag_ID %in% edata$Mass_Tag_ID)

# Create the reduced emeta object
emeta <- pep_data$e_meta[idx, ]

# Keep a copy of the original f_data object.
fdata <- pep_data$f_data

# Clean up the sample names ----------------------------------------------------

the_names <- names(edata)

# Remove any unneeded information from the sample names. I am the lord of this
# script so I will determine what is useful information!!!!! If anyone disagrees
# with my decisions they can grovel at the feet of my R powers and plead their
# case. I will be a merciful lord and grant them a small portion of my time to
# explain their petty thoughts. They will then be banished to the fourth level
# of the R inferno where they can contemplate their insignificant existence.
the_names[2:31] <- the_names[2:31] %>%
  str_split( "_") %>%
  map(`[`, 3:6) %>%
  map(paste0, collapse = "_") %>%
  unlist() %>%
  str_replace("691", "AM") %>%
  str_remove("protein_")

names(edata) <- the_names
fdata$Name <- the_names[2:31]
fdata$Virus <- str_replace(fdata$Virus, "691", "AM")

# save(edata,
#      fdata,
#      emeta,
#      file = '~/pmartR/inst/testdata/little_pairdata.RData')
