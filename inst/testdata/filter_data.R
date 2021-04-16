# The following script will select a subset of the peptide data from the
# pmartRdata package. The subsetted data frames will be used for unit testing.
# The e_data and e_meta data frames will have a mix of rows that fall above and
# below the default threshold for the various filtering methods.
#
# The purpose of this script is to demonstrate how the e_data, f_data, and
# e_meta data frames were created. It should NOT be rerun because the data sets
# in pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Load the pepData object ------------------------------------------------------

# Load the peptide data from pmartRdata.
pdata <- pmartRdata::pep_object

# Forge data for molecular filter testing --------------------------------------

# Filter by molecule.
pdata_mol <- molecule_filter(pdata)

# Determine which columns contain a 1 or 2 (these will be filtered out later).
leq <- which(pdata_mol$Num_Observations <= 2)

# Obtain the indices that have values greater than 2.
geq <- which(pdata_mol$Num_Observations >= 3)

set.seed(40)

# Randomly select a proportion of data that will be filtered.
prop_filter <- rbeta(1, shape1 = 1, shape2 = 1)

# compute the number of observations that will come from either leq or geq.
filtered <- floor(150 * prop_filter)
non_filtered <- 150 - filtered

# Subset the pep_edata and pep_emeta objects.
edata_mol <- pmartRdata::pep_edata[c(leq[1:filtered],
                                     geq[1:non_filtered]), ]
emeta_mol <- pmartRdata::pep_emeta[c(leq[1:filtered],
                                     geq[1:non_filtered]), ]

# Import the original f_data data frame because it is not affected by subsetting
# e_data and e_meta.
fdata_mol <- pmartRdata::pep_fdata
