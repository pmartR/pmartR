# The following script will select a subset of the peptide data from the
# pmartRdata package. The e_data and e_meta data frames will have some rows with
# a sum (row sum of non-missing values) less than the filter threshold in order
# to test the molecule_filter functions.

# Load the peptide data from pmartRdata.
pdata <- pmartRdata::pep_object

# Filter by molecule.
pdata_mol <- molecule_filter(pdata)

# Determine which columns contain a 1 or 2 (these will be filtered out later).
leq_2 <- which(pdata_mol$Num_Observations <= 2)

# Obtain the indices that have values greater than 2.
geq_3 <- which(pdata_mol$Num_Observations >= 3)

set.seed(40)

# Randomly select a proportion of data that will be filtered.
prop_filter <- rbeta(1, shape1 = 1, shape2 = 1)

# compute the number of observations that will come from either leq_2 or geq_3.
filtered <- floor(150 * prop_filter)
non_filtered <- 150 - filtered

# Subset the pep_edata and pep_emeta objects.
edata_mol <- pmartRdata::pep_edata[c(leq_2[1:filtered],
                                     geq_3[1:non_filtered]), ]
emeta_mol <- pmartRdata::pep_emeta[c(leq_2[1:filtered],
                                     geq_3[1:non_filtered]), ]

# Import the original f_data data frame because it is not affected by subsetting
# e_data and e_meta.
fdata_mol <- pmartRdata::pep_fdata

# Save the reduced e_data and e_meta data frames along with the original f_data
# data frame.
save(edata_mol,
     emeta_mol,
     fdata_mol,
     file = '~/pmartR/inst/testdata/filter_data_mol.RData')
