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

# Copy part of the pepData object for proteomics filter testing ----------------

# Extract the entire e_data, f_data, and e_meta objects from pmartRdata.
edata_pro <- pmartRdata::pep_edata[1:500, ]
fdata_pro <- pmartRdata::pep_fdata
emeta_pro <- pmartRdata::pep_emeta[1:500, ]

# Fabricate a pepData object with the reduced data set.
pdata_pro <- as.pepData(e_data = edata_pro,
                        f_data = fdata_pro,
                        e_meta = emeta_pro,
                        edata_cname = "Mass_Tag_ID",
                        fdata_cname = "SampleID",
                        emeta_cname = "Protein")

# Create the proteomics filter standards to compare test output to.
pfStandard <- proteomics_filter(pdata_pro)

# Generate the filtered omicsData standards to compare test output to.
afStandard <- applyFilt(filter_object = pfStandard,
                        omicsData = pdata_pro,
                        min_num_peps = 2)

# Create a reduced data frame for the imdanova filter --------------------------

# Without singleton groups ---------------

# Load the reduced peptide data frame.
load(system.file('testdata',
                 'little_pdata.RData',
                 package = 'pmartR'))

# Forge a pepData object.
pdata <- as.pepData(e_data = edata,
                    f_data = fdata,
                    e_meta = emeta,
                    edata_cname = "Mass_Tag_ID",
                    fdata_cname = "SampleID",
                    emeta_cname = "Protein")

# Run the group_designation function on pdata.
pdata <- group_designation(omicsData = pdata,
                           main_effects = "Condition")

# Run imdanova_filter on pdata with group information.
ifStandard <- imdanova_filter(omicsData = pdata)

# Filter the reduced pepData object using anova. af: anova filter.
afStandard <- applyFilt(filter_object = ifStandard,
                        omicsData = pdata,
                        min_nonmiss_anova = 3, 
                        min_nonmiss_gtest = NULL, 
                        remove_singleton_groups = FALSE)

# Filter the reduced pepData object using gtest. gf: gtest filter.
gfStandard <- applyFilt(filter_object = ifStandard,
                        omicsData = pdata,
                        min_nonmiss_anova = NULL, 
                        min_nonmiss_gtest = 3, 
                        remove_singleton_groups = FALSE)

# Filter the reduced pepData object using gtest. bf: both (anova, gtest) filter.
bfStandard <- applyFilt(filter_object = ifStandard,
                        omicsData = pdata,
                        min_nonmiss_anova = 3, 
                        min_nonmiss_gtest = 3, 
                        remove_singleton_groups = FALSE)

# With singleton groups ---------------

# Forge a pepData object with only one sample from Mock. This will create a
# singleton group when group designating the data. sg: singleton group.
pdata_sg <- as.pepData(e_data = edata[, 1:11],
                       f_data = fdata[1:10, ],
                       e_meta = emeta,
                       edata_cname = "Mass_Tag_ID",
                       fdata_cname = "SampleID",
                       emeta_cname = "Protein")

# Run the group_designation function on pdata_sg.
pdata_sg <- group_designation(omicsData = pdata_sg,
                              main_effects = "Condition")

# Run imdanova_filter on pdata_sg with group information.
ifStandard_sg <- imdanova_filter(omicsData = pdata_sg)

# Save data to the testdata directory ------------------------------------------

# Create a path to the folder where the data will be saved.
sPath <- file.path('/Users/mart077/OneDrive - PNNL/Documents/multi_probe',
                   '/pmartR/inst/testdata/imdanova_standards.RData')

# Save the data extracted from pmartRdata.
save(ifStandard,
     ifStandard_sg,
     afStandard,
     gfStandard,
     bfStandard,
     file = sPath)
