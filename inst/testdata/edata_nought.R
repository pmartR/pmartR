# The data sets created in this file were generated with the old edata_replace
# function. These data sets cannot be recreated with this script. The purpose of
# this script is to show how the data sets were initially created.

# Load: isobaricpepData --------------------------------------------------------

load(system.file('testdata',
                 'little_isodata.RData',
                 package = 'pmartR'))

# Replace NAs: isobaricpepData -------------------------------------------------

# Fabricate an isobaricpepData object.
isodata <- as.isobaricpepData(e_data = edata,
                              f_data = fdata,
                              e_meta = emeta,
                              edata_cname = 'Peptide',
                              fdata_cname = 'Sample',
                              emeta_cname = 'Protein')

# Use edata_replace to change NAs to 0s with the isobaric peptide test data.
isodata2 <- edata_replace(omicsData = isodata,
                          x = NA,
                          y = 0)

# Extract edata (with zeros) from the isodata2 object.
iso_nil <- isodata2$e_data

# Ensure the number of zeros is correct. This should be 400.
sum(iso_nil == 0)

# Load: lipidData --------------------------------------------------------------

load(system.file('testdata',
                 'lipidData.RData',
                 package = 'pmartR'))

# Replace NAs: lipidData -------------------------------------------------------

# Produce a lipidData object.
ldata <- as.lipidData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'LipidCommonName',
                      fdata_cname = 'Sample_Name',
                      emeta_cname = 'LipidClass')

# Use edata_replace to change NAs to 0s with the lipid test data.
ldata2 <- edata_replace(omicsData = ldata,
                        x = NA,
                        y = 0)

# Extract edata (with zeros) from the ldata2 object.
lip_nil <- ldata2$e_data

# Ensure the number of zeros is correct. This should be 884.
sum(lip_nil == 0)

# Load: metabData --------------------------------------------------------------

load(system.file('testdata',
                 'metaboliteData.RData',
                 package = 'pmartR'))

# Replace NAs: metabData -------------------------------------------------------

# Forge a metabData object.
mdata <- as.metabData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Metabolite',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'MClass')

# Use edata_replace to change NAs to 0s with the metabolite test data.
mdata2 <- edata_replace(omicsData = mdata,
                        x = NA,
                        y = 0)

# Extract edata (with zeros) from the mdata2 object.
met_nil <- mdata2$e_data

# Ensure the number of zeros is correct. This should be 148.
sum(met_nil == 0)

# Load: nmrData ----------------------------------------------------------------

load(system.file('testdata',
                 'nmrData.RData',
                 package = 'pmartR'))

# Alter the original nmr data frame to contain zeros ---------------------------

set.seed(22)

# Select row indices.
rows <- sample(1:38, 338, replace = TRUE)

# Select column indices.
columns <- sample(2:42, 338, replace = TRUE)

# Loop through each row column pair.
for (e in 1:338) {
  
  # Replace current value in edata with zero.
  edata[rows[[e]], columns[[e]]] <- 0
  
}

# Create the nmr edata object with zeros.
nmr_nil <- edata

# Ensure the number of zeros is correct. This should be 302.
sum(nmr_nil == 0)

# Load: pepData ----------------------------------------------------------------

load(system.file('testdata',
                 'little_pdata.RData',
                 package = 'pmartR'))

# Replace NAs: pepData ---------------------------------------------------------

# Construct a pepData object.
pdata <- as.pepData(e_data = edata,
                    f_data = fdata,
                    e_meta = emeta,
                    edata_cname = 'Mass_Tag_ID',
                    fdata_cname = 'SampleID',
                    emeta_cname = 'Protein')

# Use edata_replace to change NAs to 0s with the peptide test data.
pdata2 <- edata_replace(omicsData = pdata,
                        x = NA,
                        y = 0)

# Extract edata (with zeros) from the pdata2 object.
pep_nil <- pdata2$e_data

# Ensure the number of zeros is correct. This should be 341.
sum(pep_nil == 0)

# Load: proData ----------------------------------------------------------------

load(system.file('testdata',
                 'little_prdata.RData',
                 package = 'pmartR'))

# Replace NAs: proData ---------------------------------------------------------

# Generate a proData object.
prdata <- as.proData(e_data = edata,
                     f_data = fdata,
                     e_meta = emeta,
                     edata_cname = 'Reference',
                     fdata_cname = 'SampleID',
                     emeta_cname = 'PClass')

# Use edata_replace to change NAs to 0s with the protein test data.
prdata2 <- edata_replace(omicsData = prdata,
                         x = NA,
                         y = 0)

# Extract edata (with zeros) from the prdata2 object.
pro_nil <- prdata2$e_data

# Ensure the number of zeros is correct. This should be 234.
sum(pro_nil == 0)

# Save the edata object with zeros in the data frame.
# save(iso_nil,
#      lip_nil,
#      met_nil,
#      nmr_nil,
#      pep_nil,
#      pro_nil,
#      file = '~/pmartR/inst/testdata/edata_nought.RData')
