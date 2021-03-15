# The purpose of this script is to demonstrate how the transfigured pepData
# objects were created. It should NOT be rerun because the data sets in
# pmartRdata will change over time and this will lead to errors in the unit
# tests.

# Load: pepData ----------------------------------------------------------------

load(system.file('testdata',
                 'little_pdata.RData',
                 package = 'pmartR'))

# Transform the pepData objects ------------------------------------------------

# Use all possible transform options to create additional pepData objects.

# Construct a pepData object.
pdata <- as.pepData(e_data = edata,
                    f_data = fdata,
                    e_meta = emeta,
                    edata_cname = 'Mass_Tag_ID',
                    fdata_cname = 'SampleID',
                    emeta_cname = 'Protein')

# Natural log transfiguration.
lpdata <- edata_transform(pdata,
                          "log")

# Log base 2 transfiguration.
l2pdata <- edata_transform(pdata,
                           "log2")

# Log base 10 transfiguration.
l10pdata <- edata_transform(pdata,
                            "log10")

# Save the transfigured pepData objects.
# save(pdata,
#      lpdata,
#      l2pdata,
#      l10pdata,
#      file = '~/pmartR/inst/testdata/transmuted_pdata.RData')
