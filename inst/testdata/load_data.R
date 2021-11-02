#### Helper script to load objects for testing.

# pepData (pdata)
load(system.file('testdata',
                 'little_pdata.RData',
                 package = 'pmartR'))

pdata <- as.pepData(e_data = edata,
                    f_data = fdata,
                    e_meta = emeta,
                    edata_cname = 'Mass_Tag_ID',
                    fdata_cname = 'SampleID',
                    emeta_cname = 'Protein')

# proData (prdata)
load(system.file('testdata',
                 'little_prdata.RData',
                 package = 'pmartR'))

prdata <- as.proData(e_data = edata,
                     f_data = fdata,
                     e_meta = emeta,
                     edata_cname = 'Reference',
                     fdata_cname = 'SampleID',
                     emeta_cname = 'PClass')

# metabData (mdata)
load(system.file('testdata',
                 'metaboliteData.RData',
                 package = 'pmartR'))

mdata <- as.metabData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Metabolite',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'MClass')

# lipidData (ldata)
load(system.file('testdata',
                 'lipidData.RData',
                 package = 'pmartR'))

ldata <- as.lipidData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'LipidCommonName',
                      fdata_cname = 'Sample_Name',
                      emeta_cname = 'LipidClass')

load(system.file('testdata',
                 'nmrData.RData',
                 package = 'pmartR'))

# nmrData (nmrdata)
nmrdata <- as.nmrData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Metabolite',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'nmrClass')

# isobaricpepData (isodata)
load(system.file('testdata',
                 'little_isodata.RData',
                 package = 'pmartR'))

isodata <- as.isobaricpepData(e_data = edata,
                              f_data = fdata,
                              e_meta = emeta,
                              edata_cname = 'Peptide',
                              fdata_cname = 'Sample',
                              emeta_cname = 'Protein')
