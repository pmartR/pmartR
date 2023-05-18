#### Helper script to load objects for testing plots. 

# Run SPANS now so that it doesn't have to be done at the time of the tests,
# since we aren't testing SPANS, we're testing the plot function.

load(system.file('testdata', 'little_pdata.RData', package = 'pmartR'))
pep_object <- as.pepData(
  e_data = edata,
  f_data = fdata,
  e_meta = emeta,
  edata_cname = 'Mass_Tag_ID',
  fdata_cname = 'SampleID',
  emeta_cname = 'Protein'
)
mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
mypep <- group_designation(omicsData = mypep, main_effects = "Condition")
pep_spans_result <- spans_procedure(omicsData = mypep)

# Save the result --------------------------------------------------------------

# save(
# 	pep_spans_result,
# 	file = '~/pmartR/inst/testdata/plot_objects.RData'
# )