library(vdiffr)
context('plot functions')

test_that('plot functions are producing desired output',{
  # some tests inexplicably fail on github actions, skip them if this is true
  IS_CI = Sys.getenv("CI") == "true"  

  set.seed(31415926)
  
  ## Create each of the omicsData objects that we'll use throughout ------------
  
  # pepData
  load(system.file('testdata', 'little_pdata.RData', package = 'pmartR'))
  pep_object <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Mass_Tag_ID',
    fdata_cname = 'SampleID',
    emeta_cname = 'Protein'
  )
  
  # isobaricpepData
  load(system.file('testdata', 'little_isodata.RData', package = 'pmartR'))
  isobaric_object <- as.isobaricpepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Peptide',
    fdata_cname = 'Sample',
    emeta_cname = 'Protein'
  )
  
  # proData
  load(system.file('testdata', 'little_prdata.RData', package = 'pmartR'))
  pro_object <- as.proData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Reference',
    fdata_cname = 'SampleID',
    emeta_cname = 'PClass'
  )
  
  # lipidData
  load(system.file('testdata', 'lipidData.RData', package = "pmartR"))
  lipid_pos_object <- as.lipidData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'LipidCommonName',
    fdata_cname = 'Sample_Name',
    emeta_cname = 'LipidClass'
  )
  
  # metabData
  load(system.file('testdata', 'metaboliteData.RData', package = 'pmartR'))
  metab_object <- as.metabData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Metabolite',
    fdata_cname = 'SampleID',
    emeta_cname = 'MClass'
  )
  
  # nmrData
  load(system.file('testdata', 'nmrData.RData', package = "pmartR"))
  nmr_identified_object <- as.nmrData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Metabolite',
    fdata_cname = 'SampleID',
    emeta_cname = 'nmrClass'
  )
  
  # seqData
  load(system.file('testdata', 'little_seqdata.RData', package = 'pmartR'))
  rnaseq_object <- as.seqData(
    e_data = edata,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  )
  
  ## Load a pre-calculated SPANS result ----------------------------------------
  
  load(system.file('testdata', 'plot_objects.RData', package = "pmartR"))
  
  ## Test plot.dataRes ---------------------------------------------------------
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
  result <- edata_summary(omicsData = mylipid, by = "molecule", groupvar = "Condition")
  expect_doppelganger("plot.dataRes", plot(result))
  expect_doppelganger("plot.dataRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger("plot.dataRes (bw_theme)", plot(result, bw_theme = TRUE))
  
  ## Test plot.naRes -----------------------------------------------------------
  
  mylipid <- group_designation(omicsData = lipid_pos_object, main_effects = "Condition")
  result <- missingval_result(omicsData = mylipid)
  expect_doppelganger("plot.naRes (bar)", 
                      plot(naRes_obj = result, 
                           omicsData = mylipid,
                           plot_type = "bar", 
                           x_lab_angle = 50, 
                           order_by = "Condition", 
                           color_by = "Condition")
  )
  expect_doppelganger("plot.naRes (scatter)", 
                      plot(naRes_obj = result,
                           omicsData = mylipid,
                           plot_type = "scatter",
                           x_lab_angle = 50,
                           color_by = "Condition")
  )
  
  ## Test plot.nmrnormRes ------------------------------------------------------
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  mynmrnorm <- normalize_nmr(omicsData = mynmr,
                             apply_norm = FALSE,
                             metabolite_name = "unkm1.53")
  expect_doppelganger("plot.nmrnormRes", plot(mynmrnorm))
  
  mynmrnorm2 <- normalize_nmr(omicsData = mynmr,
                              apply_norm = FALSE,
                              sample_property_cname = "Concentration")
  expect_doppelganger("plot.nmrnormRes (2)", plot(mynmrnorm))
  expect_doppelganger("plot.nmrnormRes (color_by)", plot(mynmrnorm, nmrData=mynmr, color_by="Time"))
  
  ## Test plot.SPANSRES --------------------------------------------------------
  
  expect_doppelganger("plot.SPANSRes", plot(pep_spans_result))
  expect_doppelganger("plot.SPANSRes (bw_theme)", plot(pep_spans_result, bw_theme = TRUE))
  expect_doppelganger("plot.SPANSRes (color_high color_low)", plot(pep_spans_result, color_high = "#00FFFF", color_low = "#FF0000"))
  
  ## Test plot.isobaricnormRes -------------------------------------------------
  
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  result <- normalize_isobaric(
    myiso, 
    exp_cname = "Set", 
    apply_norm = FALSE,
    refpool_cname = "Reference", 
    refpool_notation = "Yes"
  )
  
  result_obj <- normalize_isobaric(
    myiso, 
    exp_cname = "Set", 
    apply_norm = TRUE,
    refpool_cname = "Reference", 
    refpool_notation = "Yes"
  )
  
  result_obj_norm <- normalize_global(
    result_obj,
    norm_fn = "mean",
    subset_fn = "all",
    apply_norm = TRUE
  )
  
  expect_doppelganger("plot.isobaricnormRes", plot(result))
  expect_doppelganger("plot.isobaricnormRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger("plot.isobaricnormRes (bw_theme)", plot(result, bw_theme = FALSE))
  
  expect_doppelganger("plot.isobaricnormRes (global normalized)", plot(result_obj_norm, bw_theme = FALSE))
  
  ## Test plot.corRes ----------------------------------------------------------
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab, main_effects = "Condition")
  my_correlation <- cor_result(omicsData = mymetab)
  expect_doppelganger("plot.corRes", plot(my_correlation, omicsData = mymetab, order_by = "Condition"))
  
  ## Test plot.dimRes ----------------------------------------------------------
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale="log2")
  mylipid <- group_designation(omicsData = mylipid, main_effects = "Condition")
  pca_lipids <- dim_reduction(omicsData = mylipid)
  expect_doppelganger("plot.dimRes", plot(pca_lipids))
  
  ## Test plot.moleculeFilt ----------------------------------------------------
  
  molfilt <- molecule_filter(omicsData = pep_object)
  expect_doppelganger("plot.moleculeFilt", plot(molfilt, min_num = 5))
  expect_doppelganger("plot.moleculeFilt (cumulative)", plot(molfilt, min_num = 3, cumulative = FALSE))
  
  ## Test plot.imdanovaFilt ----------------------------------------------------
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Condition")
  to_filter <- imdanova_filter(omicsData = mypep)
  expect_doppelganger("plot.imdanovaFilt", plot(to_filter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3))
  
  ## Test plot.proteomicsFilt --------------------------------------------------
  
  my_filter <- proteomics_filter(omicsData = pep_object)
  expect_doppelganger("plot.proteomicsFilt", plot(my_filter, min_num_peps = 3))
  expect_doppelganger("plot.proteomicsFilt (redundancy)", plot(my_filter, plot_type = "redundancy"))
  
  ## Test plot.rmdFilt ---------------------------------------------------------
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab, main_effects = "Condition")
  rmd_results <- rmd_filter(omicsData = mymetab, metrics=c("MAD", "Skewness", "Correlation"))
  expect_doppelganger("plot.rmdFilt", plot(rmd_results, pvalue_threshold = 0.01, order_by = "Condition"))
  
  ## Test plot.cvFilt ----------------------------------------------------------
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Condition")
  cvfilt <- cv_filter(omicsData = mypep)
  expect_doppelganger("plot.cvFilt", plot(cvfilt, cv_threshold = 20))
  expect_doppelganger("plot.cvFilt (log_scale)", plot(cvfilt, cv_threshold = 10, log_scale = FALSE))
  
  ## Test plot.normRes ---------------------------------------------------------
  
  mymetab <- edata_transform(omicsData = metab_object,
                             data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab,
                               main_effects = "Condition")
  norm_object <- normalize_global(omicsData = mymetab,
                                  subset_fn = "all",
                                  norm_fn = "median")
  expect_doppelganger("plot.normRes", plot(norm_object, order_by = "Condition", color_by = "Condition"))
  
  ## Test plot.statRres
  
  mypro <- edata_transform(omicsData = pro_object, data_scale = "log2")
  mypro <- group_designation(omicsData = mypro,
                             main_effects = "Condition")
  imdanova_Filt <- imdanova_filter(omicsData = mypro)
  mypro <- applyFilt(filter_object = imdanova_Filt,
                     omicsData = mypro,
                     min_nonmiss_anova=2)
  anova_res <- imd_anova(omicsData = mypro, test_method = 'anova')
  expect_doppelganger("plot.statRes (anova)", plot(anova_res))

  if(!IS_CI) {
    expect_doppelganger("plot.statRes (anova volcano)", plot(anova_res, plot_type = "volcano"))
  }

  imd_res <- imd_anova(omicsData = mypro, test_method = 'gtest')
  expect_doppelganger("plot.statRes (gtest)", plot(imd_res))
  imd_anova_res <- imd_anova(omicsData = mypro,
                             test_method = 'comb',
                             pval_adjust_a_multcomp ='bon',
                             pval_adjust_g_multcomp = 'bon')
  expect_doppelganger("plot.statRes (combined)", plot(imd_anova_res, bw_theme = TRUE))

  if(!IS_CI) {
    expect_doppelganger("plot.statRes (combined volcano)", plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE))
  }
  ## Test plot.totalcountFilt --------------------------------------------------
  
  seqfilt <- total_count_filter(omicsData = rnaseq_object)
  if(!IS_CI) {
    expect_doppelganger("plot.totalCountFilt", plot(seqfilt, min_count = 5))
  }
  
  ## Test plot.RNAFilt ---------------------------------------------------------
  
  seqfilt <- RNA_filter(omicsData = rnaseq_object)
  if(!IS_CI) {
    expect_doppelganger("plot.RNAFilt", plot(seqfilt))
  }
  
  ## Test plot.(omicsData_type) ------------------------------------------------
  
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  expect_doppelganger("plot.isobaricpepData", plot(myiso))
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
  expect_doppelganger("plot.lipidData", plot(omicsData = mylipid, order_by = "Condition", color_by = "Condition"))
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  expect_doppelganger("plot.metabData", plot(omicsData = mymetab, order_by = "Condition", color_by = "Condition"))
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  expect_doppelganger("plot.nmrData", plot(omicsData = mynmr, order_by = "Condition", color_by = "Condition"))
  
  mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
  expect_doppelganger("plot.pepData", plot(omicsData = mypep, order_by = "Condition", color_by = "Condition"))
  
  mypro <- edata_transform(omicsData = pro_object, data_scale = "log2")
  expect_doppelganger("plot.proData", plot(omicsData = pro_object, order_by = "Condition", color_by = "Condition"))
  
  myseq <- group_designation(omicsData = rnaseq_object, main_effects = "Tissue")
  expect_doppelganger("plot.seqData", plot(omicsData = rnaseq_object, transformation = "lcpm"))
  
})