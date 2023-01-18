context('plot functions')

test_that('plot functions are producing desired output',{
  library(pmartRdata)
  library(vdiffr)
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
  result <- edata_summary(omicsData = mylipid, by = "molecule", groupvar = "Virus")
  expect_doppelganger("plot.dataRes", plot(result))
  expect_doppelganger("plot.dataRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger("plot.dataRes (bw_theme)", plot(result, bw_theme = FALSE))
  
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  result <- normalize_isobaric(myiso, exp_cname = "Plex", 
                                   apply_norm = FALSE,
                                   refpool_cname = "Virus", 
                                   refpool_notation = "Pool")
  expect_doppelganger("plot.isobaricnormRes", plot(result))
  expect_doppelganger("plot.isobaricnormRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger("plot.isobaricnormRes (bw_theme)", plot(result, bw_theme = FALSE))
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  mynmrnorm <- normalize_nmr(omicsData = mynmr,
                             apply_norm = FALSE,
                             metabolite_name = "unkm1.53")
  expect_doppelganger("plot.nmrnormRes", plot(mynmrnorm))
  expect_doppelganger("plot.nmrnormRes", plot(mynmrnorm))
  mynmrnorm2 <- normalize_nmr(omicsData = mynmr,
                              apply_norm = FALSE,
                              sample_property_cname = "Concentration")
  expect_doppelganger("plot.nmrnormRes (2)", plot(mynmrnorm))
  expect_doppelganger("plot.nmrnormRes (color_by)", plot(mynmrnorm, nmrData=mynmr, color_by="Time"))
  
  mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
  mypep <- group_designation(omicsData = mypep, main_effects = "Phenotype")
  result <- spans_procedure(omicsData = mypep)
  expect_doppelganger("plot.SPANSRes", plot(result))
  expect_doppelganger("plot.SPANSRes (bw_theme)", plot(result, bw_theme = FALSE))
  expect_doppelganger("plot.SPANSRes (color_high color_low)", plot(result, color_high = "#00FFFF", color_low = "#FF0000"))
  
  mylipid <- group_designation(omicsData = lipid_neg_object, main_effects = "Virus")
  mylipid <- group_designation(omicsData = lipid_neg_object, main_effects = "Virus")
  result <- missingval_result(omicsData = mylipid)
  expect_doppelganger("plot.naRes (bar)", 
                      plot(naRes_obj = result, 
                           omicsData = mylipid,
                           plot_type = "bar", 
                           x_lab_angle = 50, 
                           order_by = "Virus", 
                           color_by = "Virus")
  )
  expect_doppelganger("plot.naRes (scatter)", 
                      plot(naRes_obj = result,
                           omicsData = mylipid,
                           plot_type = "scatter",
                           x_lab_angle = 50,
                           color_by = "Virus")
  )
  
  result <- missingval_result(omicsData = rnaseq_object)
  expect_doppelganger("plot.naRes (bar 2)", 
                      plot(naRes_obj = result,
                           omicsData = rnaseq_object,
                           plot_type = "bar")
  )
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab, main_effects = "Phenotype")
  my_correlation <- cor_result(omicsData = mymetab)
  expect_doppelganger("plot.corRes", plot(my_correlation, omicsData = mymetab, order_by = "Phenotype"))
  
  myseq_correlation <- cor_result(omicsData = rnaseq_object)
  expect_doppelganger("plot.corRes (2)", plot(myseq_correlation))
  
  set.seed(31415926)
  # dim_reduction does not produce a consistent result, need to store the result somewhere else
  mylipid <- edata_transform(omicsData = lipid_neg_object, data_scale="log2")
  mylipid <- group_designation(omicsData = mylipid, main_effects = "Virus")
  pca_lipids <- dim_reduction(omicsData = mylipid)
  expect_doppelganger("plot.dimRes", plot(pca_lipids))
  
  myseq <- group_designation(omicsData = rnaseq_object, main_effects = "Virus")
  pca_seq <- dim_reduction(omicsData = myseq)
  expect_doppelganger("plot.dimRes (2)", plot(pca_seq))
  
  data(pep_object)
  molfilt <- molecule_filter(omicsData = pep_object)
  expect_doppelganger("plot.moleculeFilt", plot(molfilt, min_num = 5))
  expect_doppelganger("plot.moleculeFilt (cumulative)", plot(molfilt, min_num = 3, cumulative = FALSE))
  
  seqfilt <- total_count_filter(omicsData = rnaseq_object)
  expect_doppelganger("plot.totalCountFilt", plot(seqfilt, min_count = 5))
  
  seqfilt <- RNA_filter(omicsData = rnaseq_object)
  expect_doppelganger("plot.RNAFilt", plot(seqfilt))
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
  to_filter <- imdanova_filter(omicsData = mypep)
  expect_doppelganger("plot.imdanovaFilt", plot(to_filter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3))
  
  my_filter <- proteomics_filter(omicsData = pep_object)
  expect_doppelganger("plot.proteomicsFilt", plot(my_filter, min_num_peps = 3))
  expect_doppelganger("plot.proteomicsFilt (redundancy)", plot(my_filter, plot_type = "redundancy"))
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab, main_effects = "Phenotype")
  rmd_results <- rmd_filter(omicsData = mymetab, metrics=c("MAD", "Skewness", "Correlation"))
  expect_doppelganger("plot.rmdFilt", plot(rmd_results, pvalue_threshold = 0.0001, order_by = "Phenotype"))
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
  cvfilt <- cv_filter(omicsData = mypep)
  expect_doppelganger("plot.cvFilt", plot(cvfilt, cv_threshold = 20))
  expect_doppelganger("plot.cvFilt (log_scale)", plot(cvfilt, cv_threshold = 10, log_scale = FALSE))
  
  mymetab <- edata_transform(omicsData = metab_object,
                                  data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab,
                                    main_effects = "Phenotype")
  norm_object <- normalize_global(omicsData = mymetab,
                                  subset_fn = "all",
                                  norm_fn = "median")
  expect_doppelganger("plot.normRes", plot(norm_object, order_by = "Phenotype", color_by = "Phenotype"))
  
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  expect_doppelganger("plot.isobaricpepData", plot(myiso))
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
  expect_doppelganger("plot.lipidData", plot(omicsData = mylipid, order_by = "Virus", color_by = "Virus"))
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  expect_doppelganger("plot.metabData", plot(omicsData = mymetab, order_by = "Phenotype", color_by = "Phenotype"))
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  expect_doppelganger("plot.nmrData", plot(omicsData = mynmr))
  
  expect_doppelganger("plot.seqData", plot(omicsData = rnaseq_object, transformation = "lcpm"))
  
  mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
  expect_doppelganger("plot.pepData", plot(omicsData = mypep, order_by = "Phenotype", color_by = "Phenotype"))
  
  expect_doppelganger("plot.proData", plot(omicsData = pro_object, order_by = "Phenotype", color_by = "Phenotype"))
  
  mypro <- group_designation(omicsData = pro_object,
                             main_effects = c("Phenotype"))
  imdanova_Filt <- imdanova_filter(omicsData = mypro)
  mypro <- applyFilt(filter_object = imdanova_Filt,
                         omicsData = mypro,
                         min_nonmiss_anova=2)
  anova_res <- imd_anova(omicsData = mypro, test_method = 'anova')
  expect_doppelganger("plot.statRes (anova)", plot(anova_res))
  expect_doppelganger("plot.statRes (anova volcano)", plot(anova_res, plot_type = "volcano"))
  imd_res <- imd_anova(omicsData = mypro, test_method = 'gtest')
  expect_doppelganger("plot.statRes (gtest)", plot(imd_res))
  imd_anova_res <- imd_anova(omicsData = mypro,
                             test_method = 'comb',
                             pval_adjust_a ='bon',
                             pval_adjust_g = 'bon')
  expect_doppelganger("plot.statRes (combined)", plot(imd_anova_res, bw_theme = TRUE))
  expect_doppelganger("plot.statRes (combined volcano)", plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE))
})