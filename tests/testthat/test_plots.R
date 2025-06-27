context('plot functions')

expect_doppelganger_ci <- function(...) {
  is_ci = Sys.getenv("CI") == "true"
  
  if (is_ci) {
    invisible()
  } else {
    vdiffr::expect_doppelganger(...)
  }
}

test_that('plot helper functions are working good!', {
  load(system.file('testdata', 'little_pdata.RData', package = 'pmartR'))
  
  # need a test fixture for this...
  edata['Mass Tag ID'] <- edata$Mass_Tag_ID
  edata$Mass_Tag_ID <- NULL
  
  emeta['Mass Tag ID'] <- emeta$Mass_Tag_ID
  emeta['Protein Space'] <- emeta$Protein
  emeta$Mass_Tag_ID <- NULL
  emeta$Protein <- NULL
  
  pep_object <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Mass Tag ID',
    fdata_cname = 'SampleID',
    emeta_cname = 'Protein Space'
  )
  
  mypep<- edata_transform(omicsData = pep_object, data_scale = "log2")
  mypep <- group_designation(omicsData = mypep,
                            main_effects = "Condition")
  imdanova_Filt <- imdanova_filter(omicsData = mypep)
  mypep <- applyFilt(filter_object = imdanova_Filt,
                     omicsData = mypep,
                     min_nonmiss_anova=2)
  anova_res <- imd_anova(omicsData = mypep, test_method = 'combined')
  
  volcano_df <- pmartR:::make_volcano_plot_df(anova_res)
  
  expect_equal(dim(volcano_df), c(dim(mypep$e_data)[1]*2, 10))
})

test_that('plot functions are producing desired output',{
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("vdiffr")
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
  expect_doppelganger_ci("plot.dataRes", plot(result))
  expect_doppelganger_ci("plot.dataRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger_ci("plot.dataRes (bw_theme)", plot(result, bw_theme = TRUE))
  
  ## Test plot.naRes -----------------------------------------------------------
  mylipid <- group_designation(omicsData = lipid_pos_object, main_effects = "Condition")
  result <- missingval_result(omicsData = mylipid)
  expect_doppelganger_ci("plot.naRes (bar, group color)", 
                      plot(result, 
                           omicsData = mylipid,
                           plot_type = "bar", 
                           x_lab_angle = 50, 
                           order_by = "Condition", 
                           color_by = "Group")
  )
  expect_doppelganger_ci("plot.naRes (bar, group order)", 
                      plot(result, 
                           omicsData = mylipid,
                           plot_type = "bar", 
                           x_lab_angle = 50, 
                           order_by = "Group", 
                           color_by = "Condition")
  )
  expect_doppelganger_ci("plot.naRes (scatter)", 
                      plot(result,
                           omicsData = mylipid,
                           plot_type = "scatter",
                           x_lab_angle = 50,
                           groups = T)
  )
  
  ## Test plot.nmrnormRes ------------------------------------------------------
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  mynmrnorm <- normalize_nmr(omicsData = mynmr,
                             apply_norm = FALSE,
                             metabolite_name = "unkm1.53")
  expect_doppelganger_ci("plot.nmrnormRes", plot(mynmrnorm))
  
  mynmrnorm2 <- normalize_nmr(omicsData = mynmr,
                              apply_norm = FALSE,
                              sample_property_cname = "Concentration")
  expect_doppelganger_ci("plot.nmrnormRes (2)", plot(mynmrnorm))
  expect_doppelganger_ci("plot.nmrnormRes (color_by)", plot(mynmrnorm, nmrData=mynmr, color_by="Time"))
  
  ## Test plot.SPANSRES --------------------------------------------------------
  
  expect_doppelganger_ci("plot.SPANSRes", 
                         plot(pep_spans_result))
  expect_doppelganger_ci("plot.SPANSRes (bw_theme)", 
                         plot(pep_spans_result, bw_theme = TRUE))
  expect_doppelganger_ci("plot.SPANSRes (color_high color_low)", 
                         plot(pep_spans_result, color_high = "#00FFFF", color_low = "#FF0000"))
  expect_doppelganger_ci("plot.SPANSRes (N biomolecule bar)", 
                         plot(pep_spans_result, Npep_bar = T))
  
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
  
  expect_doppelganger_ci("plot.isobaricnormRes", plot(result))
  expect_doppelganger_ci("plot.isobaricnormRes (palette)", plot(result, palette = "YlOrRd"))
  expect_doppelganger_ci("plot.isobaricnormRes (bw_theme)", plot(result, bw_theme = FALSE))
  
  expect_doppelganger_ci("plot.isobaricnormRes (global normalized)", plot(result_obj_norm, bw_theme = FALSE))
  
  ## Test plot.corRes ----------------------------------------------------------
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab, main_effects = "Condition")
  my_correlation <- cor_result(omicsData = mymetab)
  expect_doppelganger_ci("plot.corRes", plot(my_correlation, omicsData = mymetab, order_by = "Condition"))
  
  ## Test plot.dimRes ----------------------------------------------------------
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale="log2")
  mylipid <- group_designation(omicsData = mylipid, main_effects = "Condition")
  pca_lipids <- dim_reduction(omicsData = mylipid)
  expect_doppelganger_ci("plot.dimRes", plot(pca_lipids))
  
  ## Test plot.moleculeFilt ----------------------------------------------------
  
  molfilt <- molecule_filter(omicsData = pep_object)
  expect_doppelganger_ci("plot.moleculeFilt", plot(molfilt, min_num = 5))
  expect_doppelganger_ci("plot.moleculeFilt (cumulative)", plot(molfilt, min_num = 3, cumulative = FALSE))
  
  ## Test plot.imdanovaFilt ----------------------------------------------------
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Condition")
  to_filter <- imdanova_filter(omicsData = mypep)
  expect_doppelganger_ci("plot.imdanovaFilt", plot(to_filter, min_nonmiss_anova = 2, min_nonmiss_gtest = 3))
  
  ## Test plot.proteomicsFilt --------------------------------------------------
  
  my_filter <- proteomics_filter(omicsData = pep_object)
  expect_doppelganger_ci("plot.proteomicsFilt", plot(my_filter, min_num_peps = 3))
  expect_doppelganger_ci("plot.proteomicsFilt (redundancy)", plot(my_filter, plot_type = "redundancy"))
  
  ## Test plot.rmdFilt ---------------------------------------------------------
  
  mypep <- edata_transform(omicsData = pmartRdata::pep_object, data_scale = "log2")
  mypep$f_data$pair_id <- c(rep(1:4, 2), rep(5:8, 2), rep(9:12, 2))
  mypep$f_data$Third_pheno <- 1:2
  
  ## Multiple groups, w/ pairing
  
  mypep <- group_designation(omicsData = mypep, main_effects = c("Phenotype", "Third_pheno"), 
                             pair_id = "pair_id", pair_group = "SecondPhenotype", pair_denom = "B")
  groups_pair <- rmd_filter(omicsData = mypep)
  
  ## Multiple groups
  
  mypep <- group_designation(omicsData = mypep, main_effects = c("Phenotype", "Third_pheno"))
  groups <- rmd_filter(omicsData = mypep)
  
  ## One group, w/ pairing
  
  mypep <- group_designation(omicsData = mypep, main_effects = c("Phenotype"), 
                             pair_id = "pair_id", pair_group = "SecondPhenotype", pair_denom = "B")
  group_pair <- rmd_filter(omicsData = mypep)
  
  ## One group
  
  mypep <- group_designation(omicsData = mypep, main_effects = c("Phenotype"))
  group <- rmd_filter(omicsData = mypep)
  
  ## Only pairing
  
  mypep <- group_designation(omicsData = mypep,
                             pair_id = "pair_id", pair_group = "SecondPhenotype", pair_denom = "B")
  pair <- rmd_filter(omicsData = mypep)
  
  ## silly names
  legend_lab <- c("Hey", "There")
  
  rme_plots <- purrr::map2(list(groups_pair, groups, group_pair, group, pair),
                           list("gsp", "gs", "gp", "g", "p"),
                           function(rmd_results, label){
    
    suppressWarnings({
      suppressMessages({
        expect_doppelganger_ci(paste0("plot.rmdFilt.",label), 
                               plot(rmd_results))
        expect_doppelganger_ci(paste0("plot.rmdFilt.color.",label), 
                               plot(rmd_results, color_by = "SecondPhenotype"))
        expect_doppelganger_ci(paste0("plot.rmdFilt.hist.",label),
                               plot(rmd_results, hist = T))
        expect_doppelganger_ci(paste0("plot.rmdFilt.color.hist.",label),
                               plot(rmd_results, hist = T, color_by = "SecondPhenotype"))
        expect_doppelganger_ci(paste0("plot.rmdFilt.pval.",label),
                               plot(rmd_results, pvalue_threshold = 0.05))
        expect_doppelganger_ci(paste0("plot.rmdFilt.sample",label), 
                               plot(rmd_results, order_by = "SecondPhenotype", sampleID = "Sample_47_Phenotype3_A"))
        
        expect_doppelganger_ci(paste0("plot.rmdFilt.legend.",label),
                               plot(rmd_results, legend_lab = legend_lab))
        expect_doppelganger_ci(paste0("plot.rmdFilt.color.legend.",label),
                               plot(rmd_results, color_by = "SecondPhenotype", legend_lab = legend_lab))
        expect_doppelganger_ci(paste0("plot.rmdFilt.hist.legend.",label),
                               plot(rmd_results, hist = T, legend_lab = legend_lab))
        expect_doppelganger_ci(paste0("plot.rmdFilt.color.hist.legend.",label),
                               plot(rmd_results, hist = T, color_by = "SecondPhenotype", legend_lab = legend_lab))
        expect_doppelganger_ci(paste0("plot.rmdFilt.pval.legend.",label),
                               plot(rmd_results, pvalue_threshold = 0.05, legend_lab = legend_lab))
        expect_doppelganger_ci(paste0("plot.rmdFilt.sample.legend.",label),
                               plot(rmd_results, order_by = "SecondPhenotype", sampleID = "Sample_47_Phenotype3_A", legend_lab = legend_lab))
      })
    })
  })
  
  ## Test plot.cvFilt ----------------------------------------------------------
  
  mypep <- group_designation(omicsData = pep_object, main_effects = "Condition")
  cvfilt <- cv_filter(omicsData = mypep)
  expect_doppelganger_ci("plot.cvFilt", plot(cvfilt, cv_threshold = 20))
  expect_doppelganger_ci("plot.cvFilt (log_scale)", plot(cvfilt, cv_threshold = 10, log_scale = FALSE))
  
  ## Test plot.normRes ---------------------------------------------------------
  
  mymetab <- edata_transform(omicsData = metab_object,
                             data_scale = "log2")
  mymetab <- group_designation(omicsData = mymetab,
                               main_effects = "Condition")
  norm_object <- normalize_global(omicsData = mymetab,
                                  subset_fn = "all",
                                  norm_fn = "median")
  expect_doppelganger_ci("plot.normRes", plot(norm_object, order_by = "Condition", color_by = "Condition"))
  
  ## Test plot.statRres
  
  mypro <- edata_transform(omicsData = pro_object, data_scale = "log2")
  mypro <- group_designation(omicsData = mypro,
                             main_effects = "Condition")
  imdanova_Filt <- imdanova_filter(omicsData = mypro)
  mypro <- applyFilt(filter_object = imdanova_Filt,
                     omicsData = mypro,
                     min_nonmiss_anova=2)
  anova_res <- imd_anova(omicsData = mypro, test_method = 'anova')
  expect_doppelganger_ci("plot.statRes (anova)", plot(anova_res))

  expect_doppelganger_ci("plot.statRes (anova volcano)", plot(anova_res, plot_type = "volcano"))

  imd_res <- imd_anova(omicsData = mypro, test_method = 'gtest')
  expect_doppelganger_ci("plot.statRes (gtest)", plot(imd_res))
  imd_anova_res <- imd_anova(omicsData = mypro,
                             test_method = 'comb',
                             pval_adjust_a_multcomp ='bon',
                             pval_adjust_g_multcomp = 'bon')
  expect_doppelganger_ci("plot.statRes (combined)", plot(imd_anova_res, bw_theme = TRUE))

  expect_doppelganger_ci("plot.statRes (combined volcano)", plot(imd_anova_res, plot_type = "volcano", bw_theme = TRUE))
  ## Test plot.totalcountFilt --------------------------------------------------
  
  seqfilt <- total_count_filter(omicsData = rnaseq_object)
  expect_doppelganger_ci("plot.totalCountFilt", plot(seqfilt, min_count = 5))
  
  ## Test plot.RNAFilt ---------------------------------------------------------
  
  seqfilt <- RNA_filter(omicsData = rnaseq_object)
  expect_doppelganger_ci("plot.RNAFilt", plot(seqfilt))
  
  ## Test plot.(omicsData_type) ------------------------------------------------
  
  myiso <- edata_transform(omicsData = isobaric_object, data_scale = "log2")
  expect_doppelganger_ci("plot.isobaricpepData", plot(myiso))
  
  mylipid <- edata_transform(omicsData = lipid_pos_object, data_scale = "log2")
  expect_doppelganger_ci("plot.lipidData", plot(mylipid, order_by = "Condition", color_by = "Condition"))
  
  mymetab <- edata_transform(omicsData = metab_object, data_scale = "log2")
  expect_doppelganger_ci("plot.metabData", plot(mymetab, order_by = "Condition", color_by = "Condition"))
  
  mynmr <- edata_transform(omicsData = nmr_identified_object, data_scale = "log2")
  expect_doppelganger_ci("plot.nmrData", plot(mynmr, order_by = "Condition", color_by = "Condition"))
  
  mypep <- edata_transform(omicsData = pep_object, data_scale = "log2")
  expect_doppelganger_ci("plot.pepData", plot(mypep, order_by = "Condition", color_by = "Condition"))
  
  mypro <- edata_transform(omicsData = pro_object, data_scale = "log2")
  expect_doppelganger_ci("plot.proData", plot(pro_object, order_by = "Condition", color_by = "Condition"))
  
  myseq <- group_designation(omicsData = rnaseq_object, main_effects = "Tissue")
  expect_doppelganger_ci("plot.seqData", plot(rnaseq_object, transformation = "lcpm"))
  
})