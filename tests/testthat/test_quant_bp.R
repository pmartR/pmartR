context('quantitation: Reverend Bayes')

test_that('bpquant produces correct isoforms', {
  # Load and preprocess data ---------------------------------------------------

  # Load peptide data.
  load(system.file('testdata',
    'little_pdata.RData',
    package = 'pmartR'
  ))

  # Construct a pepData object.
  pdata <- as.pepData(
    e_data = edata,
    f_data = fdata,
    e_meta = emeta,
    edata_cname = 'Mass_Tag_ID',
    fdata_cname = 'SampleID',
    emeta_cname = 'Protein'
  )

  # Log transmute the peptide data.
  pdata <- edata_transform(pdata, "log")

  # Add some groups.
  pdata <- group_designation(
    omicsData = pdata,
    main_effects = "Condition"
  )

  # Apply an IMD-ANOVA filter.
  ifilter <- imdanova_filter(omicsData = pdata)
  pdata <- applyFilt(
    filter_object = ifilter,
    omicsData = pdata,
    min_nonmiss_anova = 2
  )

  # Run some statisiticalness on the filtered data.
  inova <- imd_anova(
    omicsData = pdata,
    test_method = 'comb',
    pval_adjust_a_multcomp = 'bon',
    pval_adjust_g_multcomp = 'bon'
  )

  # Quantitate with Reverend Bayes ---------------------------------------------

  # Run bpquant pdata with the IMD-ANOVA output as the statRes object.
  bayes <- bpquant(statRes = inova, pepData = pdata, parallel = FALSE)

  # Vector of unique proteins from the reduced pdata object. Will be used to
  # test the isoformRes_subset attribute from the bpquant function.
  prots <- c(
    "ALBU_HUMAN", "ALBU_HUMAN", "ALBU_HUMAN", "ALBU_HUMAN",
    "RL40_HUMAN", "RL40_HUMAN", "CYC_HUMAN", "G3P_HUMAN", "G3P_HUMAN",
    "G3P_HUMAN", "G3P_HUMAN", "G3P_HUMAN", "G3P_HUMAN", "G3P_HUMAN",
    "G3P_HUMAN", "PYGB_HUMAN", "PYGB_HUMAN", "PYGB_HUMAN",
    "PYGB_HUMAN", "TRFE_HUMAN", "CATA_HUMAN", "HBB_HUMAN",
    "SUCA_HUMAN", "ATPA_HUMAN", "ATPA_HUMAN", "TRFL_HUMAN",
    "AACT_HUMAN", "AACT_HUMAN", "AACT_HUMAN", "AACT_HUMAN",
    "AACT_HUMAN", "HBA_HUMAN", "P20D2_HUMAN", "VTNC_HUMAN",
    "K1C19_HUMAN", "HSP72_HUMAN", "SYEP_HUMAN", "UBP7_HUMAN",
    "TMOD3_HUMAN", "PGAM1_HUMAN", "EFTU_HUMAN", "EIF3H_HUMAN",
    "CH60_HUMAN", "CH60_HUMAN", "CH60_HUMAN", "PDIA3_HUMAN",
    "PDIA3_HUMAN", "RA1L2_HUMAN", "RBBP7_HUMAN", "K2C7_HUMAN",
    "K2C7_HUMAN", "K2C7_HUMAN", "K1C18_HUMAN", "SFPQ_HUMAN",
    "ATPB_HUMAN", "RS17_HUMAN", "ANXA2_HUMAN", "ANXA2_HUMAN",
    "ANXA2_HUMAN", "PPIA_HUMAN", "PPIA_HUMAN", "K1C17_HUMAN",
    "XRCC5_HUMAN", "H2B1B_HUMAN", "KPYM_HUMAN", "KPYM_HUMAN",
    "HNRPM_HUMAN", "ACTN1_HUMAN", "ACTB_HUMAN", "ACTB_HUMAN",
    "LMNA_HUMAN", "PHB_HUMAN", "NONO_HUMAN", "ROA2_HUMAN",
    "ROA2_HUMAN", "GRP78_HUMAN", "GRP78_HUMAN", "GRP78_HUMAN",
    "EF1A1_HUMAN", "H90B3_HUMAN", "TRAP1_HUMAN", "GBLP_HUMAN",
    "IF4A1_HUMAN", "ROA0_HUMAN", "HNRPC_HUMAN", "AT1A1_HUMAN",
    "SND1_HUMAN", "HNRPU_HUMAN", "RS7_HUMAN", "H2A1_HUMAN",
    "COR1C_HUMAN", "NPM_HUMAN", "HNRPF_HUMAN", "RTCB_HUMAN",
    "HNRPK_HUMAN", "H2AY_HUMAN", "EF1G_HUMAN", "4F2_HUMAN"
  )

  # Vector of peptide IDs from the reduced pdata object. This will be used to
  # test the isoformRes_subset attribute.
  pepe_ids <- as.character(c(
    1104, 1237, 1768, 4198254, 1246, 1253, 1110, 1214,
    1406, 7601, 10343, 15747, 21470, 24307, 6948822,
    1649, 11055, 12404, 21149, 11956, 12894, 13878,
    66123, 216191, 6948877, 4207823, 6639885, 6679059,
    6709059, 6719683, 6915905, 6643916, 6701524,
    6706932, 6753571, 6948901, 6781312, 6793445,
    6854457, 6859822, 6880071, 6907124, 6948817,
    6948829, 6948884, 6948819, 6948849, 6948821,
    6948823, 6948824, 6948886, 6948902, 6948837,
    6948826, 6948827, 6948830, 6948831, 6948867,
    6948908, 6948832, 6948850, 6948833, 6948834,
    6948836, 6948839, 6948866, 6948894, 6948841,
    6948842, 6948882, 6948843, 6948844, 6948913,
    6948847, 6948863, 6948851, 6948891, 6948899,
    6948852, 6948854, 6948855, 6948856, 6948865,
    6948878, 6948880, 6948887, 6948895, 6948896,
    6948897, 6948898, 6948900, 6948903, 6948904,
    6948905, 6948906, 6948907, 6948911, 6948915
  ))

  # Sleuth around the isoformRes object.
  expect_s3_class(bayes, "isoformRes")
  expect_equal(length(bayes), 66)
  expect_equal(
    dim(attr(bayes, "isoformRes_subset")),
    c(98, 3)
  )
  expect_equal(
    bayes[[1]],
    data.frame(check.names = FALSE, 
      Protein = rep("ALBU_HUMAN", 5),
      Mass_Tag_ID = c(
        "1047", "1104", "1237", "1768",
        "4198254"
      ),
      proteoformID = c(0, 1, 1, 1, 1)
    )
  )
  expect_equal(
    bayes[[5]],
    data.frame(check.names = FALSE, 
      Protein = rep("PYGB_HUMAN", 5),
      Mass_Tag_ID = c(
        "1649", "11055", "11078", "12404",
        "21149"
      ),
      proteoformID = c(1, 1, 0, 1, 1)
    )
  )
  expect_equal(
    bayes[[22]],
    data.frame(check.names = FALSE, 
      Protein = "EFTU_HUMAN",
      Mass_Tag_ID = c("6880071"),
      proteoformID = c(1)
    )
  )
  expect_equal(
    bayes[[36]],
    data.frame(check.names = FALSE, 
      Protein = "XRCC5_HUMAN",
      Mass_Tag_ID = c("6948834"),
      proteoformID = c(1)
    )
  )
  expect_equal(
    bayes[[63]],
    data.frame(check.names = FALSE, 
      Protein = "HNRPK_HUMAN",
      Mass_Tag_ID = c("6948906"),
      proteoformID = c(1)
    )
  )
  expect_equal(
    attr(bayes, "isoformRes_subset")$Protein,
    prots
  )
  expect_equal(
    attr(bayes, "isoformRes_subset")$Protein_Isoform,
    prots
  )
  expect_equal(
    attr(bayes, "isoformRes_subset")$Mass_Tag_ID,
    pepe_ids
  )

  # Test internals of Rev Bayes quantitation -----------------------------------

  # Run bpquant_mod on the flags for ALBU_HUMAN.
  bayes_mod <- pmartR:::bpquant_mod(
    protein_sig = data.frame(check.names = FALSE, flags1 = c(-1, 1, 1, 1, 1)),
    pi_not = 0.9,
    max_proteoforms = 5
  )

  # Go over output with a fine-tooth comb.
  expect_equal(
    round(bayes_mod$post_prob, 5),
    round(
      c(4.827795e-04, 8.469562e-01, 8.691286e-05, 1.524741e-01),
      5
    )
  )
  expect_equal(
    bayes_mod$peptide_idx,
    c(0, 1, 1, 1, 1)
  )
  expect_equal(
    bayes_mod$unique_sigs$flags1,
    c(-1, 1)
  )
  expect_equal(
    bayes_mod$num_proteoforms,
    1
  )
  expect_equal(
    bayes_mod$proteoform_configs,
    matrix(
      c(
        0, 0,
        1, 0,
        0, 1,
        1, 1
      ),
      nrow = 4,
      byrow = TRUE
    )
  )

  # Run isoformRes_func on the output for ALBU_HUMAN.
  iso_fun <- pmartR:::isoformRes_func(
    df = data.frame(check.names = FALSE, 
      Protein = rep("ALBU_HUMAN", 4),
      Mass_Tag_ID = c("1104", "1237", "1768", "4198254"),
      proteoformID = c(1, 1, 1, 1)
    ),
    emeta_cname = "Protein",
    edata_cname = "Mass_Tag_ID"
  )

  # Sniff around the output for isoformRes_func.
  expect_identical(
    iso_fun,
    data.frame(check.names = FALSE, 
      Protein = rep("ALBU_HUMAN", 4),
      Protein_Isoform = rep("ALBU_HUMAN", 4),
      Mass_Tag_ID = c("1104", "1237", "1768", "4198254")
    )
  )

  # Run bpquant_mod on the flags for 6PGL_HUMAN.
  bayes_mod <- pmartR:::bpquant_mod(
    protein_sig = data.frame(check.names = FALSE, flags1 = c(-1, -1, -1, 1, 1)),
    pi_not = 0.9,
    max_proteoforms = 5
  )

  # Go over output with a fine-tooth comb.
  expect_equal(
    round(bayes_mod$post_prob, 5),
    round(
      c(0.006577728, 0.298581994, 0.014977305, 0.679862973),
      5
    )
  )
  expect_equal(
    bayes_mod$peptide_idx,
    c(1, 1, 1, 2, 2)
  )
  expect_equal(
    bayes_mod$unique_sigs$flags1,
    c(-1, 1)
  )
  expect_equal(
    bayes_mod$num_proteoforms,
    2
  )
  expect_equal(
    bayes_mod$proteoform_configs,
    matrix(
      c(
        0, 0,
        1, 0,
        0, 1,
        1, 1
      ),
      nrow = 4,
      byrow = TRUE
    )
  )

  # Run isoformRes_func on the output for 6PGL_HUMAN.
  iso_fun <- pmartR:::isoformRes_func(
    df = data.frame(check.names = FALSE, 
      Protein = rep("6PGL_HUMAN", 5),
      Mass_Tag_ID = c(
        "8622908", "8655070", "9513231", "34862026",
        "65465565"
      ),
      proteoformID = c(1, 1, 1, 2, 2)
    ),
    emeta_cname = "Protein",
    edata_cname = "Mass_Tag_ID"
  )

  # Sniff around the output for isoformRes_func.
  expect_identical(
    iso_fun,
    data.frame(check.names = FALSE, 
      Protein = rep("6PGL_HUMAN", 5),
      Protein_Isoform = c(
        "6PGL_HUMAN;1", "6PGL_HUMAN;1",
        "6PGL_HUMAN;1", "6PGL_HUMAN;2",
        "6PGL_HUMAN;2"
      ),
      Mass_Tag_ID = c(
        "8622908", "8655070", "9513231", "34862026",
        "65465565"
      )
    )
  )
})
