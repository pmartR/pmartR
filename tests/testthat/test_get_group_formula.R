context('class: seqData')


### TO do: add emeta

test_that('get_group_formula returns the correct data frame and results', {
  # Load the reduced peptide data frames ---------------------------------------

  load(system.file('testdata',
    'little_seqdata.RData',
    package = 'pmartR'
  ))

  # Run as.seqData with agreeable data frames ----------------------------------

  # Construct a seqData object with the edata, fdata, and emeta data frames.
  seqdata <- as.seqData(
    e_data = edata,
    f_data = fdata,
    edata_cname = 'ID_REF',
    fdata_cname = 'Samples'
  )


  expect_error(get_group_formula(seqdata), "group_designation has not been run")

  ## Always factor levels
  seqdata$f_data$Pair <- paste0("Pair_", rep(1:20, 2))
  seqdata$f_data$Pair_group <- paste0("Potato_", c(rep(1, 20), rep(2, 20)))

  ## Alt with numeric main effect
  seqdata$f_data$group2 <- rnorm(40)

  seqdata$f_data$covar <- rep(rnorm(20), 2) ## numeric, redundant with pair id
  seqdata$f_data$covar2 <- rep(c(1, 1, 1, 1, rnorm(16)), 2) ## numeric, not redundant with pair id
  seqdata$f_data$covar4 <- as.numeric(as.factor(seqdata$f_data$Treatment)) ## factor, redundant with group


  ## Factor main effects ## works
  seqdata_grp1 <- group_designation(seqdata, main_effects = c("Tissue", "Treatment"))

  ## numeric main effects ## -> this is automatically read as a factor level
  seqdata_grp2 <- suppressWarnings(group_designation(seqdata, main_effects = c("group2")))

  ## Factor main effects and covar
  seqdata_grp3 <- group_designation(seqdata,
    main_effects = c("Tissue", "Treatment"),
    covariates = "covar"
  )

  seqdata_grp4 <- group_designation(seqdata,
    main_effects = c("Treatment"),
    covariates = "covar4"
  )

  ## Treatment and paired groupp ## works
  seqdata_grp5 <- group_designation(seqdata,
    main_effects = c("Treatment"),
    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1"
  )

  ## Test rejection of redundant -- not considered redundant since pair_group is changed for processing
  seqdata_grp6 <- group_designation(seqdata,
    main_effects = c("Treatment"), covariates = "covar",
    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1"
  )

  ## Non-redundant ## works
  seqdata_grp7 <- group_designation(seqdata,
    main_effects = c("Treatment"), covariates = "covar2",
    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1"
  )

  ## pairs only
  seqdata_grp8 <- group_designation(seqdata,
    pair_id = "Pair", pair_group = "Pair_group", pair_denom = "Potato_1"
  )

  group_list <- list(
    seqdata_grp1, ## factor main effects only
    seqdata_grp2, ## numeric main effects -- valid?
    seqdata_grp3, ## factor main effects w/ covariate redundant with pair id
    seqdata_grp4, ## factor main effects w/ covariate redundant with group
    seqdata_grp5, ## factor main effect w/ pair
    seqdata_grp6, ## factor main effect w/ pair and covariate
    seqdata_grp7, ## factor main effect w/ pair and numeric covariate
    seqdata_grp8 ## pairs only
  )

  expect_warning(
    get_group_formula(seqdata_grp4),
    "At least 1 detected covariate is confounded with Group"
  )

  res_get_group_formula <- suppressWarnings(purrr::map(group_list, get_group_formula))


  ## Expected dimensions df
  df_close_list <- list(
    get_group_DF(seqdata_grp1),
    get_group_DF(seqdata_grp2),
    cbind(get_group_DF(seqdata_grp3), seqdata_grp3$f_data[8]),
    get_group_DF(seqdata_grp4),
    cbind(get_group_DF(seqdata_grp5), seqdata_grp5$f_data[c("Pair", "Pair_group")]), ## refactored for levels
    cbind(get_group_DF(seqdata_grp6), seqdata_grp6$f_data[c("Pair", "Pair_group", "covar")]), ## refactored for levels
    cbind(get_group_DF(seqdata_grp7), seqdata_grp7$f_data[c("Pair", "Pair_group", "covar2")]),
    cbind(get_group_DF(seqdata_grp8), seqdata_grp8$f_data[c("Pair", "Pair_group")])
  )

  ## Compatible grouping for seqData stats
  Expected_form <- list(
    "~0 + Group", ## factor main effects only
    "~0 + Group", ## numeric main effects -- valid?
    "~0 + covar + Group", ## factor main effects w/ covariate
    "~0 + Group", ## factor main effects w/ covariate redundant with group
    "~0 + Pair  + Pair_group  + Group", # pair and main effect
    "~0 + Pair  + Pair_group  + covar + Group", # pair and main effect and covar
    "~0 + Pair  + Pair_group  + covar2 + Group", # pair and main effect and covar
    "~0 + Pair + Pair_group" ## pair only
  )

  ## formulas are identical
  expect_identical(purrr::map(res_get_group_formula, 2), Expected_form)

  dfs_only <- purrr::map(res_get_group_formula, 1)

  ## dimensions are identical
  expect_identical(lapply(dfs_only, dim), lapply(df_close_list, dim))


  ## Columns are identical
  res <- purrr::map2(dfs_only, df_close_list, function(df, df_close) {
    purrr::map(colnames(dfs_only), function(nm) {
      if (nm != "Pair") { ## pair is refactored
        expect_identical(dfs_only[[nm]], df_close[[nm]])
      } else {
        if ("Group" %in% colnames(df_close)) {
          grouping_info <- dplyr::arrange(
            df_close,
            Group,
            "Pair_group", "Pair"
          )

          all_mult_levels <- table(apply(grouping_info[c("Group", "Pair_group")],
            1, paste,
            collapse = ""
          ))

          grouping_info[["Pair"]] <- unlist(purrr::map(all_mult_levels, function(el) {
            paste0("new_pair_name", 1:el)
          }))

          expect_identical(dfs_only[[nm]], grouping_info[[nm]])
        } else expect_identical(dfs_only[[nm]], df_close[[nm]])
      }
    })
  })
})
