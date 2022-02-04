source(system.file('testdata', 'load_data.R', package = 'pmartR'))

f_meta <- data.frame(
  "Proteins" = c(paste0("Mock", 1:3), paste0("Infection", c(1:7)), NA,  "Infection9"),
  "Lipids" = c(paste0("Mock", 1:3), paste0("Infection", c(1:4)), NA, paste0("Infection", c(6:9))),
  "Metabolites" = c(paste0("Mock", 1:3), paste0("Infection", c(1:9))),
  "Condition" = c(rep("A", 3), rep("B", 9))
)

# fmeta where the Lipids column is bad
bad_fmeta1 <- data.frame(
  "Proteins" = c(paste0("Mock", 1:3), paste0("Infection", c(1:7)), NA,  "Infection9"),
  "Lipids" = c(paste0("Mock", 1:3), paste0("Infection", c(1:3)), NA, NA, paste0("Infection", c(6:9))),
  "Metabolites" = c(paste0("Mock", 1:3), paste0("Infection", c(1:9))),
  "Condition" = c(rep("A", 3), rep("B", 9))
)

lipid_object_proc <- edata_transform(ldata, "log2")
lipid_object_proc <- normalize_global(lipid_object_proc, "all", "median", apply_norm = T)

metab_object_proc <- edata_transform(mdata, "log2")
metab_object_proc <- normalize_global(lipid_object_proc, "all", "median", apply_norm = T)

pro_object_proc <- edata_transform(prdata, "log2")
pro_object_proc <- normalize_global(pro_object_proc, "all", "median", apply_norm = T)

pro_grouped = group_designation(pro_object_proc, main_effects = "Condition")
metab_grouped = group_designation(metab_object_proc, main_effects = "Condition")
lipid_grouped = group_designation(lipid_object_proc, main_effects = "Condition")

mdata_man_fmeta1 = as.multiData(lipid_object_proc, metab_object_proc, pro_object_proc, f_meta = f_meta)
mdata_man_fmeta2 = as.multiData(metab_object_proc, pro_object_proc, lipid_object_proc, f_meta = f_meta)
mdata_auto_fmeta = as.multiData(metab_object_proc, pro_object_proc, lipid_object_proc, auto_fmeta = T)
mdata_auto_fmeta_sinter = as.multiData(metab_object_proc, lipid_object_proc, pro_object_proc, auto_fmeta = T, sample_intersect = T)
mdata_auto_fmeta_noarr = as.multiData(metab_object_proc, lipid_object_proc, pro_object_proc, auto_fmeta = T, match_samples = F)
mdata_noarr_sint = as.multiData(metab_object_proc, lipid_object_proc, pro_object_proc, auto_fmeta = T, match_samples = F, sample_intersect = T)
mdata_grouped = as.multiData(metab_grouped, lipid_grouped, pro_grouped, auto_fmeta = T, sample_intersect = T, keep_sample_info = T)

obj_list = list(mdata_man_fmeta1 , mdata_man_fmeta2, mdata_auto_fmeta, mdata_auto_fmeta_sinter, mdata_grouped)

test_that("Bad input throws error", {
  # objects with different log2/normalization
  expect_error(as.multiData(ldata, prdata))
  expect_error(as.multiData(edata_transform(ldata, "log2"), prdata))
  expect_error(as.multiData(normalize_global(mdata, "all", "median", apply_norm = T), pro_object_proc))
  
  # not both grouped
  expect_error(as.multiData(group_designation(lipid_object_proc, main_effects = "Condition"), pro_object_proc))
  
  # f_meta not provided and f_datas are not valid f_metas
  expect_error(as.multiData(lipid_object_proc, pro_object_proc))
  expect_error(as.multiData(metab_object_proc, pro_object_proc, lipid_object_proc))
  
  # isobaric object is not reference normalized
  expect_error(as.multiData(ldata, isodata))
})

test_that("multiData attributes are aligned with <object>$omicsData", {
  for(obj in obj_list) {
    matches = pmartR:::fmeta_matches(obj$omicsData, obj$f_meta)
    for(i in 1:length(matches)) {
      fdata_cname = get_fdata_cname(obj$omicsData[[i]])
      # The i-th fmeta_samp_cname should be a match for the i-th object's sample ID's
      expect_true(attributes(obj)$fmeta_samp_cname[i] %in% matches[[i]])
    }
  }
})

test_that("multiData options produce expected output", {
  na_inds_counts <- function(x) {
    na_rows = which(rowSums(is.na(x$f_meta[,attr(x, "fmeta_samp_cname")])) > 0) 
    row_unique_counts = sapply(1:(nrow(x$f_meta) - length(na_rows)), function(i) {
      length(unique(unlist(x$f_meta[i,]))) == 1
    })
    return(list("nas" = na_rows, "runq" = row_unique_counts))
  }
  # auto fmeta should match samples by default, others should have NA's
  res <- na_inds_counts(mdata_auto_fmeta)
  expect_true(length(res[["nas"]]) > 0)
  expect_true(all(res[["runq"]]))
  
  # taking the sample intersect should remove samples Infection 8 and Infection 5
  expect_true(nrow(mdata_auto_fmeta_sinter$f_meta) < nrow(mdata_auto_fmeta$f_meta))
  res <- na_inds_counts(mdata_auto_fmeta_sinter)

  expect_true(length(res[["nas"]]) == 0)
  expect_true(all(res[["runq"]]))
  
  # if match_samples = False, then there should be mis-aligned f_meta
  res <- na_inds_counts(mdata_auto_fmeta_noarr)
  expect_true(length(res[["nas"]]) == 0)
  expect_true(!all(res[["runq"]]))
  
  # sample intersect should remove bad alignment resulting from match_samples = False
  res <- na_inds_counts(mdata_noarr_sint)
  expect_true(length(res[["nas"]]) == 0)
  expect_true(all(res[["runq"]]))
  
  # TODO:  Group assignments are aligned
})