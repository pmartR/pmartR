context("data loads and objects created successfully")
library(testthat)
library(pmartRdata)
library(pmartR)

# instantiate all objects...
pep_edata <- pmartRdata::pep_edata
pep_fdata <- pmartRdata::pep_fdata
pep_emeta <- pmartRdata::pep_emeta
lipid_edata <- pmartRdata::lipid_edata
lipid_fdata <- pmartRdata::lipid_fdata
metab_edata <- pmartRdata::metab_edata
metab_fdata <- pmartRdata::metab_fdata
pro_edata <- pmartRdata::pro_edata
pro_fdata <- pmartRdata::pro_fdata
techrep_edata <- pmartRdata::techrep_edata
techrep_fdata <- pmartRdata::techrep_fdata

#....and store in a list to iterate
data_list <- list(pep_edata = pep_edata, pep_fdata = pep_fdata, pep_emeta = pep_emeta, lipid_edata = lipid_edata, lipid_fdata = lipid_fdata,
                  metab_edata = metab_edata, metab_fdata = metab_fdata, pro_edata = pro_edata, pro_fdata = pro_fdata, techrep_edata = techrep_edata, techrep_fdata = techrep_fdata)
iter_names <- c("pep", "pro", "lipid", "metab", "techrep")

#### Most of these are pretty simple and stable, but tests can easily be added to the loops in either block to test ALL objects. ####

test_that("bad input in data object creation throw error", {
  for(class in iter_names){
    edata_cname <- colnames(data_list[[paste0(class, "_edata")]])[1]
    fdata_cname <- colnames(data_list[[paste0(class, "_fdata")]])[1]
    emeta_cname <- if(class == "pep") colnames(data_list[[paste0(class, "_emeta")]])[1] else NULL
    techrep_cname = if(class == "techrep") "TECH_REP" else NULL
    e_meta = if(class == "pep") data_list[[paste0(class, "_emeta")]] else NULL
    edata <- data_list[[paste0(class, "_edata")]]
    fdata <- data_list[[paste0(class, "_fdata")]]
    
    fn = switch(class, "pep" = as.pepData, "pro" = as.proData, "lipid" = as.lipidData, "metab" = as.metabData, "techrep" = as.pepData)
    # expect error when 1 column is randomly added to edata or 1 row is randomly dropped from fdata
    expect_error(fn(e_data = cbind(edata, edata[,sample(2:ncol(edata))]), f_data = fdata, e_meta = e_meta,
                    edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = emeta_cname, techrep_cname = techrep_cname))
    expect_error(fn(e_data = edata, f_data = fdata[-sample(1:nrow(fdata), 1),], e_meta = e_meta,
                    edata_cname = edata_cname, fdata_cname = fdata_cname, emeta_cname = emeta_cname, techrep_cname = techrep_cname))
    
    # bad column names
    expect_error(fn(e_data = edata, f_data = fdata, e_meta = e_meta,
                    edata_cname = edata_cname, fdata_cname = "badname", emeta_cname = emeta_cname, techrep_cname = techrep_cname))
    expect_error(fn(e_data = edata, f_data = f_data, e_meta = e_meta,
                    edata_cname = "badname", fdata_cname = fdata_cname, emeta_cname = emeta_cname, techrep_cname = techrep_cname))
  }  
})

# correct object creation calls
pep <- as.pepData(e_data = pep_edata, f_data = pep_fdata, e_meta = pep_emeta, edata_cname = "Mass_Tag_ID", fdata_cname = "SampleID", emeta_cname = "Mass_Tag_ID")
pro <- as.proData(e_data = pro_edata, f_data = pro_fdata, edata_cname = "Reference", fdata_cname = "SampleID")
lipid <- as.lipidData(e_data = lipid_edata, f_data = lipid_fdata, edata_cname = "LipidCommonName", fdata_cname = "Sample_Name")
metab <- as.metabData(e_data = metab_edata, f_data = metab_fdata, edata_cname = "Metabolite" , fdata_cname = "SampleID")

object_list = list(pep, pro, lipid, metab)

test_that("integrity of as.omicsData output", {
  for(el in object_list){
    expect_equal(colnames(el$e_data[-which(colnames(el$e_data) == get_edata_cname(el))]), as.character(el$f_data[,get_fdata_cname(el)]))
    expect_true(inherits(attr(el, "data_info")$norm_info, "list"))
    expect_true(inherits(attr(el, "filters"), "list"))
    expect_true(attr(el, "data_info")$norm_info$is_normalized == FALSE)
  }
})
