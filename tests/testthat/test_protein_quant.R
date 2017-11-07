#protein_quant testthat 
library(testthat)
library(pmartRqc)
library(pmartRdata)

data("pep_object")
data("lipid_object")
mypepData = edata_transform(pep_object, "log2")

edata_cname = attr(mypepData, "cnames")$edata_cname
#for some reason attr(pep_object, "cnames")$emeta_cname is NULL
attr(mypepData, "cnames")$emeta_cname<- "Protein"

#zrollup_result = protein_quant(mypepData, method = 'zrollup', combine_fn = 'median', single_pep = TRUE, single_observation = TRUE)

vector = c(1,2,3)

#hardcoding results for zrollup of peptides that map to a single protein
#first we subset the peptides that map to the chosen protein
unique_proteins = as.character(unique(mypepData$e_meta$Protein))
single_protein = unique_proteins[[2]]
single_res_zrollup = as.matrix(zrollup_result$e_data[which(zrollup_result$e_data$Protein == single_protein), -1])

#now lets see what mass_tag_ids are mapped to single_protein in emeta
inds = which(mypepData$e_meta$Protein == single_protein)
mass_tags = mypepData$e_meta$Mass_Tag_ID[inds]

pep_inds = which(mypepData$e_data$Mass_Tag_ID %in% mass_tags)

pep_subset = mypepData$e_data[pep_inds,]
pep_subset = pep_subset[,-which(names(pep_subset) == edata_cname)]

pep_mat = as.matrix(pep_subset)
num_peps = nrow(pep_mat)

#now we begin the zrollup 
mds<- vector(mode = "numeric", length = num_peps)
for(i in 1:num_peps){
  temp = median(pep_mat[i,], na.rm = TRUE)
  mds[i]<- temp
}

sds<- vector(mode = "numeric", length = num_peps)
for(i in 1:num_peps){
  temp = sd(pep_mat[i,], na.rm = TRUE)
  sds[i]<- temp
}

medians_mat = matrix(mds, nrow = num_peps, ncol = ncol(pep_mat), byrow = F)
standiv_mat = matrix(sds, nrow = num_peps, ncol = ncol(pep_mat), byrow = F)

peptides_scaled = (pep_mat - medians_mat)/standiv_mat

res = matrix(NA, nrow = 1, ncol =  ncol(pep_mat))
for(i in 1:ncol(peptides_scaled)){
  res[1,i]<- median(peptides_scaled[,i], na.rm = TRUE)
}

hardcoded_zrollup = res

context("output tests for protein_quant(method = 'zrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(zrollup_result), equals("proData"))
  expect_that(length(zrollup_result), equals(3))
  expect_that(length(attr(zrollup_result, "filters")), equals(2))
})

context("input tests for protein_quant(method = 'zrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', single_pep = FALSE, single_observation = FALSE), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', isoformRes = vector, single_pep = TRUE, single_observation = TRUE, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'zrollup', single_pep = TRUE, single_observation = TRUE, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'zrollup', single_pep = TRUE, single_observation = TRUE), throws_error())
  })


context("tests using hard coded results of protein_quant(method = 'zrollup')")

test_that("output of protein_quant() matches zrollup hard coded results",{
  expect_that(as.vector(single_res_zrollup), equals(as.vector(hardcoded_zrollup)))
  expect_that(ncol(single_res_zrollup), equals(ncol(hardcoded_zrollup)))
})

##################### tests for protein_quant(method = 'qrollup') #########################
#qrollup_result = protein_quant(mypepData, method = 'qrollup', qrollup_thresh = 2, combine_fn = 'mean')
qrollup_thresh = 2

#hardcoding results for qrollup of peptides that map to a single protein
#first we subset the peptides that map to the chosen protein
unique_proteins = as.character(unique(mypepData$e_meta$Protein))
single_protein = unique_proteins[[2]]
single_res_qrollup = as.matrix(qrollup_result$e_data[which(qrollup_result$e_data$Protein == single_protein), -1])

#now lets see what mass_tag_ids are mapped to single_protein in emeta
inds = which(mypepData$e_meta$Protein == single_protein)
mass_tags = mypepData$e_meta$Mass_Tag_ID[inds]

pep_inds = which(mypepData$e_data$Mass_Tag_ID %in% mass_tags)

pep_subset = mypepData$e_data[pep_inds,]
pep_subset = pep_subset[,-which(names(pep_subset) == edata_cname)]

pep_mat = as.matrix(pep_subset)
num_peps = nrow(pep_mat)

#Subset peptides whose abundance is >= to qrollup_thresh 
non_na_cnt<- vector(mode = "numeric", length = num_peps)
for(i in 1:num_peps){
  count = sum(!is.na(pep_mat[i,]))
  non_na_cnt[i]<- count
}

new_subset = pep_mat[which(non_na_cnt >= 2), ]

#Set protein abundance as the mean of peptide abundances from qrollup_thresh_subset
res = matrix(NA, nrow = 1, ncol =  ncol(new_subset))
for(i in 1:ncol(new_subset)){
  res[1,i]<- mean(pep_mat[,i], na.rm = TRUE)
}

hardcoded_qrollup = res


context("output tests for protein_quant(method = 'qrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(qrollup_result), equals("proData"))
  expect_that(length(qrollup_result), equals(3))
})

context("input tests for protein_quant(method = 'qrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = "abcdef", combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = 2, isoformRes = vector, combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', qrollup_thresh = 2, combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'qrollup', combine_fn = "mean"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'qrollup', qrollup_thresh = 2, combine_fn = "mean"), throws_error())
})


context("tests using hard coded results of protein_quant(method = 'qrollup')")

test_that("output of protein_quant() matches qrollup hard coded results",{
  expect_that(as.vector(single_res_qrollup), equals(as.vector(hardcoded_qrollup)))
  expect_that(ncol(single_res_qrollup), equals(ncol(hardcoded_qrollup)))
})

##################### tests for protein_quant(method = 'qrollup') #########################
#rrollup_result = protein_quant(mypepData, method = 'rrollup', combine_fn = 'median')

#hardcoding results for qrollup of peptides that map to a single protein
#first we subset the peptides that map to the chosen protein
unique_proteins = as.character(unique(mypepData$e_meta$Protein))
single_protein = unique_proteins[[2]]
single_res_rrollup = as.matrix(rrollup_result$e_data[which(rrollup_result$e_data$Protein == single_protein), -1])

#now lets see what mass_tag_ids are mapped to single_protein in emeta
inds = which(mypepData$e_meta$Protein == single_protein)
mass_tags = mypepData$e_meta$Mass_Tag_ID[inds]

pep_inds = which(mypepData$e_data$Mass_Tag_ID %in% mass_tags)

pep_subset = mypepData$e_data[pep_inds,]
pep_subset = pep_subset[,-which(names(pep_subset) == edata_cname)]

pep_mat = as.matrix(pep_subset)
num_peps = nrow(pep_mat)

#Select Reference Peptide -- peptide with least amount of missing data 
na_cnt<- vector(mode = "numeric", length = num_peps)
for(i in 1:num_peps){
  count = sum(is.na(pep_mat[i,]))
  na_cnt[i]<- count
}

least_na = which(na_cnt == min(na_cnt))

if(length(least_na)>1){
  mds<- vector(mode = "numeric", length = length(least_na))
  for(i in least_na){
    temp = median(pep_mat[i,], na.rm = TRUE)
    mds[i]<- temp
  }

  least_na = least_na[which(mds==max(mds))]		
}

prot_val = as.vector(pep_mat[least_na,])
prot_mat = matrix(prot_val, nrow = num_peps, ncol = ncol(pep_mat), byrow=T)
temp_mat = prot_mat - pep_mat

scaling_factor<- vector(mode = "numeric", length = num_peps)
for(i in 1:num_peps){
  temp = median (temp_mat[i,], na.rm = TRUE )
  scaling_factor[i]<- temp
}

x_scaled = pep_mat + matrix(scaling_factor, nrow = num_peps, ncol = ncol(pep_mat))

protein_val<- vector(mode = "numeric", length = ncol(x_scaled))
for(i in 1:ncol(x_scaled)){
  temp = median(x_scaled[,i], na.rm = TRUE)
  protein_val[i]<- temp
}

hardcoded_rrollup = protein_val


context("output tests for protein_quant(method = 'rrollup')")

test_that("results are of appropriate class and length",{ 
  expect_that(class(rrollup_result), equals("proData"))
  expect_that(length(rrollup_result), equals(3))
})

context("input tests for protein_quant(method = 'rrollup')")

test_that("incorrect argument input throws error",{  
  expect_that(protein_quant(pepData = mypepData, method = 'rrollup', isoformRes = vector, combine_fn = "median"), throws_error())
  expect_that(protein_quant(pepData = mypepData, method = 'rrollup', combine_fn = "abcdef"), throws_error())
  expect_that(protein_quant(pepData = lipid_object, method = 'rrollup', combine_fn = "median"), throws_error())
})


context("tests using hard coded results of protein_quant(method = 'rrollup')")

test_that("output of protein_quant() matches rrollup hard coded results",{
  expect_that(as.vector(single_res_rrollup), equals(as.vector(hardcoded_rrollup)))
  expect_that(length(as.vector(single_res_rrollup)), equals(length(hardcoded_rrollup)))
})
