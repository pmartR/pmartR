#bpquant for pmartRqc

bpquant<- function(statRes, pepData, pi_not = .9, max_proteoforms = 5){
  
  #some checks
  if(pi_not < 0 | pi_not > 1) stop("pi_not must be between 0 and 1")
  if(max_proteoforms <= 0) stop("max_proteoforms must be 1 or greater")
  
  if(class(statRes)[1] != "statRes") stop("statRes must be an object of class statRes")
  if(class(pepData)[1] != "pepData") stop("pepData must be an object of class pepData")
  
  
  #pulling e_meta from pepData, aswell as cnames
  e_meta<- pepData$e_meta
  edata_cname<- attr(pepData, "cnames")$edata_cname
  emeta_cname<- attr(pepData, "cnames")$emeta_cname
  
  #pulling signatures from statRes object
  signatures<- statRes$Flags
  
  #check that rows in signatures equals rows in e_data
  if(nrow(signatures) != nrow(pepData$e_data)) stop("rows in signatures must equal rows in e_data")

  #changing entries in signatures to 1 or -1
  for(j in 2:ncol(signatures)){
    
    if(!is.null(inds<- which(signatures[, j] > 1))){
      signatures[inds, j] = 1
    }
    
    if(!is.null(inds2<- which(signatures[, j] < -1))){
      signatures[inds2, j] = -1
    }
  }
  
  #merging e_meta with signatures from statRes
  protein_sig_data<- data.table:::merge.data.table(e_meta, signatures, by = edata_cname, all.x = TRUE)
  
  protein_sig_data<- as.data.frame(protein_sig_data)

  #removing any signatures with NA's
  protein_sig_data<- na.omit(protein_sig_data)

 #pull protein column from protein_sig_data and apply unique function
  unique_proteins<- unique(protein_sig_data[[emeta_cname]])

  res<- list()
  
  library(doParallel)
  cores<- detectCores()
  cl<- makeCluster(cores)
  registerDoParallel(cl)
  
 isoformRes<- foreach(i=1:length(unique_proteins)) %dopar%{
    
    row_ind<- which(protein_sig_data[, emeta_cname] == unique_proteins[i])
    cur_protein<- protein_sig_data[row_ind, ]
    cur_protein_sigs = cur_protein[, -(1:2)]
       
    result<- bpquant_mod(protein_sig = cur_protein_sigs, pi_not = pi_not, max_proteoforms = max_proteoforms)
    peptide_id<- result$peptide_idx
    
    temp<- data.frame(Protein = as.character(cur_protein[, emeta_cname]), Mass_Tag_ID = as.character(cur_protein[, edata_cname]), proteoformID = peptide_id, stringsAsFactors = FALSE)
    names(temp)<- c(emeta_cname, edata_cname, "proteoformID")
    
    res[[i]]<- temp
    
 }

 stopCluster(cl)
 
 bound_result<- do.call(rbind, isoformRes)
 zero_inds<- which(bound_result$proteoformID == 0)
 bound_result<- bound_result[-zero_inds, ]
 attr(isoformRes, "isoformRes_subset")<- bound_result
 
 class(isoformRes)<- "isoformRes"
 return(isoformRes)
 
#do.call(rbind, r)
  
}
