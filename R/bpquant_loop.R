#' Applies bpquant_mod function to each unique protein in pepData object
#'
#' The bpquant_loop function takes in a statRes object and a pepData object.It organizes the molecule and protein data into a data frame called protein_sig_data. Next bpquant_mod is applied to each unique protein and the results are stored and returned.  
#'
#' @param statRes an object of the class 'statRes'
#' @param pepData is an omicsData object of the class 'pepData'
#' @param pi_not is a numeric value between 0 and 1 indicating the background probability/frequency of a zero signature.
#' @param max_proteoforms a numeric value, a maximum threshold for the number of possible proteoforms.
#' 
#' @return a list of data frames, one for each unique protein. The data frames have three columns, "Protein", "Mass_Tag_ID" (which is molecule ID), and "ProteoformID". The class of this list is 'isoformRes'. 
#' 
#' @details The statRes object contains the signatures data frame, the pepData object is used for its e_meta data frame. Next the signatures data frame and e_meta are merged by their edata_cname (“Peptide”) columns, this new data frame called protein_sig_data will be input to bpquant_mod in a “foreach” statement. “Foreach” will subset protein_sig_data for each unique protein and apply bpquant_mod to each subset and store the results. 
#' 
#' @examples 
#' dontrun{
#' library(pmartR)
#' library(pmartRdata)
#' 
#' mypepData <- group_designation(omicsData = pep_object, main_effects = c("Condition"))
#' mypepData = edata_transform(mypepData, "log2")
#' 
#' imdanova_Filt <- imdanova_filter(omicsData = mypepData)
#' mypepData <- applyFilt(filter_object = imdanova_Filt, omicsData = mypepData, min_nonmiss_anova=2)
#' 
#' imd_anova_res <- imd_anova(omicsData = mypepData, test_method = 'comb', pval_adjust='bon')
#' 
#' result = bpquant_loop(statRes = imd_anova_res, pepData = mypepData)
#' 
#' }
#' 
#' @rdname bpquant_loop
#' @export

bpquant_loop<- function(statRes, pepData, pi_not = .9, max_proteoforms = 5){
  
  #some checks
  if(pi_not < 0 | pi_not > 1) stop("pi_not must be between 0 and 1")
  if(max_proteoforms <= 0) stop("max_proteoforms must be 1 or greater")
  
  if(!inherits(statRes, "statRes")) stop("statRes must be an object of class statRes")
  if(!inherits(pepData, "pepData")) stop("pepData must be an object of class pepData")
  
  
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
  protein_sig_data<- merge(e_meta, signatures, by = edata_cname, all.x = TRUE)
  
  protein_sig_data<- as.data.frame(protein_sig_data)

  #removing any signatures with NA's
  protein_sig_data<- na.omit(protein_sig_data)

 #pull protein column from protein_sig_data and apply unique function
  unique_proteins<- unique(protein_sig_data[[emeta_cname]])

  res<- list()
  
  suppressMessages(suppressPackageStartupMessages({
    library(parallel)
  })
  )
  cores<- detectCores()
  cl<- makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
 isoformRes<- foreach::foreach(i=1:length(unique_proteins))%dopar%{
    
    row_ind<- which(protein_sig_data[, emeta_cname] == unique_proteins[i])
    cur_protein<- protein_sig_data[row_ind, ]
    cur_protein_sigs = cur_protein[which(names(cur_protein) %in% colnames(signatures[-which(colnames(signatures) == edata_cname)]))]
    colnames(cur_protein_sigs) = paste(rep("flags", ncol(cur_protein_sigs)), 1:ncol(cur_protein_sigs), sep = "")
       
    result<- pmartR::bpquant_mod(protein_sig = cur_protein_sigs, pi_not = pi_not, max_proteoforms = max_proteoforms)
    peptide_id<- result$peptide_idx
    
    temp<- data.frame(Protein = as.character(cur_protein[, emeta_cname]), Mass_Tag_ID = as.character(cur_protein[, edata_cname]), proteoformID = peptide_id, stringsAsFactors = FALSE)
    names(temp)<- c(emeta_cname, edata_cname, "proteoformID")
    
    res[[i]]<- temp
    
 }

 stopCluster(cl)
 
 #filter out isoformRes with max(proteoformID) of zero
 list_inds = lapply(isoformRes,function(x){max(x$proteoformID)})
 inds = which(list_inds == 0)
 isoformRes = isoformRes[-inds]
 
 isoformRes2<- lapply(isoformRes, isoformRes_func, emeta_cname = emeta_cname, edata_cname = edata_cname)
 isoformRes2<- do.call(rbind, isoformRes2) 
 
 attr(isoformRes, "isoformRes_subset")<- isoformRes2
 
 class(isoformRes)<- "isoformRes"
 return(isoformRes)
}

# function to use with lapply on isoformRes
isoformRes_func<- function(df, emeta_cname, edata_cname){
  temp<- vector("list", max(df$proteoformID))
  
  for(i in 1:max(df$proteoformID)){
    cur_subset<- df[which(df$proteoformID == i), ]
    new_df<- data.frame(cur_subset[[emeta_cname]], paste(cur_subset[[emeta_cname]], cur_subset$proteoformID, sep = ';'), as.numeric(cur_subset[[edata_cname]]))
    names(new_df)<- c(names(df)[1], "Protein_Isoform", names(df)[2])
    temp[[i]]<- new_df
  }
  
  do.call(rbind, temp)
}


