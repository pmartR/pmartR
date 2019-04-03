#' Applies zrollup function 
#' 
#' This function applies the zrollup method to a pepData object for each unique protein and returns a proData object. 
#' 
#' @param pepData an omicsData object of class 'pepData'
#' @param combine_fn logical indicating what combine_fn to use, defaults to median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in applying zrollup function. Defaults to TRUE.
#' 
#' @return an omicsData object of class 'proData'
#' 
#' @details In the zrollup method, peptides are scaled as, pep_scaled = (pep - median)/sd, and protein abundance is set as the mean of these scaled peptides.  
#'
#' @references Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins, J. N., Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical tool for quantitative analysis of -omics data}. Bioinformatics (Oxford, England), 24(13), 1556-1558. 
#'  
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = zrollup(pepData = pep_object) 
#'}
#' 
#' @rdname zrollup
#' 

zrollup<- function(pepData, combine_fn = "median", parallel = TRUE){
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  
  #check that combine_fn is one of 'mean', 'median'
  if(!(combine_fn %in% c('median', 'mean'))) stop("combine_fn has to be one of 'mean' or 'median'")
  
  pep_id = attr(pepData, "cnames")$edata_cname
  pro_id = attr(pepData, "cnames")$emeta_cname
  
  pep = data.table::data.table(pepData$e_data)
  pro = data.table::data.table(pepData$e_meta[,c(pep_id, pro_id)])
  temp = merge(x = pro, y = pep, by = pep_id, all.x = F, all.y = T)
  temp = as.data.frame(temp, check.names=check_names)[,-which(names(temp)==pep_id)]
  
  #pull protein column from temp and apply unique function
  unique_proteins<- unique(temp[[pro_id]])
  
  #assigning function to chosen_combine_fn 
  if(combine_fn == "median"){
    chosen_combine_fn<- combine_fn_median
  }else{chosen_combine_fn = combine_fn_mean}
  
  if(parallel == TRUE){
    
    final_list<- vector("list", length(unique_proteins))
    
  
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    r<-foreach::foreach(i=1:length(unique_proteins))%dopar%{
      
      row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
      current_subset<- temp[row_ind,]
      current_subset<- current_subset[,-which(names(temp) == pro_id)]
      
      #### Perform Z_Rollup ####
      ## store number of peptides ##
      num_peps = nrow(current_subset)
      
      res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
      ## Step 1: Compute mean and sd of peptides ##
      mds = apply(current_subset, 1, median, na.rm = T)
      sds = apply(current_subset, 1, sd, na.rm = T)
      
      ## Step 2: Scale peptide data as pep_scaled = (pep - median)/sd
      medians_mat = matrix(mds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)
      standiv_mat = matrix(sds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)
      
      peptides_scaled = apply((current_subset - medians_mat)/standiv_mat, 2, chosen_combine_fn)
      
      res[1,] = peptides_scaled
      res<- data.frame(res)
      names(res)<- names(current_subset)
      
      final_list[[i]]<- res
    }
    parallel::stopCluster(cl)
    
    final_result<- do.call(rbind, r)
    final_result<- cbind(unique_proteins, final_result)
    names(final_result)[1]<-pro_id
    
    samp_id = attr(pepData, "cnames")$fdata_cname
    data_scale = attr(pepData, "data_info")$data_scale
    is_normalized = attr(pepData, "data_info")$norm_info$is_normalized
    
    #subsetting pepData$e_meta by 'unique_proteins' 
    emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
    
    if(ncol(pepData$e_meta) == 2){
      e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
      names(e_meta)<-pro_id
    }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
    
    prodata = as.proData(e_data = data.frame(final_result, check.names=check_names), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, is_normalized = is_normalized, check.names = check_names)
    
    #check for isobaricpepData class
    if(inherits(pepData, "isobaricpepData")){
      #update attributes in prodata
      attr(prodata, "isobaric_info") = attr(pepData, "isobaric_info")
      attr(prodata, "isobaric_info")$norm_info$is_normalized = attr(pepData, "isobaric_info")$norm_info$is_normalized
    }
    
    #updating prodata attributes
    attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
    attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
    attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
    
    attr(prodata, "filters")<- attr(pepData, "filters")
    attr(prodata, "group_DF")<- attr(pepData, "group_DF")
    attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  }
  
  #applying zrollup without doParallel
  else{
    final_list<- vector("list", length(unique_proteins))
    
    for(i in 1:length(unique_proteins)){
      
      row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
      current_subset<- temp[row_ind,]
      current_subset<- current_subset[,-which(names(temp) == pro_id)]
      
      #### Perform Z_Rollup ####
      ## store number of peptides ##
      num_peps = nrow(current_subset)
      
      res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
      
      ## Step 1: Compute mean and sd of peptides ##
      mds = apply(current_subset, 1, median, na.rm = T)
      sds = apply(current_subset, 1, sd, na.rm = T)
      
      ## Step 2: Scale peptide data as pep_scaled = (pep - median)/sd
      medians_mat = matrix(mds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)
      standiv_mat = matrix(sds, nrow = num_peps, ncol = ncol(current_subset), byrow = F)
      
      peptides_scaled = apply((current_subset - medians_mat)/standiv_mat, 2, chosen_combine_fn)
      
      res[1,] = peptides_scaled
      res<- data.frame(res)
      names(res)<- names(current_subset)
      
      final_list[[i]]<- res
    }
    
    final_result<- do.call(rbind, final_list)
    final_result<- cbind(unique_proteins, final_result)
    names(final_result)[1]<-pro_id
    
    samp_id = attr(pepData, "cnames")$fdata_cname
    data_scale = attr(pepData, "data_info")$data_scale
    is_normalized = attr(pepData, "data_info")$norm_info$is_normalized
    
    #subsetting pepData$e_meta by 'unique_proteins' 
    emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
    
    if(ncol(pepData$e_meta) == 2){
      e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
      names(e_meta)<-pro_id
    }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
    
    prodata = as.proData(e_data = data.frame(final_result, check.names=check_names), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, is_normalized = is_normalized, check.names = check_names)
    
    #check for isobaricpepData class
    if(inherits(pepData, "isobaricpepData")){
      #update attributes in prodata
      attr(prodata, "isobaric_info") = attr(pepData, "isobaric_info")
      attr(prodata, "isobaric_info")$norm_info$is_normalized = attr(pepData, "isobaric_info")$norm_info$is_normalized
    }
    
    #updating prodata attributes
    attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
    attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
    attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
    
    attr(prodata, "filters")<- attr(pepData, "filters")
    attr(prodata, "group_DF")<- attr(pepData, "group_DF")
    attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  }
  
  return(prodata)
} 
