#' Applies qrollup function 
#' 
#' This function applies the qrollup method to a pepData object for each unique protein and returns a proData object. 
#' 
#' @param pepData an omicsData object of class 'pepData'
#' @param qrollup_thresh is a numeric value; is the peptide abundance cutoff value. 
#' @param combine_fn logical indicating what combine_fn to use, defaults to median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in applying qrollup function. Defaults to TRUE.
#' 
#' @return an omicsData object of class 'proData'
#' 
#' @details In the qrollup method, peptides are selected according to a user selected abundance cutoff value (qrollup_thresh), and protein abundance is set as the mean of these selected peptides.  
#' 
#' @references Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins, J. N., Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical tool for quantitative analysis of -omics data}. Bioinformatics (Oxford, England), 24(13), 1556-1558. 
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = qrollup(pepData = pep_object, qrollup_thresh = 2) 
#'}
#' 
#' @rdname qrollup
#' 

qrollup<- function(pepData, qrollup_thresh, combine_fn = "median", parallel = TRUE){
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  # check that a qrollup_thresh is numeric and between 0 and 1#
  if(!is.numeric(qrollup_thresh) | qrollup_thresh > 1 | qrollup_thresh < 0){
    stop("qrollup_thresh must be Numeric and between 0 and 1")
  }
  #check that combine_fn is one of 'mean', 'median'
  if(!(combine_fn %in% c('median', 'mean'))) stop("combine_fn has to be one of 'mean' or 'median'")
  
  pep_id = attr(pepData, "cnames")$edata_cname
  pro_id = attr(pepData, "cnames")$emeta_cname
  
  pep = data.table(pepData$e_data)
  pro = data.table(pepData$e_meta[,c(pep_id, pro_id)])
  temp = data.table:::merge.data.table(x = pro, y = pep, by = pep_id, all.x = F, all.y = T)
  temp = as.data.frame(temp, check.names=check_names)[,-which(names(temp)==pep_id)]
  
  #pull protein column from temp and apply unique function
  unique_proteins<- unique(temp[[pro_id]])
  
  #assigning function to chosen_combine_fn 
  if(combine_fn == "median"){
    chosen_combine_fn<- combine_fn_median
  }else{chosen_combine_fn = combine_fn_mean}
  
  if(parallel == TRUE){
    
    final_list<- vector("list", length(unique_proteins))
    
    suppressMessages(suppressPackageStartupMessages({
      library(doParallel)
    })
    )
    cores<- detectCores()
    cl<- makeCluster(cores)
    registerDoParallel(cl)
    
    r<-foreach(i=1:length(unique_proteins))%dopar%{
      
      row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
      current_subset<- temp[row_ind,]
      current_subset<- current_subset[,-which(names(temp) == pro_id)]
      
      #### Perform Q_Rollup ####
      ## store number of peptides ##
      num_peps = nrow(current_subset)
      
      res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
      ## if only 1 peptide, set the protein value to the peptide ##
      if(num_peps==1){
        protein_val = unlist(current_subset)
      }else{
        ## Step 1: Subset peptides whose abundance is >= to qrollup_thresh ##
        means = apply(current_subset,1,mean,na.rm=T)
        quantil = quantile(means, probs = qrollup_thresh, na.rm = T)

        new_subset = current_subset[which(means >= quantil), ]
        
        #after step 1 if only 1 peptide, set the protein value to the peptide
        if(nrow(new_subset) == 1){
          protein_val = unlist(new_subset)
        }else{
          ## Step 2: Set protein abundance as the mean/median of peptide abundances 
          protein_val = apply(new_subset, 2, chosen_combine_fn)
        }
      }
      
      res[1,] = protein_val
      res<- data.frame(res)
      names(res)<- names(current_subset)
      
      final_list[[i]]<- res
    }
    stopCluster(cl)
    
    final_result<- do.call(rbind, r)
    final_result<- cbind(unique_proteins, final_result)
    names(final_result)[1]<-pro_id
    
    samp_id = attr(pepData, "cnames")$fdata_cname
    data_scale = attr(pepData, "data_info")$data_scale
    data_norm = attr(pepData, "data_info")$data_norm
    
    #subsetting pepData$e_meta by 'unique_proteins' 
    emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
    
    if(ncol(pepData$e_meta) == 2){
      e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
      names(e_meta)<-pro_id
    }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
    
    prodata = as.proData(e_data = data.frame(final_result, check.names=check_names), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, data_norm = data_norm, check.names = check_names)
    
    #check for isobaricpepData class
    if(inherits(pepData, "isobaricpepData")){
      #update attributes in prodata
      attr(prodata, "isobaric_info") = attr(pepData, "isobaric_info")
      attr(prodata, "data_info")$isobaric_norm = attr(pepData, "data_info")$isobaric_norm
    }
    
    #updating prodata attributes
    attr(prodata, "data_info")$norm_info = attr(pepData, "data_info")$norm_info
    attr(prodata, "data_info")$data_types = attr(pepData, "data_info")$data_types
    attr(prodata, "data_info")$norm_method = attr(pepData, "data_info")$norm_method
    
    attr(prodata, "filters")<- attr(pepData, "filters")
    attr(prodata, "group_DF")<- attr(pepData, "group_DF")
    attr(prodata, "imdanova")<- attr(pepData, "imdanova")
  }
  
  #applying rrollup without doParallel
  
  else{
    final_list<- vector("list", length(unique_proteins))
    
    for(i in 1:length(unique_proteins)){
      
      row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
      current_subset<- temp[row_ind,]
      current_subset<- current_subset[,-which(names(temp) == pro_id)]
      
      #### Perform Q_Rollup ####
      ## store number of peptides ##
      num_peps = nrow(current_subset)
      
      res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
      ## if only 1 peptide, set the protein value to the peptide ##
      if(num_peps==1){
        protein_val = unlist(current_subset)
      }else{
        ## Step 1: Subset peptides whose abundance is >= to qrollup_thresh ##
        non_na_cnt = apply(!is.na(current_subset), 1, sum)
        new_subset = current_subset[which(non_na_cnt >= qrollup_thresh), ]
        
        #after step 1 if only 1 peptide, set the protein value to the peptide
        if(nrow(new_subset) == 1){
          protein_val = unlist(new_subset)
        }else{
          ## Step 2: Set protein abundance as the mean of peptide abundances from qrollup_thresh_subset
          protein_val = apply(new_subset, 2, chosen_combine_fn)
        }
      }
      
      res[1,] = protein_val
      res<- data.frame(res)
      names(res)<- names(current_subset)
      
      final_list[[i]]<- res
    }
    
    final_result<- do.call(rbind, final_list)
    final_result<- cbind(unique_proteins, final_result)
    names(final_result)[1]<-pro_id
    
    samp_id = attr(pepData, "cnames")$fdata_cname
    data_scale = attr(pepData, "data_info")$data_scale
    data_norm = attr(pepData, "data_info")$data_norm
    
    #subsetting pepData$e_meta by 'unique_proteins' 
    emeta_indices<- match(unique_proteins, pepData$e_meta[[pro_id]])
    
    if(ncol(pepData$e_meta) == 2){
      e_meta = as.data.frame(pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)])
      names(e_meta)<-pro_id
    }else {e_meta = pepData$e_meta[emeta_indices, -which(names(pepData$e_meta)==pep_id)]} 
    
    prodata = as.proData(e_data = data.frame(final_result, check.names=check_names), f_data = pepData$f_data, e_meta = e_meta ,edata_cname = pro_id, fdata_cname = samp_id, emeta_cname = pro_id, data_scale = data_scale, data_norm = data_norm, check.names = check_names)
    
    #check for isobaricpepData class
    if(inherits(pepData, "isobaricpepData")){
      #update attributes in prodata
      attr(prodata, "isobaric_info") = attr(pepData, "isobaric_info")
      attr(prodata, "data_info")$isobaric_norm = attr(pepData, "data_info")$isobaric_norm
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