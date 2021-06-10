#' Applies qrollup function
#'
#' This function applies the qrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#' @param qrollup_thresh numeric value between 0 and 1 inclusive. Peptides above
#'   this threshold are used to roll up to the protein level
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying qrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the qrollup method, peptides are selected according to a user
#'   selected abundance cutoff value (qrollup_thresh), and protein abundance is
#'   set as the mean of these selected peptides.
#'
#' @references Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins,
#'   J. N., Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical
#'   tool for quantitative analysis of -omics data}. Bioinformatics (Oxford,
#'   England), 24(13), 1556-1558.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = qrollup(pepData = pep_object, qrollup_thresh = 2)
#' }
#'
#' @rdname qrollup
#'   
qrollup <- function (pepData, qrollup_thresh,
                     combine_fn, parallel = TRUE) {
  
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
  }
  # check that a qrollup_thresh is numeric and between 0 and 1#
  if(!is.numeric(qrollup_thresh) || qrollup_thresh > 1 || qrollup_thresh < 0){
    stop("qrollup_thresh must be Numeric and between 0 and 1")
  }
  
  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname
  
  # Combine e_data and e_meta (just the peptide and protein ID columns) into one
  # data frame by peptide ID. This is a right join with e_meta being the data
  # frame on the right.
  temp <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
                y = pepData$e_data,
                by = pep_id,
                all.x = FALSE,
                all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    data.frame()
  
  #pull protein column from temp and apply unique function
  unique_proteins<- unique(temp[[pro_id]])
  
  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
    
  r<-foreach::foreach(i=1:length(unique_proteins))%dopar%{
    
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
      peps_used <- 1
    }else{
      ## Step 1: Subset peptides whose abundance is >= to qrollup_thresh ##
      means = apply(current_subset,1,mean,na.rm=T)
      quantil = quantile(means, probs = qrollup_thresh, na.rm = T)

      new_subset = current_subset[which(means >= quantil), ]
      peps_used <- nrow(new_subset)
      
      #after step 1 if only 1 peptide, set the protein value to the peptide
      if(nrow(new_subset) == 1){
        protein_val = unlist(new_subset)
      }else{
        ## Step 2: Set protein abundance as the mean/median of peptide abundances 
        protein_val = apply(new_subset, 2, combine_fn)
      }
    }
    
    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] <- protein_val
    
    # Convert the single row matrix to a single row data frame and rename the
    # columns to the original column names.
    res <- data.frame(res)
    names(res) <- names(current_subset)
    
    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case the list containing res and peps_used will be assigned to the ith
    # element of r.
    list(res, peps_used)
    
  }
  
  # Combine the protein abundances (or is it abundanci?).
  final_result <- data.frame(unique_proteins,
                             data.table::rbindlist(
                               lapply(r, function(x) x[[1]])
                             ))
  names(final_result)[1] <- pro_id
  
  # Combine the peptide counts with their corresponding proteins.
  temp_pepes <- data.frame(final_result[, 1],
                           n_peps_used = sapply(r, function(x) x[[2]]))
  names(temp_pepes)[1] <- pro_id
  
  # Combine the peptide counts with pepData$e_meta by protein. This is done to
  # preserve information in the rows of e_meta. For example, there could be
  # fewer unique proteins than there are rows of e_meta. If we do not combine
  # the peptide counts this way then only one row per protein will be kept and
  # the additional information will be lost.
  temp_emeta <- dplyr::left_join(x = pepData$e_meta,
                                 y = temp_pepes,
                                 by = pro_id) %>%
    # Remove peptide id column.
    dplyr::select(-rlang::sym(pep_id))
  
  # Extricate attribute info for creating the proData object.
  check_names <- attr(pepData, "check.names")
  data_scale <- attr(pepData, "data_info")$data_scale
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized
  
  # Create a proData object with the quantitated proteins.
  prodata <- as.proData(e_data = final_result,
                        f_data = pepData$f_data,
                        e_meta = temp_emeta,
                        edata_cname = pro_id,
                        fdata_cname = samp_id,
                        emeta_cname = pro_id,
                        data_scale = data_scale,
                        is_normalized = is_normalized,
                        check.names = check_names)
  
  return(prodata)
  
}
