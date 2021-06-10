#' Applies rrollup function
#'
#' This function applies the rrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#' 
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#'   
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying rrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the rrollup method, peptides are scaled based on a reference
#'   peptide and protein abundance is set as the mean of these scaled peptides.
#'
#' @references Matzke, M. M., Brown, J. N., Gritsenko, M. A., Metz, T. O.,
#'   Pounds, J. G., Rodland, K. D., ... Webb-Robertson, B.-J. (2013). \emph{A
#'   comparative analysis of computational approaches to relative protein
#'   quantification using peptide peak intensities in label-free LC-MS
#'   proteomics experiments}. Proteomics, 13(0), 493-503.
#'   
#'   Polpitiya, A. D., Qian, W.-J., Jaitly, N., Petyuk, V. A., Adkins, J. N.,
#'   Camp, D. G., ... Smith, R. D. (2008). \emph{DAnTE: a statistical tool for
#'   quantitative analysis of -omics data}. Bioinformatics (Oxford, England),
#'   24(13), 1556-1558.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pep_object)
#' result = rrollup(pepData = pep_object)
#' }
#'
#' @rdname rrollup
#' 
rrollup <- function (pepData, combine_fn, parallel = TRUE) {
  
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) {
    
    stop("pepData is not an object of the appropriate class")
    
  }
  
  # check that a protein mapping is provided #
  if(is.null(pepData$e_meta)){
    stop("A mapping to proteins must be provided in order to use the protein_filter function.")
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
  unique_proteins <- unique(temp[[pro_id]])
  
  # set up parallel backend
  if(parallel == TRUE){
    cores<- parallel::detectCores()
    cl<- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }
  
  r <- foreach::foreach(i=1:length(unique_proteins))%dopar%{
    
    row_ind<- which(temp[ ,pro_id] == unique_proteins[i])
    current_subset<- temp[row_ind,]
    current_subset<- current_subset[,-which(names(temp) == pro_id)]
    
    #### Perform R_Rollup ####
    ## store number of peptides ##
    num_peps = nrow(current_subset)
    
    # Create a matrix with one row and the same number of columns as the number
    # of samples. This will be used to create a data frame with only one row
    # after using the apply function columnwise.
    res = matrix(NA, nrow = 1, ncol =  ncol(current_subset))
    
    ## if only 1 peptide, set the protein value to the peptide ##
    if(num_peps==1){
      protein_val = unlist(current_subset)
    }else{
      ## Step 1: Select Reference Peptide -- peptide with least amount of
      ## missing data ##
      na.cnt = apply(is.na(current_subset),1,sum)
      least.na = which(na.cnt == min(na.cnt))
      
      ## If tied, select one with highest median abundance##
      if(length(least.na)>1){
        mds = apply(current_subset,1,median,na.rm=T)[least.na]
        least.na = least.na[which(mds==max(mds))]		
      }
      prot_val = unlist(current_subset[least.na,])
      
      ## Step 2: Ratio all peptides to the reference.  Since the data is on the
      ## log scale, this is the difference ##
      scaling_factor = apply(rep(as.numeric(prot_val),
                                 each = nrow(current_subset)) - current_subset,
                             1,
                             median,
                             na.rm=T)
      
      ## Step 3: Use the median of the ratio as a scaling factor for each
      ## peptide ##
      x_scaled = current_subset + rep(scaling_factor, ncol(current_subset))
      
      ## Step 4: Set Abundance as Median Peptide Abundance ##
      protein_val = apply(x_scaled, 2, combine_fn)
      
    }
    
    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] = protein_val
    
    # Convert the single row matrix to a single row data frame and rename the
    # columns to the original column names.
    res <- data.frame(res)
    names(res) <- names(current_subset)
    
    # Using the foreach function with %dopar% will assign the last element
    # within the curly brackets to the object when foreach is called. In this
    # case res will be assigned to the ith element of r.
    res
    
  }
  
  # Combine the protein abundances (or is it abundanci?).
  final_result <- data.frame(unique_proteins,
                             data.table::rbindlist(r))
  names(final_result)[1] <- pro_id
  
  # Extricate attribute info for creating the proData object.
  check_names <- attr(pepData, "check.names")
  data_scale <- attr(pepData, "data_info")$data_scale
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized
  
  # Create a proData object with the quantitated proteins.
  prodata <- as.proData(e_data = final_result,
                        f_data = pepData$f_data,
                        e_meta = dplyr::select(pepData$e_meta,
                                               -rlang::sym(pep_id)),
                        edata_cname = pro_id,
                        fdata_cname = samp_id,
                        emeta_cname = pro_id,
                        data_scale = data_scale,
                        is_normalized = is_normalized,
                        check.names = check_names)
  
  return (prodata)
  
}
