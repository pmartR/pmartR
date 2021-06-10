#' Applies zrollup function
#'
#' This function applies the zrollup method to a pepData object for each unique
#' protein and returns a proData object.
#'
#' @param pepData an omicsData object of class 'pepData'
#' @param combine_fn logical indicating what combine_fn to use, defaults to
#'   median, other option is mean
#' @param parallel logical indicating whether or not to use "doParallel" loop in
#'   applying zrollup function. Defaults to TRUE.
#'
#' @return an omicsData object of class 'proData'
#'
#' @details In the zrollup method, peptides are scaled as, pep_scaled = (pep -
#'   median)/sd, and protein abundance is set as the mean of these scaled
#'   peptides.
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
#' result = zrollup(pepData = pep_object)
#' }
#'
#' @rdname zrollup
#'   
zrollup <- function (pepData, combine_fn, parallel = TRUE) {
  
  check_names = getchecknames(pepData)
  
  # check that pepData is of appropraite class #
  if(!inherits(pepData, "pepData")) stop("pepData is not an object of the appropriate class")
  
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
    
    peptides_scaled = apply((current_subset - medians_mat)/standiv_mat,
                            2, combine_fn)
    
    # Convert protein_val to a matrix. This needs to be done because a vector
    # cannot be converted to a single row data frame.
    res[1,] = peptides_scaled
    
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
