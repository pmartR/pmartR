#' Protein Quantitation using Mean or Median Peptide Abundances
#' 
#' This function takes in a pepData object and returns a proData object
#' 
#' @param pepData omicsData object of class 'pepData'
#' 
#' @param combine_fn A character string that can either be 'mean' or 'median'.
#' 
#' @return An omicsData object of class 'proData'
#' 
pquant <- function (pepData,
                    combine_fn) {
  
  # check that pepData is of appropraite class #
  if (!inherits(pepData, "pepData")) {
    
    stop ("pepData is not an object of the appropriate class")
    
  }
  
  # check that a protein mapping is provided #
  if (is.null(pepData$e_meta)) {
    
    stop (paste("A mapping to proteins must be provided in order to use the",
                "protein_filter function.",
                sep = " "))
    
  }
  
  # Fish out the e_data, f_data, and e_meta column names corresponding to the
  # peptide, sample, and protein IDs.
  pep_id <- attr(pepData, "cnames")$edata_cname
  samp_id = attr(pepData, "cnames")$fdata_cname
  pro_id <- attr(pepData, "cnames")$emeta_cname

  # Quantitate the heck out of the peptides!
  res <- merge(x = pepData$e_meta[, c(pep_id, pro_id)],
               y = pepData$e_data,
               by = pep_id,
               all.x = FALSE,
               all.y = TRUE) %>%
    dplyr::select(-rlang::sym(pep_id)) %>%
    dplyr::group_by(!!rlang::sym(pro_id)) %>%
    dplyr::mutate(dplyr::across(.cols = -dplyr::any_of(pro_id),
                                .fns = combine_fn)) %>%
    data.frame()
  
  
  # Extricate attribute info for creating the proData object.
  check_names <- attr(pepData, "check.names")
  data_scale <- attr(pepData, "data_info")$data_scale
  is_normalized <- attr(pepData, "data_info")$norm_info$is_normalized

  # Create a proData object with the quantitated proteins.
  prodata <- as.proData(e_data = res,
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
