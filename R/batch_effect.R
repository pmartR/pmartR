#' ComBat adjusted object
#' 
#' This function returns a pmart object that has been undergone for ComBat
#' batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param use_groups logical argument indicating if the group designation should
#'   be applied to the data for batch correction. Defaults to FALSE. If TRUE,
#'   ComBat runs with both a batch effect and a group designation
#'
#' @return Object of same class as omicsData that has been undergone
#'   ComBat batch effect correction
#' 
#' @example
#' \dontrun{
#' library(pmartRdata)
#' data("pep_object")
#' batch_correct(pep_object)
#' }
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_combat <- function(omicsData,use_groups = FALSE){
  
  # run checks -----------------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    stop (paste("omicsData must have batch_id attribute for batch correction",
                sep = ' '))
  }
  
  # check that omicsData has group info if use_groups is true
  if (use_groups == TRUE & is.null(attr(omicsData,"group_DF"))) {
    stop (paste("omicsData must have group designation if use_groups = TRUE"))
  }
  
  # important information for downstream analysis ------------------------------
  
  # get the variables we need to create pmart object
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdat = omicsData$f_data
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  emet =  NULL
  emet_cname = NULL
  if(!is.null(omicsData$e_meta)){
    emet = omicsData$e_meta
    emet_cname = pmartR::get_emeta_cname(omicsData)
  }
  
  
  # make the data available for combat -----------------------------------------
  # pull out the edata
  edata <- omicsData$e_data
  # find which column id is the edata cname
  id_col <- which(colnames(omicsData$e_data) == edat_cname)
  # add rownames
  rownames(edata) <- edata[,id_col]
  # remove the first column now that we have the molecules stored under rownames
  edata <- edata[,-id_col]
  # save the column ordering of edata
  edatOrder <- colnames(edata)
  # convert to matrix
  edataMatrix <- as.matrix(edata)
  
  # obtain vector of the batch
  batchVector <- attributes(attr(omicsData,"group_DF"))$batch_id
  # make sure that the batch vector matches the order of the sample ID
  batchVector <- batchVector %>%
    dplyr::arrange(match(edat_cname,edatOrder))
  batchVector <- batchVector[,2]
  
  # order may need to be reordered
  # determine if we need group or not
  
  if(use_groups == TRUE) {
    # make sure the group vector matches the order of the dataframe
    groupDat <- attr(omicsData,"group_DF") %>%
      dplyr::arrange(match(edat_cname,edatOrder))
    mod_group <- model.matrix(~as.factor(groupDat$Group))
  } else {
    mod_group <- NULL
  }
  
  # run ComBat
  combat_obj <- sva::ComBat(dat = edataMatrix, batch = batchVector, mod = mod_group)
  
  # convert back to pmart object -----------------------------------------------
  
  # convert back to pmart object
  combat_obj <- data.frame(combat_obj)
  combat_obj <- tibble::rownames_to_column(combat_obj, var = pmartR::get_edata_cname(omicsData))
  
  # create the pmart object
  if(attr(omicsData,"class") == "pepData"){
    pmartObj = pmartR::as.pepData(e_data = combat_obj,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(attr(omicsData,"class") == "proData"){
    pmartObj = pmartR::as.proData(e_data = combat_obj,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(attr(omicsData,"class") == "metabData"){
    pmartObj = pmartR::as.metabData(e_data = combat_obj,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(attr(omicsData,"class") == "lipidData"){
    pmartObj = pmartR::as.lipidData(e_data = combat_obj,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(attr(omicsData,"class") == "nmrData"){
    pmartObj = pmartR::as.nmrData(e_data = combat_obj,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  get_data_info(pmartObj)

  # might want to track if batch correction was done
  
  # Update the data_info attribute.
  attr(pmartObj, 'data_info') <- set_data_info(
    e_data = omicsData$e_data,
    edata_cname = pmartR::get_edata_cname(omicsData),
    data_scale_orig = pmartR::get_data_scale_orig(omicsData),
    data_scale = pmartR::get_data_scale(omicsData),
    data_types = pmartR::get_data_info(omicsData)$data_types,
    norm_info = pmartR::get_data_info(omicsData)$norm_info,
    is_normalized = pmartR::get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = pmartR::get_data_info(omicsData)$batch_info,
    is_bc = pmartR::get_data_info(omicsData)$batch_info$is_bc
  )
  
  # Add the group information to the group_DF attribute in the omicsData object.
  attr(pmartObj, "group_DF") = attr(omicsData,"group_DF")

  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "combat",
    params = list()
  )
  
  # Update the meta_info attribute.
  attr(pmartObj, 'meta_info') <- pmartR:::set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = get_emeta_cname(omicsData)
  )
  
  return(pmartObj)
}
