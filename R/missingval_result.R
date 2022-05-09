#' Creates an object of class naRes
#'
#' This function takes in an omicsData object, outputs a list of two data
#' frames, one contains NA by sample, the second contains NA by molecule
#'
#' @param omicsData an object of class "pepData", "proData", "metabData",
#'   "lipidData", or "nmrData", created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @return outputs a list of two data frames, one contains NA by sample (missing
#'   values per sample), the second contains NA by molecule (missing values per
#'   molecule). The output is assigned class 'naRes'.
#'
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data("lipid_object")
#' data("metab_object")
#'
#' result = missingval_result(lipid_object)
#' result2 = missingval_result(metab_object)
#'
#' }
#'
#' @rdname missingval_result
#'
#' @export
#' 
missingval_result<- function(omicsData){
  
  # Add a check 
  
  #check for correct class
  if(!inherits(omicsData, c("pepData", "proData", "lipidData",
                            "metabData", "nmrData", "seqData"))) {
    
    stop (paste("omicsData must have class of the following, 'pepData',",
                "'proData', 'lipidData', 'metabData', 'nmrData', 'seqData'",
                sep = " "))
    
  }
  
  
  # pulling cname attr
  # edata_cname<- attr(omicsData, "cnames")$edata_cname
  edata_cname<- get_edata_cname(omicsData)
  edata_cname_id<- which(names(omicsData$e_data) == edata_cname)
  # fdata_cname<- attr(omicsData, "cnames")$fdata_cname
  fdata_cname<- get_fdata_cname(omicsData)
  
  # Count the number of NA or zeros values per column.
  if(inherits(omicsData, "seqData")){
    res_per_col<- colSums((omicsData$e_data[, -edata_cname_id]) == 0)
    res_by_sample<- data.frame(
      "sample_names" = names(omicsData$e_data[, -edata_cname_id]),
      "num_zeros" = as.numeric(res_per_col)
    )
  } else {
    res_per_col<- colSums(is.na(omicsData$e_data[, -edata_cname_id]))
    res_by_sample<- data.frame(
      "sample_names" = names(omicsData$e_data[, -edata_cname_id]),
      "num_NA" = as.numeric(res_per_col)
    )
  }
  
  names(res_by_sample)[1] <- fdata_cname
  
  # Merge res_by_sample with f_data to get additional columns of f_data. For
  # example, the Group and VizSampNames columns. Group is used to color the plot
  # and VizSampNames is used to display shorter sample names.
  res_by_sample <- merge(res_by_sample, omicsData$f_data, by = fdata_cname)
  
  # Check if the group designation function has been run. The group_DF info will
  # be used to add the "Group" column to res_by_sample. This column may contain
  # the same data as another column in f_data but it will have a different name
  # from the f_data column.
  if (!is.null(attr(omicsData, "group_DF"))) {

    res_by_sample <- merge(res_by_sample, attr(omicsData, "group_DF"))

  }
  
  # Count the number of NA values per row.
  
  if(inherits(omicsData, "seqData")){
    
    res_per_row <- rowSums(omicsData$e_data[, -edata_cname_id] == 0)
    
    res_by_molecule <- data.frame("molecule"= omicsData$e_data[, edata_cname_id],
                                 "num_zeros"= as.numeric(res_per_row))
    names(res_by_molecule)[1] <- edata_cname
    
    result<- list("zeros.by.sample" = res_by_sample,
                  "zeros.by.molecule" = res_by_molecule)
    class(result) <- "naRes"
    
    attr(result, "cnames") <- list("edata_cname" = edata_cname,
                                   "fdata_cname" = fdata_cname)
    
  } else {
    
    res_per_row <- rowSums(is.na(omicsData$e_data[, -edata_cname_id]))
    res_by_molecule <- data.frame("molecule"= omicsData$e_data[, edata_cname_id],
                                 "num_NA"= as.numeric(res_per_row))
    names(res_by_molecule)[1] <- edata_cname
    
    result<- list("na.by.sample" = res_by_sample,
                  "na.by.molecule" = res_by_molecule)
    class(result) <- "naRes"
    
    attr(result, "cnames") <- list("edata_cname" = edata_cname,
                                   "fdata_cname" = fdata_cname)
  }

  return (result)  
  
}
