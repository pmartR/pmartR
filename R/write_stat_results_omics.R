#' Creates a list of three sheets, Normalized Data, DA Test Results, and Only DA biomolecules - for OMICS project
#
#' 
#' @param omicsData an object of one of the classes "pepData", "proData", "metabData", or "lipidData", usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, or \code{\link{lipidData}}, respectively.
#' @param statResData an object of class statRes, created by \code{\link{imd_anova}}
#' @param refCondition character string, what the reference condition looks like (e.g. "MOCK", "Mock", "mock", etc.)
#' @param filePath default is NULL, if NULL then workbook will not be written out, if it is a character string specifying the file path, file name and .xlsx extension
#' 
#' @return a list of three data frames Normalized_Data, DA_Test_Results, and DA_"data type"_Only
#' 
#' @author Natalie Heller
#' 
#' 
write_stat_results_omics <- function(omicsData, statResData, refCondition, filePath = NULL){
  
  if(!require(openxlsx)){
    stop("openxlsx package is required")
  }
  
  ## INITIAL CHECKS ##
  ## Make sure both ojects are present
  if(is.null(omicsData)){
    stop("'omicsData' must be provided")
  }
  if(is.null(statResData)){
    stop("'statResData' must be provided")
  }
  ## Make sure both are the correct class
  if(!(class(omicsData) %in% c("pepData", "proData", "metabData", "lipidData"))){
    stop("omicsData is not of the appopriate class")
  }
  if(!(class(statResData) %in% c("statRes"))){
    stop("statResData is not of the appropriate class")
  }
  ## END OF CHECKS ##
  
  ## First get the Normalized Data tab with metadata
  if("e_meta" %in% attributes(omicsData)$names){
    Normalized_Data <- left_join(omicsData$e_meta, omicsData$e_data, by = intersect(names(omicsData$e_meta), names(omicsData$e_data)))
  }else{
    Normalized_Data <- omicsData$e_data
  }
  
  
  ## Get the query condition(s)
  ii <- unique(as.character(attributes(statResData)$group_DF$VIRUS))
  q_cond <- ii[which(ii != refCondition)]
  
  pval <- statResData$P_values
  flag <- statResData$Flags
  fc <- statResData$Full_results[, c(which(attributes(omicsData)$cname$edata_cname == colnames(statResData$Full_results)), grep("fold_change", tolower(colnames(statResData$Full_results))))]
  
  # Rename columns of pval
  a1 <- which(attributes(omicsData)$cname$edata_cname == colnames(pval))
  colnames(pval)[-a1] <- gsub(pattern = "pvals_", replacement = "", x = tolower(colnames(pval)[-a1]))
  colnames(pval)[-a1] <- paste(q_cond, colnames(pval)[-a1], "Pval", sep = "_")
  
  # Renames the columns of flag
  a2 <- which(attributes(omicsData)$cname$edata_cname == colnames(flag))
  colnames(flag)[-a2] <- gsub(pattern = "flags_", replacement = "", x = tolower(colnames(flag)[-a2]))
  colnames(flag)[-a2] <- paste(q_cond, colnames(flag)[-a1], "Flag", sep = "_")
  
  a3 <- which(attributes(omicsData)$cname$edata_cname == colnames(fc))
  colnames(fc)[-a3] <- gsub(pattern = "Fold_change_", replacement = "", x = colnames(fc)[-a3])
  colnames(fc)[-a3] <- paste(colnames(fc)[-a3], "Log2FC", sep = "_")
  
  DA_Test_Results <- merge(x = pval, y = fc, by = intersect(colnames(pval), colnames(fc)))
  DA_Test_Results <- merge(x = DA_Test_Results, y = flag, by = intersect(colnames(DA_Test_Results), colnames(flag)))
  
  ## Now get the DA_Molecule_Only tab
  tmp <- flag
  tmp$TMP <- apply(flag[, -a1], 1, function(x){abs(sum(x))})
  indx <- which(tmp$TMP != 0)
  
  DA_Only <- DA_Test_Results[indx, ]
  
  outputData <- list()
  outputData$Normalized_Data <- Normalized_Data
  outputData$DA_Test_Results <- DA_Test_Results
  
  if(class(omicsData) == "metabData"){
    outputData$DA_Metabolites_Only <- DA_Only
  }else{
    if(class(omicsData) == "pepData"){
      outputData$DA_Peptides_Only <- DA_Only
    }else{
      if(class(omicsData) == "proData"){
        outputData$DA_Proteins_Only <- DA_Only
      }else{
        outputData$DA_Lipids_Only <- DA_Only
      }
    }
  }
  if(!is.null(filePath)){
    write.xlsx(x = outputData, file = filePath)
  }
}