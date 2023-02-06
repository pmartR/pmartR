#' Summary of pmartR Analysis Functions
#'
#' Provide basic summaries for results objects from the pmartR package.
#'
#' @return a summary table or list for the pmartR results object
#'
#' @examples
#' library(pmartRdata)
#' mypep <- group_designation(omicsData = pep_object, main_effects = "Phenotype")
#' mypep <- edata_transform(omicsData = mypep, data_scale = "log2")
#' 
#' norm_result <- normalize_global(omicsData = mypep, norm_fn = "median", subset_fn = "all")
#' summary(norm_result)
#' 
#' \dontrun{
#' spans_results <- spans_procedure(omicsData = mypep)
#' summary(spans_results)
#' }
#' 
#' dim_results <- dim_reduction(omicsData = mypep)
#' summary(dim_results)
#' 
#' cor_results <- cor_result(omicsData = mypep)
#' summary(cor_results) 
#'
#' @author Lisa Bramer, Kelly Stratton, Thomas Johansen
#'
#' @export

#' @export
#' @rdname summary-pmartR-results
#' @name summary-pmartR-results
#' @param omicsNorm object of class normRes, produced by calling
#'   \code{\link{normalize_global} with option apply_norm=FALSE}
summary.normRes <- function(omicsNorm) {
  # summary for the normalization object (if normalize=FALSE when calling normalization_calc)
  
  # get values #
  res <- list(subset_fn = paste0(toupper(substr(omicsNorm$subset_fn, 1, 1)), substr(omicsNorm$subset_fn, 2, nchar(omicsNorm$subset_fn))),
              norm_fn = paste0(toupper(substr(omicsNorm$norm_fn, 1, 1)), substr(omicsNorm$norm_fn, 2, nchar(omicsNorm$norm_fn))),
              bt = ifelse(is.null(omicsNorm$parameters$backtransform$location), FALSE, TRUE),
              n_features_calc = omicsNorm$n_features_calc,
              prop_features_calc = omicsNorm$prop_features_calc)
  
  # construct output #
  res <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
  catmat <- data.frame(unlist(res, use.names = FALSE))
  
  # assemble text strings #
  colnames(catmat) <- NULL
  rownames(catmat) <- c("Subset of data used for normalization calculation ", "Normalization function used for normalization calculation ", "Backtransform requested ", "Number of biomolecules in subset ", "Proportion of biomolecules in subset ")
  
  cat("\nSummary of 'normRes' Object\n---------------------------")
  cat(capture.output(catmat), sep = "\n")
  cat("\n")
  
  return(invisible(res)) # should this be returning "catmat" instead?
}

#' @export
#' @rdname summary-pmartR-results
#' @name summary-pmartR-results
#' @param SPANSres_obj object of class SPANSRes created by calling
#'   \code{\link{spans_procedure}} on a grouped pepData or proData object
summary.SPANSRes <- function(SPANSRes_obj){
  
  spanscores <- sort(unique(SPANSRes_obj$SPANS_score), decreasing = TRUE)
  SPANSRes_obj <- SPANSRes_obj %>% 
    dplyr::mutate(rank = dplyr::dense_rank(dplyr::desc(SPANS_score)))
  
  cat("\nSummary of spans procedure\n")
  
  cat(paste0("\nHighest ranked method(s)\n"))
  cat(capture.output(head(SPANSRes_obj)), sep="\n")
  
  cat(paste0("\nNumber of input methods:  ", nrow(SPANSRes_obj)))
  cat(paste0("\nNumber of methods scored:  ", sum(as.logical(SPANSRes_obj$passed_selection))))
  cat(paste0("\nAverage molecules used in normalization:  ", round(mean(as.numeric(SPANSRes_obj$mols_used_in_norm), na.rm = TRUE))))
  
}

#' @export
#' @rdname summary-pmartR-results
#' @name summary-pmartR-results
#' @param dimRes_object object of class dimRes, which is a list containing
#'   sample identifiers and the principle components scores
summary.dimRes <- function(dimRes_object){
  
  
  # get R^2 values #
  r2 <- rep(NA, length(dimRes_object) -1)
  for(i in 1:length(r2)){
    r2[i] <- round(attr(dimRes_object, "R2")[i],3)
  }
  
  dim_summary = data.frame(R_squared = r2)
  row.names(dim_summary) <- attributes(dimRes_object)$names[-1]
  
  message("Summary of 'dimRes' Object\n\n")
  message(capture.output(dim_summary), sep="\n")
  
  return(invisible(dim_summary))
}

#' @export
#' @rdname summary-pmartR-results
#' @name summary-pmartR-results
#' @param corRes_object object of class corRes
summary.corRes <- function(corRes_object){
  
  cor_summary <- summary(corRes_object[upper.tri(corRes_object)])
  
  cat("Summary of Correlation Matrix (Upper Triangle)\n\n")
  cat(capture.output(cor_summary), sep="\n")
  
  return(invisible(cor_summary))
}