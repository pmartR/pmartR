#' Reduce Dimension of Data for Exploratory Data Analysis
#'
#' Calculate principal components using projection pursuit estimation, which implements an expectation-maximization (EM) estimation algorithm when data is missing.
#'
#' @param omicsData an object of the class 'pepdata', 'prodata', 'metabData', 'lipidData', 'nmrData' usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, respectively.
#' @param k integer number of principal components to return. Defaults to 2.
#'
#' @return a data.frame with first \code{k} principal component scores, sample identifiers, and group membership for each sample (if group designation was previously run on the data).
#'
#' @references Redestig, H., Stacklies, W., Scholz, M., Selbig, J., & Walther, D. (2007). \emph{pcaMethods - a bioconductor package providing PCA methods for incomplete data}. Bioinformatics. 23(9): 1164-7.
#'
#' @details Any biomoleculs seen in only one sample or with a variance less than 1E-6 across all samples are not included in the PCA calculations. This function leverages code from \code{\link[pcaMethods]{pca}}.
#' 
#' @examples 
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object <- edata_transform(omicsData = lipid_object, data_scale="log2")
#' lipid_object <- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' pca_lipids <- dim_reduction(omicsData = lipid_object)
#' plot(pca_lipids)
#' summary(pca_lipids)
#' }
#' 
#' @export
#' @rdname dim_reduction
#' @name dim_reduction
#'
dim_reduction <- function(omicsData, k = 2){
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData","proData","metabData", "lipidData", "nmrData"))) stop("omicsR_data must be an object of class 'pepdata','prodata', 'metabData', 'lipidData', or 'nrmData'.")

  # check that group designation has been run #
  if(!("group_DF" %in% names(attributes(omicsData)))) warning("group_designation has not been run on this data and may limit plotting options")

  # data should be log transformed #
  if(!attr(omicsData, "data_info")$data_scale %in% c("log2", "log10", "log")){
    warning("omicsData$e_data should be log transformed prior to calling dim_reduction. See documentation for edata_transform function for more information.")
  }

  samp_id = attr(omicsData, "cnames")$fdata_cname
  pep_id = attr(omicsData, "cnames")$edata_cname

  temp_data = omicsData$e_data[, -which(names(omicsData$e_data) == pep_id)]

  ## check for samples seen in only one sample or no samples and remove ##
  minsamps = which(apply(!is.na(temp_data), 1, sum) < 2)
  if(length(minsamps) > 0){
    temp_data = temp_data[-minsamps,]
  }

  ## check for near zero variance features and remove ##
  minvars = which(apply(temp_data, 1, var, na.rm = T) < 0.000001)
  if(length(minvars) > 0){
    temp_data[-minvars, ]
  }

  pca_res = pcaMethods::pca(object = as.matrix(t(temp_data)), method = "ppca", scale = "vector", nPcs = k)
  pca_ests = pca_res@scores[,1:k]

  temp_res = data.frame(SampleID = names(temp_data), pca_ests)

  class(temp_res) <- "dimRes"

  if(!is.null(attr(omicsData, "group_DF"))){
  attr(temp_res, "group_DF") <- attr(omicsData, "group_DF")
  }else{attr(temp_res, "group_DF") <- NULL}

  attr(temp_res, "R2") <- pca_res@R2

  return(temp_res)
}
