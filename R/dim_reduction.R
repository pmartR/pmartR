#' Reduce Dimension of Data for Exploratory Data Analysis
#'
#' For data types other than seqData, this function calculates principal
#' components using projection pursuit estimation, which implements an
#' expectation-maximization (EM) estimation algorithm when data is missing. For
#' seqData counts, a generalized version of principal components analysis for
#' non-normally distributed data is calculated under the assumption of a
#' negative binomial distribution with global dispersion.
#'
#' @param omicsData an object of the class 'pepdata', 'prodata', 'metabData',
#'   'lipidData', 'nmrData', or 'seqData', created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, or
#'   \code{\link{as.seqData}}, respectively.
#' @param k integer number of principal components to return. Defaults to 2.
#'
#' @return a data.frame with first \code{k} principal component scores, sample
#'   identifiers, and group membership for each sample (if group designation was
#'   previously run on the data). The object is of class dimRes (dimension
#'   reduction Result).
#'
#' @references Redestig H, Stacklies W, Scholz M, Selbig J, & Walther
#'   D (2007). \emph{pcaMethods - a bioconductor package providing PCA methods
#'   for incomplete data}. Bioinformatics. 23(9): 1164-7.
#'   
#'   Townes FW, Hicks SC, Aryee MJ, Irizarry RA (2019). \emph{Feature selection
#'   and dimension reduction for single-cell RNA-seq based on a multinomial
#'   model.} Genome Biol. 20, 1–16.
#'   
#'   Huang H, Wang Y, Rudin C, Browne EP (2022). \emph{Towards a comprehensive
#'   evaluation of dimension reduction methods for transcriptomic data
#'   visualization.} Communications Biology 5, 719.
#'   
#' @details Any biomolecules seen in only one sample or with a variance less
#'   than 1E-6 across all samples are not included in the PCA calculations. This
#'   function leverages code from \code{\link[pcaMethods]{pca}} and
#'   \code{\link[glmpca]{glmpca}} .
#' 
#' @examples 
#' library(pmartRdata)
#' 
#' mylipid <- edata_transform(omicsData = lipid_neg_object, data_scale="log2")
#' mylipid <- group_designation(omicsData = mylipid, main_effects = "Virus")
#' pca_lipids <- dim_reduction(omicsData = mylipid)
#' 
#' myseq <- group_designation(omicsData = rnaseq_object, main_effects = "Virus")
#' pca_seq <- dim_reduction(omicsData = myseq)
#' 
#' @export
#'
dim_reduction <- function (omicsData, k = 2){
  
  # check that omicsData is of appropriate class #
  if(!inherits(omicsData, c("pepData", "proData", "metabData",
                            "lipidData", "nmrData", "seqData")))
    stop(paste("omicsData must be an object of class 'pepdata','prodata',",
               "'metabData', 'lipidData', 'nrmData', or 'seqData'.",
               sep = " "))

  # check that group designation has been run #
  if(!("group_DF" %in% names(attributes(omicsData))))
    warning(paste("group_designation has not been run on this data and may",
                  "limit plotting options",
                  sep = " "))

  # data should be log transformed #
  if (get_data_scale(omicsData) == "abundance") {
    stop ("omicsData must be log transformed prior to calling dim_reduction.")
  }
  
  # Check the data scale. Data must be on one of the log scales.
  if (inherits(omicsData, "seqData") && 
      get_data_scale(omicsData) != "counts" ) {
    
    # Welcome to the pit of despair!
    stop ("seqData must be untransformed prior to calling dim_reduction.")
    
  }

  samp_id = attr(omicsData, "cnames")$fdata_cname
  pep_id = attr(omicsData, "cnames")$edata_cname

  temp_data = omicsData$e_data[, -which(names(omicsData$e_data) == pep_id)]

  ## check for samples seen in only one sample or no samples and remove ##
  minsamps = which(rowSums(!is.na(temp_data)) < 2)
  if(length(minsamps) > 0){
    temp_data = temp_data[-minsamps,]
  }

  ## check for near zero variance features and remove ##
  minvars = which(apply(temp_data, 1, var, na.rm = T) < 0.000001)
  if(length(minvars) > 0){
    temp_data = temp_data[-minvars, ]
  }
  
  if(inherits(omicsData, "seqData")){
    # GLM-PCA
    # https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pca-plot-using-generalized-pca
    
    ## Supposed to be performed on raw counts
    ## Temp data looks like an e_data w/o the identifier column
    pca_res <- glmpca::glmpca(temp_data, L=k, fam = "nb")
    pca_ests <- pca_res$factors ### Plot from here for PC1 and PC2
    colnames(pca_ests) <- paste0("PC", 1:ncol(pca_ests))

  } else {
    
    pca_res = pcaMethods::pca(object = as.matrix(t(temp_data)),
                              method = "ppca",
                              scale = "vector",
                              nPcs = k)
    pca_ests = pca_res@scores[,1:k]
    
  }

  temp_res = data.frame(SampleID = names(temp_data), pca_ests)

  class(temp_res) <- "dimRes"

  attr(temp_res, "group_DF") <- get_group_DF(omicsData)

  
  if(!inherits(pca_res, "glmpca")){
    attr(temp_res, "R2") <- pca_res@R2
  }

  return(temp_res)
  
}
