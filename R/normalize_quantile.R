#' Quantile Normalization
#'
#' Perform quantile normalization
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', 'nmrData', created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively. 
#'
#' @details Quantile normalization is an algorithm for normalizing a set of data
#'   vectors by giving them the same distribution. It is applied to data on the
#'   abundance scale (e.g. not a log scale). It is often used for microarray
#'   data.
#'
#' @return The normalized data is returned in an object of the appropriate S3
#'   class (e.g. pepData), on the same scale as omicsData (e.g. if omicsData
#'   contains log2 transformed data, the normalization will be performed on the
#'   non-log2 scale and then re-scaled after normalization to be returned on the
#'   log2 scale).
#'   
#' @details The method is implemented as described in Bolstad et al. (2003).
#'
#' @examples
#' library(pmartRdata)
#' myfilt <- molecule_filter(omicsData = metab_object)
#' # quantile normalization requires complete data
#' # summary(filter_object = myfilt, min_num = 50)
#' mymetab <- applyFilt(filter_object = myfilt, omicsData = metab_object, min_num = 50)
#' norm_data <- normalize_quantile(omicsData = mymetab)
#'
#' @author Kelly Stratton
#'
#' @references Bolstad, B. M., Irizarry, R. A., Åstrand, M., & Speed, T. P.
#'   (2003). A comparison of normalization methods for high density
#'   oligonucleotide array data based on variance and bias. Bioinformatics,
#'   19(2), 185-193.
#'
#' @export
#'
normalize_quantile <- function (omicsData) {
  
  ## initial checks ##
  
  # check that omicsData is of the appropriate class
  if (!inherits(omicsData, c("proData", "pepData", "lipidData",
                             "metabData", "nmrData"))) {
    
    # Stop the user with an error for sullying the fine reputation of pmart.
    # RESPECT THE STANDARDS!
    stop ("omicsData must be of class 'pepData', 'proData', 'lipidData', 'metabData' or 'nmrData'")
    
  }
  
  # give message if proportion of missing data is > 0
  if (attributes(omicsData)$data_info$prop_missing > 0) {
    
    # Stop the user for their barbaric use of pmart tools on the data.
    stop (paste("The proportion of missing data is ",
                round(attributes(omicsData)$data_info$prop_missing, 2),
                ". Quantile normalization only works with complete data. ",
                "Consider using SPANS to choose an appropriate normalization ",
                "method for a dataset that includes missing values.",
                sep = ""))
    
  }
  
  # Fish out the data scale of the input data. This will be used to convert the
  # data back to the correct scale if the data are input on the log scale.
  input_data_scale <- get_data_scale(omicsData)
  
  # data should be on raw scale #
  if (get_data_scale(omicsData) %in% c("log2", "log10", "log")) {
    omicsData <- edata_transform(omicsData, "abundance")
  }
  
  ## end of initial checks ##
  
  # Pluck out the index of the column containing the metabolite IDs.
  id_col <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Rearrange the names of edata (if the metabolite ID column is not the first
  # colum) because after normalization the metabolite ID column will be the
  # first column of e_data.
  edata_cnames <- c(names(omicsData$e_data[id_col]),
                    names(omicsData$e_data[-id_col]))
  
  # Extract the column containing the metabolite IDs. This will be used to
  # reassemble e_data after normalization.
  edata_idcol <- omicsData$e_data[, id_col]
  
  # 1. Transpose and sort/order each column of e_data in ascending order.
  x_sort <- apply(t(omicsData$e_data[, -id_col]), 2,
                  function(c) sort(c, na.last = FALSE))
  x_ord <- apply(t(omicsData$e_data[, -id_col]), 2,
                 function(c) order(c, na.last = FALSE))
  
  # 2. Compute the mean for each sorted row. These values will replace the
  # original counts according to their order within each colum.
  row_means <- rowMeans(x_sort, na.rm = TRUE)
  
  # 3. Replace, within every row, each abundance value with its corresponding
  # mean. For example, in row one the smallest abundance value will be replaced
  # with the smallest mean, the second smallest abundance value will be replaced
  # with the second smallest mean, and so on. This will be done for each row of
  # x_sort.
  x_norm <- omicsData$e_data[, -id_col]
  for (i in 1:ncol(x_sort)) {
    
    # Replace abundance values with their corresponding means. The order of the
    # means corresponds with the original order of the abundance values. For
    # example, if the smallest abundance value in the first column of e_data
    # occurs in row 5 then the smallest mean for this column will also be in row
    # five (since we are working on the transposed data the column/row indices
    # are switched).
    x_norm[i, x_ord[, i]] <- row_means
    
  }
  
  # Reassemble e_data with the normalized data. NOTE: The ID column will now be
  # the first column of e_data if it wasn't already.
  omicsData$e_data <- data.frame(edata_idcol, x_norm)
  names(omicsData$e_data) <- edata_cnames
  
  
  # Update the norm_info list within the data_info attribute.
  attributes(omicsData)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "quantile" # new attribute as of 12/21/17
  )
  
  # Check if the input data was on a log scale.
  if (input_data_scale != "abundance") {
    
    # Transfigure the data back to the correct log scale.
    omicsData <- edata_transform(omicsData, input_data_scale)
    
  }
  
  return(omicsData)
  
}
