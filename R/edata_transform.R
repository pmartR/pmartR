#' Apply a Transformation to the Data
#'
#' This function applies a transformation to the e_data element of omicsData
#'
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', or 'seqData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, \code{\link{as.nmrData}}, \code{\link{as.seqData}}, respectively.
#' @param data_scale a character string indicating the type of transformation to be applied to the data. Valid values for 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData': 'log2', 'log', 'log10', or 'abundance'. Valid values for 'seqData': 'upper', 'median', 'lcpm'. A value of 'abundance' indicates the data has previously undergone one of the log transformations and should be transformed back to raw values with no transformation applied. For 'seqData', 'lcpm' transforms by log2 counts per million, 'upper' transforms by the upper quartile of non-zero counts, and 'median' transforms by the median of non-zero counts. 
#' 
#' @details This function is intended to be used before analysis of the data begins. Data are typically analyzed on a log scale.
#'
#' @return data object of the same class as omicsData
#'
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(metab_object)
#' metab_object2 <- edata_transform(omicsData = metab_object, data_scale="log2")
#' attr(metab_object2, "data_info")$data_scale
#'}
#' @author Kelly Stratton, Natalie Heller
#'
#' @export
#' 
edata_transform <- function (omicsData, data_scale) {
  
  # Initial checks -------------------------------------------------------------

  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData",
                             "lipidData", "nmrData"))) {
    
    # Throw an error that the input for omicsData is not the appropriate class.
    stop(paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
               "'lipidData', or 'nmrData'",
               sep = ' '))
    
  } 
  
  # check that data_scale is one of the acceptable options #
  if (!(data_scale %in% c('log2', 'log10', 'log', 'abundance'))) {
    
    # Tell the user that the input to data_scale is an abomination!
    stop (paste(data_scale, "is not a valid option for 'data_scale'.",
                "See details of as.pepData for specifics.",
                sep=" "))
    
  }

  # Check to make sure the data isn't already on the scale input by the user.
  if(get_data_scale(omicsData) == data_scale) {
    
    # Stop all further calculations with an error message.
    stop(paste("Data is already on",
               data_scale,
               "scale.",
               sep = " "))
    
  }
  
  # Perform the actual transmogrification --------------------------------------
  
  # Fish out the column index where edata_cname occurs.
  iCol <- which(names(omicsData$e_data) == get_edata_cname(omicsData))
  
  # Extract the data_scale from the omics data object.
  scale <- get_data_scale(omicsData)
  
  # Execute the transmogrification given the current scale and the input scale.
  switch(scale,
         
         # Transmogrify the data from abundance to something else.
         'abundance' = {
           
           # Find input data scale and "make that change".
           if (data_scale == "log"){
             
             # Natural logify the data.
             omicsData$e_data[, -iCol] <- log(omicsData$e_data[, -iCol])
             
           } else if (data_scale == "log2") {
             
             # Log base 2ify the data.
             omicsData$e_data[, -iCol] <- log2(omicsData$e_data[, -iCol])
             
           } else if (data_scale == "log10") {
             
             # Log base 10ify the data.
             omicsData$e_data[, -iCol] <- log10(omicsData$e_data[, -iCol])
             
           }
           
         },
         
         # Mutate the data from log to another scale.
         'log' = {
           
           # Find input data scale and "make that change".
           if (data_scale == "abundance"){
             
             # Natural logify the data.
             omicsData$e_data[, -iCol] <- exp(omicsData$e_data[, -iCol])
             
           } else if (data_scale == "log2") {
             
             # Log base 2ify the data.
             omicsData$e_data[, -iCol] <- log2(exp(omicsData$e_data[, -iCol]))
             
           } else if (data_scale == "log10") {
             
             # Log base 10ify the data.
             omicsData$e_data[, -iCol] <- log10(exp(omicsData$e_data[, -iCol]))
             
           }
           
         },
         
         # Recast the data from the log2 scale to another scale.
         'log2' = {
           
           # Find input data scale and "make that change".
           if (data_scale == "abundance"){
             
             # Natural logify the data.
             omicsData$e_data[, -iCol] <- 2^(omicsData$e_data[, -iCol])
             
           } else if (data_scale == "log") {
             
             # Log base 2ify the data.
             omicsData$e_data[, -iCol] <- log(2^(omicsData$e_data[, -iCol]))
             
           } else if (data_scale == "log10") {
             
             # Log base 10ify the data.
             omicsData$e_data[, -iCol] <- log10(2^(omicsData$e_data[, -iCol]))
             
           }
           
         },
         
         # Change the data from the log10 scale to a different one.
         'log10' = {
           
           # Find input data scale and "make that change".
           if (data_scale == "abundance"){
             
             # Natural logify the data.
             omicsData$e_data[, -iCol] <- 10^(omicsData$e_data[, -iCol])
             
           } else if (data_scale == "log") {
             
             # Log base 2ify the data.
             omicsData$e_data[, -iCol] <- log(10^(omicsData$e_data[, -iCol]))
             
           } else if (data_scale == "log2") {
             
             # Log base 10ify the data.
             omicsData$e_data[, -iCol] <- log2(10^(omicsData$e_data[, -iCol]))
             
           }
           
         }
         
         )
  
  # Update data_scale in the data_info attribute.
  attr(omicsData, 'data_info')$data_scale <- data_scale
  
  # Return the transmogrified omics object along with its attributes (some of
  # them updated and others left alone).
  return (omicsData)

}
