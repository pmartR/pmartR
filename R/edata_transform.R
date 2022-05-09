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
           
         },
         
         ## seqData only
         'counts' = {
           
           temp_data <- omicsData$e_data[, -iCol]
           
           if (data_scale == 'lcpm'){
             
             ## log cpm, limma voom method and a similar method used for visualizations in edgeR
             
             ## EdgeR
             # First scales the prior.count/pseudo-count and adds 2x the scaled prior count to the libsize
             # prior.count.scaled <- lib.size/mean(lib.size)*prior.count
             # lib.size <- lib.size+2*prior.count.scaled
             # lib.size <- 1e-6*lib.size
             # Calculates log2 log2(t( (t(x)+prior.count.scaled) / lib.size ))
             
             # sum library size
             samp_sum <- apply(temp_data, 
                               2, 
                               sum,
                               na.rm = TRUE) + 1
             
             # divide adjusted (ensure non-zero) counts by library size
             div_sum <- sweep((temp_data + .5), 2, samp_sum, `/`)
             
             # Apply per million multiplier and log2
             omicsData$e_data[, -iCol] <- log2(div_sum * 10^6)
             
           } else if (data_scale == 'upper'){
             
             warning("Zeros will be regarded as NA for 'upper' transformation")
             
             temp_data[temp_data == 0] <- NA
             
             # Grab non-zero upper quantile of data
             samp_upper <- apply(temp_data, 
                                 2, 
                                 quantile,
                                 na.rm = TRUE,
                                 probs = .75)
             
             g.q <- quantile(unlist(temp_data), probs = .75, na.rm=TRUE)
             
             # Divide each count by the upper quantile in respective columns
             div_75 <- sweep(temp_data, 2, samp_upper, `/`)
             
             # Set new data
             omicsData$e_data[, -iCol] <- div_75*g.q # back transform
             
           } else if(data_scale == 'median'){
             
             warning("Zeros will be regarded as NA for 'median' transformation")
             
             temp_data[temp_data == 0] <- NA
             
             # Grab non-zero median of data
             samp_med <- apply(temp_data, 
                               2, 
                               median,
                               na.rm = TRUE)
             
             # Divide each count by the upper quantile in respective columns
             div_med <- sweep(temp_data, 2, samp_med, `/`)
             
             g.q <- median(unlist(temp_data), na.rm=TRUE)
             
             # Set new data
             omicsData$e_data[, -iCol] <- div_med*g.q
             
           } 
           # else if(data_scale == "vst"){
           # 
           #   ### Requires Grouping information ###
           #   if(is.null(get_group_DF(omicsData))) stop(
           #     "iqlr requires group_designation to calculate within group variance")
           # 
           #   edata_deseq <- DESeqDataSetFromMatrix(
           #     omicsData$e_data[colnames(omicsData$e_data) != get_edata_cname(omicsData)],
           #     colData = get_group_DF(omicsData)[-1],
           #     design = ~Tissue + Treatment ## Figure out design designation
           #   )
           # 
           #  # DEseq method
           # res <- DESeq2::vst(edata_deseq)
           # omicsData$e_data[, -iCol] <- SummarizedExperiment::assay(res)
           # 
           # } else if(data_scale == "rlog"){
           # 
           #   ### Requires Grouping information ###
           #   if(is.null(get_group_DF(omicsData))) stop(
           #     "iqlr requires group_designation to calculate within group variance")
           #   
           #   edata_deseq <- DESeqDataSetFromMatrix(
           #     omicsData$e_data[colnames(omicsData$e_data) != get_edata_cname(omicsData)],
           #     colData = get_group_DF(omicsData)[-1],
           #     design = ~Tissue + Treatment ## Figure out design designation
           #   )
           # 
           #   # DEseq method
           #   res <- DESeq2::rlog(edata_deseq)
           #   omicsData$e_data[, -iCol] <- SummarizedExperiment::assay(res)
           # 
           # } else if(data_scale == "clr"){
           #   # Aldex2 method, centered log ratio
           #   # per https://github.com/ggloor/ALDEx2_dev/blob/master/ALDEx2/R/clr_function.r
           #   # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
           #   # i.e., do a centered logratio transformation as per Aitchison
           #   
           #   # sum library size
           #   samp_mean <- apply(temp_data, 
           #                     2, 
           #                     mean,
           #                     na.rm = TRUE)
           #   
           #   # divide adjusted (ensure non-zero) counts by library mean
           #   div_mean <- sweep((temp_data + .5), 2, samp_mean, `/`)
           #   
           #   # Apply per million multiplier and log2
           #   omicsData$e_data[, -iCol] <- log2(div_mean)
           #   
           #   
           # } else if(data_scale == "iqlr"){
           #   # Aldex2 method, inter quartile log ratio
           #   # https://www.bioconductor.org/packages/devel/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html
           #   
           #   ### Requires Grouping information ###
           #   if(is.null(get_group_DF(omicsData))) stop(
           #     "iqlr requires group_designation to calculate within group variance")
           #   
           #   # sum library size
           #   samp_mean <- apply(temp_data, 
           #                      2, 
           #                      mean,
           #                      na.rm = TRUE)
           #   
           #   # divide adjusted (ensure non-zero) counts by library mean
           #   div_mean <- sweep((temp_data + .5), 2, samp_mean, `/`)
           #   
           #   # Get log2 feature variance by group
           #   grps <- get_group_DF(omicsData)[1:2]
           #   values <- cbind(omicsData$e_data[iCol], log2(div_mean))
           #   
           #   grouped_data <- dplyr::left_join(
           #     reshape2::melt(values, id.var = get_edata_cname(omicsData)), 
           #     grps, 
           #     by = c("variable" = get_fdata_cname(omicsData)))
           #   
           #   group_cols <- vars(one_of(c("Group", get_edata_cname(omicsData))))
           #   
           #   all_vars <- grouped_data %>% 
           #     dplyr::group_by_at(group_cols) %>%
           #     summarise(variance = var(value, na.rm = T))
           #   
           #   # Use only rows from inner quantile
           #    cutoff25 <- quantile(all_vars$variance, .25)
           #    cutoff75 <- quantile(all_vars$variance, .75)
           #   
           #    all_vars <- ungroup(all_vars) %>% 
           #      group_by_at(vars(one_of(get_edata_cname(omicsData))))
           #    
           #    keep_rows <- all_vars %>% summarise(
           #      all_within = all(Reduce("&", list(variance > cutoff25,
           #                                    variance < cutoff75)))
           #    )
           #   
           #    ## Modify e_data for retained features
           #    keep_mols_edata <- which(omicsData$e_data[[iCol]] %in% 
           #                         keep_rows$ID_REF[keep_rows$all_within])
           #    keep_mols_emeta <- which(omicsData$e_meta[[iCol]] %in% 
           #                               keep_rows$ID_REF[keep_rows$all_within])
           #    
           #    
           #    omicsData$e_data <- omicsData$e_data[keep_mols_edata,]
           #    omicsData$e_meta <- omicsData$e_meta[keep_mols_emeta,]
           #    
           #    attr(omicsData, "data_info") <- pmartR:::set_data_info(
           #      e_data = omicsData$e_data,
           #      edata_cname = get_edata_cname(omicsData),
           #      data_scale_orig =get_data_scale_orig(omicsData),
           #      data_scale = get_data_scale(omicsData),
           #      data_types = c('RNA transcripts', 'as.seqData'),
           #      norm_info = get_data_info(omicsData)$norm_info,
           #      is_normalized = get_data_norm(omicsData)
           #      )
           #    
           #    # set meta data attributes #
           #    attr(omicsData, "meta_info") <- pmartR:::set_meta_info(
           #      e_meta = omicsData$e_meta,
           #      emeta_cname = get_emeta_cname(omicsData))
           #    
           #    ## Re-apply clr
           #    temp_data <- omicsData$e_data[, -iCol]
           #    
           #    # sum library size
           #    samp_mean <- apply(temp_data, 
           #                       2, 
           #                       mean,
           #                       na.rm = TRUE)
           #    
           #    # divide adjusted (ensure non-zero) counts by library mean
           #    div_mean <- sweep((temp_data + .5), 2, samp_mean, `/`)
           #    
           #    # Apply per million multiplier and log2
           #    omicsData$e_data[, -iCol] <- log2(div_mean)
           #   
           # }
         }
         
         )
  
  # Update data_scale in the data_info attribute.
  attr(omicsData, 'data_info')$data_scale <- data_scale
  
  # Return the transmogrified omics object along with its attributes (some of
  # them updated and others left alone).
  return (omicsData)

}
