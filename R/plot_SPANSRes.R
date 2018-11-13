#' Plots an object of class SPANSRes
#' 
#' For plotting an S3 object of type 'SPANSRes'
#'
#' @param SPANSRes_obj an object of the class 'SPANSRes', usually created by \code{\link{spans_procedure()}}.
#' @return plots a plotly object
#' 
#' @rdname plot-SPANSRes
#' @export

plot.SPANSRes <- function(SPANSRes_obj){
  SPANSRes_obj <- SPANSRes_obj %>% 
    dplyr::mutate(ss_par = paste(subset_method, parameters, sep = " | "), SPANS_score = as.numeric(SPANS_score)) %>%
    dplyr::left_join(attr(SPANSRes_obj, "method_selection_pvals"), by = c("subset_method", "normalization_method", "parameters"))
  
  
  p <- plotly::plot_ly(SPANSRes_obj, x = ~normalization_method, y = ~ss_par, z = ~SPANS_score,
                       hoverinfo = 'text', text = ~paste('</br> F(-Log10(HSmPV)):', F_log_HSmPV,
                                                         '</br> F(Log10(NSmPV)): ', F_log_NSmPV,
                                                         '</br> Scale p-value: ', scale_p_value,
                                                         '</br> Location p-value', location_p_value)) %>%
    plotly::add_heatmap() %>% 
    plotly::colorbar(title = "SPANS score") %>% 
    plotly::layout(plot_bgcolor = 'black', xaxis = list(title = "Normalization Method"), yaxis = list(title = "Subset Method"))
  
  return(p)
}