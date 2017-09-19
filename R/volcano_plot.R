#' missingval_volcanoplot
#' 
#' fold change vs t-Test p-value 
#' 
#'@param comparison can either be a character string name of the comparison to plot, or an integer index refering to the comparisons attribute vector
#' 
#'@rdname missingval_volcanoplot
#'
#' \tabular{ll}{
#' \code{x_lab} \tab character string to be used for x-axis label. Defaults to NULL \cr
#' \code{y_lab} \tab character string to be used for y-axis label. Defaults to NULL \cr
#' \code{title_plot} \tab character string to be used for the plot title. Defaults to NULL. \cr
#' \code{title_size} \tab integer value specifying the font size for the plot title. Default is 14. \cr
#' \code{x_lab_size} \tab integer value indicating the font size for the x-axis. Defaults to 11. \cr
#' \code{y_lab_size} \tab integer value indicating the font size for the y-axis. Defaults to 11. \cr
#' \code{bw_theme} \tab logical indicator of whether to use the "theme_bw". Defaults to FALSE, in which case the ggplot2 default theme is used. \cr
#' }
#'
#'@export
#'

missingval_volcanoplot<- function(statRes, comparison, x_lab = NULL, ...) {
  .missingval_volcanoplot(statRes, comparison, x_lab, ...)
}

.missingval_volcanoplot<- function(statRes, comparison, x_lab = NULL, y_lab = NULL, title_plot = NULL, title_size = 14, x_lab_size = 11, y_lab_size = 11, bw_theme = FALSE){

#check that staRes object is of 'statRes' class
if(class(statRes) != "statRes") stop("object must be of class 'statRes'")

#check that comparison is of correct class, either integer or character
if(!(class(comparison)) %in% c("character", "numeric")) stop("comparison must be either of 'character' or 'numeric' class")  

#plot label checks
if(!is.null(title_plot)) {
  if(!is.character(title_plot)) stop("title_plot must be a character vector")
}
if(!is.null(x_lab)) {
  if(!is.character(x_lab)) stop("x_lab must be a character vector")
}
if(!is.null(y_lab)) {
  if(!is.character(y_lab)) stop("y_lab must be a character vector")
}  
  
#make sure comparison is in attr(statRes, "comparisons")  
if(class(comparison) == "character"){
  if(!(comparison %in% attr(statRes, "comparisons"))) stop(paste(comparison, "was not found in statRes attributes", sep = " ")) 
}

#check that when comparison is numeric it is a valid index  
if(class(comparison) == "numeric"){
  if(!(comparison %in% c(1:length(attr(statRes, "comparisons"))))) stop("comparison must be valid index of comparisons attr of statRes object")
  
  #now we take comparison and re-assign it a character string instead of a numeric
  comparison = attr(statRes, "comparisons")[comparison]
}    

#lets check the len of comparison
if(length(comparison) == 1){
  
  #forming pvalue column name based on comparison
  pval_col_name = paste("P_value_T", comparison, sep = "_")
  
  #lets check if names(statRes$Full_results) contains pval_col_name 
  if(!(pval_col_name %in% colnames(statRes$Full_results))) stop(paste(pval_col_name, "was not found in statRes object", sep = " "))
  
  #forming fold change column name based on comparison
  fold_change_col_name = paste("Fold_change", comparison, sep ="_") 
  
  #lets check if names(statRes$Full_results) contains fold_change_col_name 
  if(!(fold_change_col_name %in% colnames(statRes$Full_results))) stop(paste(fold_change_col_name, "was not found in statRes object", sep = " "))
  
  #now lets look for pval_col_name in statRes$Full_results
  pval_ind = which(names(statRes$Full_results) %in% pval_col_name)
  
  #looking for fold_change_col_name in statRes$Full_results
  fold_change_ind = which(names(statRes$Full_results) %in% fold_change_col_name)

  #now lets subset the pvalue col and fold change col from Full_results data frame 
  pvalue_data = statRes$Full_results[,pval_ind]
  fold_change_data = statRes$Full_results[,fold_change_ind]
  
  #cbind data
  plotdata = cbind(pvalue_data, fold_change_data)
  plotdata = as.data.frame(plotdata)

  #now we can apply -log10() to the pvalue_data
  plotdata$pvalue_data = -log10(plotdata$pvalue_data)
  #take log2 of fold chang aswell
  
  # make labels #
  xlabel <- ifelse(is.null(x_lab), "Fold Change", x_lab)
  ylabel <- ifelse(is.null(y_lab), "-log10 t-Test P-value", y_lab)
  plot_title <- ifelse(is.null(title_plot), paste("Volcano Plot", comparison, sep = " "), title_plot)
  
  if(bw_theme == FALSE){
    p <- ggplot(plotdata, aes(x = fold_change_data, y = pvalue_data)) + geom_point(color = "blue") +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size))
  }
  else{
    p <- ggplot(plotdata, aes(x = fold_change_data, y = pvalue_data)) + geom_point(color = "blue") +
      ggplot2::xlab(xlabel) +
      ggplot2::ylab(ylabel) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size), axis.title.x = ggplot2::element_text(size = x_lab_size), axis.title.y = ggplot2::element_text(size = y_lab_size)) +
      ggplot2::theme_bw()
  }
    
  return(p)
}
  
}




