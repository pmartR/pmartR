#volcano plot function, takes in a statRes object and comparison_ind(is an index or indices of 'attr(statRes, "comparisons")' vector)

missingval_volcanoplot<- function(statRes, comparison_ind){

#extract comparison attr from statRes object
comparisons = attr(statRes, "comparisons")

if(!(comparison_ind %in% c(1:length(comparisons)))) stop("comparison_ind must be valid index of comparisons attr of statRes object")

#lets check the len of compare
if(length(comparison_ind) == 1){
  
  #extract comparisons that correspond to comparison_ind
  compare = comparisons[comparison_ind]
  
  #forming pvalue column name based on compare
  pval_col_name = paste("P_value_T", compare, sep = "_")
  
  #lets check if names(statRes$Full_results) contains pval_col_name 
  if(!(pval_col_name %in% colnames(statRes$Full_results))) stop(paste(pval_col_name, "was not found in statRes object", sep = " "))
  
  #forming fold change column name based on compare
  fold_change_col_name = paste("Fold_change", compare, sep ="_") 
  
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
  
  p <- ggplot(plotdata, aes(x = fold_change_data, y = pvalue_data)) + geom_point(color = "blue")
  return(p)
}


  
}




