#' Test normalization
#' 
#' Tests normalization results of a normRes object
#' 
#' @param normRes_object an object of class normRes, produced by normalize function
#' @param parametric logical indicating what test to use (anova test or Mann-Whitney U test), defaults to TRUE(anova test is used).
#' 
#' @return an object of class normRes, but with an additional list item 'test'
#' 
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data(lipid_object)
#' lipid_object<- edata_transform(omicsData = lipid_object, data_scale="log2")
#' lipid_object<- group_designation(omicsData = lipid_object, main_effects = "Condition")
#' norm_object<- normalize(omicsData = lipid_object, subset_fn = "all", norm_fn = "median")
#' 
#' test_result<- test_normalization(normRes_object = norm_object, parametric = TRUE)
#'}
#'
#'  @export  
#'  

test_normalization<- function(normRes_object, parametric = TRUE){
  #pulling omicsData from normRes_object
  omicsData = attr(normRes_object, "omicsData")
  group_DF = attr(omicsData, "group_DF")
  
  #case where scale is NULL and location is non-NULL
  if(is.null(normRes_object$parameters$normalization$scale) & !is.null(normRes_object$parameters$normalization$location)){
    location = normRes_object$parameters$normalization$location
    location = data.frame("Sample_Name" = names(location), "location" = location, stringsAsFactors = FALSE)
    row.names(location) = NULL
    
    my_data = merge(location, group_DF, by = "Sample_Name")
    
    if(parametric == TRUE){
      res_anova = aov(location ~ Group, data = my_data)
      sum_anova = summary(res_anova)
      
      result = data.frame("Test" = "anova", "Test Statistic" = sum_anova[[1]]$`F value`[1], "P-value" = sum_anova[[1]]$`Pr(>F)`[1], check.names = F, stringsAsFactors = F)
      
    }else{
      res_wilcox = wilcox.test(location ~ Group, data = my_data)
      result = data.frame("Test" = "Mann Whittney", "Test Statistic" = res_wilcox$statistic, "P-value" = res_wilcox$p.value, check.names = F, stringsAsFactors = F)
    }
    
    normRes_object$test = list("scale" = NULL, "location" = result) 
  }
  
  #case where location is NULL and scale is non-NULL
  else if(is.null(normRes_object$parameters$normalization$location) & !is.null(normRes_object$parameters$normalization$scale)){
    scale = normRes_object$parameters$normalization$scale
    scale = data.frame("Sample_Name" = names(scale), "scale" = scale, stringsAsFactors = FALSE)
    row.names(scale) = NULL
    
    my_data = merge(scale, group_DF, by = "Sample_Name")
    
    if(parametric == TRUE){
      res_anova = aov(scale ~ Group, data = my_data)
      sum_anova = summary(res_anova)
      
      result = data.frame("Test" = "anova", "Test Statistic" = sum_anova[[1]]$`F value`[1], "P-value" = sum_anova[[1]]$`Pr(>F)`[1], check.names = F, stringsAsFactors = F)
      
    }else{
      res_wilcox = wilcox.test(scale ~ Group, data = my_data)
      result = data.frame("Test" = "Mann Whittney", "Test Statistic" = res_wilcox$statistic, "P-value" = res_wilcox$p.value, check.names = F, stringsAsFactors = F)
      
    }
  
    normRes_object$test = list("scale" = result, "location" = NULL) 
  }
  
  #case where both location is non-NULL and scale is non-NULL
  else if(!is.null(normRes_object$parameters$normalization$location) & !is.null(normRes_object$parameters$normalization$scale)){
    location = normRes_object$parameters$normalization$location
    location = data.frame("Sample_Name" = names(location), "location" = location, stringsAsFactors = FALSE)
    row.names(location) = NULL
    
    my_data_location = merge(location, group_DF, by = "Sample_Name")
    
    scale = normRes_object$parameters$normalization$scale
    scale = data.frame("Sample_Name" = names(scale), "scale" = scale, stringsAsFactors = FALSE)
    row.names(scale) = NULL
    
    my_data_scale = merge(scale, group_DF, by = "Sample_Name")
    
    if(parametric == TRUE){
      res_anova_location = aov(location ~ Group, data = my_data_location)
      sum_anova_location = summary(res_anova_location)
      
      result_location = data.frame("Test" = "anova", "Test Statistic" = sum_anova_location[[1]]$`F value`[1], "P-value" = sum_anova_location[[1]]$`Pr(>F)`[1], check.names = F, stringsAsFactors = F)
      
      res_anova_scale = aov(scale ~ Group, data = my_data_scale)
      sum_anova_scale = summary(res_anova_scale)
      
      result_scale = data.frame("Test" = "anova", "Test Statistic" = sum_anova_scale[[1]]$`F value`[1], "P-value" = sum_anova_scale[[1]]$`Pr(>F)`[1], check.names = F, stringsAsFactors = F)
      
    }else{
      res_wilcox_location = wilcox.test(location ~ Group, data = my_data_location)
      result_location = data.frame("Test" = "Mann Whittney", "Test Statistic" = res_wilcox_location$statistic, "P-value" = res_wilcox_location$p.value, check.names = F, stringsAsFactors = F)
      
      res_wilcox_scale = wilcox.test(scale ~ Group, data = my_data_scale)
      result_scale = data.frame("Test" = "Mann Whittney", "Test Statistic" = res_wilcox_scale$statistic, "P-value" = res_wilcox_scale$p.value, check.names = F, stringsAsFactors = F)
    }
    
    normRes_object$test = list("scale" = result_scale, "location" = result_location)
  }
  return(normRes_object)
}
