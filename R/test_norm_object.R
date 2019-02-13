#' Test the location and scale parameters from a normalization procedure
#' 
#' Computes p-values from a test of dependence between normalization parameters and group assignment of a normalized omicsData or normRes object.
#' 
#' @param norm_object An object of class 'pepData', 'proData', 'lipidData', 'metabData', 'isobaricpepData' that has had \code{normalize_global()} run on it, or a 'normRes' object
#' @param test_fn The statistical test to use, current tests are "anova" and "kw" for a Kruskal-Wallis test.
#' @param group_vec character vector of the same length as the number of normalization parameters (number of samples)
#'
#' @return A list with 2 entries containing the p_value of the test performed on the location and scale(if they exist) parameters.
#' 
#' @export
#' 

test_normRes <- function(norm_object, test_fn = "kw", group_vec = NULL){
  obj_classes <-  c("pepData", "proData", "lipidData", "metabData", "isobaricpepData")
  
  # object is of correct class and group_vec is a character vector
  if(!inherits(norm_object, c(obj_classes, "normRes"))) stop("This object is not of class 'pepData', 'proData', 'lipidData', 'metabData', 'isobaricpepData' or 'normRes'")
  if(!is.null(group_vec) & !inherits(group_vec, "character")) stop("group_vec must be a character vector")
  
  # if this is an omicsData object, test that is has been normalized and has a group DF...
  # ...if so store the location and scale vectors that should be under $data_info$norm_info and get the group vector from group DF is none was specified
  
  # if its a normRes object, check that there is a group_DF inside the 'omicsData' attribute or a group_vec was specified, then store vectors etc.
  if(inherits(norm_object, obj_classes)){
    if(is.null(attributes(norm_object)$data_info$norm_info$is_normalized)) stop("Normalization has not been run on this data")
    if(is.null(attributes(norm_object)$group_DF) & is.null(group_vec)) stop("No grouping structure present in object or specified in 'group_vec' argument")
    
    location <- attributes(norm_object)$data_info$norm_info$params$norm_location
    scale <- attributes(norm_object)$data_info$norm_info$params$norm_scale
    if(!is.null(attributes(norm_object)$group_DF) & is.null(group_vec)) group_vec <- attributes(norm_object)$group_DF$Group %>% as.character()
  }
  else if(inherits(norm_object, "normRes")){
   if(is.null(attr(attr(norm_object, "omicsData"), "group_DF")) & is.null(group_vec)) stop("No grouping structure present in object or specified in 'group_vec' argument")
   
   location <- norm_object$parameters$normalization$location
   scale <- norm_object$parameters$normalization$scale
   if(!is.null(attr(attr(norm_object, "omicsData"), "group_DF")) & is.null(group_vec)) group_vec <- attr(attr(norm_object, "omicsData"), "group_DF")$Group %>% as.character()
  }
  
  # check that the final group vector is valid
  if(length(location) != length(group_vec)) stop("group vector must have the same number of entries as the number of normalization parameters")
  if(!is.null(scale)){
    if(length(scale) != length(group_vec)) stop("group vector must have the same number of entries as the number of normalization parameters")
  }
  
  # get p values using the value and group vectors and selected test
  if(test_fn == "kw"){
    p_location <- kw_rcpp(matrix(location, nrow = 1), group_vec)
    p_scale <- if(!is.null(scale)) kw_rcpp(matrix(scale, nrow = 1), group_vec) else NULL
  }
  else if(test_fn == "anova"){
    p_location = aov(location ~ group_vec) %>% summary() %>% {.[[1]]$`Pr(>F)`[1]}
    p_scale = if(!is.null(scale)) aov(scale ~ group_vec) %>% summary() %>% {.[[1]]$`Pr(>F)`[1]} else NULL
  }
  
  res <- list(p_location = p_location, p_scale = p_scale)
  
  return(res)
  
}