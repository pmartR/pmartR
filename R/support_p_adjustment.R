#I haven't found a function that quickly gets the rank of the original p-values so do it in R for now,
#still much faster than a pure R version
ranked_holm_cpp <- function(ps){
  return(holm_cpp(ps)[rank(ps)])
}