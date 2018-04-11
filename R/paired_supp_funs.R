#Function to translate paired data into differences
build_pair_mat <- function(pair_df){
  #pair_df - data.frame with two columns: sampleIDs and PairIDs (in that order)
  
  #Creat matrix of 1, -1 that will turn single observations into paired differences
  coli_levels <- unique(pair_df[,2])
  n_levels <- length(coli_levels)
  Xmatrix <- matrix(0,nrow(pair_df),n_levels)
        
  for(j in 1:n_levels){
    if(length(which(pair_df[,2]==coli_levels[j]))>2){
      stop(paste("Only two samples can be associated with the same 'pair'.  More than two samples are associated with pair",coli_levels[j],"."))
    }
    Xmatrix[pair_df[,2]==coli_levels[j],j] <- c(1,-1)
  }

  return(Xmatrix)
}
