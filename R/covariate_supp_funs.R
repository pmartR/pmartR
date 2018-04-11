#Function to create an X matrix based on a covariate data frame
build_x_mat <- function(cov_df){
  
  if(is.null(ncol(cov_df))){
    cov_df <- matrix(cov_df,ncol=1)
  }
  
  #If all covariates are numeric, simply return the same matrix back
  if(is.numeric(cov_df)){
    return(data.matrix(cov_df))
  }
  
  #If the covariates are a mix of numeric, factors and characters, return matrix of group identifiers
  Xmatrix <- NULL
  for(i in 1:ncol(cov_df)){
    
    if(is.numeric(cov_df[,i])){
      #If column i is numeric, append it to X
      Xmatrix <- cbind(Xmatrix,cov_df[,i])
    }else{
      coli_levels <- unique(cov_df[,i])
      n_levels <- length(coli_levels)
      if(n_levels!=length(cov_df[,i])){
        Xcoli <- matrix(0,nrow(cov_df),n_levels)
        
        for(j in 1:n_levels){
          Xcoli[cov_df[,i]==coli_levels[j],j] <- 1
        }
        Xmatrix <- cbind(Xmatrix,Xcoli)
      }
    }
  }
  return(Xmatrix)
}

#---Function to remove linearly dependent columns in x corresponding to group levels----####
reduce_xmatrix <- function(x,ngroups){
  p <- ncol(x)
  orig_rank <- qr(x)$rank
  if(orig_rank<p){
    stdx <- apply(x,2,function(mat) (mat-mean(mat))/sd(mat))
    dmat <- as.matrix(dist(t(stdx)))
    dmat <- dmat[(1:ngroups),-(1:ngroups)]
    if(any(dmat==0)){
      ind_mat <- matrix(1:ngroups,ngroups,(p-ngroups))
      cols_to_remove <- ind_mat[which(dmat==0)]
      x[,cols_to_remove] <- 0
    }
  }
  return(x)
}

##------Function to construct the null space projection matrix------##

proj_mat <- function(X){
  #If the X matrix has atleast two rows, find projection matrix 
  #into the null space corresponding to X
  if(!is.null(nrow(X))){
    Imat <- diag(1,nrow(X))
    Px <- X%*%MASS::ginv(t(X)%*%X)%*%t(X)
    return(Imat-Px)
  }
  
  #If only one datapoint is left then just return 1
  return(1)
}


##-----Project each row of the data.matrix into X's null space-----###
project_to_null <- function(data_mat, Xmatrix){
  
  data_no_x <- data_mat
  
  for(i in 1:nrow(data_mat)){
    
    to_rem <- which(is.na(data_mat[i,]))
    if(length(to_rem)>0){
      roi <-  data_mat[i,-to_rem]
      IPxi <- proj_mat(Xmatrix[-to_rem,])
      data_no_x[i,-to_rem] <- IPxi%*%matrix(roi,ncol=1)
    }else{
      roi <- data_mat[i,]
      IPxi <- proj_mat(Xmatrix)
      data_no_x[i,] <- IPxi%*%matrix(roi,ncol=1)
    }

  }
  
  return(data.frame(data_no_x))
  
}