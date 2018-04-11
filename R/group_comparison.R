
##----------------------------------------------------------------##
##-------- Do the group comparisons for the anova test   ---------##
##----------------------------------------------------------------##

group_comparison_anova <- function(groupData,comparisons,anova_results_full){
  #Takes the results of anova_test() and returns group comparison p-values
  
  #groupData - data frame that assigns sample names to groups
  #comparisons - dataframe that defiens the comparsions of interest
  #anova_results - the results of the MSomicsSTAT::anova_test() function
  
  #The group means include the word "Group" in them so that is what will be passed 
  #to the group_comparison(...) and fold_change(...) functions along with estimated variance
  anova_results <- anova_results_full$Results
  means <- data.matrix(anova_results[,grep("Mean",colnames(anova_results))])
  sigma2 <- data.matrix(anova_results$Variance)
  sizes <- data.matrix(anova_results_full$Sizes)
  #rm(anova_results_full)

  #If "comparisons" is null, then do all pairwise comparisons, else use the "comparisons" df to 
  #create the appropriate C matrix
  
  ##---- I used to deal with one vs two factor differently, that code is at the end of this script -----##
  Cmat_res <- create_c_matrix(group_df = groupData, to_compare_df = comparisons)
  
  Cmat <- Cmat_res$cmat
  group_comp <- group_comparison_anova_cpp(means, sizes, sigma2, Cmat)
  
  group_comp$diff_mat <- data.frame(group_comp$diff_mat)
  colnames(group_comp$diff_mat) <- Cmat_res$names
  colnames(group_comp$diff_ses) <- Cmat_res$names
  colnames(group_comp$t_tests) <- paste(Cmat_res$names,"Test-Stat")
  group_comp$p_values <- data.frame(group_comp$p_values)
  colnames(group_comp$p_values) <- paste(Cmat_res$names,"p-value")
  return(group_comp)
}

##----------------------------------------------------------------##
##-------- Do the group comparisons for the imd G-test   ---------##
##----------------------------------------------------------------##

group_comparison_imd <- function(groupData,comparisons,observed,absent){
  #Takes the results of anova_test() and returns group comparison p-values
  
  #groupData - data frame that assigns sample names to groups
  #comparisons - dataframe that defiens the comparsions of interest
  #observed - matrix of number of observed counts
  #absent - matrix of number of observed counts
  

  #If "comparisons" is null, then do all pairwise comparisons, else use the "comparisons" df to 
  #create the appropriate C matrix
  Cmat_res <- create_c_matrix(groupData, comparisons)
  Cmat <- abs(Cmat_res$cmat)
  
  #Number of samples per group
  nK <- observed[1,]+absent[1,]
  
  moMat <- observed%*%t(Cmat)
  maMat <- absent%*%t(Cmat)
  N <- moMat[1,]+maMat[1,]
  
  #CokMat <- CakMat <- EokMat <- EakMat <- NULL
  Gstats <- matrix(0,ncol=nrow(Cmat),nrow=nrow(moMat))
  for(i in 1:ncol(moMat)){
    nk_N_ratio <- nK/N[i]*Cmat[i,]
    
    EokMat <- matrix(moMat[,i],ncol=1)%*%matrix(nk_N_ratio,nrow=1)
    EakMat <- matrix(maMat[,i],ncol=1)%*%matrix(nk_N_ratio,nrow=1)
    CokMat <- observed*matrix(Cmat[i,],nrow=nrow(observed),ncol=ncol(Cmat),byrow=TRUE)
    CakMat <- absent*matrix(Cmat[i,],nrow=nrow(observed),ncol=ncol(Cmat),byrow=TRUE)
    GmatObs <- CokMat*log(CokMat/EokMat)
    GmatAbs <- CakMat*log(CakMat/EakMat)
    GmatObs[is.nan(GmatObs)] <- GmatAbs[is.nan(GmatAbs)] <- 0
    Gstats[,i] <- 2*rowSums(GmatObs+GmatAbs)
  }
  
  colnames(Gstats) <- Cmat_res$names
  pvals <- pchisq(Gstats, 1, lower.tail = FALSE)

  #Return signs for imd_test flags, i.e., diff in proportion of observed values
  signs <- sign((observed/(observed+absent))%*%t(Cmat_res$cmat))
  
  return(list(GStats=Gstats,pvalues=pvals,signs=signs))
}

##-----------------------------------------------------------------------------------##
##-------- Function to create the C matrix that will do all the comparisons ---------##
##-----------------------------------------------------------------------------------##

create_c_matrix <- function(group_df,to_compare_df=NULL){
  #This function will create the C-matrix that defines the comparisons that are requested
  #and returns a list of names that defines which two groups are being compared for each 
  #row of the C matrix 
  
  #group_df - defines which columns (`sample_name`) belong to which groups (`Group`)
  groups_two_fac <- groups <- as.character(unique(group_df$Group))
  ngroups <- length(groups)
  compi <- NULL
  
  if(is.null(to_compare_df)){
    #####
    #Create the Cmatrix that corresponds to all pairwise comparisons
    #With two factors, all pairwise includes the higher order terms
    #THIS COULD PROBABLY BE MADE MUCH FASTER IN C++
    n_comparisons <- choose(ngroups,2)
    Cmat <- matrix(0,n_comparisons,ngroups)
    row <- 1
    for(i in 1:(ngroups-1)){
      for(j in (i+1):ngroups){
        Cmat[row,i] <- 1
        Cmat[row,j] <- -1
        row <- row+1
        compi <- c(compi,paste0(groups[i],"_vs_",groups[j]))
      }
    }
  }else{
    n_comparisons <- nrow(to_compare_df)
    Cmat <- matrix(0,n_comparisons,ngroups)
    all_comp_levels <- as.character(unique(unlist(to_compare_df)))
    
    if(ncol(group_df)>2){#Two factor case, allow for main effect tests so grow groups to include those
      group_df[,-1] <- apply(group_df[,-1],2,as.factor)
      groups_two_fac <- as.character(unique(unlist(group_df[,-1])))
    }
    
    #Check that all of the groups in "to_compare_df" are present in "group_df"
    if(!all(all_comp_levels%in%groups_two_fac) & !all(all_comp_levels%in%groups)){
      stop("Some groups listed in the 'comparisons' argument do not appear in the 'groupData' argument.")
    }

    for(i in 1:n_comparisons){
      control_i <- grep(to_compare_df$Control[i],groups)
      test_i <- grep(to_compare_df$Test[i],groups)
      Cmat[i,control_i] <- (-1/length(control_i))
      Cmat[i,test_i] <- 1/length(test_i)
      compi <- c(compi,paste0(to_compare_df$Test[i],"_vs_",to_compare_df$Control[i]))
    }
    
  }
  
  return(list(cmat=Cmat, names=compi))
}


##----------------------------------------------------------------------------##
##------ This is how I used to deal with one and two factor comaprisons ------##
##----------------------------------------------------------------------------##
# if(ncol(groupData)==2){
#   #One factor anova
#   Cmat_res <- create_c_matrix(groupData, comparisons)
# }else{
#   #Two factor anova, create C matrix for each factor
#   groupData_fac1 <- groupData[,c(1,3)]
#   groupData_fac2 <- groupData[,c(1,4)]
#   colnames(groupData_fac1)[2] <- colnames(groupData_fac2)[2] <- c("Group")
#   Cmat_res1 <- create_c_matrix(groupData_fac1, comparisons)
#   Cmat_res2 <- create_c_matrix(groupData_fac2, comparisons)
#   Cmat_res <- list(cmat=Cmat_res1$cmat,names=c(Cmat_res1$names,Cmat_res2$names))
#   Cmat_res$cmat <- cbind(Cmat_res$cmat,matrix(0,nrow=nrow(Cmat_res1$cmat),ncol=ncol(Cmat_res2$cmat)))
#   Cmat_res$cmat <- rbind(Cmat_res$cmat,cbind(matrix(0,nrow=nrow(Cmat_res2$cmat),ncol=ncol(Cmat_res1$cmat)),Cmat_res2$cmat))
# }