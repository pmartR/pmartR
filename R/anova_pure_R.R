####
#These two functions perform the ANOVA part of IMD ANOVA.  I include these only
#to be able to determine how much slower a pure R implementation compares to that of 
#an Rcpp implementation

anova_test_R <- function(omicsData, groupData, return_sizes = FALSE){
  
  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("At least two groups are necessary to perform an ANOVA.")
  }
  
  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData$sample_name)]
  
  
  gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))
  results <- data.frame(t(apply(data,1,mean_sse_fun,gp=groupData$Group)))
  #The C++ function returns a list so we need to make it into a data.frame
  #results <- data.frame(raw_results$group_means, raw_results$Sigma2, raw_results$Fstats, raw_results$pvalue)
  
  #Rename the columns to match group names
  group_names <- paste("Group",as.character(unique(groupData$Group)))
  #colnames(results) <- c(group_names,"Variance","F-Statistic","p-value")
  
  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  results <- cbind(omicsData$e_data[edatacname],results)
  
  if(return_sizes){
    sizes_df <- raw_results$group_sizes
    colnames(sizes_df) <- group_names
    return(list(Results=results,Sizes=sizes_df))
  }
  
  return(results)
}


mean_sse_fun <- function(x,gp){
  #Compute groups means, ANOVA F-stat and p-value   
  groups <- unique(gp)
  n <- sum(!is.na(x))
  ngroups <- length(groups)
  means <- rep(NA,ngroups)
  SSB <- SSE <- 0
  
  non_missing_groups <- length(unique(gp[!is.na(x)]))
  
  for(i in 1:ngroups){
    means[i] <- mean(x[gp==groups[i]],na.rm=T)
    if(!is.na(means[i])){
      SSE <- SSE+sum((x[gp==groups[i]]-means[i])^2,na.rm=T)
    }
  }
  SSB <- sum((x-mean(x,na.rm=TRUE))^2,na.rm=TRUE)-SSE
  
  FStat <- (SSB/(non_missing_groups-1))/(SSE/(n-non_missing_groups))
  if(!is.na(FStat)){
    pval <- pf(FStat,df1=(non_missing_groups-1),df2=(n-non_missing_groups),lower.tail = FALSE)
  }else{
    pval <- NA
  }
  
  return(c(means,FStat,pval))
}
