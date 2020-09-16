####
#These two functions perform the IMD part of IMD ANOVA.  I include these only
#to be able to determine how much slower a pure R implementation compares to that of 
#an Rcpp implementation

imd_test_R <- function(omicsData){
  # check that omicsData is of the appropriate class
  if(!inherits(omicsData, c("proData","pepData","lipidData", "metabData", "nmrData"))) stop("omicsData is not an object of appropriate class")
  
  # Check for group_DF attribute #
  if(is.null(attr(omicsData, "group_DF"))){
    stop("group_designation must be called in order to create a 'group_DF' attribute for omicsData.")
  }else{
    groupData <- attr(omicsData, "group_DF")
  }
  
  #Catch if number of groups is too small
  k <- length(unique(groupData$Group))
  if(k<2){
    stop("This test doesn't make sense for less than two groups.")
  }
  
  #Create a data matrix that only includes the samples in "groupData" and put them in the order specified by
  #"groupData"
  data <- omicsData$e_data[,as.character(groupData$sample_name)]
  
  #Number of samples associated with each group
  nK <- as.numeric(table(groupData$Group))
  
  #Find the number of times each peptide was absent for each of the experimental groups, matrix of CAk values
  gp <- factor(groupData$Group,labels=1:k,levels=unique(groupData$Group))
  #absent <- count_missing_cpp(data.matrix(data),gp)
  absent <-  t(apply(data,1,count_missing,gp=groupData$Group))
  
  
  #Number of times each peptide was observed for each of the experimental groups, matrix of COk values
  observed <- matrix(nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)
  observed <- observed-absent
  
  #compute total number of samples from which each peptide is absent (mA) or observed (mO)
  mA <- rowSums(absent)
  mO <- rowSums(observed)
  N <- mA+mO
  
  #Expected number of absences (and observations) for each group
  EAk <- (mA/N)*matrix(nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)
  EOk <- (mO/N)*matrix(nK,nrow=nrow(absent),ncol=ncol(absent),byrow=TRUE)
  
  #Compute test statistic but summing over all nonzero counts (na.rm=TRUE removes the zero counts)
  Gk <- 2*rowSums(observed*log(observed/EOk)+absent*log(absent/EAk), na.rm=TRUE)
  
  #Put the statistic, degrees of freedom and p-value into a data frame
  results <- data.frame(Statistic=Gk, DF=(k-1), p_value=pchisq(Gk,df=k-1,lower.tail = FALSE))
  
  #Pull off edata_cname and add to results df
  edatacname <- attr(omicsData,"cnames")$edata_cname
  results <- cbind(omicsData$e_data[edatacname],results)
  return(results)
}




#Counts the number of missing observations in each group for each row of x
count_missing <- function(x,gp){
  unq_gps <- unique(gp)
  res <- rep(0,length(unq_gps))
  for(i in 1:length(unq_gps)){
    res[i] <- sum(is.na(x[gp==unq_gps[i]]))
  }
  return(res)
}

#imd_test_R(omicsData = filtered_pep_data, groupData = group_df)
#microbenchmark(imd_test_R(omicsData = filtered_pep_data, groupData = group_df),imd_test(omicsData = filtered_pep_data, groupData = group_df))
