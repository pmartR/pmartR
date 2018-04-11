#include <Rcpp.h>
using namespace Rcpp;

//For each row (peptide), this function computes the group means, counts the number of non-NA
//observations, estimates sigma^2, then computes the ANOVA F-statistic and associated p-value
//
//The results of this function are fed into group_comparison() which takes the group means,
//group sizes and sigma^2 value to do whatever group comparisons the user asks for

// [[Rcpp::export]]
List gtest_cpp(NumericMatrix data, NumericVector gp) {
  //data is the n-by-p matrix of data
  //gp is a vector that identifies what group each column belongs to, 
  //  assumed to be numeric 1-m where m is the total number of groups
  int i,j,k,groupi,rowi_size;
  int n = data.nrow();  //number of rows in data matrix, i.e., peptides/proteins/...
  int p = data.ncol();  //number of samples
  int m = max(gp);      //number of groups
  int nonmissing_m; //number of groups with at least one observation
  double overall_mean = 0; //mean across groups
  
  NumericVector Gstats(n), p_value(n); //F-statistic & p-value for each row
  NumericMatrix group_sums(n,m); //matrix to save the group sums in
  NumericMatrix group_means(n,m); //matrix to save the group means in
  NumericMatrix group_sizes(n,m);
  NumericVector rowi(m), Gstat_rowi(m);
  NumericVector rowi_gsize(m); //vector to contain the number of non-na observations
  //per group
  
  //Iterate over the matrix rows (peptides, proteins,...) to get group means
  //and SSEs
  for(i=0; i<n; i++){
    
    //Reset the number of non-missing groups
    nonmissing_m = m;
    
    //Iterate over the matrix columns (samples) to get groups means for each row
    for(j=0; j<p; j++){
      
      //Get the groupi
      groupi = gp[j]-1;
      
      //Compute group sums by adding each observations that's not an NA
      if(!NumericMatrix::is_na(data(i,j))){
        group_sums(i,groupi) += data(i,j);
        group_sizes(i,groupi) += 1; 
      }
      
    }
    
    //Store total number of non-na obs per row
    rowi_gsize = group_sizes(i,_);
    rowi_size = std::accumulate(rowi_gsize.begin(),rowi_gsize.end(),0.0);
    
    //Translate the group sums into group means for each row
    for(k=0; k<m; k++){
      group_means(i,k) = group_sums(i,k)/group_sizes(i,k);
      
      //If an entire group is missing (which shouldn't happen), the number of groups needs to be decreased
      if(group_sizes(i,k)<1){
        nonmissing_m = nonmissing_m - 1;
      }
      
      //group_sizes(k) = 0;
    }
    
    //compute overall mean for row i 
    rowi = group_means(i,_);
    overall_mean = std::accumulate(rowi.begin(),rowi.end(),0.0);
    overall_mean = overall_mean/nonmissing_m;
    
    //compute g-statistic
    Gstat_rowi = log(group_means(i,_)/overall_mean);
    Gstat_rowi = group_means(i,_)*Gstat_rowi;
    Gstats(i) = 2*std::accumulate(Gstat_rowi.begin(), Gstat_rowi.end(),0.0);
    
    //Arguments passed to pchisq are: value, df, lower tail?, log scale?
    p_value(i) = R::pchisq(Gstats(i),nonmissing_m-1,false,false);
  }//end iteration over rows
  
  return List::create(Named("group_means") = group_means,
                      Named("group_sizes") = group_sizes,
                      Named("Gstat") = Gstats,
                      Named("pvalue") = p_value);
}


/*
** R
#R-code to double check c++ function
#R version of the above
gtest_fun <- function(x,gp){
#Compute groups means, ANOVA F-stat and p-value   
  groups <- levels(gp)
  n <- sum(!is.na(x))
  ngroups <- length(groups)
  means <- rep(NA,ngroups)
  
  non_missing_groups <- length(unique(gp[!is.na(x)]))
  
  for(i in 1:ngroups){
    means[i] <- mean(x[gp==groups[i]],na.rm=T)
  }
  ngroups <- length(which(!is.na(means)))
  omean <- mean(means,na.rm=T)

  Gstat <- 2*sum(means*log(means/omean))
  pval <- pchisq(Gstat,ngroups-1,lower.tail = FALSE)
  return(c(means,Gstat,pval))
}

R_comparison <- function(data,gp){
  return(t(apply(data,1,gtest_fun,gp=gp)))
}

library(mintJansson)
peptide_data <- mintR::group_designation(mintR_data = peptide_data, main_effects = c("site","treatment"))
group_df <- attr(peptide_data, "group_DF")
gp <- group_df$Group
pep_data <- data.matrix(peptide_data$e_data[c(1:10),-1])
  
cpp_res <- gtest_cpp(pep_data, gp)
R_res <- R_comparison(pep_data, gp)
  
#Check group means
round(cpp_res$group_means-R_res[,1:6],digits=13)
  
#Check f-stats
round(cpp_res$Gstat-R_res[,7],digits=10)
  
#Check p-values
round(cpp_res$pvalue-R_res[,8],digits=10)
  
library(microbenchmark)
microbenchmark(gtest_cpp(pep_data, gp),R_comparison(pep_data, gp))  
*/
