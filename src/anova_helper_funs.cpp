#include <Rcpp.h>
using namespace Rcpp;

//For each row (peptide), this function computes the group means, counts the number of non-NA
//observations, estimates sigma^2, then computes the ANOVA F-statistic and associated p-value
//
//The results of this function are fed into group_comparison() which takes the group means,
//group sizes and sigma^2 value to do whatever group comparisons the user asks for

// [[Rcpp::export]]
List anova_cpp(NumericMatrix data, NumericVector gp, int unequal_var, NumericVector df_red) {
  //data is the n-by-p matrix of data
  //gp is a vector that identifies what group each column belongs to, 
  //  assumed to be numeric 1-m where m is the total number of groups
  //unequal_var is a 0/1 depending on if variances are allowed to be unequal
  //df_red number of degrees of freedom spent elsewhere
  
  int i,j,k,groupi,rowi_size;
  int n = data.nrow();  //number of rows in data matrix, i.e., peptides/proteins/...
  int p = data.ncol();  //number of samples
  int m = max(gp);      //number of groups
  int nonmissing_m; //number of groups with at least one observation
  double overall_mean = 0; //mean across groups
  NumericVector diff_vars; //Allow for unequal variances in two group situation
  double dfi; //Satterthwaite df approximation if unequal variances
  
  double SSE,SST, SSB; //SSEs, SSTs and SSBs computed by row
  NumericVector Fstats(n), p_value(n), sigma2(n); //F-statistic & p-value for each row
  NumericMatrix group_sums(n,m); //matrix to save the group sums in
  NumericMatrix group_means(n,m); //matrix to save the group means in
  NumericMatrix group_sizes(n,m);
  NumericVector rowi(m);
  NumericVector rowi_gsize(m); //vector to contain the number of non-na observations
                                //per group
  
  //Iterate over the matrix rows (peptides, proteins,...) to get group means
  //and SSEs
  for(i=0; i<n; i++){
    
    //Be sure each SS starts out at zero
    SSE = 0;
    SSB = 0;
    SST = 0;
    
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
    rowi = group_sums(i,_);
    overall_mean = std::accumulate(rowi.begin(),rowi.end(),0.0);
    overall_mean = overall_mean/rowi_size;
    
    if(m==2 & unequal_var==1){ //If there are only two groups, allow the variances to be different, i.e., Welch's t-test
      diff_vars[0] = 0;
      diff_vars[1] = 0;
      
      for(j=0; j<p; j++){
        groupi = gp[j] - 1;
        if(!NumericMatrix::is_na(data(i,j))){
          diff_vars[groupi] += pow(data(i,j)-group_means(i,groupi),2);  
        }
      }
      diff_vars[0] /= (rowi_gsize[0]-1);
      diff_vars[1] /= (rowi_gsize[1]-1);
      sigma2(i) = diff_vars[0]/rowi_gsize[0]+diff_vars[1]/rowi_gsize[1]; //Welch's estimate of variance
      Fstats(i) = pow(group_means(i,0)-group_means(i,1),2)/sigma2(i);
      //Satterthwaite approximation to degrees of freedom
      dfi = pow(sigma2(i),2)/(pow(diff_vars[0]/rowi_gsize[0],2)/(rowi_gsize[0]-1)+pow(diff_vars[1]/rowi_gsize[1],2)/(rowi_gsize[1]-1));
      p_value(i) = R::pf(Fstats(i),1,dfi,false,false);
      sigma2(i) *= dfi; //scale by estimated df to get sample variance not mean standard error

    }else{ //If there are more than two groups, compute pooled variance assuming equal variance across groups
      //Iterate over columns (again) to get a SSE for each row
      for(j=0; j<p; j++){
        groupi = gp[j] - 1;
        if(!NumericMatrix::is_na(data(i,j))){
          SSE += pow(data(i,j)-group_means(i,groupi),2);  
          SST += pow(data(i,j)-overall_mean,2);  
        }
      }
      //Compute F-statistic
      SSB = SST-SSE;
      sigma2(i) = (SSE/(rowi_size-nonmissing_m-df_red(i)));
      Fstats(i) = (SSB/(nonmissing_m - 1))/(sigma2(i));
      //Arguments passed to pf are: value, df1, df2, lower tail?, log scale?
      p_value(i) = R::pf(Fstats(i),nonmissing_m-1,rowi_size-nonmissing_m-df_red(i),false,false);
      
    }


  }//end iteration over rows
  
    return List::create(Named("group_means") = group_means,
                      Named("group_sizes") = group_sizes,
                      Named("Sigma2") = sigma2,
                      Named("Fstats") = Fstats,
                      Named("pvalue") = p_value);
}


/*

 ** R
#R-code to double check c++ function
#R version of the above
mean_sse_fun <- function(x,gp){
  #Compute groups means, ANOVA F-stat and p-value   
  groups <- levels(gp)
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

R_comparison <- function(data,gp){
  return(t(apply(data,1,mean_sse_fun,gp=gp)))
}

library(mintJansson)
peptide_data <- mintR::group_designation(mintR_data = peptide_data, main_effects = c("site","treatment"))
group_df <- attr(peptide_data, "group_DF")
gp <- group_df$Group
pep_data <- data.matrix(peptide_data$e_data[c(1:10),-1])

cpp_res <- anova_cpp(pep_data, gp)
R_res <- R_comparison(pep_data, gp)

#Check group means
round(cpp_res$group_means-R_res[,1:6],digits=13)

#Check f-stats
round(cpp_res$Fstats-R_res[,7],digits=10)

#Check p-values
round(cpp_res$pvalue-R_res[,8],digits=10)

#library(microbenchmark)
#microbenchmark(group_means_sse(fake_data, fake_groups),R_comparison(fake_data, fake_groups))  
*/
