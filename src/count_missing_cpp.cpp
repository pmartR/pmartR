#include <Rcpp.h>
using namespace Rcpp;
//Used by imd_test in order to count the number of missing observations for each row in 
//the data matrix

// [[Rcpp::export]]
NumericMatrix count_missing_cpp(NumericMatrix data, NumericVector gp) {
  //data is the n-by-p matrix of data
  //gp is a vector that identifies what group each column belongs to, 
  //  assumed to be numeric 1-m where m is the total number of groups
  int i,j,groupi;
  int n = data.nrow();  //number of observations
  int p = data.ncol();  //number of samples
  int m = max(gp);      //number of groups
  
  NumericMatrix num_missing(n,m);
  
  printf("n is %d \n", n);
  printf("p is %d \n", p);
  
  //Iterate over the matrix rows
  for(i=0; i<n; i++){
    
    //Iterate over the matrix columns
    for(j=0; j<p; j++){
      
      //Get the groupi
      groupi = gp[j]-1;
      
      if(NumericMatrix::is_na(data(i,j))){
        num_missing(i,groupi) += 1; 
      }
      
    }
  }
  return num_missing;
}

