#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix ptukey_speed(NumericMatrix qstats, NumericVector sizes) {
  int n_samps = qstats.nrow();
  int n_compares = qstats.ncol();
  NumericMatrix new_pvals(n_samps,n_compares);
  int i,j;
  
  for(i=0; i<n_samps; i++){
    for(j=0; j<n_compares; j++){
      //definition of ptukey
      //double 	ptukey(double q, double rr, double cc, double df, int lt, int lg)
      new_pvals(i,j) = R::ptukey(fabs(qstats(i,j)),1,n_compares, sizes(i)-n_compares,false,false);
    }
  }
  
  return new_pvals;
}


/*
** R
stats <- matrix(rnorm(100,20,2),20,5)/sqrt(2)
sizes <- rpois(20,20)
qvals <- ptukey(abs(stats),5,sizes-5,1,0,0)
qvals-ptukey_speed(stats,sizes)

microbenchmark(ptukey(stats,5,sizes-4,1,0,0),ptukey_speed(stats,sizes))

*/
