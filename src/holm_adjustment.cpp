#include <Rcpp.h>
using namespace Rcpp;

// This function is used with the std::sort function and moves NA or NaN values
// to the end of the vector. These functions are used together so the p-values
// will be sorted in ascending order while ignoring missing p-values. Without
// using the withNaN function the sort function will leave missing values in
// their initial location in the vector. This will lead to incorrectly adjusted
// p-values. (VERY VERY BAD!)
bool withNaN(double d0, double d1) {
  if (std::isnan(d0)) return false;
  if (std::isnan(d1)) return true;
  return d0 < d1;
}

// [[Rcpp::export]]
NumericVector holm_cpp(NumericVector ps) {

  // Check for the finite values in the p-value input vector (ps). Using
  // is_finite will return false for +/- infinity, NA, and NaN. This vector will
  // be used to determine the actual number of tests performed.
  LogicalVector is_pval = is_finite(ps);

  int n = std::accumulate(is_pval.begin(), is_pval.end(), 0);

  NumericVector sorted_ps = clone(ps);
  NumericVector adj_ps(n);
  int keep_going = 1;

  // Sort p-values in ascending order with NA/NaN at the end of the vector.
  std::sort(sorted_ps.begin(), sorted_ps.end(), withNaN);

  adj_ps[0] = fmin(sorted_ps[0]*n,1.0);

  for(int i=1; i<n; i++){

    if(keep_going){

      adj_ps[i] = fmin(sorted_ps[i]*(n-i),1.0);

    }else{
      adj_ps[i] = fmax(adj_ps[i-1],sorted_ps[i]);
    }

    if(adj_ps[i]<=adj_ps[i-1]){
      keep_going=0;
      adj_ps[i] = adj_ps[i-1];
    }

  }
  return adj_ps;
}

/*
** R

my_holm_r <- function(ps){
  n <- length(ps)
  orig_order <- rank(ps)
  adj_ps <- pmin(sort(ps)*(n-seq(0,n-1)),1)
  diffs <- diff(adj_ps)

  if(any(diffs<0)){
    adj_ps[(min(which(diffs<0))+1):n] <- adj_ps[min(which(diffs<0))]
  }

  return(adj_ps[orig_order])
}

ps <- rbeta(100,1,3)
all(holm_cpp(ps)==sort(my_holm_r(ps)))
all(holm_cpp(ps)==sort(p.adjust(ps,method = 'holm')))
all(my_second_holm_r(ps)==p.adjust(ps,method = 'holm'))


library(microbenchmark)
microbenchmark(holm_cpp(ps),my_holm_r(ps))
microbenchmark(my_second_holm_r(ps),p.adjust(ps,method = 'holm'))


psmat <- matrix(rbeta(3000,1,3),ncol=3)
holm_r <- t(apply(psmat,1,p.adjust,method = 'holm'))
my_res <- t(apply(psmat,1,my_second_holm_r))
plot(holm_r,my_res)

*/
