#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix ptukey_speed(NumericMatrix qstats, NumericVector sizes) {
  int n_samps = qstats.nrow();
  int n_compares = qstats.ncol();
  NumericMatrix new_pvals(n_samps,n_compares);
  int i,j;

  for (i=0; i<n_samps; i++) {

    // Extract the ith row of qstats as a vector.
    NumericVector da_stats = qstats(i, _);

    // Check for the finite values in the statistics row vector (qstats[i, ]).
    // Using is_finite will return false for +/- infinity, NA, and NaN. This
    // vector will be used to determine the actual number of tests performed.
    LogicalVector isValid = is_finite(da_stats);

    int nActual = std::accumulate(isValid.begin(), isValid.end(), 0);

    for (j=0; j<n_compares; j++) {

      //definition of ptukey
      //double 	ptukey(double q, double rr, double cc, double df, int lt, int lg)
      new_pvals(i,j) = R::ptukey(fabs(qstats(i,j)),
                                 1,
                                 nActual,
                                 sizes(i) - nActual,
                                 false,
                                 false);
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
