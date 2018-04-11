#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


//------Function to construct the null space projection matrix------//
// [[Rcpp::export]]
List proj_mat_cpp(arma::mat X, int ngroups){
  //If the X matrix has atleast two rows, find projection matrix 
  //into the null space corresponding to X

  int n = X.n_rows;
  arma::mat Px;
  arma::mat Imat(n,n);
  arma::uword PxRank;
  Imat.eye();
  
  Px = pinv(X.t()*X)*X.t(); //Moore-Penrose pseudo-inverse, can always (?) be found but takes longer
  Px.head_rows(ngroups).zeros(); //Set first ngroups rows to be zero
  Px = X*Px;
  PxRank = rank(Px);
  return List::create(Named("Ipx") = Imat-Px,
                      Named("PxRank") = PxRank);
}

//-----Project each row of the data.matrix into X's null space-----//
// [[Rcpp::export]]
List project_to_null_cpp(arma::mat data_mat, arma::mat Xmatrix, int ngroups){
  
  arma::mat data_no_x;
  data_no_x = data_mat;
  int i;
  int n = data_mat.n_rows;
  List proj_res;
  arma::rowvec rowi, rowi_no_x;
  arma::colvec red_rowi;
  arma::mat ImPx, Xtemp, ImpxTemp;
  arma::uvec elems_to_keep;
  arma::uvec row_ind;
  //vector of integers that represent how many degrees of freedom lost by correcting for covariates
  arma::uvec lost_df(n);
  lost_df.fill(0);
  arma::uword full_df;
  
  
  //Compute projection matrix once using full Xmatrix
  proj_res = proj_mat_cpp(Xmatrix, ngroups);
  ImPx = as<arma::mat>(proj_res["Ipx"]);
  full_df = as<arma::uword>(proj_res["PxRank"]);
  
  for(i=0;i<n;i++){
    
    rowi = data_mat.row(i);
    
    if(rowi.has_nan()){
      
      //Find the finite values in row i and put them in a column vector called "red_rowi"
      elems_to_keep = find_finite(rowi);
      red_rowi = rowi.elem(elems_to_keep);

      //Compute the projection matrix after removing the rows with missing data 
      Xtemp = Xmatrix.rows(elems_to_keep);
      proj_res = proj_mat_cpp(Xtemp, ngroups);
      ImpxTemp = as<arma::mat>(proj_res["Ipx"]);;
      lost_df(i) = as<arma::uword>(proj_res["PxRank"]);
      
      //Project the data into the null space of Xtemp
      rowi_no_x = arma::conv_to<arma::rowvec>::from(ImpxTemp*red_rowi);
      
      //Replace the ith row, finite elements with the projected data
      row_ind = i;
      data_no_x.submat(row_ind,elems_to_keep) = rowi_no_x;
      
    }else{
      //Translate 'rowi' into a column vector, multiply by 'ImPx' then translate back into row vector
      data_no_x.row(i) = arma::conv_to<arma::rowvec>::from(ImPx*arma::conv_to<arma::colvec>::from(rowi));
      lost_df(i) = full_df;
    }
  }
  
  return List::create(Named("data_no_x") = data_no_x,
                      Named("lost_df") = lost_df);
    
}





/*
** R
source('~/Documents/MinT/MSomics/MSomicsSTAT/R/covariate_supp_funs.R')
xmat <- matrix(c(1,1,1,0,0,0,0,0,1,1),ncol=2)
data_mat <- matrix(rnorm(100),20,5)

Rversion <- proj_mat(xmat)
cppversion <- proj_mat_cpp(xmat)
Rversion-cppversion

proj_mat_cpp(matrix(1,1,1))

#microbenchmark::microbenchmark(proj_mat(xmat),proj_mat_cpp(xmat))

data_mat[2] <- NA
project_to_null_cpp(data_mat,xmat)-project_to_null(data_mat,xmat)

#microbenchmark::microbenchmark(project_to_null_cpp(data_mat,xmat),project_to_null(data_mat,xmat))

*/
