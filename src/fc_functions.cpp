#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/* Old verion of fold_change_diff that can't handle missing values (kept for posterity only)
arma::mat fold_change_diff(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change using differencing
  //means - matrix of group means
  // C - matrix that defines the group comparisons you want to make
  
  int n = data.n_rows;
  int p = data.n_cols;
  int num_comparisons = C.n_rows;
  arma::mat fc_diff(n,num_comparisons);
  arma::colvec rowi_means(p);
  int i;
  
  for(i=0;i<n;i++){
    rowi_means = arma::conv_to<arma::colvec>::from(data.row(i));
    //if(rowi_means.has_nan()){
    //  for(j=0;j<p;j++){
    //    if(NumericVector::is_na(rowi_means(j))){
    //      rowi_means(j) = 1e6;
    //    }
    //  }
    //}
    //Rcout << "The value is " << rowi_means << std::endl;
    fc_diff.row(i) = arma::conv_to<arma::rowvec>::from(C*rowi_means);
  }
  
  return fc_diff;
}*/

// [[Rcpp::export]]
arma::mat fold_change_diff(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change using differencing
  //means - matrix of group means
  // C - matrix that defines the group comparisons you want to make
  
  int n = data.n_rows;
  int p = data.n_cols;
  int num_comparisons = C.n_rows;
  arma::mat fc_diff(n,num_comparisons);
  arma::colvec rowi_means(p);
  arma::uvec found_finite;
  arma::uvec zero_C;
  arma::colvec not_nan;
  arma::mat good_C;
  arma::mat bad_C;
  arma::rowvec temp_fc_diff(p);
  arma::colvec rsums;
  
  for(int i=0;i<n;i++){
    rowi_means = arma::conv_to<arma::colvec>::from(data.row(i));
    
    //Treat rows with NaN as special cases
    if(rowi_means.has_nan()){
      //Which rows of C do not involve the NaN element(s)?
      bad_C = C.cols(find_nonfinite(rowi_means)); //Submatrix of C involving the NaN element
      
      temp_fc_diff.fill(arma::datum::nan); //Fill the vector with NAs
      
      if(bad_C.n_cols<=(p-2)){ //Only proceed if there are at least two non-NaN means, otherwise return all NaNs
        //Absolute row sum of bad_C to see where the differences we can compute are
        rsums = arma::sum(abs(bad_C),1);
        
        zero_C = arma::find(rsums<0.01);    //Which rows of bad_C have zero elements? i.e., should compute differences of
        
        //non-nan elements of rowi_means
        found_finite = find_finite(rowi_means);
        not_nan = rowi_means.rows(found_finite); //rowi_means is a column vector so take it's rows
        good_C = C.cols(found_finite);          //C is a matrix so take the columns I need
        good_C = good_C.rows(zero_C);
        
        
        temp_fc_diff.cols(zero_C) = arma::conv_to<arma::rowvec>::from(good_C*not_nan);
      }
      
      fc_diff.row(i) = temp_fc_diff;
    }else{
      fc_diff.row(i) = arma::conv_to<arma::rowvec>::from(C*rowi_means);
    }
  }
  return fc_diff;
}


// [[Rcpp::export]]
arma::mat fold_change_ratio(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change ratios
  //means - matrix of group means
  // C - matrix that defines the group comparisons you want to make
  
  arma::mat fc_mat;
  arma::mat log_data;
  arma::mat fc_diff;
  //ratio is same as difference on log scale, so take log, then difference, then exponentiate

  log_data = log(data);
  fc_diff = fold_change_diff(log_data, C);
  fc_mat = exp(fc_diff);
    
  return fc_mat;
}

// [[Rcpp::export]]
arma::mat fold_change_logbase2(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change ratios
  //means - matrix of group means on the log base 2 scale
  // C - matrix that defines the group comparisons you want to make
  
  arma::mat fc_mat;
  arma::mat fc_diff;
  
  //ratio on original scale (x/y) is exp{ln(2)[log(x,base2)-log(y,base2)]}
  fc_diff = log(2)*fold_change_diff(data, C);
  fc_mat = exp(fc_diff);
  
  return fc_mat;
}


// [[Rcpp::export]]
arma::mat fold_change_diff_na_okay(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change using differencing
  //To avoid returning NAs even if all groups of interest are not NA, it will do elementwise 
  //multiplication instead of matrix multiplication
  
  //means - matrix of group means
  // C - matrix that defines the group comparisons you want to make
  
  int n = data.n_rows;
  int p = data.n_cols;
  int num_comparisons = C.n_rows;
  arma::mat fc_diff(n,num_comparisons);
  arma::colvec rowi_means(p);
  int i,j,k;
  arma::uvec of_int;
  fc_diff.zeros();
  
  for(i=0;i<n;i++){
    rowi_means = arma::conv_to<arma::colvec>::from(data.row(i));
    if(rowi_means.has_nan()){
      for(j=0;j<num_comparisons;j++){
          of_int = find(C.row(j));
          for(k=0;k<of_int.size();k++){
            fc_diff(i,j) += C(j,of_int[k])*rowi_means(of_int[k]);
          }
        }
    }else{
      fc_diff.row(i) = arma::conv_to<arma::rowvec>::from(C*rowi_means);
    }
  }
  
  return fc_diff;
}

// The R code below will make sure I'm doing it righ in C++,

/*
** R

fake_data <- matrix(c(abs(rnorm(20)),rpois(20,2),rgamma(20,2,5)),nrow=10)
fake_data[c(3,6,40,30,52,59)] <- NA
means <- matrix(c(rowMeans(fake_data[,c(1,2)],na.rm=T),rowMeans(fake_data[,c(3,4)],na.rm=T),rowMeans(fake_data[,c(5,6)],na.rm=T)),ncol=3)
meansb2 <- log(means,base=2)
R_res <- matrix(0,nrow(fake_data),3)
Cmat <- matrix(c(1,-1,0,0,1,-1,1,0,-1),nrow=3,byrow=T)

for(i in 1:nrow(means)){
  R_res[i,1] <- means[i,1]/means[i,2]
  R_res[i,2] <- means[i,2]/means[i,3]
  R_res[i,3] <- means[i,1]/means[i,3]
}

cpp_res <- fold_change_ratio(means,Cmat)
cpp_res_b2 <- fold_change_logbase2(meansb2,Cmat)
cpp_res
cpp_res_b2
R_res
*/
