#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//This is a copy of `fold_change_diff` that will be called by `group_comparison`
// [[Rcpp::export]]
arma::mat fold_change_diff_copy(arma::mat data, arma::mat C)  {
  //Given the group means, and wanted group comparisons, compute the fold change using differencing
  //means - matrix of group means
  // C - matrix that defines the group comparisons you want to make
  
  int n = data.n_rows;
  int p = data.n_cols;
  int num_comparisons = C.n_rows;
  arma::mat fc_diff(n,num_comparisons);
  arma::colvec rowi_means(p);
  for(int i=0;i<n;i++){
    rowi_means = arma::conv_to<arma::colvec>::from(data.row(i));
    fc_diff.row(i) = arma::conv_to<arma::rowvec>::from(C*rowi_means);
  }
  return fc_diff;
}

// [[Rcpp::export]]
List group_comparison_anova_cpp(arma::mat means, arma::mat sizes, arma::vec sigma2, arma::mat C) {
  //Given the group means, group sizes and esimtaes of sigma^2, do all the 
  //group comparisons requested.  Returns estimated difference, standard error, t-statistic
  //and p-value
  //means - matrix of group means
  //sizes- matrix of group sizes
  //sigma2 - vector of sigma^2 estimates
  // C - matrix that defines the group comparisons you want to make
  
  int i, j, k;
  int num_comparisons = C.n_rows; //Number of comparisons to be made
  int n_groups = means.n_cols; //Number of groups
  int n = means.n_rows; //Number of rows (peptides)
  int N = 0; //Total number of observations for rowi;
  arma::mat XpXInv(n_groups, n_groups); //Matrix used to compute SEs: (X'X)^{-1}
  arma::colvec rowi_means(n_groups), rowi_sizes(n_groups); //Storage for each row of means and sample sizes

  //arma::mat diff_mat(n,num_comparisons); //Matrix of comparison means
  arma::mat diff_mat;
  
  arma::mat diff_ses(n,num_comparisons); //Matrix of comparison ses
  arma::mat t_tests(n,num_comparisons); //Matrix of t-test statistics
  arma::mat p_values(n,num_comparisons); //Matrix of comparison p-values
  p_values.fill(1.0);
  arma::mat ses_mat(num_comparisons, num_comparisons); //Storage for the sqrt of the comparison var-cov matrix
  
  diff_mat.zeros();
  p_values.zeros();
  
  //Compute the differences between means using the fold_change_diff function
  diff_mat = fold_change_diff_copy(means, C);
  
  for(i=0;i<n;i++){
    XpXInv.zeros(); //Fill with zeros
    
    //Create (X'X)^{-1} for each row
    for(j=0; j<n_groups; j++){
      if(sizes(i,j)>0){
        XpXInv(j,j) = 1/sizes(i,j);
        N += sizes(i,j);
      }
    }
    
    ses_mat = sqrt(C*XpXInv*C.t()*sigma2(i)); // Compute the square root of the variance-covariance matrix
    diff_ses.row(i) = arma::conv_to<arma::rowvec>::from(ses_mat.diag()); //Take the diagonal elements as Vars for each comparison
    t_tests.row(i) = diff_mat.row(i)/diff_ses.row(i);
    
    //Peel off sizes for row i and compute DFs
    rowi_sizes = arma::conv_to<arma::colvec>::from(sizes.row(i));
    //Rcpp::Rcout << "df_vec = " << df_vec << std::endl;
    for(k=0; k<num_comparisons; k++){
      
      //Only compute p-values if the degrees of freedom are atleast 3
      if(arma::is_finite(t_tests(i,k)) && (N-n_groups)>2){
        p_values(i,k) = 2*R::pt(fabs(t_tests(i,k)),(N-n_groups),false,false);
      }else{
        p_values(i,k) = arma::datum::nan;
      }
    }
    N = 0; //reset row size counter
  }
  
  return List::create(Named("diff_mat") = diff_mat,
                      Named("diff_ses") = diff_ses,
                      Named("t_tests") = t_tests,
                      Named("p_values") = p_values);
  
}

/*
** R
Rcpp::sourceCpp('mintR/inst/cpp/anova_helper_funs.cpp')
fake_data <- matrix(c(rnorm(20),rpois(20,2),rgamma(20,2,5)),nrow=10)
fake_data[c(3,6,40,30,52,59)] <- NA
fake_groups <- c(1,1,2,2,3,3)
R_res <- matrix(0,nrow(fake_data),4)
Cmat <- matrix(c(1,-1,0,0,1,-1),nrow=2,byrow=T)

cpp_res_p1 <- anova_cpp(fake_data,fake_groups)
for(i in 1:nrow(fake_data)){
  lmres <- lm(fake_data[i,]~as.factor(fake_groups))
  R_res[i,c(1:3)] <- c(mean(fake_data[i,c(1,2)],na.rm=T),mean(fake_data[i,c(3,4)],na.rm=T),mean(fake_data[i,c(5,6)],na.rm=T))
  R_res[i,4] <- summary(lmres)$sigma^2
}

R_res2 <- t(Cmat%*%t(R_res[,1:3]))

final_res <- group_comparison(cpp_res_p1$group_means,cpp_res_p1$group_sizes,cpp_res_p1$Sigma2,Cmat)

final_res$diff_mat-R_res2

*/
