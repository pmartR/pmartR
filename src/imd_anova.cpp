#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// ANOVA FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

// [[Rcpp::export]]
List two_factor_anova_cpp(arma::mat y,
                          arma::mat X_full,
                          arma::mat X_red,
                          NumericVector red_df,
                          arma::colvec group_ids){

  int i,j;
  int n = y.n_rows;
  int p_red = X_red.n_cols;
  int p_full = X_full.n_cols;
  arma::rowvec yrowi(y.n_cols);
  arma::uvec to_remove;
  arma::mat X_red_nona, X_full_nona, invXtX_red, invXtX_full, XFinal;
  arma::mat PxRed, PxFull;
  arma::rowvec yrowi_nona;
  int num_to_remove;
  int df_red, df_full;
  double sigma2_red, sigma2_full;
  NumericVector sig_est(n), pval(n), Fstat(n),  atl_one;
  arma::mat diag_mat;
  arma::colvec par_ests(p_full), par_ests_temp,group_ids_nona(y.n_cols), group_ids_nona_unq;
  arma::mat parmat(n,p_full), group_sizes(n,p_full);
  arma::rowvec gsizes(p_full), csums(p_full);
  arma::uvec zero_cols, non_zero_cols;
  arma::uvec missing_gp,size_j;

  //Loop over rows in y
  for(i=0; i<n; i++){
    yrowi = y.row(i);
    to_remove = arma::find_nonfinite(yrowi);
    //Rcout << "Indices to remove " << to_remove << std::endl;

    yrowi_nona = yrowi;
    X_red_nona = X_red;
    X_full_nona = X_full;
    num_to_remove = to_remove.size();
    group_ids_nona = group_ids-1; //Make group_ids zero indexed

    //Remove NAs if applicable
    for(j=(num_to_remove-1);j>=0;j--){
      yrowi_nona.shed_col(to_remove[j]);
      X_red_nona.shed_row(to_remove[j]);
      X_full_nona.shed_row(to_remove[j]);
      group_ids_nona.shed_row(to_remove[j]);
    }

    //Remove completely empty columns
    csums = sum(X_full_nona,0);
    zero_cols = arma::find(csums==0);
    non_zero_cols = arma::find(csums);

    df_red = X_red_nona.n_rows-rank(X_red_nona)-red_df(i); //Subtract off df spent elsewhere, e.g., on covariates
    df_full = X_full_nona.n_rows-rank(X_full_nona)-red_df(i);

    PxRed = X_red_nona*pinv(X_red_nona.t()*X_red_nona)*X_red_nona.t();
    diag_mat.resize(size(PxRed));
    diag_mat.eye();

    sigma2_red = arma::conv_to<double>::from(yrowi_nona*(diag_mat-PxRed)*yrowi_nona.t()/(df_red));

    if((df_red-df_full)<=0){

      //If interaction can't be estimated, automatically select smaller model
      pval(i) = 1;
      Fstat(i) = 0;
    }else{
      PxFull = X_full_nona*pinv(X_full_nona.t()*X_full_nona)*X_full_nona.t();
      diag_mat.resize(size(PxFull));
      diag_mat.eye();
      sigma2_full = arma::conv_to<double>::from(yrowi_nona*(diag_mat-PxFull)*yrowi_nona.t()/(df_full));

      Fstat(i) = (sigma2_red*df_red-sigma2_full*df_full)/sigma2_full;
      pval(i) = R::pf(Fstat(i),df_red-df_full,df_full,0,0);
    }

    if(pval(i)<0.05){
      //Reject null hypothesis that reduced model is good enough, use full model
      XFinal = X_full_nona;
      sig_est(i) = sigma2_full;
    }else{
      XFinal = X_full_nona;

      XFinal.cols(p_red,p_full-1).fill(0.0); //Zero out the interaction terms if they're insignificant
      sig_est(i) = sigma2_red;
    }

    //"Parameter estiamtes" are the group means: Xbeta_hat=X(XpX)^{-1}XpY
    par_ests_temp = XFinal*pinv(XFinal.t()*XFinal)*XFinal.t()*yrowi_nona.t();

    //Find groups that had at least one non-missing value
    group_ids_nona_unq=arma::unique(group_ids_nona);

    //Fill in the par_ests vector, put NaN if group was missing or the average effect in groups with no missing data
    if(group_ids_nona_unq.n_elem<p_full){
      par_ests.zeros();
      //cntr = 0;

      for(j=0;j<p_full;j++){
        missing_gp = find(group_ids_nona_unq==j);
        if(missing_gp.n_elem>0){
          //Rcpp::Rcout <<"group_ids_nona \n"<< find(group_ids_nona==j)<<"\n";
          par_ests(j) = mean(par_ests_temp(find(group_ids_nona==j)));
          //cntr++;
        }else{
          par_ests(j) = arma::datum::nan;
        }
      }

    }else{
      for(j=0;j<p_full;j++){
        par_ests(j) = mean(par_ests_temp(find(group_ids_nona==j)));
      }
    }



    parmat.row(i) = par_ests.t();

    //Compute group sizes after accounting for NaNs
    for(j=0;j<p_full;j++){
      size_j = find(group_ids_nona==j);
      gsizes(j) = size_j.n_elem;
      //Rcpp::Rcout <<"size_j \n"<<size_j;
    }

    group_sizes.row(i) = gsizes; //For now don't return the interaction groups

  }
  return List::create(Named("par_estimates") = parmat,
                      Named("group_sizes") = group_sizes,
                      Named("Sigma2") = sig_est,
                      Named("Fstats") = Fstat,
                      Named("pvalue") = pval);
}

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
  arma::uvec found_finite;
  arma::uvec zero_C;
  arma::colvec not_nan;
  arma::mat good_C;
  arma::mat bad_C;
  arma::rowvec temp_fc_diff(num_comparisons);
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

    // Create a counter that keeps track of the number of groups with all
    // missing data. This value will be used in calculating the degrees of
    // freedom in the R::pt function. For example, if there are three groups but
    // one of them has all missing values the degrees of freedom for the number
    // of groups should be 3 - 1 = 2.
    int n_na_grps = 0;

    //Create (X'X)^{-1} for each row
    for(j=0; j<n_groups; j++){
      if(sizes(i,j)>0){
        XpXInv(j,j) = 1/sizes(i,j);
        N += sizes(i,j);
      } else {

        // Add one to n_na_grps because group j had a count of 0. Which means
        // all samples in the group had missing values.
        n_na_grps++;

      }

    }

    // Compute the square root of the variance-covariance matrix
    ses_mat = sqrt(C*XpXInv*C.t()*sigma2(i));
    // Take the diagonal elements as Vars for each comparison
    diff_ses.row(i) = arma::conv_to<arma::rowvec>::from(ses_mat.diag());
    t_tests.row(i) = diff_mat.row(i)/diff_ses.row(i);

    //Peel off sizes for row i and compute DFs
    rowi_sizes = arma::conv_to<arma::colvec>::from(sizes.row(i));
    //Rcpp::Rcout << "df_vec = " << df_vec << std::endl;
    for(k=0; k<num_comparisons; k++){

      //Only compute p-values if the degrees of freedom are atleast 3
      if (arma::is_finite(t_tests(i,k)) && (N - (n_groups - n_na_grps)) > 0) {
        p_values(i,k) = 2*R::pt(fabs(t_tests(i,k)),
                 (N - (n_groups - n_na_grps)),
                 false,
                 false);
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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// P-VALUE FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// FOLD CHANGE FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  arma::rowvec temp_fc_diff(num_comparisons);
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

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// COVARIATE FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  // Moore-Penrose pseudo-inverse, can always (?) be found but takes longer
  Px = pinv(X.t()*X)*X.t();

  // Set first ngroups rows to be zero
  Px.head_rows(ngroups).zeros();
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
  // vector of integers that represent how many degrees of freedom lost by
  // correcting for covariates
  arma::uvec lost_df(n);
  lost_df.fill(0);

  for(i=0;i<n;i++){

    rowi = data_mat.row(i);

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

  }

  return List::create(Named("data_no_x") = data_no_x,
                      Named("lost_df") = lost_df);

}
