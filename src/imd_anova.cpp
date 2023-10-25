#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

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

      if(static_cast<int>(bad_C.n_cols)<=(p-2)){ //Only proceed if there are at least two non-NaN means, otherwise return all NaNs
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
arma::rowvec fold_change_diff_row(arma::rowvec means, arma::mat C) {
  // Given the vector of means, and wanted group comparisons, compute the fold change using differencing
  // means - row vector of group means
  // C - matrix that defines the group comparisons you want to make
  int p = means.n_cols;
  int num_comparisons = C.n_rows;
  arma::rowvec fc_diff(num_comparisons);
  arma::uvec found_finite;
  arma::uvec zero_C;
  arma::colvec not_nan;
  arma::mat good_C;
  arma::mat bad_C;
  arma::colvec rsums;
  
  //Treat rows with NaN as special cases
  if(means.has_nan()){
    //Which rows of C do not involve the NaN element(s)?
    bad_C = C.cols(find_nonfinite(means)); //Submatrix of C involving the NaN element
    
    fc_diff.fill(arma::datum::nan); //Fill the vector with NAs
    
    if(static_cast<int>(bad_C.n_cols)<=(p-2)){ //Only proceed if there are at least two non-NaN means, otherwise return all NaNs
      //Absolute row sum of bad_C to see where the differences we can compute are
      rsums = arma::sum(abs(bad_C),1);
      
      zero_C = arma::find(rsums<0.01);    //Which rows of bad_C have zero elements? i.e., should compute differences of
      
      //non-nan elements of rowi_means
      found_finite = find_finite(means);
      not_nan = means.cols(found_finite).t(); 
      good_C = C.cols(found_finite); //C is a matrix so take the columns I need
      good_C = good_C.rows(zero_C);
      
      
      fc_diff.cols(zero_C) = arma::conv_to<arma::rowvec>::from(good_C*not_nan);
      
    }

  }else{
    fc_diff = arma::conv_to<arma::rowvec>::from(C*means.t());
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
  arma::uvec of_int;
  fc_diff.zeros();

  for(int i=0;i<n;i++){
    rowi_means = arma::conv_to<arma::colvec>::from(data.row(i));
    if(rowi_means.has_nan()){
      for(int j=0;j<num_comparisons;j++){
        of_int = find(C.row(j));
        for(int k=0;k<of_int.size();k++){
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
// ANOVA FUNCTIONS
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//For each row (peptide), this function computes the group means, counts the number of non-NA
//observations, estimates sigma^2, then computes the ANOVA F-statistic and associated p-value
//
//The results of this function are fed into group_comparison() which takes the group means,
//group sizes and sigma^2 value to do whatever group comparisons the user asks for

// [[Rcpp::export]]
List anova_cpp(
  arma::mat data, 
  NumericVector gp, 
  int unequal_var, 
  arma::mat X, 
  arma::mat Beta,
  arma::mat pred_grid,
  arma::uvec continuous_covar_inds,
  int n_covar_levels
  ) {
  //data is the n-by-p matrix of data
  //gp is a vector that identifies what group each column belongs to,
  //  assumed to be numeric 1-m where m is the total number of groups
  //unequal_var is a 0/1 depending on if variances are allowed to be unequal
  //df_red number of degrees of freedom spent elsewhere

  int n = data.n_rows;  //number of rows in data matrix, i.e., peptides/proteins/...
  int p = data.n_cols;  //number of samples
  int m = max(gp);      //number of groups
  double overall_mean = 0; //mean across groups
  NumericVector diff_vars(n); //Allow for unequal variances in two group situation

  NumericVector Fstats(n), p_value(n), sigma2(n); //F-statistic & p-value for each row
  NumericMatrix group_sums(n,m); //matrix to save the group sums in
  NumericMatrix group_means(n,m); //matrix to save the group means in
  arma::mat adj_group_means(n,m); // matrix of means adjusted for covariates
  NumericMatrix group_sizes(n,m);
  NumericVector rowi(m);
  NumericVector rowi_gsize(m); //vector to contain the number of non-na observations
  //per group

  arma::mat covar_vals;
  if (X.n_cols > m) {
    covar_vals = X.cols(m, X.n_cols - 1);
  }

  arma::rowvec Xrowi(p);
  arma::uvec elems_to_keep;
  arma::uvec dof(n);
  
  //Iterate over the matrix rows (peptides, proteins,...) to get group means
  //and SSEs
  for(int i=0; i<n; i++){

    //Be sure each SS starts out at zero
    double SSE = 0;
    double SSB = 0;
    double SST = 0;

    //Reset the number of non-missing groups
    int missing_m = 0;

    Xrowi = data.row(i);
    //Find the finite values in row i and put them in a column vector called "red_rowi"
    elems_to_keep = find_finite(Xrowi);

    int row_rank =  rank(X.rows(elems_to_keep)); // determines degrees of freedom, accounts for things like missing groups etc.
    arma::rowvec beta = Beta.row(i);

    // Average numeric covariates to get the adjusted means (lsmeans)
    if(continuous_covar_inds.size() > 0) {
      arma::mat X_nona = X.rows(elems_to_keep);
      // Replace each column indexed by covar_inds with it's mean
      for (int k = 0; k < continuous_covar_inds.size(); k++) {
        double colmean = arma::mean(X_nona.col(continuous_covar_inds[k] - 1));
        arma::colvec mean_vec = arma::ones(pred_grid.n_rows) * colmean;
        pred_grid.col(continuous_covar_inds[k] - 1) = mean_vec;
      }
    }

    // Compute the marginal means by average over covariates
    arma::colvec grid_preds = pred_grid * beta.t();

    if (n_covar_levels > 0) {
      arma::rowvec marginal_means(m);
      for (int k = 0; k < m; k++) {
        marginal_means(k) = arma::mean(grid_preds.rows(k * n_covar_levels, (k + 1) * n_covar_levels - 1));
      }

      adj_group_means.row(i) = marginal_means;
    } else {
      adj_group_means.row(i) = grid_preds.t();
    }

    //Iterate over the matrix columns (samples) to get groups means for each row
    for(int j=0; j<p; j++){

      //Get the groupi
      int groupi = gp[j]-1;

      //Compute group sums by adding each observations that's not an NA
      if(!NumericMatrix::is_na(data(i,j))){
        group_sums(i,groupi) += data(i,j);
        group_sizes(i,groupi) += 1;
      }

    }

    //Store total number of non-na obs per row
    rowi_gsize = group_sizes(i,_);
    int rowi_size = std::accumulate(rowi_gsize.begin(),rowi_gsize.end(),0.0);

    //Translate the group sums into group means for each row
    for(int k=0; k<m; k++){
      group_means(i,k) = group_sums(i,k)/group_sizes(i,k);

      //If an entire group is missing (which shouldn't happen), the number of groups needs to be decreased
      if(group_sizes(i,k)<1){
        missing_m = missing_m + 1;
        adj_group_means(i,k) = arma::datum::nan;
      }

      //group_sizes(k) = 0;
    }

    //compute overall mean for row i
    rowi = group_sums(i,_);
    overall_mean = std::accumulate(rowi.begin(),rowi.end(),0.0);
    overall_mean = overall_mean/rowi_size;

    if(m==2 && unequal_var==1){ //If there are only two groups, allow the variances to be different, i.e., Welch's t-test
      diff_vars[0] = 0;
      diff_vars[1] = 0;

      for(int j=0; j<p; j++){
        int groupi = gp[j] - 1;
        if(!NumericMatrix::is_na(data(i,j))){
          diff_vars[groupi] += pow(data(i,j)-group_means(i,groupi),2);
        }
      }
      diff_vars[0] /= (rowi_gsize[0]-1);
      diff_vars[1] /= (rowi_gsize[1]-1);
      sigma2(i) = diff_vars[0]/rowi_gsize[0]+diff_vars[1]/rowi_gsize[1]; //Welch's estimate of variance
      Fstats(i) = pow(group_means(i,0)-group_means(i,1),2)/sigma2(i);
      //Satterthwaite approximation to degrees of freedom
      double dfi = pow(sigma2(i),2)/(pow(diff_vars[0]/rowi_gsize[0],2)/(rowi_gsize[0]-1)+pow(diff_vars[1]/rowi_gsize[1],2)/(rowi_gsize[1]-1));
      p_value(i) = R::pf(Fstats(i),1,dfi,false,false);
      sigma2(i) *= dfi; //scale by estimated df to get sample variance not mean standard error
      dof(i) = dfi;
    }else{ //If there are more than two groups, compute pooled variance assuming equal variance across groups
      //Iterate over columns (again) to get a SSE for each row
      for(int j=0; j<p; j++){
        if(!NumericMatrix::is_na(data(i,j))){
          SST += pow(data(i,j)-overall_mean,2);
        }
      }
    
      arma::colvec residual_vec(p);
      residual_vec = arma::conv_to<arma::colvec>::from(data.row(i)).rows(elems_to_keep);
      residual_vec = residual_vec - X.rows(elems_to_keep)*Beta.row(i).t();
      SSE = arma::as_scalar(residual_vec.t()*residual_vec);

      //Compute F-statistic
      SSB = SST-SSE;
      sigma2(i) = (SSE/(rowi_size-row_rank));

      Fstats(i) = (SSB/(row_rank - 1))/(sigma2(i));
      //Arguments passed to pf are: value, df1, df2, lower tail?, log scale?
      p_value(i) = R::pf(Fstats(i),row_rank - 1,rowi_size-row_rank,false,false);
      dof(i) = rowi_size-row_rank;
    }

  }//end iteration over rows

  return List::create(Named("group_means") = group_means,
                      Named("adj_group_means") = adj_group_means,
                      Named("group_sizes") = group_sizes,
                      Named("Sigma2") = sigma2,
                      Named("Fstats") = Fstats,
                      Named("pvalue") = p_value,
                      Named("dof") = dof);
}

// [[Rcpp::export]]
List two_factor_anova_cpp(arma::mat y,
                          arma::mat X_full,
                          arma::mat X_red,
                          arma::colvec group_ids,
                          arma::mat beta_to_mu_full,
                          arma::mat beta_to_mu_red,
                          arma::uvec continuous_covar_inds,
                          int n_covar_levels
                          ){

  int i,j;
  int n = y.n_rows;
  int p_red = X_red.n_cols;
  int p_full = X_full.n_cols;
  arma::rowvec yrowi(y.n_cols);
  arma::uvec to_remove;
  arma::mat X_red_nona, X_full_nona, XFinal;
  arma::mat PxRed, PxFull;
  arma::rowvec yrowi_nona;
  int num_to_remove;

  int df_red, df_full;
  arma::uvec dof(n);
  arma::uvec whichX(n);
  
  double sigma2_red, sigma2_full;
  NumericVector sig_est(n), pval(n), Fstat(n);
  arma::mat diag_mat;
  arma::colvec par_ests(p_full),group_ids_nona(y.n_cols), group_ids_nona_unq;
  arma::mat parmat(n,p_full);
  arma::rowvec csums(p_full);
  arma::uvec zero_cols, non_zero_cols;
  arma::uvec missing_gp,size_j;

  // For a given biomolecule, the number of groups may be less
  arma::colvec unq_group_ids = arma::unique(group_ids);
  int n_groups = unq_group_ids.n_elem;
  arma::mat lsmeans(n,n_groups);
  arma::rowvec gsizes(n_groups);
  arma::mat group_sizes(n,n_groups);

  //Loop over rows in y
  for(i=0; i<n; i++){
    yrowi = y.row(i);
    to_remove = arma::find_nonfinite(yrowi);

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

    df_red = X_red_nona.n_rows-rank(X_red_nona);
    df_full = X_full_nona.n_rows-rank(X_full_nona);

    PxRed = pinv(X_red_nona.t()*X_red_nona)*X_red_nona.t()*yrowi_nona.t();
    
    arma::colvec residual_vec = yrowi_nona.t();
    residual_vec = residual_vec - X_red_nona*PxRed;

    sigma2_red = arma::conv_to<double>::from(residual_vec.t()*residual_vec/df_red);

    if((df_red-df_full)<=0){
      //If interaction can't be estimated, automatically select smaller model
      pval(i) = 1;
      Fstat(i) = 0;
      dof(i) = df_red;
    }else{
      PxFull = X_full_nona*pinv(X_full_nona.t()*X_full_nona)*X_full_nona.t();
      diag_mat.resize(size(PxFull));
      diag_mat.eye();
      sigma2_full = arma::conv_to<double>::from(yrowi_nona*(diag_mat-PxFull)*yrowi_nona.t()/(df_full));

      PxFull = pinv(X_full_nona.t()*X_full_nona)*X_full_nona.t()*yrowi_nona.t();
    
      arma::colvec residual_vec = yrowi_nona.t();
      residual_vec = residual_vec - X_full_nona*PxFull;

      sigma2_full = arma::conv_to<double>::from(residual_vec.t()*residual_vec/df_full);

      Fstat(i) = (sigma2_red*df_red-sigma2_full*df_full)/sigma2_full;
      pval(i) = R::pf(Fstat(i),df_red-df_full,df_full,0,0);
    }

    if(pval(i)<0.05){
      //Reject null hypothesis that reduced model is good enough, use full model
      XFinal = X_full_nona;
      sig_est(i) = sigma2_full;
      dof(i) = df_full;
      whichX(i) = 1; // 1 indicates full model
    }else{
      XFinal = X_full_nona;
      XFinal.cols(p_red,p_full-1).fill(0.0); //Zero out the interaction terms if they're insignificant
      sig_est(i) = sigma2_red;
      dof(i) = df_red;
      whichX(i) = 0;
    }

    // Parameter estimates using the full X
    arma::colvec beta = pinv(XFinal.t()*XFinal)*XFinal.t()*yrowi_nona.t();

    //Find groups that had at least one non-missing value
    group_ids_nona_unq=arma::unique(group_ids_nona);

    arma::mat beta_to_mu;
    if (whichX[i] == 0) {
      beta_to_mu = beta_to_mu_red;
    } else {
      beta_to_mu = beta_to_mu_full;
    }
    
    // Average numeric covariates in preparation for model predictions.
    if(continuous_covar_inds.size() > 0) {
      // Replace each column indexed by covar_inds with it's mean
      for (int k = 0; k < continuous_covar_inds.size(); k++) {
        double colmean = arma::mean(XFinal.col(continuous_covar_inds[k] - 1));
        arma::colvec mean_vec = arma::ones(beta_to_mu.n_rows) * colmean;
        beta_to_mu.col(continuous_covar_inds[k] - 1) = mean_vec;
      }
    }

    // Predicted values over all levels of the main effects and categorical covariates.  Continuous covariates are averaged.
    arma::colvec grid_preds = beta_to_mu * beta.head_rows(beta_to_mu.n_cols);

    // Compute the lsmeans by averaging over covariates
    if (n_covar_levels > 1) {
      arma::rowvec marginal_means(n_groups);
      int n_per_covar_group = grid_preds.n_elem / n_covar_levels;

      for (int k = 0; k < n_groups; k++) {
        marginal_means(k) = arma::mean(grid_preds.rows(k * n_per_covar_group, (k + 1) * n_per_covar_group - 1));
      }

      lsmeans.row(i) = marginal_means;
    } else {
      lsmeans.row(i) = grid_preds.t();
    }
    
    parmat.row(i) = beta.t();

    //Compute group sizes after accounting for NaNs
    for(j=0;j<n_groups;j++){
      size_j = find(group_ids_nona==j);
      gsizes(j) = size_j.n_elem;
    }
    
    group_sizes.row(i) = gsizes; //For now don't return the interaction groups

  }
  return List::create(Named("lsmeans") = lsmeans,
                      Named("par_estimates") = parmat,
                      Named("group_sizes") = group_sizes,
                      Named("Sigma2") = sig_est,
                      Named("Fstats") = Fstat,
                      Named("pvalue") = pval,
                      Named("dof") = dof,
                      Named("which_X") = whichX);
}

// [[Rcpp::export]]
List group_comparison_anova_cpp(arma::mat means, arma::mat data, arma::mat sizes, arma::mat which_xmatrix, arma::mat Xfull, arma::mat Xred, arma::mat Cfull, arma::mat Cred, arma::mat Cmu) {
  //Given the raw data, group sizes, design matrix consisting of group levels, and comparison matrix return the group comparisons requested.  Returns the estimated difference, standard errors, t-statistics, and p-values.
  //Returns estimated difference, standard error, t-statistic, p-values
  //data - raw data matrix of n_biomolecules by n_samples
  //sizes - matrix of group sizes n_biomolecules x n_groups
  //C - matrix that defines the group comparisons you want to make

  int num_comparisons = Cred.n_rows; //Number of comparisons to be made
  int n_groups = sizes.n_cols; //Number of groups 
  int n = sizes.n_rows; //Number of rows (peptides)
  int N = 0; //Total number of observations for rowi;
  arma::colvec rowi_means(n_groups), rowi_sizes(n_groups); //Storage for each row of means and sample sizes

  //arma::mat diff_mat(n,num_comparisons); //Matrix of comparison means
  arma::mat diff_mat(n, num_comparisons);

  arma::mat diff_ses(n,num_comparisons); //Matrix of comparison ses
  arma::mat t_tests(n,num_comparisons); //Matrix of t-test statistics
  arma::mat p_values(n,num_comparisons); //Matrix of comparison p-values
  p_values.fill(1.0);
  arma::mat ses_mat(num_comparisons, num_comparisons); //Storage for the sqrt of the comparison var-cov matrix

  arma::rowvec rowi;
  arma::uvec elems_to_keep;

  diff_mat.zeros();
  p_values.zeros();

  for(int i=0;i<n;i++){
    arma::mat X;
    arma::mat C;

    if (which_xmatrix(i) == 0) {
      X = Xred;
      C = Cred;
    } else {
      X = Xfull;
      C = Cfull;
    }
    rowi = data.row(i);
    elems_to_keep = find_finite(rowi);
    arma::mat XpXInv = pinv(X.rows(elems_to_keep).t() * X.rows(elems_to_keep));

    // keep track of the number of missing groups
    for(int j=0; j<n_groups; j++){
      if(sizes(i,j)>0){
        N += sizes(i,j);
      } 
    }

    int dof = N - rank(X.rows(elems_to_keep));

    arma::uvec zerosize_idx = find(sizes.row(i).cols(0, n_groups - 1) == 0);
    arma::rowvec row_mean = means.row(i); 
    if (zerosize_idx.size() > 0) {
      row_mean.cols(zerosize_idx).fill(arma::datum::nan);
    }

    // Take the group means using the group means contrast matrix Cmu.  This needs to be done since the Betas to not take into account missing groups.
    arma::rowvec mean_diffs = fold_change_diff_row(row_mean, Cmu.cols(0, n_groups - 1));
    diff_mat.row(i) = mean_diffs;
    
    arma::colvec residual_vec(X.n_rows);
    residual_vec = arma::conv_to<arma::colvec>::from(data.row(i)).rows(elems_to_keep);

    arma::colvec Beta = XpXInv*X.rows(elems_to_keep).t()*residual_vec;

    residual_vec = residual_vec - X.rows(elems_to_keep)*Beta;
    double SSE = arma::as_scalar(residual_vec.t()*residual_vec);
    double sigma_hat = SSE/dof;

    // Compute the square root of the variance-covariance matrix
    ses_mat = sqrt(C*XpXInv*C.t()*sigma_hat);
    // Take the diagonal elements as Vars for each comparison
    diff_ses.row(i) = arma::conv_to<arma::rowvec>::from(ses_mat.diag());
    t_tests.row(i) = mean_diffs/diff_ses.row(i);

    for(int k=0; k<num_comparisons; k++){

      //Only compute p-values if the degrees of freedom are atleast 3
      if (arma::is_finite(t_tests(i,k)) && dof > 0) {
        p_values(i,k) = 2*R::pt(fabs(t_tests(i,k)),
                 dof,
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

  // Sort p-values in ascending order with NA/NaN at the end of the vector.
  std::sort(sorted_ps.begin(), sorted_ps.end(), withNaN);

  adj_ps[0] = fmin(sorted_ps[0]*n, 1.0);

  for(int i=1; i<n; i++){

    adj_ps[i] = fmin(sorted_ps[i]*(n - i), 1.0);

    if(adj_ps[i]<=adj_ps[i - 1]){

      adj_ps[i] = adj_ps[i - 1];
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

//-----Compute the matrix of coefficients for the simple additive model-----//
// [[Rcpp::export]]
arma::mat compute_betas(arma::mat data_mat, arma::mat Xmatrix){
  int n = data_mat.n_rows;  //number of biomolecules
  arma::rowvec rowi;
  arma::colvec red_rowi, beta;
  arma::mat effects, Xtemp;
  arma::mat B = arma::zeros(n, Xmatrix.n_cols);
  arma::uvec elems_to_keep;
  arma::uvec row_ind;

  for(int i=0;i<n;i++){

    rowi = data_mat.row(i);

    //Find the finite values in row i and put them in a column vector called "red_rowi"
    elems_to_keep = find_finite(rowi);
    red_rowi = rowi.elem(elems_to_keep);

    //Compute the projection matrix after removing the rows with missing data
    Xtemp = Xmatrix.rows(elems_to_keep);
    arma::mat P = pinv(Xtemp.t()*Xtemp)*Xtemp.t();
    beta = P*red_rowi;
    B.row(i) = beta.t();
  }
    
  return B;

}
