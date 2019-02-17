#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List two_factor_anova_cpp(arma::mat y, arma::mat X_full, arma::mat X_red, NumericVector red_df, arma::colvec group_ids){
  
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
    //if(zero_cols.n_rows>0){
    //  for(j=(zero_cols.n_rows-1);j>=0;j--){
    //    X_full_nona.shed_col(zero_cols[j]);
    //  }  
    //}
    
    //while(rank(X_full_nona)<X_full_nona.n_cols){
    //  X_full_nona.shed_col((X_full_nona.n_cols-1));
    //}
    
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
      //Rcpp::Rcout << "Dim XFinal: " << XFinal.n_cols << "; p_red "<<p_red<<"; p_full "<<p_full << std::endl;
      XFinal.cols(p_red,p_full-1).fill(0.0); //Zero out the interaction terms if they're insignificant
      sig_est(i) = sigma2_red;
    }
    //Rcpp::Rcout <<"Xfinal\n"<<XFinal;  
    par_ests_temp = XFinal*pinv(XFinal.t()*XFinal)*XFinal.t()*yrowi_nona.t();
    //Rcpp::Rcout <<"par_ests_temp \n"<<par_ests_temp<<"\n";  
    group_ids_nona_unq=arma::unique(group_ids_nona);
    //Rcpp::Rcout <<"size_na_gps \n"<< group_ids_nona;  
    
    
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
    
    //Rcpp::Rcout <<"par_ests \n"<<par_ests;  
    
    
    parmat.row(i) = par_ests.t();
    //Rcpp::Rcout <<"group_ids_nona \n"<<group_ids_nona;  
    //Sum columns of to get number of observations in each group
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


/*
** R
nrows_to_kp <- nrow(data)
ytst <- matrix(data.matrix(unname(data[1:nrows_to_kp,])),nrow=nrows_to_kp)
grp_nums <- as.numeric(as.factor(groupData$Group))
Xfull_nona <- Xfull[-which(is.na(ytst)),]
Xred_nona <- Xred[-which(is.na(ytst)),]

all_dat <- cbind(groupData,Y=t(ytst))
(ad_sum <- all_dat%>%group_by(Media,Location)%>%summarize(Size=length(which(!is.na(Y))),Mean=mean(Y,na.rm=T)))
colSums(Xfull_nona)
colSums(Xred_nona)
tfan <- two_factor_anova_cpp(y = ytst,X_red = Xred,X_full = Xfull,red_df,group_ids = grp_nums)
cbind(ad_sum,Tst=t(tfan$par_estimates),Size=t(tfan$group_sizes))

*/
