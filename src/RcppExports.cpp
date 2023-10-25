// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pooled_cv_rcpp
std::list<double> pooled_cv_rcpp(arma::mat mtr, std::vector<std::string> group);
RcppExport SEXP _pmartR_pooled_cv_rcpp(SEXP mtrSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mtr(mtrSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(pooled_cv_rcpp(mtr, group));
    return rcpp_result_gen;
END_RCPP
}
// unpooled_cv_rcpp
std::list<double> unpooled_cv_rcpp(NumericMatrix mtr);
RcppExport SEXP _pmartR_unpooled_cv_rcpp(SEXP mtrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mtr(mtrSEXP);
    rcpp_result_gen = Rcpp::wrap(unpooled_cv_rcpp(mtr));
    return rcpp_result_gen;
END_RCPP
}
// count_missing_cpp
NumericMatrix count_missing_cpp(NumericMatrix data, NumericVector gp);
RcppExport SEXP _pmartR_count_missing_cpp(SEXP dataSEXP, SEXP gpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gp(gpSEXP);
    rcpp_result_gen = Rcpp::wrap(count_missing_cpp(data, gp));
    return rcpp_result_gen;
END_RCPP
}
// fold_change_diff
arma::mat fold_change_diff(arma::mat data, arma::mat C);
RcppExport SEXP _pmartR_fold_change_diff(SEXP dataSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(fold_change_diff(data, C));
    return rcpp_result_gen;
END_RCPP
}
// fold_change_diff_row
arma::rowvec fold_change_diff_row(arma::rowvec means, arma::mat C);
RcppExport SEXP _pmartR_fold_change_diff_row(SEXP meansSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type means(meansSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(fold_change_diff_row(means, C));
    return rcpp_result_gen;
END_RCPP
}
// fold_change_ratio
arma::mat fold_change_ratio(arma::mat data, arma::mat C);
RcppExport SEXP _pmartR_fold_change_ratio(SEXP dataSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(fold_change_ratio(data, C));
    return rcpp_result_gen;
END_RCPP
}
// fold_change_diff_na_okay
arma::mat fold_change_diff_na_okay(arma::mat data, arma::mat C);
RcppExport SEXP _pmartR_fold_change_diff_na_okay(SEXP dataSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(fold_change_diff_na_okay(data, C));
    return rcpp_result_gen;
END_RCPP
}
// anova_cpp
List anova_cpp(arma::mat data, NumericVector gp, int unequal_var, arma::mat X, arma::mat Beta, arma::mat pred_grid, arma::uvec continuous_covar_inds, int n_covar_levels);
RcppExport SEXP _pmartR_anova_cpp(SEXP dataSEXP, SEXP gpSEXP, SEXP unequal_varSEXP, SEXP XSEXP, SEXP BetaSEXP, SEXP pred_gridSEXP, SEXP continuous_covar_indsSEXP, SEXP n_covar_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gp(gpSEXP);
    Rcpp::traits::input_parameter< int >::type unequal_var(unequal_varSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Beta(BetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pred_grid(pred_gridSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type continuous_covar_inds(continuous_covar_indsSEXP);
    Rcpp::traits::input_parameter< int >::type n_covar_levels(n_covar_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(anova_cpp(data, gp, unequal_var, X, Beta, pred_grid, continuous_covar_inds, n_covar_levels));
    return rcpp_result_gen;
END_RCPP
}
// two_factor_anova_cpp
List two_factor_anova_cpp(arma::mat y, arma::mat X_full, arma::mat X_red, arma::colvec group_ids, arma::mat pred_grid_full, arma::mat pred_grid_red, arma::uvec continuous_covar_inds, int n_covar_levels);
RcppExport SEXP _pmartR_two_factor_anova_cpp(SEXP ySEXP, SEXP X_fullSEXP, SEXP X_redSEXP, SEXP group_idsSEXP, SEXP pred_grid_fullSEXP, SEXP pred_grid_redSEXP, SEXP continuous_covar_indsSEXP, SEXP n_covar_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_full(X_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_red(X_redSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type group_ids(group_idsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pred_grid_full(pred_grid_fullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pred_grid_red(pred_grid_redSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type continuous_covar_inds(continuous_covar_indsSEXP);
    Rcpp::traits::input_parameter< int >::type n_covar_levels(n_covar_levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(two_factor_anova_cpp(y, X_full, X_red, group_ids, pred_grid_full, pred_grid_red, continuous_covar_inds, n_covar_levels));
    return rcpp_result_gen;
END_RCPP
}
// group_comparison_anova_cpp
List group_comparison_anova_cpp(arma::mat means, arma::mat data, arma::mat sizes, arma::mat which_xmatrix, arma::mat Xfull, arma::mat Xred, arma::mat Cfull, arma::mat Cred, arma::mat Cmu);
RcppExport SEXP _pmartR_group_comparison_anova_cpp(SEXP meansSEXP, SEXP dataSEXP, SEXP sizesSEXP, SEXP which_xmatrixSEXP, SEXP XfullSEXP, SEXP XredSEXP, SEXP CfullSEXP, SEXP CredSEXP, SEXP CmuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type means(meansSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sizes(sizesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type which_xmatrix(which_xmatrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xfull(XfullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xred(XredSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cfull(CfullSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cred(CredSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Cmu(CmuSEXP);
    rcpp_result_gen = Rcpp::wrap(group_comparison_anova_cpp(means, data, sizes, which_xmatrix, Xfull, Xred, Cfull, Cred, Cmu));
    return rcpp_result_gen;
END_RCPP
}
// holm_cpp
NumericVector holm_cpp(NumericVector ps);
RcppExport SEXP _pmartR_holm_cpp(SEXP psSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ps(psSEXP);
    rcpp_result_gen = Rcpp::wrap(holm_cpp(ps));
    return rcpp_result_gen;
END_RCPP
}
// ptukey_speed
NumericMatrix ptukey_speed(NumericMatrix qstats, NumericVector sizes);
RcppExport SEXP _pmartR_ptukey_speed(SEXP qstatsSEXP, SEXP sizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type qstats(qstatsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sizes(sizesSEXP);
    rcpp_result_gen = Rcpp::wrap(ptukey_speed(qstats, sizes));
    return rcpp_result_gen;
END_RCPP
}
// compute_betas
arma::mat compute_betas(arma::mat data_mat, arma::mat Xmatrix);
RcppExport SEXP _pmartR_compute_betas(SEXP data_matSEXP, SEXP XmatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data_mat(data_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmatrix(XmatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_betas(data_mat, Xmatrix));
    return rcpp_result_gen;
END_RCPP
}
// kw_rcpp
std::list<double> kw_rcpp(arma::mat mtr, std::vector<std::string> group);
RcppExport SEXP _pmartR_kw_rcpp(SEXP mtrSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mtr(mtrSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(kw_rcpp(mtr, group));
    return rcpp_result_gen;
END_RCPP
}
// nonmissing_per_grp
arma::Mat<int> nonmissing_per_grp(arma::mat mtr, std::vector<std::string> group);
RcppExport SEXP _pmartR_nonmissing_per_grp(SEXP mtrSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mtr(mtrSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(nonmissing_per_grp(mtr, group));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pmartR_pooled_cv_rcpp", (DL_FUNC) &_pmartR_pooled_cv_rcpp, 2},
    {"_pmartR_unpooled_cv_rcpp", (DL_FUNC) &_pmartR_unpooled_cv_rcpp, 1},
    {"_pmartR_count_missing_cpp", (DL_FUNC) &_pmartR_count_missing_cpp, 2},
    {"_pmartR_fold_change_diff", (DL_FUNC) &_pmartR_fold_change_diff, 2},
    {"_pmartR_fold_change_diff_row", (DL_FUNC) &_pmartR_fold_change_diff_row, 2},
    {"_pmartR_fold_change_ratio", (DL_FUNC) &_pmartR_fold_change_ratio, 2},
    {"_pmartR_fold_change_diff_na_okay", (DL_FUNC) &_pmartR_fold_change_diff_na_okay, 2},
    {"_pmartR_anova_cpp", (DL_FUNC) &_pmartR_anova_cpp, 8},
    {"_pmartR_two_factor_anova_cpp", (DL_FUNC) &_pmartR_two_factor_anova_cpp, 8},
    {"_pmartR_group_comparison_anova_cpp", (DL_FUNC) &_pmartR_group_comparison_anova_cpp, 9},
    {"_pmartR_holm_cpp", (DL_FUNC) &_pmartR_holm_cpp, 1},
    {"_pmartR_ptukey_speed", (DL_FUNC) &_pmartR_ptukey_speed, 2},
    {"_pmartR_compute_betas", (DL_FUNC) &_pmartR_compute_betas, 2},
    {"_pmartR_kw_rcpp", (DL_FUNC) &_pmartR_kw_rcpp, 2},
    {"_pmartR_nonmissing_per_grp", (DL_FUNC) &_pmartR_nonmissing_per_grp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_pmartR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
