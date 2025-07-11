// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rinvgamma_rcpp
double rinvgamma_rcpp(const double& a, const double& b);
RcppExport SEXP _BayCausal_rinvgamma_rcpp(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvgamma_rcpp(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rmvn_rcpp
arma::mat rmvn_rcpp(const int& n, arma::vec& mean, arma::mat& sigma);
RcppExport SEXP _BayCausal_rmvn_rcpp(SEXP nSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvn_rcpp(n, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// rinvgaussian_rcpp
double rinvgaussian_rcpp(const double& mu, const double& lambda);
RcppExport SEXP _BayCausal_rinvgaussian_rcpp(SEXP muSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvgaussian_rcpp(mu, lambda));
    return rcpp_result_gen;
END_RCPP
}
// compute_Y_tilde_q
arma::vec compute_Y_tilde_q(const int& q, const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const int& Q, const int& S, const int& P, const arma::mat& Y, const arma::mat& X, const bool& minus_mu, const bool& minus_BY, const bool& minus_AX, const bool& minus_LC);
RcppExport SEXP _BayCausal_compute_Y_tilde_q(SEXP qSEXP, SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP QSEXP, SEXP SSEXP, SEXP PSEXP, SEXP YSEXP, SEXP XSEXP, SEXP minus_muSEXP, SEXP minus_BYSEXP, SEXP minus_AXSEXP, SEXP minus_LCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const int& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const int& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const int& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const bool& >::type minus_mu(minus_muSEXP);
    Rcpp::traits::input_parameter< const bool& >::type minus_BY(minus_BYSEXP);
    Rcpp::traits::input_parameter< const bool& >::type minus_AX(minus_AXSEXP);
    Rcpp::traits::input_parameter< const bool& >::type minus_LC(minus_LCSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Y_tilde_q(q, B, A, L, C, mu, Q, S, P, Y, X, minus_mu, minus_BY, minus_AX, minus_LC));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglik_B_rcpp
double compute_loglik_B_rcpp(const int& q, const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_compute_loglik_B_rcpp(SEXP qSEXP, SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglik_B_rcpp(q, B, A, L, C, mu, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// compute_marginal_loglik_rcpp
double compute_marginal_loglik_rcpp(const arma::mat& sparsity_matrix, const int& q, const arma::mat& C, const arma::mat& tau, const arma::mat& L_0, const arma::mat& Y, const arma::mat& B, const arma::mat& X, const arma::mat& A, const arma::vec& mu, const double& a_sigma, const double& b_sigma);
RcppExport SEXP _BayCausal_compute_marginal_loglik_rcpp(SEXP sparsity_matrixSEXP, SEXP qSEXP, SEXP CSEXP, SEXP tauSEXP, SEXP L_0SEXP, SEXP YSEXP, SEXP BSEXP, SEXP XSEXP, SEXP ASEXP, SEXP muSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sparsity_matrix(sparsity_matrixSEXP);
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L_0(L_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type b_sigma(b_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_marginal_loglik_rcpp(sparsity_matrix, q, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglik_rjmcmc_cpp
double compute_loglik_rjmcmc_cpp(int& q, const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_compute_loglik_rjmcmc_cpp(SEXP qSEXP, SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglik_rjmcmc_cpp(q, B, A, L, C, mu, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_mu_rcpp
arma::mat update_mu_rcpp(const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::mat& tau, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_update_mu_rcpp(SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP tauSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(update_mu_rcpp(B, A, L, C, tau, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_A_rcpp
arma::mat update_A_rcpp(const arma::mat& gamma_alpha, const arma::mat& nu_alpha, const arma::mat& B, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const arma::mat& tau, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_update_A_rcpp(SEXP gamma_alphaSEXP, SEXP nu_alphaSEXP, SEXP BSEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP tauSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma_alpha(gamma_alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nu_alpha(nu_alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(update_A_rcpp(gamma_alpha, nu_alpha, B, L, C, mu, tau, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_L_free_rcpp
arma::mat update_L_free_rcpp(const arma::vec& mu, const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& tau, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_update_L_free_rcpp(SEXP muSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP tauSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(update_L_free_rcpp(mu, A, B, C, tau, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_C_rcpp
arma::mat update_C_rcpp(const arma::vec& mu, const arma::mat& A, const arma::mat& B, const arma::mat& L, const arma::mat& tau, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_update_C_rcpp(SEXP muSEXP, SEXP ASEXP, SEXP BSEXP, SEXP LSEXP, SEXP tauSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(update_C_rcpp(mu, A, B, L, tau, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_tau_rcpp
arma::mat update_tau_rcpp(const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const arma::vec& sigma2, const arma::mat& Y, const arma::mat& X);
RcppExport SEXP _BayCausal_update_tau_rcpp(SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP sigma2SEXP, SEXP YSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(update_tau_rcpp(B, A, L, C, mu, sigma2, Y, X));
    return rcpp_result_gen;
END_RCPP
}
// update_sigma2_rcpp
arma::vec update_sigma2_rcpp(const arma::mat& B, const arma::mat& A, const arma::mat& L, const arma::mat& C, const arma::vec& mu, const arma::mat& tau, const arma::mat& Y, const arma::mat& X, const double& a_sigma, const double& b_sigma);
RcppExport SEXP _BayCausal_update_sigma2_rcpp(SEXP BSEXP, SEXP ASEXP, SEXP LSEXP, SEXP CSEXP, SEXP muSEXP, SEXP tauSEXP, SEXP YSEXP, SEXP XSEXP, SEXP a_sigmaSEXP, SEXP b_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type a_sigma(a_sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type b_sigma(b_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_sigma2_rcpp(B, A, L, C, mu, tau, Y, X, a_sigma, b_sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayCausal_rinvgamma_rcpp", (DL_FUNC) &_BayCausal_rinvgamma_rcpp, 2},
    {"_BayCausal_rmvn_rcpp", (DL_FUNC) &_BayCausal_rmvn_rcpp, 3},
    {"_BayCausal_rinvgaussian_rcpp", (DL_FUNC) &_BayCausal_rinvgaussian_rcpp, 2},
    {"_BayCausal_compute_Y_tilde_q", (DL_FUNC) &_BayCausal_compute_Y_tilde_q, 15},
    {"_BayCausal_compute_loglik_B_rcpp", (DL_FUNC) &_BayCausal_compute_loglik_B_rcpp, 9},
    {"_BayCausal_compute_marginal_loglik_rcpp", (DL_FUNC) &_BayCausal_compute_marginal_loglik_rcpp, 12},
    {"_BayCausal_compute_loglik_rjmcmc_cpp", (DL_FUNC) &_BayCausal_compute_loglik_rjmcmc_cpp, 9},
    {"_BayCausal_update_mu_rcpp", (DL_FUNC) &_BayCausal_update_mu_rcpp, 8},
    {"_BayCausal_update_A_rcpp", (DL_FUNC) &_BayCausal_update_A_rcpp, 10},
    {"_BayCausal_update_L_free_rcpp", (DL_FUNC) &_BayCausal_update_L_free_rcpp, 8},
    {"_BayCausal_update_C_rcpp", (DL_FUNC) &_BayCausal_update_C_rcpp, 8},
    {"_BayCausal_update_tau_rcpp", (DL_FUNC) &_BayCausal_update_tau_rcpp, 8},
    {"_BayCausal_update_sigma2_rcpp", (DL_FUNC) &_BayCausal_update_sigma2_rcpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayCausal(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
