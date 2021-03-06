// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// approxB
Rcpp::NumericMatrix approxB(Rcpp::NumericVector y, Rcpp::IntegerMatrix d_id, Rcpp::NumericMatrix pi_mat, bool directed);
RcppExport SEXP _NetMix_approxB(SEXP ySEXP, SEXP d_idSEXP, SEXP pi_matSEXP, SEXP directedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type d_id(d_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pi_mat(pi_matSEXP);
    Rcpp::traits::input_parameter< bool >::type directed(directedSEXP);
    rcpp_result_gen = Rcpp::wrap(approxB(y, d_id, pi_mat, directed));
    return rcpp_result_gen;
END_RCPP
}
// getZ
Rcpp::IntegerMatrix getZ(Rcpp::NumericMatrix pi_mat);
RcppExport SEXP _NetMix_getZ(SEXP pi_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pi_mat(pi_matSEXP);
    rcpp_result_gen = Rcpp::wrap(getZ(pi_mat));
    return rcpp_result_gen;
END_RCPP
}
// alphaLBound
double alphaLBound(arma::vec par, arma::uvec tot_nodes, arma::umat c_t, arma::mat x_t, arma::umat s_mat, arma::uvec t_id, arma::cube var_beta, arma::cube mu_beta);
RcppExport SEXP _NetMix_alphaLBound(SEXP parSEXP, SEXP tot_nodesSEXP, SEXP c_tSEXP, SEXP x_tSEXP, SEXP s_matSEXP, SEXP t_idSEXP, SEXP var_betaSEXP, SEXP mu_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tot_nodes(tot_nodesSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type c_t(c_tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_t(x_tSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type s_mat(s_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t_id(t_idSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type var_beta(var_betaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mu_beta(mu_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(alphaLBound(par, tot_nodes, c_t, x_t, s_mat, t_id, var_beta, mu_beta));
    return rcpp_result_gen;
END_RCPP
}
// alphaGrad
arma::vec alphaGrad(arma::vec par, arma::uvec tot_nodes, arma::umat c_t, arma::mat x_t, arma::umat s_mat, arma::uvec t_id, arma::cube var_beta, arma::cube mu_beta);
RcppExport SEXP _NetMix_alphaGrad(SEXP parSEXP, SEXP tot_nodesSEXP, SEXP c_tSEXP, SEXP x_tSEXP, SEXP s_matSEXP, SEXP t_idSEXP, SEXP var_betaSEXP, SEXP mu_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tot_nodes(tot_nodesSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type c_t(c_tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_t(x_tSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type s_mat(s_matSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t_id(t_idSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type var_beta(var_betaSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type mu_beta(mu_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(alphaGrad(par, tot_nodes, c_t, x_t, s_mat, t_id, var_beta, mu_beta));
    return rcpp_result_gen;
END_RCPP
}
// mmsbm_fit
Rcpp::List mmsbm_fit(const arma::mat& z_t, const arma::mat& x_t, const arma::vec& y, const arma::uvec& time_id_dyad, const arma::uvec& time_id_node, const arma::uvec& nodes_per_period, const arma::umat& node_id_dyad, const arma::field<arma::uvec>& node_id_period, const arma::mat& mu_b, const arma::mat& var_b, const arma::cube& mu_beta, const arma::cube& var_beta, const arma::vec& mu_gamma, const arma::vec& var_gamma, const arma::mat& pi_init, arma::mat& kappa_init_t, arma::mat& b_init_t, arma::cube& beta_init_r, arma::vec& gamma_init_r, Rcpp::List& control);
RcppExport SEXP _NetMix_mmsbm_fit(SEXP z_tSEXP, SEXP x_tSEXP, SEXP ySEXP, SEXP time_id_dyadSEXP, SEXP time_id_nodeSEXP, SEXP nodes_per_periodSEXP, SEXP node_id_dyadSEXP, SEXP node_id_periodSEXP, SEXP mu_bSEXP, SEXP var_bSEXP, SEXP mu_betaSEXP, SEXP var_betaSEXP, SEXP mu_gammaSEXP, SEXP var_gammaSEXP, SEXP pi_initSEXP, SEXP kappa_init_tSEXP, SEXP b_init_tSEXP, SEXP beta_init_rSEXP, SEXP gamma_init_rSEXP, SEXP controlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type z_t(z_tSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x_t(x_tSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time_id_dyad(time_id_dyadSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time_id_node(time_id_nodeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type nodes_per_period(nodes_per_periodSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type node_id_dyad(node_id_dyadSEXP);
    Rcpp::traits::input_parameter< const arma::field<arma::uvec>& >::type node_id_period(node_id_periodSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mu_b(mu_bSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type var_b(var_bSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type mu_beta(mu_betaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type var_beta(var_betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_gamma(mu_gammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type var_gamma(var_gammaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type pi_init(pi_initSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type kappa_init_t(kappa_init_tSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_init_t(b_init_tSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type beta_init_r(beta_init_rSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type gamma_init_r(gamma_init_rSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type control(controlSEXP);
    rcpp_result_gen = Rcpp::wrap(mmsbm_fit(z_t, x_t, y, time_id_dyad, time_id_node, nodes_per_period, node_id_dyad, node_id_period, mu_b, var_b, mu_beta, var_beta, mu_gamma, var_gamma, pi_init, kappa_init_t, b_init_t, beta_init_r, gamma_init_r, control));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NetMix_approxB", (DL_FUNC) &_NetMix_approxB, 4},
    {"_NetMix_getZ", (DL_FUNC) &_NetMix_getZ, 1},
    {"_NetMix_alphaLBound", (DL_FUNC) &_NetMix_alphaLBound, 8},
    {"_NetMix_alphaGrad", (DL_FUNC) &_NetMix_alphaGrad, 8},
    {"_NetMix_mmsbm_fit", (DL_FUNC) &_NetMix_mmsbm_fit, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_NetMix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
