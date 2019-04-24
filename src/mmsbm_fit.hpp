#include "MMModelClass.hpp"


Rcpp::List mmsbm_fit(const Rcpp::NumericMatrix& z_t,
                     const Rcpp::NumericMatrix& x_t,
                     const Rcpp::IntegerVector& y,
                     const Rcpp::IntegerVector& time_id_dyad,
                     const Rcpp::IntegerVector& time_id_node,
                     const Rcpp::IntegerVector& nodes_per_period,
                     const Rcpp::IntegerMatrix& node_id_dyad,
                     const Rcpp::NumericMatrix& mu_b,
                     const Rcpp::NumericMatrix& var_b,
                     const Rcpp::NumericMatrix& phi_init,
                     Rcpp::NumericMatrix& kappa_init_t,
                     Rcpp::NumericMatrix& b_init_t,
                     Rcpp::NumericVector& beta_init,
                     Rcpp::NumericVector& gamma_init,
                     Rcpp::List control
		     );
