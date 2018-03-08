#include "MMModelClass.hpp"


Rcpp::List mmsbm_fit(const Rcpp::NumericMatrix z_t,
		     const Rcpp::NumericMatrix x_t,
		     const Rcpp::IntegerVector y,
		     Rcpp::List control
		     );
