#include <Rcpp.h>
using Rcpp::NumericMatrix;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::CharacterVector;
using Rcpp::CharacterMatrix;
using R::runif;

NumericMatrix approxB(IntegerVector y,
                      IntegerMatrix d_id,
                      NumericMatrix pi_mat);

