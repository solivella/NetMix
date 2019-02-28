#include "approxB.hpp"


// [[Rcpp::export]]
NumericMatrix approxB(NumericVector y,
                      IntegerMatrix d_id,
                      NumericMatrix pi_mat)
{
  int N_BLK = pi_mat.nrow();
  int N_DYAD = d_id.nrow();
  NumericMatrix den(N_BLK, N_BLK), num(N_BLK, N_BLK), B_t(N_BLK, N_BLK);
  //NumericMatrix prob(N_BLK, N_BLK);
  int s, r;
  double prob_temp;
  //double probsum, rd;
  for(int d = 0; d < N_DYAD; ++d){
    s = d_id(d, 0);
    r = d_id(d, 1);
    //probsum = 0.0;
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        prob_temp = pi_mat(g, s) * pi_mat(h, r);
        num(h, g) += y[d] * prob_temp;
        den(h, g) += prob_temp;
        //probsum += prob(h, g);
      }
    }
    //rd = runif(0,1);
    // for(int g = 0; g < N_BLK; ++g){
    //   for(int h = 0; h < N_BLK; ++h){
    //     prob_temp = prob(h, g) / probsum;
    //     if(rd < prob_temp){
    //       num(h, g) += y[d];
    //       den(h, g)++;
    //       break;
    //     } 
    //     rd -= prob_temp;
    //   }
    // }
  }
  std::transform(num.begin(), num.end(),
                 den.begin(),
                 B_t.begin(),
                 std::divides<double>());
  return B_t;
}
