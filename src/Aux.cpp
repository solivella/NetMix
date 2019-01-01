#include "Aux.hpp"
#include "MMModelClass.hpp"

//[[Rcpp::export]]
Rcpp::IntegerMatrix getZ(Rcpp::NumericMatrix pmat)
{
  int NROW = pmat.nrow();
  int NCOL = pmat.ncol();
  int mflag, bloc;
  double u, acc;
  Rcpp::NumericVector cprob(NROW); 
  Rcpp::IntegerMatrix res(NROW, NCOL);
  for(int i = 0; i < NCOL; ++i){
    u = R::runif(0, 1);
    acc = 0.0;
    for(int j = 0; j < NROW; ++j){
      acc += pmat(j, i);
      cprob[j] = acc;
    }
    bloc = findInterval(&(cprob[0]), NROW, u, FALSE, FALSE, 0, &mflag);
    res(bloc, i) = 1;
  }
  return(res);
}

//[[Rcpp::export]]
Rcpp::NumericMatrix approxBdyad(Rcpp::NumericVector y,
                                Rcpp::IntegerMatrix node_id,
                                Rcpp::NumericMatrix pi_mat,
                                bool directed)
{
  int N_BLK = pi_mat.nrow();
  int N_DYAD = y.size();
  int s, r;//, z_1, z_2, check;
  double prob_temp;//u_1, u_2;
  //Rcpp::NumericVector cs(N_BLK);
  Rcpp::NumericMatrix den(N_BLK, N_BLK), num(N_BLK, N_BLK);

  for(int dyad = 0; dyad < N_DYAD; ++dyad){
    s = node_id(dyad, 0);
    r = node_id(dyad, 1);
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        if(!directed){
          if(h < g) {
            continue;
          }
        }
        prob_temp = pi_mat(g, s) * pi_mat(h, r);
        num(h, g) += y(dyad) * prob_temp;
        den(h, g) += prob_temp;
        if(!directed){
          num(g, h) += y(dyad) * prob_temp;
          den(g, h) += prob_temp;
        }
      }
    }
  }
  std::transform(den.begin(), den.end(), den.begin(),
                 [](double& x){return x < 1e-12 ? 1e-12 : x;});
  
  Rcpp::NumericMatrix B_t(N_BLK, N_BLK);
  std::transform(num.begin(), num.end(),
                 den.begin(),
                 B_t.begin(),
                 std::divides<double>());
  return B_t;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix approxB(Rcpp::NumericMatrix y,
                      Rcpp::NumericMatrix pi_mat,
                      bool directed)
{
  int N_BLK = pi_mat.nrow();
  int N_NODE = y.ncol();
  Rcpp::NumericMatrix den(N_BLK, N_BLK), num(N_BLK, N_BLK);
  Rcpp::NumericMatrix B_t(N_BLK, N_BLK);
  Rcpp::NumericVector cs(N_BLK);
  int z_1, z_2, check;
  double u_1, u_2, prob_temp;
  for(int col = 0; col < N_NODE; ++col){
    for(int row = 0; row < N_NODE; ++row){
      if(!directed){
        if(row < col){
          continue;
        }
      }
      u_1 = R::runif(0, 1);
      u_2 = R::runif(0, 1);
      cs = Rcpp::cumsum(pi_mat.column(row)).get();
      z_1 = findInterval(&cs[0], N_BLK, u_1, FALSE, FALSE, 0, &check);
      cs = Rcpp::cumsum(pi_mat.column(col)).get();
      z_2 = findInterval(&cs[0], N_BLK, u_2, FALSE, FALSE, 0, &check);
      //probsum = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        for(int h = 0; h < N_BLK; ++h){
          if(!directed){
            if(h < g) {
              continue;
            }
          }
          prob_temp =  (g == z_1) * (h == z_2);
          num(h, g) += y(row,col) * prob_temp;
          den(h, g) += prob_temp;
          if(!directed){
            num(g, h) += y(row,col) * prob_temp;
            den(g, h) += prob_temp;
          }
        }
      }
    }
  }
  std::transform(den.begin(), den.end(), den.begin(),
                 [](double& x){return fabs(x) < 1e-10 ? 1e-10 : x;});
  
  std::transform(num.begin(), num.end(),
                 den.begin(),
                 B_t.begin(),
                 std::divides<double>());
  return B_t;
}


double digamma_approx(double x)
{
  double p;
  x=x+6;
  p=1/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
    0.008333333333333)*p-0.083333333333333)*p;
  p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
  return p;
}

double lgamma_approx(double x)
{
  double z=1/(x*x);
  
  x=x+6;
  z=(((-0.000595238095238*z+0.000793650793651)
        *z-0.002777777777778)*z+0.083333333333333)/x;
  z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
  log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  return z;
}

double lgammaDiff(double alpha, double C) {
  return lgamma(alpha + C) - lgamma(alpha);
}
double digammaDiff(double alpha, double C) {
  return R::digamma(alpha + C) - R::digamma(alpha);
}


double logSumExp(const std::vector<double>& invec)
{
  double offset = *std::max_element(invec.begin(), invec.end());
  double res = std::accumulate(invec.begin(), invec.end(),
                               0.0,
                               [&offset](double a, const double b){
                                 return a + exp(b - offset);
                               });
  return offset + log(res);
}


template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}



