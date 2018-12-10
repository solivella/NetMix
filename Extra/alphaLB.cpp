#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

template<typename T>
class Array //traverse in order of indeces (i.e. i fastest)
{
public:
  template <typename Source>
  Array(std::initializer_list<int> dim, const Source& source)
    : dims(dim),
      data(source.begin(), source.end())
  {
  }
  Array(std::initializer_list<int> dim, T val)
    : dims(dim),
      data(std::accumulate(dims.begin(),dims.end(), 1, std::multiplies<int>()), val)
  {
  }
  //typedef T* iterator;
  typename std::vector<T>::iterator begin(){
    return data.begin();
  }
  typename std::vector<T>::iterator end(){
    return data.end();
  }
  //1d
  T& operator[](int i){
    return (data[i]);
  }
  const T& operator[](int i) const {
    return (data[i]);
  }
  //2d
  T& operator()(int i, int j){//j varies most slowly.
    return (data[i + dims[0] * j]);
  }
  const T& operator()(int i, int j) const {
    return (data[i + dims[0] * j]);
  }
  //3d
  T& operator()(int i, int j, int k){
    return (data[i + dims[0] * (j + dims[1] * k)]);
  }
  const T& operator()(int i, int j, int k) const {
    return (data[i + dims[0] * (j + dims[1] * k)]);
  }
  
private:
  std::vector<int> dims;
  std::vector<T> data;
};



double lgammaDiff(double alpha, double C) {
  return lgamma(alpha + C) - lgamma(alpha);
}

double digammaDiff(double alpha, double C) {
  return R::digamma(alpha + C) - R::digamma(alpha);
}

// [[Rcpp::export]]
double alphaLB(NumericVector alpha_par,
               const int N_STATE,
               const int N_BLK,
               const int N_MONAD_PRED,
               const int N_NODE,
               const int N_TIME, 
               const NumericMatrix x_t_r,
               const NumericMatrix e_c_t,
               const NumericVector time_id_node,
               const NumericMatrix kappa_t_r,
               const double var_beta
)
{
  
  Array<double> x_t({N_MONAD_PRED, N_NODE}, x_t_r);
  Array<double> beta({N_MONAD_PRED, N_BLK, N_STATE}, alpha_par);
  Array<double> alpha({N_BLK, N_NODE, N_STATE}, 0.0);
  Array<double> kappa_t({N_TIME, N_STATE}, kappa_t_r);
  Array<double> sum_c({N_NODE}, 2.0*N_NODE);
  /**
   ALPHA COMPUTATION
   computeAlpha();
  */

  
  double linpred, row_sum;
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      row_sum = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(int x = 0; x < N_MONAD_PRED; ++x){
          linpred += x_t(x, p) * beta(x, g, m);
        }
        linpred = exp(linpred);
        row_sum += linpred;
        alpha(g, p, m) = linpred;
      }
    }
  }
  
  double res = 0.0, alpha_row = 0.0, alpha_val = 0.0;
//#pragma omp parallel for firstprivate(alpha_row, alpha_val) reduction(+:res)
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      alpha_row = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        alpha_val = alpha(g, p, m);
        alpha_row += alpha_val;
        res += lgamma(alpha_val + e_c_t(g, p)) - lgamma(alpha_val);
      }
      res += lgamma(alpha_row) - lgamma(alpha_row + sum_c[p]);
      res *= kappa_t(m, time_id_node[p]);
    }
    
    //Prior for beta
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x){
        res -= 0.5 * pow(beta(x, g, m), 2.0) / var_beta;
      }
    }
  }
  
  
  res *= -1; //VMMIN minimizes.
  return res;
}

// [[Rcpp::export]]
NumericVector alphaGr(NumericVector alpha_par,
                      int N_PAR,
                      const int N_STATE,
                      const int N_BLK,
                      const int N_MONAD_PRED,
                      const int N_NODE,
                      const int N_TIME, 
                      const NumericMatrix x_t_r,
                      const NumericMatrix e_c_t,
                      const NumericVector time_id_node,
                      const NumericMatrix kappa_t_r,
                      const double var_beta)
{
  Array<double> x_t({N_MONAD_PRED, N_NODE}, x_t_r);
  Array<double> beta({N_MONAD_PRED, N_BLK, N_STATE}, alpha_par);
  Array<double> alpha({N_BLK, N_NODE, N_STATE}, 0.0);
  Array<double> kappa_t({N_TIME, N_STATE}, kappa_t_r);
  Array<double> sum_c({N_NODE}, 2.0*N_NODE);
  /**
  ALPHA COMPUTATION
  computeAlpha();
  */
  
 
 double linpred, row_sum;
 for(int m = 0; m < N_STATE; ++m){
   for(int p = 0; p < N_NODE; ++p){
     row_sum = 0.0;
     for(int g = 0; g < N_BLK; ++g){
       linpred = 0.0;
       for(int x = 0; x < N_MONAD_PRED; ++x){
         linpred += x_t(x, p) * beta(x, g, m);
       }
       linpred = exp(linpred);
       row_sum += linpred;
       alpha(g, p, m) = linpred;
     }
   }
 }
  
  NumericVector gr(N_PAR);
 double res,  alpha_row;
 for(int m = 0; m < N_STATE; ++m){
   for(int g = 0; g < N_BLK; ++g){
     for(int x = 0; x < N_MONAD_PRED; ++x){
       res = 0.0;
       for(int p = 0; p < N_NODE; ++p){
         alpha_row = 0.0;
         for(int h = 0; h < N_BLK; ++h){
           alpha_row += alpha(h, p, m);
         }
         res += R::digamma(alpha_row) - R::digamma(alpha_row + sum_c[p]);
         res += R::digamma(alpha(g, p, m) + e_c_t(g, p)) - R::digamma(alpha(g, p, m));
         res *= kappa_t(m,  time_id_node[p]) * alpha(g, p, m) * x_t(x, p);
       }
       gr[x + N_MONAD_PRED * (g + N_BLK * m)] = -(res - beta(x, g, m) / var_beta);
     }
   }
 }
 return(gr);
}



/*** R
setwd("~/Dropbox/GitHub/NetMixRoot/NetMix/")
source("Extra/NetGenerator.R")
library(numDeriv)
set.seed(831213)
net1 <-  NetSim(BLK = 3
                  ,NODE = 100
                  ,STATE = 1
                  ,TIME = 1 
                  ,DIRECTED = TRUE
                  ,N_PRED=0
                  ,B_t = matrix(c(5, rep(-5, 3), 5, rep(-5, 3), 5), 
                                ncol=3)
                  ,beta_arr = array(c(1.25, 1.25, 1.25)*-1,
                                    c(1, 3, 1))
                  ,gamma_vec = c(0))
C_t = t(net1$pi_vecs[[1]])*200
beta_rand <- rnorm(3)
N_MONAD_PRED <- 1
N_BLK <- 3
N_STATE <- 1
grad(alphaLB, beta_rand,
        N_STATE = N_STATE,
        N_BLK = N_BLK,
        N_MONAD_PRED = N_MONAD_PRED,
        N_NODE = 100,
        N_TIME = 1, 
        x_t_r = t(net1$X[[1]]),
        e_c_t = C_t,
        time_id_node = rep(0, 100),
        kappa_t_r = matrix(1, ncol = 1, nrow = 1),
        var_beta = 1)

alphaGr(beta_rand,
        N_MONAD_PRED * N_BLK * N_STATE,
        N_STATE,
        N_BLK,
        N_MONAD_PRED,
        N_NODE = 100,
        N_TIME = 1, 
        x_t_r = t(net1$X[[1]]),
        e_c_t = C_t,
        time_id_node = rep(0, 100),
        kappa_t_r = matrix(1, ncol = 1, nrow = 1),
        var_beta = 1)
*/
