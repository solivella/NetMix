#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


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



double lgammaDiff2(double alpha, double C) {
  return lgamma(alpha + C) - lgamma(alpha);
}

// [[Rcpp::export]]
double alphaLB(NumericMatrix alpha_term,
               NumericMatrix alpha_par,
               double xi_param,
               const int N_STATE,
               const int N_BLK,
               const int N_MONAD_PRED,
               const int N_NODE,
               const NumericMatrix x_t,
               NumericMatrix beta,
               NumericMatrix e_xb,
               NumericMatrix sum_e_xb,
               NumericMatrix e_c_t,
               const NumericMatrix time_id_node,
               NumericMatrix kappa_t,
               NumericMatrix xi,
               std::vector<double> sum_c,
               NumericMatrix alpha,
               const double var_beta
)
{
  
  Array<double> sub;
  
  /**
   ALPHA COMPUTATION
   computeAlpha();
  */
  std::fill(alpha_term.begin(), alpha_term.end(), 0.0);
  
  xi_param = exp(alpha_par[0]);
  //Rprintf("xi_param %f here!\n", xi_param);
  std::vector<double>::iterator it = alpha_par.begin() + 1;
  for(int m = 0; m < N_STATE; ++m){
    for(int g = 1; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x, ++it){
        beta(x, g, m) = *it;
      }
    }
  }
  
  
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      double linpred, linpred_sum;
      linpred_sum = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(int x = 0; x < N_MONAD_PRED; ++x){
          //Rprintf("beta(%i, %i, %i) is %f\n",x, g, m,beta(x, g, m));
          linpred += x_t(x, p) * beta(x, g, m);
        }
        
        linpred = exp(linpred);
        
        e_xb(g, p, m) = linpred;
        alpha(g, p, m) = linpred;
        linpred_sum += linpred;
      }
      sum_e_xb(p, m) = linpred_sum;
      for(int g = 0; g < N_BLK; ++g) {
        if(!std::isfinite(alpha(g, p, m))){
          Rprintf("Linpredsum is %f, xi_param is %f\n",linpred_sum, xi_param);
          stop("alpha is nan\n");
        } 
        alpha(g, p, m) /= linpred_sum;
        alpha(g, p, m) *= xi_param;
        alpha_term(m, time_id_node[p]) += lgammaDiff2(alpha(g, p, m), e_c_t(g, p));
      }
    }
  }
  
  double res = 0.0;
  double kappa_val = 0.0;
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      kappa_val = kappa_t(m, time_id_node[p]);
      res += kappa_val * (lgammaDiff2(xi(p, m), sum_c[p]));
      for(int g = 0; g < N_BLK; ++g){
        res -=  kappa_val * lgammaDiff2(alpha(g, p, m),  e_c_t(g, p));
      }
    }
  }
  for(int m = 0; m < N_STATE; ++m){
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 1; x < N_MONAD_PRED; ++x){ // Don't regularize the intercept
        res += pow(beta(x, g, m), 2.0) / (2.0 * var_beta);
      }
    }
  }
  return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
