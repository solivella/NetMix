#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

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
               const int N_PAR,
               const int N_STATE,
               const int N_BLK,
               const int N_MONAD_PRED,
               const int N_NODE,
               const int N_TIME, 
               const NumericMatrix x_t_r,
               const NumericMatrix e_c_t,
               const NumericVector time_id_node,
               const NumericMatrix kappa_t_r
               
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
  
  
  double linpred;
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      for(int g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(int x = 0; x < N_MONAD_PRED; ++x){
          linpred += x_t(x, p) * beta(x, g, m);
        }
        linpred = exp(linpred);
        alpha(g, p, m) = linpred;
      }
    }
  }
  
  double res = 0.0, res_int, alpha_row = 0.0, alpha_val = 0.0;
  #pragma omp parallel for firstprivate(alpha_row, alpha_val, res_int) reduction(+:res)
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      alpha_row = 0.0;
      res_int = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        alpha_val = alpha(g, p, m);
        alpha_row += alpha_val;
        res_int += lgamma(alpha_val + e_c_t(g, p)) - lgamma(alpha_val);
      }
      res_int += lgamma(alpha_row) - lgamma(alpha_row + sum_c[p]);
      res += kappa_t(m, time_id_node[p]) * res_int;
    }
    
    //Prior for beta 
    for(int g = 0; g < N_BLK; ++g){
      res -= log(1.0 + pow(beta(0, g, m),2.0)/100);
      for(int x = 1; x < N_MONAD_PRED; ++x){
          res -= log(1.0 + pow(beta(x, g, m), 2.0)/1.5625);
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
                      const NumericMatrix kappa_t_r
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
  
  NumericVector gr(N_PAR);
  double res=0.0,  alpha_row = 0.0, prior_gr;
  for(int m = 0; m < N_STATE; ++m){
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x){
        res = 0.0;
#pragma omp parallel for firstprivate(alpha_row) reduction(+:res)
        for(int p = 0; p < N_NODE; ++p){
          alpha_row = 0.0;
          for(int h = 0; h < N_BLK; ++h){
            alpha_row += alpha(h, p, m);
          }
          res += (R::digamma(alpha_row) - R::digamma(alpha_row + sum_c[p])
                    + R::digamma(alpha(g, p, m) + e_c_t(g, p)) - R::digamma(alpha(g, p, m)))
            * kappa_t(m,  time_id_node[p]) * alpha(g, p, m) * x_t(x, p);
        }
        if(x!=0){
          prior_gr = 2 * beta(x, g, m) / (pow(beta(x, g, m), 2.0) + 1.5625);
        } else {
          prior_gr =  2 * (beta(x, g, m)) / (pow(beta(x, g, m), 2.0) + 100);
        }
        
        gr[x + N_MONAD_PRED * (g + N_BLK * m)] = -(res - prior_gr);
      }
    }
  }
  return(gr);
}



/*** R
setwd("~/Dropbox/GitHub/NetMixRoot/NetMix/")
source("Extra/NetGenerator.R")
library(numDeriv)
library(plyr)
library(doParallel)
cl <- makePSOCKcluster(22)
registerDoParallel(cl)
#set.seed(831213)

N_BLK <- 3
N_STATE <- 1
N_MONAD_PRED <- 4
N <- 300
#deriv.diffs <- rep(NA, 500)

N_sim <- 1
diffs <- array(NA, c(N_sim, N_MONAD_PRED*N_BLK))
true_beta <- array(c(3.05, rnorm(3, 0, 3),
                     3.05, rnorm(3, 0, 3),
                     3.05, rnorm(3, 0, 3)),
                   c(4, 3, 1))
Z <- list(cbind(runif(N^2), runif(N^2), runif(N^2)))  
X <- cbind(runif(N), runif(N), runif(N))
X_s <- list(cbind(1,scale(X)))
  
system.time(net2_list <-  llply(1:N_sim, 
                        function(x){
                          return(NetSim(BLK = 3
                                 ,NODE = N
                                 ,STATE = 1
                                 ,TIME = 1 
                                 ,DIRECTED = TRUE
                                 ,N_PRED = 3
                                 ,B_t = diag(4.5, 3,3) - 1.5
                                 ,beta_arr = true_beta
                                 ,gamma_vec = c(1.05, -1.05, 1.05)
                                 ,X = X_s
                                 ,Z = Z))}
                                 ,.parallel = TRUE
                        ,.paropts = list(.export=c("NetSim","true_beta","N", "X_s", "Z"))
                        ))
stopCluster(cl)
for(i in 1:N_sim){
  C_t <-  t(net2_list[[i]]$pi_vecs[[1]])*(N*2)
  #  beta_rand <- rnorm(N_MONAD_PRED*N_BLK)
  # (d1 <- grad(alphaLB, beta_rand,
  #            N_PAR =  N_MONAD_PRED * N_BLK * N_STATE,
  #            N_STATE = N_STATE,
  #            N_BLK = N_BLK,
  #            N_MONAD_PRED = N_MONAD_PRED,
  #            N_NODE = N,
  #            N_TIME = 1,
  #            x_t_r = t(net2_list[[1]]$X[[1]]),
  #            e_c_t = C_t,
  #            time_id_node = rep(0, N),
  #            kappa_t_r = matrix(1, ncol = 1, nrow = 1)))
  # 
  # (d2 <- alphaGr(beta_rand,
  #               N_MONAD_PRED * N_BLK * N_STATE,
  #               N_STATE,
  #               N_BLK,
  #               N_MONAD_PRED,
  #               N_NODE = N,
  #               N_TIME = 1,
  #               x_t_r = t(net2_list[[1]]$X[[1]]),
  #               e_c_t = C_t,
  #               time_id_node = rep(0, N),
  #               kappa_t_r = matrix(1, ncol = 1, nrow = 1)))
  #deriv.diffs[i] <- mean((d1 - d2) / d1)
  
  bfgs_res <- optim(#c(net2_list[[1]]$beta_arr[,,1]),
                    # c(lm.fit(net2_list[[1]]$X[[1]], log(t((C_t+1)/(N*2))))$coef), 
                     #rnorm(N_MONAD_PRED*N_BLK),
                     rep(0, N_MONAD_PRED*N_BLK),
                     alphaLB,
                     gr = alphaGr,
                     N_STATE = N_STATE,
                     N_PAR = N_MONAD_PRED * N_BLK * N_STATE,
                     N_BLK = N_BLK,
                     N_MONAD_PRED = N_MONAD_PRED,
                     N_NODE = N,
                     N_TIME = 1, 
                     x_t_r = t(net2_list[[i]]$X[[1]]),
                     e_c_t = C_t,
                     time_id_node = rep(0, N),
                     kappa_t_r = matrix(1, ncol = 1, nrow = 1),
                     lower=-10,
                      upper=10,
                     control=list(maxit=10000, trace=6)
                     ,method = "L-BFGS-B"
  )
  array(bfgs_res$par, c(4, 3))
  net2_list[[i]]$beta_arr
  
  diffs[i,] <- c((net2_list[[i]]$beta_arr[,,1]) - array(bfgs_res$par, c(4, 3)))
}
pdf("~/Desktop/Sims.pdf")
par(mfrow=c(3,4))
for(i in 1:ncol(diffs)){
  hist(diffs[,i], main="")
}
dev.off()

#hist(deriv.diffs, breaks=500)
*/