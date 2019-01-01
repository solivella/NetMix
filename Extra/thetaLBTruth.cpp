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
double thetaLB(NumericVector theta_par,
                const int N_DYAD_PRED,
                const int N_B_PAR,
                const int N_BLK,
                const int N_DYAD,
                const NumericMatrix mu_b_t,
                const NumericMatrix z_t,
                const NumericVector y,
                const IntegerMatrix node_id_dyad,
                const NumericMatrix phi_init,
                const bool directed)
{
  int p, q;
  NumericMatrix b_t(N_BLK, N_BLK);
  NumericVector gamma(N_DYAD_PRED);
  Array<double> send_phi({N_BLK, N_DYAD}, 0.0);
  Array<double> rec_phi({N_BLK, N_DYAD}, 0.0);
  for(int d = 0; d < N_DYAD; ++d){
    p = node_id_dyad(d, 0);
    q = node_id_dyad(d, 1);
    for(int g = 0; g < N_BLK; ++g){
      send_phi(g, d) = phi_init(g, p);
      rec_phi(g, d) = phi_init(g, q);
    }
  }
  Array<int> par_ind({N_BLK, N_BLK}, 0);
  int ind = 0;
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      if(directed){
        par_ind(h, g) = ind++;
      } else {
        if(h >= g){
          par_ind(h, g) = ind++;
        } else {
          par_ind(h, g) = par_ind(g, h);
        }
      }
    }
  }
  Array<double> theta({N_BLK, N_BLK, N_DYAD}, 0.0);
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      b_t(g, h) = theta_par[par_ind(h, g)];
    }
  }
  if(N_DYAD_PRED > 0){
    for(int z = 0; z < N_DYAD_PRED; ++z){
      gamma[z] = theta_par[N_B_PAR + z];
    }
  }
  
#pragma omp parallel for
  for(int d = 0; d < N_DYAD; ++d){
    double linpred = 0.0;
    if(N_DYAD_PRED > 0){
      for(int z = 0; z < N_DYAD_PRED; ++z){
        linpred -= z_t(z, d) * gamma[z];
      }
    }
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        theta(h, g, d) = 1./(1 + exp(linpred - b_t(h, g)));
      }
    }
  }
  
  
  /////////////////////////
  double res = 0.0;
#pragma omp parallel for collapse(3) reduction(+:res)
  for(int d = 0; d < N_DYAD; ++d){
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        res += send_phi(g, d) * rec_phi(h, d) * R::dbinom(y[d], 1, theta(h, g, d), 1);//;(y[d] * log(theta(h, g, d)) 
                                                   //+ (1.0 - y[d]) * log(1.0 - theta(h, g, d)));
      }
    }
  }
  //Rprintf("res now 2: %f\n", res);
  
  
  
  //Prior for gamma
  for(int z = 0; z < N_DYAD_PRED; ++z){
    res -= log(1.0 + pow(gamma[z], 2.0)/1.5625);
  }
  
  //Rprintf("res now 3: %f\n", res);
  
  
  //Prior for B
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      res -= log(1.0 + pow(b_t(h, g) - mu_b_t(h, g), 2.0)/100);
    }
  }
  //Rprintf("res now 4: %f\n", res);
  return -res/N_DYAD;
}



// [[Rcpp::export]]
NumericVector thetaGr(NumericVector theta_par,
                      const int N_DYAD_PRED,
                      const int N_B_PAR,
                      const int N_BLK,
                      const int N_DYAD,
                      const NumericMatrix mu_b_t,
                      const NumericMatrix z_t,
                      const NumericVector y,
                      const IntegerMatrix node_id_dyad,
                      const NumericMatrix phi_init,
                      const bool directed)
{
  int p, q;
  NumericMatrix b_t(N_BLK, N_BLK);
  NumericVector gamma(N_DYAD_PRED);
  Array<double> send_phi({N_BLK, N_DYAD}, 0.0);
  Array<double> rec_phi({N_BLK, N_DYAD}, 0.0);
  for(int d = 0; d < N_DYAD; ++d){
    p = node_id_dyad(d, 0);
    q = node_id_dyad(d, 1);
    for(int g = 0; g < N_BLK; ++g){
      send_phi(g, d) = phi_init(g, p);
      rec_phi(g, d) = phi_init(g, q);
    }
  }
  Array<int> par_ind({N_BLK, N_BLK}, 0);
  int ind = 0;
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      if(directed){
        par_ind(h, g) = ind++;
      } else {
        if(h >= g){
          par_ind(h, g) = ind++;
        } else {
          par_ind(h, g) = par_ind(g, h);
        }
      }
    }
  }
  Array<double> theta({N_BLK, N_BLK, N_DYAD}, 0.0);
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      b_t(g, h) = theta_par[par_ind(h, g)];
    }
  }
  if(N_DYAD_PRED > 0){
    for(int z = 0; z < N_DYAD_PRED; ++z){
      gamma[z] = theta_par[N_B_PAR + z];
    }
  }
  
#pragma omp parallel for
  for(int d = 0; d < N_DYAD; ++d){
    double linpred = 0.0;
    if(N_DYAD_PRED > 0){
      for(int z = 0; z < N_DYAD_PRED; ++z){
        linpred -= z_t(z, d) * gamma[z];
      }
    }
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        theta(h, g, d) = 1./(1 + exp(linpred - b_t(h, g)));
      }
    }
  }
  
  int N_PAR = N_DYAD_PRED + N_B_PAR;
  NumericVector gr(N_PAR);
  
  
  ////////////////////
  double res2, res = 0.0;
  for(int i = 0; i < N_PAR; ++i){
    gr[i] = 0.0;
  }
  
  
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      res = 0.0;
      if(!directed & (h < g)){
        continue;
      }
      for(int d = 0; d < N_DYAD; ++d){
        res2 = send_phi(g, d) * rec_phi(h, d) * (y[d] - theta(h, g, d));
        res += res2;
        for(int z = 0; z < N_DYAD_PRED; ++z){
          gr[N_B_PAR + z] -= res2 * z_t(z,d);
        }
      }
      for(int z = 0; z < N_DYAD_PRED; ++z){
        gr[N_B_PAR + z] += 2 * gamma[z] / (pow(gamma[z], 2.0) + 1.5625);
      }
      gr[par_ind(g, h)] = -res - 2 * (mu_b_t(g, h) - b_t(g, h))
        / (pow(mu_b_t(g, h), 2.0) + 100 - 2 * mu_b_t(g, h) + pow(b_t(g, h), 2.0));
    }
  }
  for(int i = 0; i < N_PAR; ++i){
    gr[i] /= N_DYAD;
  }
  return(gr);
}



/*** R
setwd("~/Dropbox/GitHub/NetMixRoot/NetMix/")
source("Extra/NetGenerator.R")
library(numDeriv)
library(plyr)
library(doParallel)
# cl <- makePSOCKcluster(22)
# registerDoParallel(cl)
set.seed(831213)
  
N_BLK <- 3
N <-100
N_DYAD <- N^2
N_DYAD_PRED <- 3
N_B_PAR <- N_BLK^2
N_sim <- 1
diffs <- array(NA, c(N_sim, N_DYAD_PRED + N_B_PAR))
true_beta <- array(c(0.05, rnorm(3, 0, 3),
                     0.05, rnorm(3, 0, 3),
                     0.05, rnorm(3, 0, 3)),
                   c(4, 3, 1))
true_gamma <- c(1.05, -1.05, 1.05)
Z <- list(cbind(runif(N^2), runif(N^2), runif(N^2)))  
X <- cbind(runif(N), runif(N), runif(N))
X_s <- list(cbind(1,scale(X)))
nodes_in_dyads <- cbind(rep(1:N, each = N), rep(1:N, times=N)) - 1
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
                                                ,gamma_vec = true_gamma
                                                ,X = X_s
                                                ,Z = Z))}
                                #,.parallel = TRUE
                                #,.paropts = list(.export=c("NetSim","true_beta","true_gamma","N", "X_s", "Z"))
))
#stopCluster(cl)
for(i in 1:N_sim){
  Z_t <-  t(Z[[1]])
  theta_rand <- rnorm(N_DYAD_PRED + N_B_PAR)
  (d1 <- grad(thetaLB, theta_rand,
              N_DYAD_PRED = net2_list[[i]]$DYAD_PRED,
              N_B_PAR = net2_list[[i]]$BLK^2,
              N_BLK = net2_list[[i]]$BLK,
              N_DYAD = net2_list[[i]]$NODE^2,
              mu_b_t = diag(10,N_BLK) - 5,  
              z_t = t(net2_list[[i]]$Z[[1]]),
              y = net2_list[[i]]$Y[[1]],
              node_id_dyad = nodes_in_dyads,
              phi_init = t(net2_list[[i]]$pi_vecs[[1]]),
              directed = net2_list[[i]]$DIRECTED))
  (d2 <- thetaGr(theta_rand,
                 N_DYAD_PRED = net2_list[[i]]$DYAD_PRED,
                 N_B_PAR = net2_list[[i]]$BLK^2,
                 N_BLK = net2_list[[i]]$BLK,
                 N_DYAD = net2_list[[i]]$NODE^2,
                 mu_b_t = diag(10,N_BLK) - 5,  
                 z_t = t(net2_list[[i]]$Z[[1]]),
                 y = net2_list[[i]]$Y[[1]],
                 node_id_dyad = nodes_in_dyads,
                 phi_init = t(net2_list[[i]]$pi_vecs[[1]]),
                 directed = net2_list[[i]]$DIRECTED))
  #deriv.diffs[i] <- mean((d1 - d2) / d1)
  # 
    bfgs_res <- optim(#c(t(net2_list[[i]]$B), net2_list[[i]]$gamma),
                      #c(glm.fit(net2_list[[i]]$Z[[1]], log(t((C_t+1)/(N*2))))$coef),
                      #rnorm( N_DYAD_PRED + N_B_PAR),
                      rep(0,  N_DYAD_PRED + N_B_PAR),
                      thetaLB,
                      gr = thetaGr,
                      N_DYAD_PRED = net2_list[[i]]$DYAD_PRED,
                      N_B_PAR = net2_list[[i]]$BLK^2,
                      N_BLK = net2_list[[i]]$BLK,
                      N_DYAD = net2_list[[i]]$NODE^2,
                      mu_b_t = diag(4.5,N_BLK) - 1.5,
                      z_t = t(net2_list[[i]]$Z[[1]]),
                      y = net2_list[[i]]$Y[[1]],
                      node_id_dyad = nodes_in_dyads,
                      phi_init = t(net2_list[[i]]$pi_vecs[[1]]),
                      directed = net2_list[[i]]$DIRECTED,
                      lower=-10,
                      upper=10,
                      control=list(maxit=10000, trace=6)
    ,method = "L-BFGS-B"
    )
    
     array(bfgs_res$par[1:N_BLK^2], c(N_BLK, N_BLK))
     net2_list[[i]]$B
     bfgs_res$par[-(1:N_BLK^2)]
     net2_list[[i]]$gamma_mat
  #   
}
  # pdf("~/Desktop/Sims.pdf")
  #   par(mfrow=c(3,4))
  #   for(i in 1:ncol(diffs)){
  #     hist(diffs[,i], main="")
  #   }
  #   dev.off()
      
#hist(deriv.diffs, breaks=500)
      */
