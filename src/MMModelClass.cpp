#include "MMModelClass.h"


using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::wrap;
using Rcpp::stop;
using Rcpp::sum;
using Rcpp::as;
using R::digamma;

/**
 CONSTRUCTOR
 */

MMModel::MMModel(const NumericMatrix& z_t,
                 const NumericMatrix& x_t,
                 const IntegerVector& y,
                 const int N_THREAD,
                 const IntegerVector& time_id_dyad,
                 const IntegerVector& time_id_node,
                 const IntegerVector& nodes_per_period,
                 const IntegerMatrix& node_id_dyad,
                 const NumericMatrix& mu_b,
                 const NumericMatrix& var_b,
                 const NumericMatrix& phi_init,
                 NumericMatrix& kappa_init_t,
                 NumericMatrix& b_init_t,
                 NumericVector& beta_init_r,
                 NumericVector& gamma_init_r,
                 List& control
)
  :
  N_NODE(sum(nodes_per_period)),
  N_DYAD(y.length()),
  N_BLK(as<int>(control["blocks"])),
  N_STATE(as<int>(control["states"])),
  N_TIME(as<int>(control["times"])),
  N_MONAD_PRED(x_t.nrow()),
  N_DYAD_PRED((z_t.nrow() == 1) && (sd(z_t) == 0.0) ? 0 : z_t.nrow()),
  N_B_PAR(as<int>(control["directed"]) ? N_BLK * N_BLK : N_BLK * (1 + N_BLK) / 2),
  OPT_ITER(as<int>(control["opt_iter"])),
  N_THREAD(N_THREAD),
  eta(as<double>(control["eta"])),
  var_gamma(as<double>(control["var_gamma"])),
  var_beta(as<double>(control["var_beta"])),
  fminAlpha(0.0),
  fminTheta(0.0),
  fncountAlpha(0),
  fncountTheta(0),
  grcountAlpha(0),
  grcountTheta(0),
  m_failAlpha(0),
  m_failTheta(0),
  verbose(as<bool>(control["verbose"])),
  directed(as<bool>(control["directed"])),
  y({N_DYAD}, y),
  time_id_dyad({N_DYAD}, time_id_dyad),
  time_id_node({N_NODE},time_id_node),
  n_nodes_time({N_TIME}, nodes_per_period),
  tot_nodes(N_NODE, 0.0),
  maskalpha(N_MONAD_PRED * N_BLK * N_STATE, 1),
  masktheta(N_B_PAR + N_DYAD_PRED, 1),
  theta_par({N_B_PAR + N_DYAD_PRED}, 0.0),
  e_wm({N_STATE}, 0.0),
  gamma({N_DYAD_PRED}, gamma_init_r),
  gamma_init({N_DYAD_PRED}, gamma_init_r),
  node_id_dyad({N_DYAD, 2}, node_id_dyad),
  par_ind({N_BLK, N_BLK}, 0),
  x_t({N_MONAD_PRED, N_NODE}, x_t),
  z_t({N_DYAD_PRED, N_DYAD}, z_t),
  mu_b_t({N_BLK, N_BLK},mu_b),
  var_b_t({N_BLK, N_BLK}, var_b),
  kappa_t({N_STATE, N_TIME}, kappa_init_t),
  b_t({N_BLK, N_BLK}, b_init_t),
  alpha_term({N_STATE, N_TIME}, 0.0),
  send_phi({N_BLK, N_DYAD}, 0.0),
  rec_phi({N_BLK, N_DYAD}, 0.0),
  e_wmn_t({N_STATE, N_STATE}, 0.0),
  e_c_t({N_BLK, N_NODE}, 0.0),
  alpha({N_BLK, N_NODE, N_STATE}, 0.0),
  theta({N_BLK, N_BLK, N_DYAD}, 0.0),
  beta({N_MONAD_PRED, N_BLK, N_STATE}, beta_init_r),
  beta_init({N_MONAD_PRED, N_BLK, N_STATE}, beta_init_r),
  new_e_c_t(N_THREAD, Array<double>({N_BLK, N_NODE}, 0.0))
{
  //Assign initial values to  W parameters
  for(int t = 1; t < N_TIME; ++t){
    for(int m = 0; m < N_STATE; ++m){
      e_wm[m] += kappa_t(m, t - 1);
      for(int n = 0; n < N_STATE; ++n){
        e_wmn_t(n, m) += kappa_t(m, t - 1) * kappa_t(n, t);
      }
    }
  }
  
  
  //Assign initial values to Phi and C
  int p, q;
  for(int d = 0; d < N_DYAD; ++d){
    p = node_id_dyad(d, 0);
    q = node_id_dyad(d, 1);
    tot_nodes[p]++;
    tot_nodes[q]++;
    for(int g = 0; g < N_BLK; ++g){
      send_phi(g, d) = phi_init(g, p);
      rec_phi(g, d) = phi_init(g, q);
      e_c_t(g, p) += send_phi(g, d);
      e_c_t(g, q) += rec_phi(g, d);
    }
  }
  
  //Create matrix of theta parameter indeces
  //(for undirected networks should force
  //symmetric blockmodel)
  int ind = 0;
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      if(directed){
        par_ind(h, g) = ind;
        ++ind;
      } else {
        if(h >= g){
          par_ind(h, g) = ind;
          ++ind;
        } else {
          par_ind(h, g) = par_ind(g, h);
        }
      }
    }
  }
  
  //Assign theta pars (which include B and gamma pars)
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      theta_par[par_ind(h, g)] = b_t(h, g);
    }
  }
  
  
  if(N_DYAD_PRED > 0)
    std::copy(gamma.begin(), gamma.end(), theta_par.begin() + N_B_PAR);
  
  
  //Assign initial values to alpha and theta
  computeAlpha();
  computeTheta();
  
  
}



/**
 ALPHA LOWER BOUND
 */

double MMModel::alphaLB()
{
  computeAlpha();
  double res = 0.0, res_int = 0.0, alpha_row = 0.0, alpha_val = 0.0;
#ifdef _OPENMP
#pragma omp parallel for firstprivate(alpha_row, alpha_val, res_int) reduction(+:res)
#endif
  for(int m = 0; m < N_STATE; ++m){
    for(int p = 0; p < N_NODE; ++p){
      alpha_row = 0.0;
      res_int = 0.0;
      for(int g = 0; g < N_BLK; ++g){
        alpha_val = alpha(g, p, m);
        alpha_row += alpha_val;
        res_int += (lgamma(alpha_val + e_c_t(g, p)) - lgamma(alpha_val));
        
      }
      res_int += (lgamma(alpha_row) - lgamma(alpha_row + tot_nodes[p]));

      res += res_int * kappa_t(m, time_id_node[p]);
    }
    
    //Prior for beta
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x){
        res -= 0.5 * pow(beta(x, g, m), 2.0) / var_beta;
      }
    }
  }
  
  return -res/N_NODE;
}

/**
 ALPHA GRADIENT
 */



void MMModel::alphaGr(int N_PAR, double *gr)
{
  //computeAlpha();
  double res=0.0, alpha_row=0.0, prior_gr=0.0;
  
  for(int m = 0; m < N_STATE; ++m){
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x){
        res = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for firstprivate(alpha_row) reduction(+:res)
  #endif
          for(int p = 0; p < N_NODE; ++p){
            alpha_row = 0.0;
            for(int h = 0; h < N_BLK; ++h){
              alpha_row += alpha(h, p, m);
            }
            res += (R::digamma(alpha_row) - R::digamma(alpha_row + tot_nodes[p])
                      + R::digamma(alpha(g, p, m) + e_c_t(g, p)) - R::digamma(alpha(g, p, m)))
              * kappa_t(m,  time_id_node[p]) * alpha(g, p, m) * x_t(x, p);
          }
          prior_gr = beta(x, g, m) / var_beta;
          gr[x + N_MONAD_PRED * (g + N_BLK * m)] = -(res - prior_gr);
        }
      }
    }
  for(int i = 0; i < N_PAR; ++i){
    gr[i] /= N_NODE;
  }



  
}



/**
 ALPHA COMPUTATION
 */

void MMModel::computeAlpha()
{
  std::fill(alpha_term.begin(), alpha_term.end(), 0.0);
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
        alpha_term(m, time_id_node[p]) += lgamma(linpred + e_c_t(g, p)) - lgamma(linpred);
      }
      alpha_term(m, time_id_node[p]) += lgamma(row_sum) - lgamma(row_sum + tot_nodes[p]);
      
    }
  }
}


/**
 THETA LB
 */

double MMModel::thetaLB(bool entropy = false)
{
  computeTheta();
  
  double res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif
  for(int d = 0; d < N_DYAD; ++d){
    for(int g = 0; g < N_BLK; ++g){
      if(entropy)
        res -= send_phi(g, d) * log(send_phi(g, d)) +
          rec_phi(g, d) * log(rec_phi(g, d));
      for(int h = 0; h < N_BLK; ++h){
        res += send_phi(g, d) * rec_phi(h, d)
        * (y[d] ? log(theta(h, g, d))
             : log(1.0 - theta(h, g, d)));
        
      }
    }
  }
  
  
  //Prior for gamma
  for(int z = 0; z < N_DYAD_PRED; ++z){
    res -= 0.5 * pow(gamma[z], 2.0) / var_gamma;
  }
  
  //Prior for B
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      res -= 0.5 * (pow(b_t(h, g) - mu_b_t(h, g), 2.0) / var_b_t(h, g));
    }
  }
  
  res *= -1; //VMMIN minimizes.
  return res / N_DYAD;
}


/**
 GRADIENT FOR THETA
 */
void MMModel::thetaGr(int N_PAR, double *gr)
{
  double res_local, res = 0.0;
  for(int i = 0; i < N_PAR; ++i)
    gr[i] = 0.0;
  
  int npar;
  for(int d = 0; d < N_DYAD; ++d){
    res = 0.0;
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        res_local = send_phi(g, d) * rec_phi(h, d) * (y[d] - theta(h, g, d));
        res += res_local;
        if((h < g) && !directed){
          continue;
        }
        npar = par_ind(h, g);
        gr[npar] -= res_local;
      }
    }
    if(N_DYAD_PRED > 0){
      for(int z = 0; z < N_DYAD_PRED; ++z){
        gr[N_B_PAR + z] -= res * z_t(z, d);
      }
    }
  }
  for(int z = 0; z < N_DYAD_PRED; ++z){
    gr[N_B_PAR + z] += gamma[z] / var_gamma;
  }
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      if((h < g) && !directed){
        continue;
      }
      npar = par_ind(h, g);
      gr[npar] += (b_t(h, g) - mu_b_t(h, g)) / var_b_t(h, g);
    }
  }
  for(int i = 0; i < N_PAR; ++i)
    gr[i] /= N_DYAD;
}


/**
 COMPUTE THETA
 */

void MMModel::computeTheta()
{
  for(int g = 0; g < N_BLK; ++g){
    for(int h = 0; h < N_BLK; ++h){
      b_t(h, g) = theta_par[par_ind(h, g)];
    }
  }
  double linpred;
  for(int d = 0; d < N_DYAD; ++d){
    linpred = 0.0;
    for(int z = 0; z < N_DYAD_PRED; ++z){
      gamma[z] = theta_par[N_B_PAR + z];
      linpred -= z_t(z, d) * gamma[z];
    }
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        theta(h, g, d) = 1./(1 + exp(linpred - b_t(h, g)));
      }
    }
  }
}

double MMModel::cLL()
{
  double res = lgamma(double(N_STATE) * eta)- lgamma(eta);
  res -= thetaLB(true);
  res -= alphaLB();
  for(int t = 0; t < N_TIME; ++t){
    for(int m = 0; m < N_STATE; ++m){
      res -= lgamma(double(N_STATE) * eta + e_wm[m]);
      for(int n = 0; n < N_STATE; ++n){
        res += log(eta + e_wmn_t(n, m));
      }
      //Entropy for kappa
      res -= kappa_t(m, t) * log(kappa_t(m,t) + 1e-8);
    }
  }
  return res;
}

/**
 BFGS OPTIMIZATION
 */

void MMModel::optim_ours(bool alpha)
{
  if(alpha){
    //std::copy(beta_init.begin(), beta_init.end(), beta.begin());
    vmmin_ours(N_MONAD_PRED * N_BLK * N_STATE, &beta[0], &fminAlpha, alphaLBW, alphaGrW, OPT_ITER, 0,
               &maskalpha[0], -1.0e+35, 1.0e-6, 1, this, &fncountAlpha, &grcountAlpha, &m_failAlpha);
  } else {
    //std::copy(gamma_init.begin(), gamma_init.end(), theta_par.begin() + N_B_PAR);
    vmmin_ours(N_B_PAR + N_DYAD_PRED, &theta_par[0], &fminTheta, thetaLBW, thetaGrW, OPT_ITER, 0,
               &masktheta[0], -1.0e+35, 1.0e-6, 1, this, &fncountTheta, &grcountTheta, &m_failTheta);
  }
}

/**
 WRAPPERS OF OPTIMIZATION FNs AND
 GRADIENTS
 */

double MMModel::thetaLBW(int n, double *par, void *ex)
{
  return(static_cast<MMModel*>(ex)->thetaLB());
}
void MMModel::thetaGrW(int n, double *par, double *gr, void *ex)
{
  static_cast<MMModel*>(ex)->thetaGr(n, gr);
}
double MMModel::alphaLBW(int n, double *par, void *ex)
{
  return(static_cast<MMModel*>(ex)->alphaLB());
}
void MMModel::alphaGrW(int n, double *par, double *gr, void *ex)
{
  static_cast<MMModel*>(ex)->alphaGr(n, gr);
}



/**
 VARIATIONAL UPDATE FOR KAPPA
 */


void MMModel::updateKappa()
{
  Rcpp::checkUserInterrupt();
  std::vector<double> kappa_vec(N_STATE);
  double res, log_denom;
  for(int t = 1; t < N_TIME; ++t){
    for(int m = 0; m < N_STATE; ++m){
      res = 0.0;
      if(t < (N_TIME - 1)){
        e_wm[m] -= kappa_t(m, t);
      }
      res -= log(double(N_STATE) * eta + e_wm[m]);
      
      if(t < (N_TIME - 1)){
        e_wmn_t(m, m) -= kappa_t(m, t) * (kappa_t(m, t + 1) + kappa_t(m, t - 1));
        res += kappa_t(m, t + 1) * kappa_t(m, t - 1) * log(eta + e_wmn_t(m, m) + 1);
        
        res += (kappa_t(m, t - 1) * (1 - kappa_t(m, t + 1))
                  + kappa_t(m, t + 1)) * log(eta + e_wmn_t(m, m));
        
        for(int n = 0; n < N_STATE; ++n){
          if(m != n){
            e_wmn_t(n, m) -= kappa_t(m, t) * kappa_t(n, t + 1);
            res += kappa_t(n, t + 1) * log(eta + e_wmn_t(n, m));
            
            e_wmn_t(m, n) -= kappa_t(m, t) * kappa_t(n, t - 1);
            res += kappa_t(n, t - 1) * log(eta + e_wmn_t(m, n));
          }
        }
      } else { // t = T
        for(int n = 0; n < N_STATE; ++n){
          e_wmn_t(m, n) -= kappa_t(m, t) * kappa_t(n, t - 1);
          res += kappa_t(n, t - 1) * log(eta + e_wmn_t(m, n));
        }
      }
      res += alpha_term(m, t);
      kappa_vec[m] = res;
    }
    log_denom = logSumExp(kappa_vec);
    for(int m = 0; m < N_STATE; ++m){
      kappa_t(m, t) = exp(kappa_vec[m] - log_denom);
      if(!std::isfinite(kappa_t(m, t))){
        stop("Kappa value became NAN.");
      }
      if(t < (N_TIME - 1)){
        e_wm[m] += kappa_t(m, t);
        e_wmn_t(m, m) += kappa_t(m, t) * (kappa_t(m, t + 1) + kappa_t(m, t - 1));
        for(int n = 0; n < N_STATE; ++n){
          if(m != n){
            e_wmn_t(n, m) += kappa_t(m, t) * kappa_t(n, t + 1);
            e_wmn_t(m, n) += kappa_t(m, t) * kappa_t(n, t - 1);
          }
        }
      } else { // t = T
        for(int n = 0; n < N_STATE; ++n){
          e_wmn_t(m, n) += kappa_t(m, t) * kappa_t(n, t - 1);
        }
      }
    }
  }
}


/**
 VARIATIONAL UPDATE FOR PHI
 */

void MMModel::updatePhiInternal(int dyad, int rec,
                                double *phi,
                                double *phi_o,
                                double *new_c,
                                int *err
)
{
  
  int t = time_id_dyad[dyad];
  double edge = y[dyad];
  int incr1 = rec ? 1 : N_BLK;
  int incr2 = rec ? N_BLK : 1;
  int node = node_id_dyad(dyad, rec);
  double *theta_temp = &theta(0, 0, dyad);
  double *te;
  
  double total = 0.0, res;
  for(int g = 0; g < N_BLK; ++g, theta_temp+=incr1){
    res = 0.0;
    for(int m = 0; m < N_STATE; ++m){
      res += kappa_t(m, t) * log(alpha(g, node, m) + e_c_t(g, node) - phi[g]);
    }
    
    te = theta_temp;
    for(int h = 0; h < N_BLK; ++h, te+=incr2){
      res += phi_o[h] * (edge * log(*te) + (1.0 - edge) * log(1.0 - *te));
    }
    phi[g] = exp(res);
    if(ISNAN(phi[g])){
      *err = 1;
    }
    total += phi[g];
  }
  
  //Normalize phi to sum to 1
  //and store new value in c
  for(int g = 0; g < N_BLK; ++g){
    phi[g] /= total;
    new_c[g] += phi[g];      
  }
}


void MMModel::updatePhi()
{
  for(int thread = 0; thread < N_THREAD; ++thread){
    std::fill(new_e_c_t[thread].begin(), new_e_c_t[thread].end(), 0.0);
  }
  int err = 0;
#ifdef _OPENMP  
#pragma omp parallel for
#endif
  for(int d = 0; d < N_DYAD; ++d){
    int thread = 0;
#ifdef _OPENMP
    thread = omp_get_thread_num();
#endif
    updatePhiInternal(d,
                      0,
                      &(send_phi(0, d)),
                      &(rec_phi(0, d)),
                      &(new_e_c_t[thread](0, node_id_dyad(d, 0))),
                      &err
    );
    updatePhiInternal(d,
                      1,
                      &(rec_phi(0, d)),
                      &(send_phi(0, d)),
                      &(new_e_c_t[thread](0, node_id_dyad(d, 1))),
                      &err
    );
  }
  
  if(err){
    stop("Phi value became NaN");
  }
  
  
  std::fill(e_c_t.begin(), e_c_t.end(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for(int p = 0; p < N_NODE; ++p){
    for(int g = 0; g < N_BLK; ++g){
      for(int i = 0; i < N_THREAD; ++i){
        e_c_t(g, p) += (new_e_c_t[i])(g, p);
      }
    }
  }
}

/**
 * Convergence checkers
 */
bool MMModel::maxDiffCheck(Array<double>& current,
                        Rcpp::NumericVector& old,
                        double tol)
{
  Rcpp::NumericVector diff(old.size());
  std::transform(current.begin(), current.end(),
                 old.begin(), diff.begin(),
                 [&](double x, double y){return fabs(x - y);});
  double max_val = *std::max_element(diff.begin(), diff.end());
  return max_val > tol ? false : true;
}

bool MMModel::checkConv(char p,
                        Rcpp::NumericVector& old,
                        double tol)
{
  bool conv;
  switch(p) {
  case 'b':
    conv = maxDiffCheck(b_t, old, tol);
    break;
  case 'g':
    conv = maxDiffCheck(gamma, old, tol);
    break;
  case 'e':
    conv = maxDiffCheck(beta, old, tol);
    break;
  case 'c':
    conv = maxDiffCheck(e_c_t, old, tol);
    break;
  }
  return conv;
}

/**
 GETTER FUNCTIONS
 */

NumericMatrix MMModel::getB()
{
  NumericMatrix res(N_BLK, N_BLK);
  std::copy(b_t.begin(), b_t.end(), res.begin());
  return res;
}

void MMModel::getB(NumericVector& res)
{
  std::copy(b_t.begin(), b_t.end(), res.begin());
}


List MMModel::getBeta()
{
  std::vector< NumericMatrix > res(N_STATE);
  for(int m = 0; m < N_STATE; ++m){
    NumericMatrix mat_res(N_MONAD_PRED, N_BLK);
    for(int g = 0; g < N_BLK; ++g){
      for(int x = 0; x < N_MONAD_PRED; ++x){
        mat_res(x, g) = beta(x, g, m);
      }
    }
    res[m] = mat_res;
  }
  return wrap(res);
}

void MMModel::getBeta(NumericVector& res)
{
  std::copy(beta.begin(), beta.end(), res.begin());
}


NumericMatrix MMModel::getC()
{
  NumericMatrix res(N_BLK, N_NODE);
  double row_total;
  for(int p = 0; p < N_NODE; ++p){
    row_total = 0.0;
    for(int g = 0; g < N_BLK; ++g){
      row_total += e_c_t(g, p);
    }
    for(int g = 0; g < N_BLK; ++g){
      res(g, p) = e_c_t(g, p) / row_total;
    }
  }
  return res;
}

void MMModel::getC(NumericVector& res)
{
  std::copy(e_c_t.begin(), e_c_t.end(), res.begin());
}

NumericVector MMModel::getGamma()
{
  NumericVector res(N_DYAD_PRED);
  std::copy(gamma.begin(), gamma.end(), res.begin());
  return res;
}

void  MMModel::getGamma(NumericVector& res)
{
  std::copy(gamma.begin(), gamma.end(), res.begin());
}

NumericMatrix MMModel::getKappa()
{
  NumericMatrix res(N_STATE, N_TIME);
  std::copy(kappa_t.begin(), kappa_t.end(), res.begin());
  return res;
}


NumericMatrix MMModel::getWmn()
{
  NumericMatrix res(N_STATE, N_STATE);
  if((N_TIME > 1) & (N_STATE > 1)){
    double row_total;
    for(int c = 0; c < N_STATE; ++c){
      row_total = 0.0;
      for(int r = 0; r < N_STATE; ++r){
        row_total += e_wmn_t(r, c);
      }
      for(int r = 0; r < N_STATE; ++r){
        res(r, c) = e_wmn_t(r, c) / row_total;
      }
    }
  }
  return res;
}
