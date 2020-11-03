#include "MMModelClass.h"


/**
 CONSTRUCTOR
 */

MMModel::MMModel(const arma::mat& z_t,
                 //const arma::mat& z_t_ho,
                 const arma::mat& x_t,
                 const arma::vec& y,
                 //const arma::vec& y_ho,
                 const arma::uvec& time_id_dyad,
                 const arma::uvec& time_id_node,
                 const arma::uvec& nodes_per_period,
                 const arma::umat& node_id_dyad,
                 //const arma::umat& node_id_dyad_ho,
                 const arma::field<arma::uvec>& node_id_period,
                 const arma::mat& mu_b,
                 const arma::mat& var_b,
                 const arma::cube& mu_beta,
                 const arma::cube& var_beta,
                 const arma::vec& mu_gamma,
                 const arma::vec& var_gamma,
                 const arma::mat& pi_init,
                 arma::mat& kappa_init_t,
                 arma::mat& b_init_t,
                 arma::cube& beta_init_r,
                 arma::vec& gamma_init_r,
                 //double sparsity,
                 Rcpp::List& control)
  :
  N_NODE(sum(nodes_per_period)),
  N_DYAD(y.n_elem),
  N_BLK(control["blocks"]),
  N_STATE(control["states"]),
  N_TIME(control["times"]),
  N_MONAD_PRED(x_t.n_rows),
  N_DYAD_PRED(arma::any(z_t.row(0)) ? z_t.n_rows : 0),
  N_B_PAR(Rcpp::as<bool>(control["directed"]) ? N_BLK * N_BLK : N_BLK * (1 + N_BLK) / 2),
  OPT_ITER(control["opt_iter"]),
  N_NODE_BATCH(arma::sum(Rcpp::as<arma::uvec>(control["batch_size"]))),
  N_THREAD(control["threads"]),
  //N_DYAD_HO(y_ho.n_elem),
  eta(Rcpp::as<double>(control["eta"])),
  forget_rate(Rcpp::as<double>(control["forget_rate"])),
  delay(Rcpp::as<double>(control["delay"])),
  //sparsity(sparsity),
  var_gamma(var_gamma),
  mu_gamma(mu_gamma),
  var_beta(var_beta),
  mu_beta(mu_beta),
  fminAlpha(0.0),
  fminTheta(0.0),
  reweightFactor(1.0),
  step_size(1.0),
  fncountAlpha(0),
  fncountTheta(0),
  grcountAlpha(0),
  grcountTheta(0),
  m_failAlpha(0),
  m_failTheta(0),
  verbose(Rcpp::as<bool>(control["verbose"])),
  directed(Rcpp::as<bool>(control["directed"])),
  y(y),
  //y_ho(y_ho),
  time_id_dyad(time_id_dyad),
  time_id_node(time_id_node),
  n_nodes_time(nodes_per_period),
  n_nodes_batch(Rcpp::as<arma::uvec>(control["batch_size"])),
  node_est(Rcpp::as<arma::uvec>(control["node_est"])),
  tot_nodes(N_NODE, arma::fill::zeros),
  node_in_batch(N_NODE, arma::fill::ones),
  dyad_in_batch(N_DYAD, arma::fill::ones),
  node_batch(N_NODE_BATCH, arma::fill::zeros),
  maskalpha(N_MONAD_PRED * N_BLK * N_STATE, 1),
  masktheta(N_B_PAR + N_DYAD_PRED, 1),
  node_id_period(node_id_period),
  theta_par(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  thetaold(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  e_wm(N_STATE, arma::fill::zeros),
  alpha_gr(N_MONAD_PRED * N_BLK * N_STATE, arma::fill::zeros),
  theta_gr(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  gamma(gamma_init_r),
  gamma_init(gamma_init_r),
  node_id_dyad(node_id_dyad),
  //node_id_dyad_ho(node_id_dyad_ho),
  par_ind(N_BLK, N_BLK, arma::fill::zeros),
  x_t(x_t),
  z_t(z_t),
  //z_t_ho(z_t_ho),
  mu_b_t(mu_b),
  var_b_t(var_b),
  kappa_t(kappa_init_t),
  b_t(b_init_t),
  alpha_term(N_STATE, N_TIME, arma::fill::zeros),
  send_phi(N_BLK, N_DYAD, arma::fill::zeros),
  rec_phi(N_BLK, N_DYAD, arma::fill::zeros),
  e_wmn_t(N_STATE, N_STATE, arma::fill::zeros),
  e_c_t(N_BLK, N_NODE, arma::fill::zeros),
  alpha(N_BLK, N_NODE, N_STATE, arma::fill::zeros),
  theta(N_BLK, N_BLK, N_DYAD, arma::fill::zeros),
  beta(beta_init_r),
  betaold(beta_init_r),
  beta_init(beta_init_r)
  //new_e_c_t(N_THREAD, Array<double>({N_BLK, N_NODE}, 0.0)
  //new_e_c_t(N_BLK, N_NODE, arma::fill::zeros)
{
  //Set number of parallel threads
#ifdef _OPENMP
  omp_set_num_threads(N_THREAD);
#endif
  

  //Assign initial values to W parameters
  for(arma::uword t = 1; t < N_TIME; ++t){
    for(arma::uword m = 0; m < N_STATE; ++m){
      e_wm[m] += kappa_t(m, t - 1);
      for(arma::uword n = 0; n < N_STATE; ++n){
        e_wmn_t(n, m) += kappa_t(m, t - 1) * kappa_t(n, t);
      }
    }
  }
  
  //Assign initial values to Phi and C
  arma::uword p, q;
  for(arma::uword d = 0; d < N_DYAD; ++d){
    p = node_id_dyad(d, 0);
    q = node_id_dyad(d, 1);
    tot_nodes[p]++;
    tot_nodes[q]++;
    for(arma::uword g = 0; g < N_BLK; ++g){
      send_phi(g, d) = pi_init(g, p);
      rec_phi(g, d) = pi_init(g, q);
      e_c_t(g, p) += send_phi(g, d);
      e_c_t(g, q) += rec_phi(g, d);
    }
  }
  //Create matrix of theta parameter indeces
  //(for undirected networks should force
  //symmetric blockmodel)
  arma::uword ind = 0;
  for(arma::uword g = 0; g < N_BLK; ++g){
    for(arma::uword h = 0; h < N_BLK; ++h){
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
  for(arma::uword g = 0; g < N_BLK; ++g){
    for(arma::uword h = 0; h < N_BLK; ++h){
      theta_par[par_ind(h, g)] = b_t(h, g);
    }
  }
  
  if(N_DYAD_PRED > 0)
    std::copy(gamma.begin(), gamma.end(), theta_par.begin() + N_B_PAR);
  
  //Assign initial values to alpha and theta
  computeAlpha();
  computeTheta();
  
  
  //Define iterator at the end of param. objects
  beta_end = beta.end();
  theta_par_end = theta_par.end();
  
  
}
MMModel::~MMModel()
{
  
}


/**
 ALPHA LOWER BOUND
 */

double MMModel::alphaLB(bool all)
{
  computeAlpha(all);
  double res = 0.0, res_int = 0.0, alpha_row = 0.0, alpha_val = 0.0;
  for(arma::uword m = 0; m < N_STATE; ++m){
#pragma omp parallel for firstprivate(alpha_val, alpha_row, res_int) reduction(+: res)
    for(arma::uword p = 0; p < N_NODE; ++p){
      if((node_in_batch[p] == 1) || all){
      alpha_row = 0.0;
      res_int = 0.0;
      for(arma::uword g = 0; g < N_BLK; ++g){
        alpha_val = alpha(g, p, m);
        alpha_row += alpha_val;
        res_int += (lgamma(alpha_val + e_c_t(g, p)) - lgamma(alpha_val));
      }
      res_int += (lgamma(alpha_row) - lgamma(alpha_row + tot_nodes[p]));
      res += res_int * kappa_t(m, time_id_node[p]);
      }
    }
    res *=  all ? 1.0 : (1. * N_NODE)/N_NODE_BATCH;
    //Prior for beta

    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
        res -= pow(beta(x, g, m) - mu_beta(x, g, m), 2.0) / var_beta(x, g, m);
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
  double res=0.0, alpha_row=0.0, prior_gr=0.0;
  arma::uword U_NPAR = N_PAR;
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
        res = 0.0;
#pragma omp parallel for firstprivate(alpha_row) reduction(+: res)
        for(arma::uword p = 0; p < N_NODE; ++p){
          if(node_in_batch[p] == 1) {
            alpha_row = 0.0;
            for(arma::uword h = 0; h < N_BLK; ++h){
              alpha_row += alpha(h, p, m);
            }
            res += (R::digamma(alpha_row) - R::digamma(alpha_row + tot_nodes[p])
                      + R::digamma(alpha(g, p, m) + e_c_t(g, p)) - R::digamma(alpha(g, p, m)))
              * kappa_t(m,  time_id_node[p]) * alpha(g, p, m) * x_t(x, p);
          }
        }
        res *= (1. * N_NODE) / N_NODE_BATCH;
        prior_gr = 2.0*(beta(x, g, m) - mu_beta(x, g, m)) / var_beta(x, g, m);
        gr[x + N_MONAD_PRED * (g + N_BLK * m)] = -(res - prior_gr);
      }
    }
  }
  for(arma::uword i = 0; i < U_NPAR; ++i)
  gr[i] /= N_NODE;
  
}



/**
 ALPHA COMPUTATION
 */

void MMModel::computeAlpha(bool all)
{
  std::fill(alpha_term.begin(), alpha_term.end(), 0.0); 
  double linpred, row_sum, correct_fact;
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword p = 0; p < N_NODE; ++p){
      if((node_in_batch[p] == 1) || all){
        correct_fact = all ? 1.0: (1. * N_NODE)/N_NODE_BATCH;//((1. * tot_nodes[p]) /  n_nodes_batch[time_id_node[p]]);
        //Rprintf("Correct fact = %f\n", correct_fact);
      row_sum = 0.0;
      for(arma::uword g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
          linpred += x_t(x, p) * beta(x, g, m);
        }
        linpred = exp(linpred);
        row_sum += linpred;
        alpha(g, p, m) = linpred;
        alpha_term(m, time_id_node[p]) +=  correct_fact *  (lgamma(linpred + e_c_t(g, p)) - lgamma(linpred));
      }
      alpha_term(m, time_id_node[p]) += correct_fact *  (lgamma(row_sum) - lgamma(row_sum + tot_nodes[p]));
    }
  }
  }
}

/**
 THETA LB
 */

double MMModel::thetaLB(bool entropy, bool all)
{
  computeTheta();

  double res = 0.0;
#pragma omp parallel for reduction(+: res)
  for(arma::uword d = 0; d < N_DYAD; ++d){
    if((dyad_in_batch[d] == 1) || all){
    for(arma::uword g = 0; g < N_BLK; ++g){
      if(entropy){
        res -= send_phi(g, d) * log(send_phi(g, d))
        + rec_phi(g, d) * log(rec_phi(g, d));
      }
      for(arma::uword h = 0; h < N_BLK; ++h){
        res += send_phi(g, d) * rec_phi(h, d)
        * (y[d] * log(theta(h, g, d))
             + (1.0 - y[d]) * log(1.0 - theta(h, g, d)));
      }
    }
    }
  }
  res *= all ? 1.0 : reweightFactor;

  //Prior for gamma
  for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
    res -= pow(gamma[z] - mu_gamma[z], 2.0) / var_gamma[z];
  }

  //Prior for B
  for(arma::uword g = 0; g < N_BLK; ++g){
    for(arma::uword h = 0; h < N_BLK; ++h){
      res -= (pow(b_t(h, g) - mu_b_t(h, g), 2.0) / var_b_t(h, g));
    }
  }

  return -res/N_DYAD;
}


/**
 GRADIENT FOR THETA
 */
void MMModel::thetaGr(int N_PAR, double *gr)
{

  arma::uword U_NPAR = N_PAR;
  double res_local, res = 0.0;
  for(arma::uword i = 0; i < U_NPAR; ++i){
    gr[i] = 0.0;
  }
  
  arma::uword npar = 0;
    for(arma::uword d = 0; d < N_DYAD; ++d){
      if(dyad_in_batch[d] == 1){
      res = 0.0;
      for(arma::uword g = 0; g < N_BLK; ++g){
        for(arma::uword h = 0; h < N_BLK; ++h){
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
        for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
          gr[N_B_PAR + z] -= res * z_t(z, d);
        }
      }
    }
  }
  for(arma::uword i = 0; i < U_NPAR; ++i){
    gr[i] *= reweightFactor; //for stochastic VI
  }
  for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
    gr[N_B_PAR + z] += 2.0*(gamma[z] - mu_gamma[z]) / var_gamma[z];
  }
  for(arma::uword g = 0; g < N_BLK; ++g){
    for(arma::uword h = 0; h < N_BLK; ++h){
      if((h < g) && !directed){
        continue;
      }
      npar = par_ind(h, g);
      gr[npar] += 2.0*(b_t(h, g) - mu_b_t(h, g)) / var_b_t(h, g);
    }
  }
  for(arma::uword i = 0; i < U_NPAR; ++i)
    gr[i] /= N_DYAD;
}


/**
 COMPUTE THETA
 */

void MMModel::computeTheta(bool all)
{
  for(arma::uword g = 0; g < N_BLK; ++g){
    for(arma::uword h = 0; h < N_BLK; ++h){
      b_t(h, g) = theta_par[par_ind(h, g)];
    }
  }
  double linpred;
  for(arma::uword d = 0; d < N_DYAD; ++d){
    if((dyad_in_batch[d] == 1) || all){
    linpred = 0.0;
    for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
      gamma[z] = theta_par[N_B_PAR + z];
      linpred -= z_t(z, d) * gamma[z];
    }
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword h = 0; h < N_BLK; ++h){
        theta(h, g, d) = 1./(1 + exp(linpred - b_t(h, g)));
      }
    }
  }
  }
}

/**
 OPTIMIZATION
 */

void MMModel::optim_ours(bool alpha)
{
  if(alpha){
    // alphaGr(N_MONAD_PRED * N_BLK * N_STATE,
    //         &alpha_gr[0], true);
    // for(std::pair<arma::cube::iterator,
    //     arma::vec::iterator> iter(beta.begin(),
    //                               alpha_gr.begin());
    //     iter.first != beta_end; 
    //     ++iter.first, ++iter.second)
    // {
    //   *iter.first -= step_size * *iter.second;  
    // }
    // computeAlpha(true);
    int npar = N_MONAD_PRED * N_BLK * N_STATE;
    betaold = beta;
    std::copy(beta_init.begin(), beta_init.end(), beta.begin());
    //beta.zeros();
    vmmin_ours(npar, &beta[0], &fminAlpha, alphaLBW, alphaGrW, OPT_ITER, 0,
                &maskalpha[0], -1.0e+35, 1.0e-6, 1, this, &fncountAlpha, &grcountAlpha, &m_failAlpha);
    
    for(arma::uword i = 0; i < npar; ++i){
      beta[i] = (1 - step_size) * betaold[i] + step_size * beta[i];
    }
    
    //Rcpp::Rcout << beta << std::endl;
  } else {
    // thetaGr(N_B_PAR + N_DYAD_PRED,
    //         &theta_gr[0], true);
    // for(arma::vec::iterator theta_val = theta_par.begin(),
    //     gr_val = theta_gr.begin();
    //     theta_val != theta_par_end; ++theta_val, ++gr_val)
    // {
    //   *theta_val -= step_size * *gr_val;  
    // }
    // computeTheta(true);
    int npar = N_B_PAR + N_DYAD_PRED;
    thetaold = theta_par;
    //theta_par.zeros();
    std::copy(gamma_init.begin(), gamma_init.end(), theta_par.begin() + N_B_PAR);
    vmmin_ours(npar, &theta_par[0], &fminTheta, thetaLBW, thetaGrW, OPT_ITER, 0,
                &masktheta[0], -1.0e+35, 1.0e-6, 1, this, &fncountTheta, &grcountTheta, &m_failTheta);
    
    for(arma::uword i = 0; i < npar; ++i){
      theta_par[i] = (1.0 - step_size) * thetaold[i] + step_size * theta_par[i];
    }
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

double MMModel::LB()
{
  double res = lgamma(double(N_STATE) * eta) - lgamma(eta);
  res -= thetaLB(true, true);
  res -= alphaLB(true);
  for(arma::uword t = 0; t < N_TIME; ++t){
    for(arma::uword m = 0; m < N_STATE; ++m){
      res -= lgamma(double(N_STATE) * eta + e_wm[m]);
      for(arma::uword n = 0; n < N_STATE; ++n){
        res += log(eta + e_wmn_t(n, m));
      }
      //Entropy for kappa
      res -= kappa_t(m, t) * log(kappa_t(m,t) + 1e-8);
    }
  }
  return -res;
}

double MMModel::LL()
{
  double ll = 0.0;
#pragma omp parallel for reduction(-:ll)
  for(arma::uword d = 0; d < N_DYAD; ++d){
    if(dyad_in_batch[d]==0){
    arma::uword p, q;
    p = node_id_dyad(d, 0);
    q = node_id_dyad(d, 1);
    ll -= log(arma::as_scalar(1.0 + exp(-(getPostMM(p).t() * b_t * getPostMM(q)
                                        + gamma.t() * z_t.col(d)))));
    }
    }
  return(ll);
}



/**
 VARIATIONAL UPDATE FOR KAPPA
 */


void MMModel::updateKappa()
{
  Rcpp::checkUserInterrupt();
  arma::vec kappa_vec(N_STATE, arma::fill::zeros);
  double res, log_denom;
  for(arma::uword t = 1; t < N_TIME; ++t){
    for(arma::uword m = 0; m < N_STATE; ++m){
      res = 0.0;
      if(t < (N_TIME - 1)){
        e_wm[m] -= kappa_t(m, t);
      }
      res -= log(double(N_STATE) * eta + std::max(e_wm[m], 0.0));
      
      if(t < (N_TIME - 1)){
        e_wmn_t(m, m) -= kappa_t(m, t) * (kappa_t(m, t + 1) + kappa_t(m, t - 1));
        res += kappa_t(m, t + 1) * kappa_t(m, t - 1) * log(eta + std::max(e_wmn_t(m, m),0.0) + 1);
        
        res += (kappa_t(m, t - 1) * (1 - kappa_t(m, t + 1))
                  + kappa_t(m, t + 1)) * log(eta + std::max(e_wmn_t(m, m), 0.0));
        
        for(arma::uword n = 0; n < N_STATE; ++n){
          if(m != n){
            e_wmn_t(n, m) -= kappa_t(m, t) * kappa_t(n, t + 1);
            res += kappa_t(n, t + 1) * log(eta + std::max(e_wmn_t(n, m),0.0));
            
            e_wmn_t(m, n) -= kappa_t(m, t) * kappa_t(n, t - 1);
            res += kappa_t(n, t - 1) * log(eta + std::max(e_wmn_t(m, n),0.0));
          }
        }
      } else { // t = T
        for(arma::uword n = 0; n < N_STATE; ++n){
          e_wmn_t(m, n) -= kappa_t(m, t) * kappa_t(n, t - 1);
          res += kappa_t(n, t - 1) * log(eta + std::max(e_wmn_t(m, n), 0.0));
        }
      }
      //Rprintf("res before: %f, alphaterm(%i, %i): %f\n", res, m, t, alpha_term(m,t));
      res += alpha_term(m, t);///n_nodes_time[t];
      //Rprintf("\tres after: %f\n", res);
      kappa_vec[m] = res;
    }
    log_denom = logSumExp(kappa_vec);
    for(arma::uword m = 0; m < N_STATE; ++m){
      kappa_t(m, t) = exp(kappa_vec[m] - log_denom);
      if(!std::isfinite(kappa_t(m, t))){
        Rcpp::stop("Kappa value became NaN.");
      }
      if(t < (N_TIME - 1)){
        e_wm[m] += kappa_t(m, t);
        e_wmn_t(m, m) += kappa_t(m, t) * (kappa_t(m, t + 1) + kappa_t(m, t - 1));
        for(arma::uword n = 0; n < N_STATE; ++n){
          if(m != n){
            e_wmn_t(n, m) += kappa_t(m, t) * kappa_t(n, t + 1);
            e_wmn_t(m, n) += kappa_t(m, t) * kappa_t(n, t - 1);
          }
        }
      } else { // t = T
        for(arma::uword n = 0; n < N_STATE; ++n){
          e_wmn_t(m, n) += kappa_t(m, t) * kappa_t(n, t - 1);
        }
      }
    }
  }
}


/**
 VARIATIONAL UPDATE FOR PHI
 */

void MMModel::updatePhiInternal(arma::uword dyad,
                                arma::uword rec,
                                double *phi,
                                double *phi_o,
                                double *new_c,
                                arma::uword *err
)
{
  
  arma::uword t = time_id_dyad[dyad];
  double edge = y[dyad];
  arma::uword incr1 = rec ? 1 : N_BLK;
  arma::uword incr2 = rec ? N_BLK : 1;
  arma::uword node = node_id_dyad(dyad, rec);
  double *theta_temp = &theta(0, 0, dyad);
  double *te;

  
  double total = 0.0, res;
  for(arma::uword g = 0; g < N_BLK; ++g, theta_temp+=incr1){
    new_c[g] -= phi[g];
    res = 0.0;
    for(arma::uword m = 0; m < N_STATE; ++m){
      // Rprintf("kappa = %f, alpha = %f, new_c = %f, phi = %f\n", kappa_t(m, t),
      //         alpha(g, node, m), new_c[g], phi[g]);
      res += kappa_t(m, t) * log(alpha(g, node, m) + std::max(new_c[g], 0.0));
    }
    //Rprintf("res before = %f\n",res);
    te = theta_temp;
    for(arma::uword h = 0; h < N_BLK; ++h, te+=incr2){
      //Rprintf("theta = %f, phi_o = %f\n",*te, phi_o[h]);
      res += phi_o[h] * (edge * log(*te) + (1.0 - edge) * log(1.0 - *te));
    }
    //Rprintf("res after = %f\n",res);
    phi[g] = exp(res);
    if(!std::isfinite(phi[g])){
      (*err)++;
    }
    total += phi[g];
  }
  
  //Normalize phi to sum to 1
  //and store new value in c
  for(arma::uword g = 0; g < N_BLK; ++g){
    phi[g] /= total;
    new_c[g] += phi[g];
  }
}


void MMModel::updatePhi()
{
  //Rcpp::Rcout << arma::mean(arma::sum(old_e)) << std::endl;
  arma::uword err = 0;
  // Update dyads with sampled nodes
  for(arma::uword d = 0; d < N_DYAD; ++d){
    Rcpp::checkUserInterrupt();
      if(node_est[node_id_dyad(d, 0)]) {
        updatePhiInternal(d,
                          0,
                          &(send_phi(0, d)),
                          &(rec_phi(0, d)),
                          &(e_c_t(0, node_id_dyad(d, 0))),
                          &err
        );
      }
      if(node_est[node_id_dyad(d, 1)]) {
        updatePhiInternal(d,
                          1,
                          &(rec_phi(0, d)),
                          &(send_phi(0, d)),
                          &(e_c_t(0, node_id_dyad(d, 1))),
                          &err
        );
      } 
  }
  //Rcpp::stop("stop here.");
  

  if(err){
    Rcpp::stop("Phi value became NaN.");
  }
  
}

void MMModel::sampleDyads(arma::uword iter)
{
  // Sample nodes for stochastic variational update
  if(N_TIME < 2){
    node_batch = arma::randperm(N_NODE, n_nodes_batch[0]);
  } else{ // N_TIME  >= 2
    arma::uvec node_batch_t;
    arma::uword ind = 0;
    for(arma::uword t = 0; t < N_TIME; ++t){
      node_batch_t = arma::randperm(n_nodes_time[t], n_nodes_batch[t]);
      for(arma::uword p = 0; p < n_nodes_batch[t]; ++p){
        node_batch[ind] = node_id_period[t][node_batch_t[p]];
        ++ind;
      }
    }
  }
  
  for(arma::uword p = 0; p < N_NODE; ++p){
    node_in_batch[p] = arma::any(node_batch == p) ? 1 : 0;
  }
  
  for(arma::uword d = 0; d < N_DYAD; ++d){
    dyad_in_batch[d] = (arma::any(node_batch == node_id_dyad(d, 0)) 
                          | arma::any(node_batch == node_id_dyad(d, 1))) ? 1 : 0;
  }
  
  reweightFactor = (1. * N_DYAD) / arma::sum(dyad_in_batch);
  step_size = 1.0 / pow(delay + iter, forget_rate);
}

/**
 * CONVERGENCE CHECKER
 */
void MMModel::convCheck(bool& conv,
                        const arma::vec& ll,
                        const double& tol)
{
  //double cor_val = arma::as_scalar(arma::cor(ll, arma::regspace(1, ll.n_elem)));
  //Rprintf("Cor is %f, std is %f\n", cor_val, arma::stddev(ll));
  //Rcpp::Rcout << ll << std::endl;
  conv = (fabs(arma::as_scalar((ll.tail(1)-ll.head(1))/ll.head(1))) < tol);
  // if((arma::stddev(ll) == 0.0) | (fabs(cor_val) <= tol)){
  //   conv = true;
  // }
}

// void MMModel::convCheck(bool& conv,
//                         const double& lb_old,
//                         const double& lb_new,
//                         const double& tol)
// {
//   if(((lb_new-lb_old)/lb_old < tol) | ((lb_new-lb_old) < 0.0)){
//     conv = false;
//   }
// }
/**
 GETTER FUNCTIONS
 */

arma::mat MMModel::getB()
{
  return(b_t);
}

void MMModel::getB(arma::mat& res)
{
  std::copy(b_t.begin(), b_t.end(), res.begin());
}


arma::cube MMModel::getBeta()
{
  return(beta);
}

void MMModel::getBeta(arma::cube& res)
{
  std::copy(beta.begin(), beta.end(), res.begin());
}

arma::vec MMModel::getPostMM(arma::uword p)
{
  arma::vec e_alpha(N_BLK);
  arma::vec res(N_BLK);
  if(node_est[p]){
  e_alpha.zeros();
  double row_total = 0.0;
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword m = 0; m < N_STATE; ++m){
        e_alpha[g] += alpha(g, p, m) * kappa_t(m, time_id_node[p]);
      }
      row_total += e_alpha[g] + e_c_t(g, p);
    }
    for(arma::uword g = 0; g < N_BLK; ++g){
      res(g) = (e_alpha[g] + e_c_t(g, p)) / row_total;
    }
  } else {
    res = e_c_t.col(p)/arma::sum(e_c_t.col(p));
  }
  return res;
}

arma::mat MMModel::getPostMM()
{
  arma::mat res(N_BLK, N_NODE);
  arma::vec e_alpha(N_BLK);
  double row_total;
  for(arma::uword p = 0; p < N_NODE; ++p){
    if(node_est[p]){
    e_alpha.zeros();
    row_total = 0.0;
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword m = 0; m < N_STATE; ++m){
        e_alpha[g] += alpha(g, p, m) * kappa_t(m, time_id_node[p]);
      }
      row_total += e_alpha[g] + e_c_t(g, p);
    }
    for(arma::uword g = 0; g < N_BLK; ++g){
      res(g, p) = (e_alpha[g] + e_c_t(g, p)) / row_total;
    }
    } else{
      res.col(p) = e_c_t.col(p)/arma::sum(e_c_t.col(p));
    }
  }
  return res;
}

arma::mat MMModel::getC()
{
  return e_c_t.t();
}


arma::mat MMModel::getPhi(bool send)
{
  if(send){
    return(send_phi);
  } else {
    return(rec_phi);
  }
}

arma::uvec MMModel::getN()
{
  return(tot_nodes);
}

arma::vec MMModel::getGamma()
{
  return(gamma);
}

void  MMModel::getGamma(arma::vec& res)
{
  std::copy(gamma.begin(), gamma.end(), res.begin());
}

arma::mat MMModel::getKappa()
{
  return(kappa_t);
}


arma::mat MMModel::getWmn()
{
  arma::mat res(N_STATE, N_STATE);
  if((N_TIME > 1) & (N_STATE > 1)){
    double row_total;
    for(arma::uword c = 0; c < N_STATE; ++c){
      row_total = 0.0;
      for(arma::uword r = 0; r < N_STATE; ++r){
        row_total += e_wmn_t(r, c);
      }
      for(arma::uword r = 0; r < N_STATE; ++r){
        res(r, c) = e_wmn_t(r, c) / row_total;
      }
    }
  }
  return res;
}
