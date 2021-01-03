#include "MMModelClassB.h"


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

MMModelB::MMModelB(const arma::mat& z_t,
                 const arma::mat& X1_t,
                 const arma::mat& X2_t,
                 const arma::vec& y,
                 const arma::uvec& time_id_dyad,
                 const arma::uvec& time_id_node1_r,
                 const arma::uvec& time_id_node2_r, 
                 const arma::uvec& nodes_per_period,
                 const arma::uvec& nodes_per_period1,
                 const arma::uvec& nodes_per_period2,
                 const arma::umat& node_id_dyad,
                 const arma::field<arma::uvec>& node_id_period1,
                 const arma::field<arma::uvec>& node_id_period2,
                 const arma::mat& mu_b,
                 const arma::mat& var_b,
                 const arma::cube& mu_beta1,
                 const arma::cube& var_beta1,
                 const arma::cube& mu_beta2,
                 const arma::cube& var_beta2,
                 const arma::vec& mu_gamma,
                 const arma::vec& var_gamma,
                 const arma::mat& phi_init1,
                 const arma::mat& phi_init2,
                 arma::mat& kappa_init_t,
                 arma::mat& b_init_t,
                 arma::cube& beta1_init_r,
                 arma::cube& beta2_init_r,
                 arma::vec& gamma_init_r,
                 Rcpp::List& control
		 )
:
  BIPART(control["bipartite"]),
  N_NODE1(sum(nodes_per_period1)),
  N_NODE2(sum(nodes_per_period2)), 
  N_DYAD(y.n_elem),
  N_BLK1(control["blocks1"]),
  N_BLK2(control["blocks2"]),
  N_STATE(control["states"]),
  N_TIME(control["times"]),
  N_MONAD_PRED1(X1_t.n_rows),
  N_MONAD_PRED2(X2_t.n_rows),
  N_DYAD_PRED(arma::any(z_t.row(0)) ? z_t.n_rows : 0),
  N_B_PAR(Rcpp::as<bool>(control["directed"]) || BIPART ? N_BLK1 * N_BLK2 : N_BLK1 * (1 + N_BLK1) / 2),
  OPT_ITER(control["opt_iter"]),
  N_NODE_BATCH1(arma::sum(Rcpp::as<arma::uvec>(control["batch_size1"]))),//**//N_NODE_BATCH1(control["batch_size1"]),
  N_NODE_BATCH2(arma::sum(Rcpp::as<arma::uvec>(control["batch_size2"]))),//**
  N_THREAD(control["threads"]),
  eta(Rcpp::as<double>(control["eta"])),
  forget_rate(Rcpp::as<double>(control["forget_rate"])),
  delay(Rcpp::as<double>(control["delay"])),
  var_gamma(var_gamma),
  mu_gamma(mu_gamma),
  var_beta1(var_beta1),
  mu_beta1(mu_beta1),
  var_beta2(var_beta2),
  mu_beta2(mu_beta2),
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
  time_id_dyad(time_id_dyad),
  time_id_node1(time_id_node1_r),
  time_id_node2(time_id_node2_r), 
  n_nodes_time1(nodes_per_period1),
  n_nodes_time2(nodes_per_period2),
  n_nodes_batch1(Rcpp::as<arma::uvec>(control["batch_size1"])),//**
  n_nodes_batch2(Rcpp::as<arma::uvec>(control["batch_size2"])),//**
  node_est1(Rcpp::as<arma::uvec>(control["node_est1"])),//**
  node_est2(Rcpp::as<arma::uvec>(control["node_est2"])),//**
  tot_nodes1(N_NODE1, arma::fill::zeros),
  tot_nodes2(N_NODE2, arma::fill::zeros),
  node_in_batch1(N_NODE1, arma::fill::ones),
  node_in_batch2(N_NODE2, arma::fill::ones),
  dyad_in_batch(N_DYAD, arma::fill::ones),
  node_batch1(N_NODE_BATCH1, arma::fill::zeros),
  node_batch2(N_NODE_BATCH2, arma::fill::zeros),
  maskalpha1(N_MONAD_PRED1 * N_BLK1 * N_STATE, 1),
  maskalpha2(N_MONAD_PRED2 * N_BLK2 * N_STATE, 1),
  masktheta(N_B_PAR + N_DYAD_PRED, 1),//mask=0 masked --don't update this parameter; mask=1 don't mask--update this parameter
  node_id_period1(node_id_period1),
  node_id_period2(node_id_period2),
  theta_par(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  thetaold(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  e_wm(N_STATE, arma::fill::zeros),
  alpha_gr1(N_MONAD_PRED1 * N_BLK1 * N_STATE, arma::fill::zeros),
  alpha_gr2(N_MONAD_PRED2 * N_BLK2 * N_STATE, arma::fill::zeros),
  theta_gr(N_B_PAR + N_DYAD_PRED, arma::fill::zeros),
  gamma(gamma_init_r),
  gamma_init(gamma_init_r),
  node_id_dyad(node_id_dyad),
  par_ind(N_BLK2, N_BLK1, arma::fill::zeros),
  x1_t(X1_t),
  x2_t(X2_t),
  z_t(z_t),
  mu_b_t(mu_b),
  var_b_t(var_b),
  kappa_t(kappa_init_t),
  b_t(b_init_t),
  alpha_term1(N_STATE, N_TIME, arma::fill::zeros),
  alpha_term2(N_STATE, N_TIME, arma::fill::zeros),
  send_phi(N_BLK1, N_DYAD, arma::fill::zeros),
  rec_phi(N_BLK2, N_DYAD, arma::fill::zeros),
  e_wmn_t(N_STATE, N_STATE, arma::fill::zeros),
  new_e_c_t1(N_BLK1, N_NODE1, arma::fill::zeros),
  new_e_c_t2(N_BLK2, N_NODE2, arma::fill::zeros),
  e_c_t1(N_BLK1, N_NODE1, arma::fill::zeros),
  e_c_t2(N_BLK2, N_NODE2, arma::fill::zeros),
  alpha1(N_BLK1, N_NODE1, N_STATE, arma::fill::zeros),
  alpha2(N_BLK2, N_NODE2, N_STATE, arma::fill::zeros),
  theta(N_BLK2, N_BLK1, N_DYAD, arma::fill::zeros),
  beta1(beta1_init_r),
  beta2(beta2_init_r),
  beta1old(beta1_init_r),
  beta2old(beta2_init_r),
  beta1_init(beta1_init_r),
  beta2_init(beta2_init_r)
{
  
  //Set number of parallel threads
#ifdef _OPENMP
  omp_set_num_threads(N_THREAD);
#endif
  
  //Assign initial values to  W parameters
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
    tot_nodes1[p]++;
    tot_nodes2[q]++;
    for(arma::uword g = 0; g < N_BLK1; ++g){
      send_phi(g, d) = phi_init1(g, p);
      e_c_t1(g, p) += send_phi(g, d);
      
      if(!BIPART){
	      rec_phi(g, d) = phi_init1(g, q);
	      e_c_t1(g, q) += rec_phi(g, d);
      }
    }
    if(BIPART){
      for(arma::uword h = 0; h < N_BLK2; ++h){
	      rec_phi(h, d) = phi_init2(h, q);
	      e_c_t2(h, q) += rec_phi(h, d);
      }
    }
  }
  //Create matrix of theta parameter indeces
  arma::uword ind = 0;
  for(arma::uword g = 0; g < N_BLK1; ++g){
    for(arma::uword h = 0; h < N_BLK2; ++h){
      if(directed || BIPART){
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
  for(arma::uword g = 0; g < N_BLK1; ++g){
    for(arma::uword h = 0; h < N_BLK2; ++h){
      theta_par[par_ind(h, g)] = b_t(h, g);//first 0 thru #items in B are B values in theta_par
    }
  }
  if(N_DYAD_PRED > 0)
    std::copy(gamma.begin(), gamma.end(), theta_par.begin() + N_B_PAR);
  
  Rprintf("Flag 1 ... \n");
  //Assign initial values to alpha and theta
  computeAlpha(N_NODE1,
	       N_BLK1,
	       N_MONAD_PRED1,
	       x1_t,
	       beta1,
	       alpha1,
         alpha_term1,
	       e_c_t1,
         time_id_node1,
         tot_nodes1,
         node_in_batch1,
         N_NODE_BATCH1);
  
  computeAlpha(N_NODE2,
		 N_BLK2,
		 N_MONAD_PRED2,
		 x2_t,
		 beta2,
		 alpha2,
     alpha_term2,
		 e_c_t2,
     time_id_node2,
     tot_nodes2,
     node_in_batch2,
     N_NODE_BATCH2);
  
  Rprintf("Flag 2 ... \n");
  computeTheta();
  //Define iterator at the end of param. objects
  beta1_end = beta1.end();
  beta2_end = beta2.end();
  theta_par_end = theta_par.end();
}

MMModelB::~MMModelB()
{
  
}


/**
   ALPHA LOWER BOUND
*/

double MMModelB::alphaLBInternal(
                        const arma::uword tmpNODE,
                        const arma::uword tmpBLK,
                        const arma::uword tmpPRED,
                        const arma::mat tmpX_t,
                        arma::cube& tmpBeta,
                        const arma::cube tmp_mu_beta,
                        const arma::cube tmp_var_beta,
                        arma::cube& tmpAlpha,
                        arma::mat& tmpAlpha_term,
                        arma::mat& tmpE_C_T,
                        const arma::uvec tmpTime_id_node,
                        arma::uvec tmpTot_nodes,
                        arma::uvec tmpNode_in_batch,
                        const arma::uword tmpNode_batch,
                        bool all
                        )  
{
 computeAlpha(
          tmpNODE,
          tmpBLK,
          tmpPRED,
          tmpX_t,
          tmpBeta,
          tmpAlpha,
          tmpAlpha_term,
          tmpE_C_T,
          tmpTime_id_node,
          tmpTot_nodes,
          tmpNode_in_batch,
          tmpNode_batch,
          all);
 double res = 0.0, res_int = 0.0, alpha_row = 0.0, alpha_val = 0.0;
 for(arma::uword m = 0; m < N_STATE; ++m){
#pragma omp parallel for firstprivate(alpha_val, alpha_row, res_int) reduction(+: res)
   for(arma::uword p = 0; p < tmpNODE; ++p){
     if((tmpNode_in_batch[p] == 1) || all){
     alpha_row = 0.0;
     res_int = 0.0;
     for(arma::uword g = 0; g < tmpBLK; ++g){
       alpha_val = tmpAlpha(g, p, m);
       alpha_row += alpha_val;
       res_int += (lgamma(alpha_val + tmpE_C_T(g, p)) - lgamma(alpha_val));
     }
     res_int += (lgamma(alpha_row) - lgamma(alpha_row + tmpTot_nodes[p] ));
     res += res_int * kappa_t(m, tmpTime_id_node[p]);
    }
   }
   res*= all ? 1.0: (1.*tmpNODE)/(1.*tmpNode_batch);
   //Prior for beta of correct family
   
   for(arma::uword g = 0; g < tmpBLK; ++g){
     for(arma::uword x = 1; x < tmpPRED; ++x){
      res -= 0.5 * pow(tmpBeta(x, g, m)- tmp_mu_beta(x, g, m), 2.0) / tmp_var_beta(x, g, m);  
     }
   }
 }
 return -res/(1.*tmpNODE);
}

double MMModelB::alphaLB(bool mode2, bool all)
{
  return mode2 ? 
  
  alphaLBInternal(
          N_NODE2,
          N_BLK2,
          N_MONAD_PRED2,
          x2_t,
          beta2,
          mu_beta2,
          var_beta2,
          alpha2,
          alpha_term2,
          e_c_t2,
          time_id_node2,
          tot_nodes2,
          node_in_batch2,
          N_NODE_BATCH2,
          all)
    :
    
    alphaLBInternal(
            N_NODE1,
            N_BLK1,
            N_MONAD_PRED1,
            x1_t,
            beta1,
            mu_beta1,
            var_beta1,
            alpha1,
            alpha_term1,
            e_c_t1,
            time_id_node1,
            tot_nodes1,
            node_in_batch1,
            N_NODE_BATCH1,
            all);
}


/**
   ALPHA GRADIENT
*/


void MMModelB::alphaGrInternal(int N_PAR, double *gr,
                      const arma::uword tmpNODE,
                      const arma::uword tmpBLK,
                      const arma::uword tmpPRED,
                      const arma::mat tmpX_t,
                      arma::cube& tmpBeta,
                      const arma::cube tmp_mu_beta,
                      const arma::cube tmp_var_beta,
                      arma::cube& tmpAlpha,
                      arma::mat& tmpAlpha_term,
                      arma::mat& tmpE_C_T,
                      const arma::uvec tmpTime_id_node,
                      arma::uvec tmpTot_nodes,
                      arma::uvec tmpNode_in_batch,
                      const arma::uword tmpNode_batch
                     )
{
 double res=0.0, alpha_row=0.0, prior_gr=0.0;
 arma::uword U_NPAR = N_PAR;
 for(arma::uword m = 0; m < N_STATE; ++m){
   for(arma::uword g = 0; g < tmpBLK; ++g){
     for(arma::uword x = 0; x < tmpPRED; ++x){
       res = 0.0;
#pragma omp parallel for firstprivate(alpha_row) reduction(+: res)
       for(arma::uword p = 0; p < tmpNODE; ++p){
         if(tmpNode_in_batch[p] == 1) {
         alpha_row = 0.0;
         for(arma::uword h = 0; h < tmpBLK; ++h){
           alpha_row += tmpAlpha(h, p, m);
         }
         res += (R::digamma(alpha_row) - R::digamma(alpha_row + tmpTot_nodes[p] )
         + R::digamma(tmpAlpha(g, p, m) + tmpE_C_T(g, p)) - R::digamma(tmpAlpha(g, p, m)))
           * kappa_t(m, tmpTime_id_node[p]) * tmpAlpha(g, p, m) * tmpX_t(x, p);
        }
       }
       res*= (1. * tmpNODE)/(1. * tmpNode_batch);
       prior_gr =  (tmpBeta(x, g, m) - tmp_mu_beta(x, g, m))/tmp_var_beta(x, g, m);//x > 0 ? (tmpBeta(x, g, m) - tmp_mu_beta(x, g, m)) / tmp_var_beta(x, g, m) : 0.0;
       gr[x + tmpPRED * (g + tmpBLK * m)] = -(res - prior_gr);//(res - prior_gr);//
     }
   }
 }
 for(arma::uword i = 0; i < U_NPAR; ++i){
   gr[i] /= (1. *tmpNODE);
 }
}

void MMModelB::alphaGr(bool mode2, int n, double* gr)
{
  mode2 ? alphaGrInternal(n,
                  gr,
                  N_NODE2,
                  N_BLK2,
                  N_MONAD_PRED2,
                  x2_t,
                  beta2,
                  mu_beta2,
                  var_beta2,
                  alpha2,
                  alpha_term2,
                  e_c_t2,
                  time_id_node2,
                  tot_nodes2,
                  node_in_batch2,
                  N_NODE_BATCH2)
    :
    
    alphaGrInternal(n,
            gr,
            N_NODE1,
            N_BLK1,
            N_MONAD_PRED1,
            x1_t,
            beta1,
            mu_beta1,
            var_beta1,
            alpha1,
            alpha_term1,
            e_c_t1,
            time_id_node1,
            tot_nodes1,
            node_in_batch1,
            N_NODE_BATCH1);
}

/**
   ALPHA COMPUTATION
*/

void MMModelB::computeAlpha(
                           const arma::uword tmpNODE,
                           const arma::uword tmpBLK,
                           const arma::uword tmpPRED,
			                     const arma::mat& tmpX_t,
			                     arma::cube& tmpBeta,
			                     arma::cube& tmpAlpha,
			                     arma::mat& tmpAlpha_term,
			                     arma::mat& tmpE_C_T,
			                     const arma::uvec tmpTime_id_node,
			                     arma::uvec tmpTot_nodes,
			                     arma::uvec tmpNode_in_batch,
			                     const arma::uword tmpNode_batch,
			                     bool all
			   )
{
std::fill(tmpAlpha_term.begin(), tmpAlpha_term.end(), 0.0);
  double linpred, row_sum, correct_fact;//**
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword p = 0; p < tmpNODE; ++p){
      if((tmpNode_in_batch[p] == 1) || all){
        correct_fact = all ? 1.0: (1. * tmpNODE)/tmpNode_batch;
      row_sum = 0.0;
      for(arma::uword g = 0; g < tmpBLK; ++g){
        linpred = 0.0;
        for(arma::uword x = 0; x < tmpPRED; ++x){
          linpred += tmpX_t(x, p) * tmpBeta(x, g, m);
        }
        linpred = exp(linpred);
        row_sum += linpred;
        tmpAlpha(g, p, m) = linpred;
        tmpAlpha_term(m, tmpTime_id_node[p]) += correct_fact * (lgamma(linpred + tmpE_C_T(g, p)) - lgamma(linpred));
      }
      tmpAlpha_term(m, tmpTime_id_node[p]) += correct_fact * (lgamma(row_sum) - lgamma(row_sum + tmpTot_nodes[p]));
     }
    }
  }
}


/**
   THETA LB
*/
double MMModelB::thetaLB(bool entropy, bool all)
{
  computeTheta();
  
  double res = 0.0;
#pragma omp parallel for reduction(+: res)
  for(arma::uword d = 0; d < N_DYAD; ++d){
    if((dyad_in_batch[d] == 1) || all){
    for(arma::uword g = 0; g < N_BLK1; ++g){
      for(arma::uword h = 0; h < N_BLK2; ++h){
        res += send_phi(g, d) * rec_phi(h, d) * (y[d] * log(theta(h, g, d)) + (1.0 - y[d]) * log(1.0 - theta(h, g, d)));
      }
    }
  }
  }
  res *= all ? 1.0 : reweightFactor;
  //Prior for gamma
  for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
    res -= 0.5*pow(gamma[z] - mu_gamma[z], 2.0) / var_gamma[z];
  }
  //Prior for B
  for(arma::uword g = 0; g < N_BLK1; ++g){
    for(arma::uword h = 0; h < N_BLK2; ++h){
      res -=0.5*(pow(b_t(h,g) - mu_b_t(h, g), 2.0) / var_b_t(h, g));
    }
  }
  return -res / N_DYAD;
}


/**
   GRADIENT FOR THETA
*/
void MMModelB::thetaGr(int N_PAR, double *gr){
  
 arma::uword U_NPAR = N_PAR;
 double res_local, res = 0.0;
 for(arma::uword i = 0; i < U_NPAR; ++i){
   gr[i] = 0.0;
 }
 arma::uword npar=0;
 for(arma::uword d = 0; d < N_DYAD; ++d){
   if(dyad_in_batch[d] == 1){
   res = 0.0;
   for(arma::uword g = 0; g < N_BLK1; ++g){
     for(arma::uword h = 0; h < N_BLK2; ++h){
       res_local = send_phi(g, d) * rec_phi(h, d) * (y[d] - theta(h, g, d));
       res += res_local;
       if(!BIPART){ if((h < g) && !directed){
          continue;
       }
       }
       npar = par_ind(h, g);
       gr[npar] -= res_local;//for B parameters, just phi*psi*(Y-theta) portion
     }
   }
     if(N_DYAD_PRED > 0){
       for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
         gr[N_B_PAR + z] -= res * z_t(z, d);//for gamma_j parameters, phi*psi(Y-theta) multiply gradient by d_pqtj  (no need for B parameters)
       }
     }
   }
 }
 for(arma::uword i = 0; i < U_NPAR; ++i){
   gr[i] *= reweightFactor; //for stochastic VI
 }
 for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
   gr[N_B_PAR + z] += (gamma[z] - mu_gamma[z]) / var_gamma[z];
 }
 for(arma::uword g = 0; g < N_BLK1; ++g){
   for(arma::uword h = 0; h < N_BLK2; ++h){
     if(!BIPART){
       if((h < g) && !directed){
       continue;
     }
     }
     npar = par_ind(h, g);
     gr[npar] += (b_t(h, g) - mu_b_t(h, g)) / var_b_t(h, g);
    }
 }
 for(arma::uword i = 0; i < U_NPAR; ++i){
   gr[i] /= N_DYAD;
 }
}


/**
   COMPUTE THETA
*/

void MMModelB::computeTheta(bool all){
  
 for(arma::uword g = 0; g < N_BLK1; ++g){
   for(arma::uword h = 0; h < N_BLK2; ++h){
      b_t(h, g) = theta_par[par_ind(h, g)];
   }
 }
 double linpred=0.0;
 for(arma::uword d = 0; d < N_DYAD; ++d){
  if((dyad_in_batch[d] == 1) || all){
   linpred = 0.0;
   for(arma::uword z = 0; z < N_DYAD_PRED; ++z){
     gamma[z] = theta_par[N_B_PAR + z];
     linpred -= z_t(z, d) * gamma[z];
   }
     for(arma::uword g = 0; g < N_BLK1; ++g){
       for(arma::uword h = 0; h < N_BLK2; ++h){
         theta(h, g, d) = 1./(1 + exp(linpred - b_t(h, g)));
       }
     }
    }
 }
}

/**
  OPTIMIZATION
*/

void MMModelB::optim_ours(bool alpha)
{
  if(alpha){
    int npar1 = N_MONAD_PRED1 * N_BLK1 * N_STATE;
    beta1old = beta1;
    std::copy(beta1_init.begin(), beta1_init.end(), beta1.begin());
    vmmin_ours(npar1, &beta1[0], &fminAlpha, alphaLBWMode1, alphaGrWMode1, OPT_ITER, 0,
               &maskalpha1[0], -1.0e+35, 1.0e-6, 1, this, &fncountAlpha, &grcountAlpha, &m_failAlpha);
    for(arma::uword i = 0; i < npar1; ++i){
      beta1[i] = (1 - step_size) * beta1old[i] + step_size * beta1[i]; 
    }
  
      if(BIPART){
      int npar2=N_MONAD_PRED2 * N_BLK2 * N_STATE;
      beta2old = beta2;
      std::copy(beta2_init.begin(), beta2_init.end(), beta2.begin());
      vmmin_ours(npar2, &beta2[0], &fminAlpha, alphaLBWMode2, alphaGrWMode2, OPT_ITER, 0,
                   &maskalpha2[0], -1.0e+35, 1.0e-6, 1, this, &fncountAlpha, &grcountAlpha, &m_failAlpha);
        for(arma::uword i = 0; i < npar2; ++i){
          beta2[i] = (1 - step_size) * beta2old[i] + step_size * beta2[i]; 
        }
      }
  } else {
    int npar = N_B_PAR + N_DYAD_PRED;
    thetaold = theta_par;
    std::copy(gamma_init.begin(), gamma_init.end(), theta_par.begin() + N_B_PAR);//copy gamma init values over to theta_par after B params
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

double MMModelB::thetaLBW(int n, double *par, void *ex)
{
 return(static_cast<MMModelB*>(ex)->thetaLB());
}
void MMModelB::thetaGrW(int n, double *par, double *gr, void *ex)
{
  static_cast<MMModelB*>(ex)->thetaGr(n, gr);
}
double MMModelB::alphaLBWMode1(int n, double *par, void *ex)
{
  return(static_cast<MMModelB*>(ex)->alphaLB(false, false)); //alphaLB(mode2,all)
}
double MMModelB::alphaLBWMode2(int n, double *par, void *ex)
{
  return(static_cast<MMModelB*>(ex)->alphaLB(true, false));
}
void MMModelB::alphaGrWMode1(int n, double *par, double *gr, void *ex)
{
  static_cast<MMModelB*>(ex)->alphaGr(false, n, gr);
}
void MMModelB::alphaGrWMode2(int n, double *par, double *gr, void *ex)
{
  static_cast<MMModelB*>(ex)->alphaGr(true, n, gr);
}

/**
 * LOG LIKELIHOOD
 */
double MMModelB::LB()
{
  double res = lgamma(double(N_STATE) * eta) - lgamma(eta);
  res -= thetaLB(true,true);
  res -= alphaLB(true,true);
  res -= alphaLB(false,true);
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
  return res;
}
double MMModelB::LL()
{
  double ll = 0.0;
#pragma omp parallel for reduction(-:ll)
  for(arma::uword d = 0; d < N_DYAD; ++d){
    if(dyad_in_batch[d]==0){
      arma::uword p, q;
      p = node_id_dyad(d, 0);
      q = node_id_dyad(d, 1);
      ll -= log(arma::as_scalar(1.0 + exp(-(getPostMM1(p).t() * b_t * getPostMM2(q)
                                              + gamma.t() * z_t.col(d)))));
    }
  }
  return(ll);
}

/**
   VARIATIONAL UPDATE FOR KAPPA
*/


void MMModelB::updateKappa()
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
     res += alpha_term1(m, t) + alpha_term2(m,t);//**two alpha terms
     kappa_vec[m] = res;
   }
   log_denom = logSumExp(kappa_vec);
   for(arma::uword m = 0; m < N_STATE; ++m){
     kappa_t(m, t) = exp(kappa_vec[m] - log_denom);
     if(!std::isfinite(kappa_t(m, t))){
       stop("Kappa value became NAN.");
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

// void MMModelB::updatePhiInternal(
//                               arma::uword dyad,
//                               arma::uword rec,
//                               double *phi,
//                               double *phi_o,
//                               double *new_c,
//                               arma::uword *err,
//                               arma::cube& theta,
//                               arma::cube& alpha,
//                               arma::mat& e_c_t,
//                               const arma::uvec& time_id_dyad,
//                               const arma::vec& y,
//                               arma::mat& kappa_t, 
//                               const arma::umat& node_id_dyad,
//                               arma::uword tmpBLK,
//                               arma::uword tmpBLK_o,
//                               arma::uword N_STATE
// 				)
// {
//   
//   arma::uword t = time_id_dyad[dyad];
//   double edge = y[dyad];
//   arma::uword incr1 = rec ? 1 : tmpBLK;
//   arma::uword incr2 = rec ? tmpBLK : 1;
//   arma::uword node = node_id_dyad(dyad, rec);
//   double *theta_temp = &theta(0, 0, dyad);
//   double *te;
//   arma::vec phi_old(tmpBLK, arma::fill::zeros);
// 
//   double total = 0.0, res;
//   for(arma::uword g = 0; g < tmpBLK; ++g, theta_temp+=incr1){
//     phi_old[g] = phi[g];
//     new_c[g] -= phi[g];
//     
//     res = 0.0;
//     for(arma::uword m = 0; m < N_STATE; ++m){
//       res += kappa_t(m, t) * log(alpha(g, node, m) + std::max(new_c[g], 0.0));
//     }
//     
//     te = theta_temp;
//     for(arma::uword h = 0; h < tmpBLK_o; ++h, te+=incr2){
//       double theta_tmp = rec ? theta(g,h,dyad): theta(h,g,dyad);
//       res += phi_o[h] * (edge * log(theta_tmp) + (1.0 - edge) * log(1.0 - theta_tmp));
//     }
//     phi[g] = exp(res);
//     if(!std::isfinite(phi[g])){
//       (*err)++;
//     }
//     total += phi[g];
//   }
//   
//   //Normalize phi to sum to 1
//   //and store new value in c
//   for(arma::uword g = 0; g < tmpBLK; ++g){
//     phi[g] /= total;
//     new_c[g] += phi[g];      
//   }
// }

//add updatePhiInternal1 and updatePhiInternal2 to check if something off with 2 structures
void MMModelB::updatePhiInternal1(
    arma::uword dyad,
    arma::uword rec,
    double *phi1,
    double *phi2,
    double *new_c1,
    arma::uword *err,
    arma::cube& theta,
    arma::cube& alpha1,
    arma::mat& e_c_t1,
    const arma::uvec& time_id_dyad,
    const arma::vec& y,
    arma::mat& kappa_t, 
    const arma::umat& node_id_dyad,
    arma::uword G,
    arma::uword H,
    arma::uword N_STATE
)
{
  arma::uword t = time_id_dyad[dyad];
  double edge = y[dyad];
  arma::uword node = node_id_dyad(dyad, 0); //rec=0
  arma::vec phi_old(N_BLK1, arma::fill::zeros);
  
  
  double total = 0.0, res;
  for(arma::uword g = 0; g < G; ++g){
    phi_old[g] = phi1[g];
    new_c1[g] -= phi1[g];
    res = 0.0;
    for(arma::uword m = 0; m < N_STATE; ++m){
      res += kappa_t(m, t) * log(alpha1(g, node, m) + std::max(new_c1[g], 0.0));// + e_c_t1(g, node) - phi1[g]);
    }
    for(arma::uword h = 0; h < H; ++h){
      double theta_tmp = theta(h,g,dyad);
      res += phi2[h] * (edge * log(theta_tmp) + (1.0 - edge) * log(1.0 - theta_tmp));
    }
    phi1[g] = exp(res);
    if(!std::isfinite(phi1[g])){
      (*err)++;
    }
    total += phi1[g];
  }
  for(arma::uword g = 0; g < G; ++g){
    phi1[g] /= total;
    new_c1[g] += phi1[g];      
  }
}

void MMModelB::updatePhiInternal2(
    arma::uword dyad,
    arma::uword rec,
    double *phi2,
    double *phi1,
    double *new_c2,
    arma::uword *err,
    arma::cube& theta,
    arma::cube& alpha2,
    arma::mat& e_c_t2,
    const arma::uvec& time_id_dyad,
    const arma::vec& y,
    arma::mat& kappa_t, 
    const arma::umat& node_id_dyad,
    arma::uword H,
    arma::uword G,
    arma::uword N_STATE
)
{
  arma::uword t = time_id_dyad[dyad];
  double edge = y[dyad];
  arma::uword node = node_id_dyad(dyad, 1);
  arma::vec phi_old(N_BLK2, arma::fill::zeros);
  
  double total = 0.0, res;
  for(arma::uword h = 0; h < H; ++h){
    phi_old[h] = phi2[h];
    new_c2[h] -= phi2[h];
    res = 0.0;
    for(arma::uword m = 0; m < N_STATE; ++m){
      res += kappa_t(m, t) * log(alpha2(h, node, m) + std::max(new_c2[h], 0.0));
    }
    
    for(arma::uword g = 0; g < G; ++g){
      double theta_tmp = theta(h,g,dyad);
      res += phi1[g] * (edge * log(theta_tmp) + (1.0 - edge) * log(1.0 - theta_tmp));
    }
    phi2[h] = exp(res);
    if(!std::isfinite(phi2[h])){
      (*err)++;
    }
    total += phi2[h];
  }
  //Normalize phi to sum to 1
  //and store new value in c
  for(arma::uword h = 0; h < H; ++h){
    phi2[h] /= total;
    new_c2[h] += phi2[h];      
  }
}



void MMModelB::updatePhi()
{

  arma::uword err = 0;
  // Update dyads with sampled nodes
  for(arma::uword d = 0; d < N_DYAD; ++d){
    Rcpp::checkUserInterrupt();
    if(node_est1[node_id_dyad(d, 0)]) {
    updatePhiInternal1(d,
                      0,
                      &(send_phi(0, d)),
                      &(rec_phi(0, d)),
                      //&(new_e_c_t1(0, node_id_dyad(d, 0))),
                      &(e_c_t1(0,node_id_dyad(d,0))),
                      &err,
                      theta,
                      alpha1,
                      e_c_t1,
                      time_id_dyad,
                      y,
                      kappa_t,
                      node_id_dyad,
                      N_BLK1,
                      N_BLK2,
                      N_STATE);
    }
    if(node_est2[node_id_dyad(d, 1)]) {
    updatePhiInternal2(d,
                       1,
                       &(rec_phi(0, d)),
                       &(send_phi(0, d)),
                       //&(new_e_c_t2(0, node_id_dyad(d, 1))),
                       &(e_c_t2(0,node_id_dyad(d,1))),
                       &err,
                       theta,
                       alpha2,
                       e_c_t2,
                       time_id_dyad,
                       y,
                       kappa_t,
                       node_id_dyad,
                       N_BLK2,
                       N_BLK1,
                       N_STATE);
    }
    //**edit out below to test out updatePhiInternal1 and 2 and separate fns
    // if(phi_order[d] >=1 ){
    // updatePhiInternal(d,
    //                   0,
    //                   &(send_phi(0, d)),
    //                   &(rec_phi(0, d)),
    //                   &(new_e_c_t1(0, node_id_dyad(d, 0))),
    //                   &err,
    //                   theta,
    //                   alpha1,
    //                   e_c_t1,
    //                   time_id_dyad,
    //                   y,
    //                   kappa_t,
    //                   node_id_dyad,
    //                   N_BLK1,
    //                   N_BLK2,
    //                   N_STATE);
    // updatePhiInternal(d,
    //                   1,
    //                   &(rec_phi(0, d)),
    //                   &(send_phi(0, d)), //using most up to date: create copy send_phi (old_send_phi) before 769, use old send_phi in 789
    //                   &(new_e_c_t2(0, node_id_dyad(d, 1))),
    //                   &err,
    //                   theta,
    //                   alpha2,
    //                   e_c_t2,
    //                   time_id_dyad,
    //                   y,
    //                   kappa_t,
    //                   node_id_dyad,
    //                   N_BLK2,
    //                   N_BLK1,
    //                   N_STATE);
    // }else{
    // updatePhiInternal(d,
    //                     1,
    //                     &(rec_phi(0, d)),
    //                     &(send_phi(0, d)), //using most up to date: create copy send_phi (old_send_phi) before 769, use old send_phi in 789
    //                     &(new_e_c_t2(0, node_id_dyad(d, 1))),
    //                     &err,
    //                     theta,
    //                     alpha2,
    //                     e_c_t2,
    //                     time_id_dyad,
    //                     y,
    //                     kappa_t,
    //                     node_id_dyad,
    //                     N_BLK2,
    //                     N_BLK1,
    //                     N_STATE);
    //   updatePhiInternal(d,
    //                     0,
    //                     &(send_phi(0, d)),
    //                     &(rec_phi(0, d)),
    //                     &(new_e_c_t1(0, node_id_dyad(d, 0))),
    //                     &err,
    //                     theta,
    //                     alpha1,
    //                     e_c_t1,
    //                     time_id_dyad,
    //                     y,
    //                     kappa_t,
    //                     node_id_dyad,
    //                     N_BLK1,
    //                     N_BLK2,
    //                     N_STATE);
    // }//**edit out above to test out updatePhiInternal1 and 2 and separate fns
  }
  
  if(err){
    Rcpp::stop("Phi value became NaN.");
  }
  
}


void MMModelB::sampleDyads(arma::uword iter)
{
  // Sample nodes for stochastic variational update
  if(N_TIME < 2){
    node_batch1 = arma::randperm(N_NODE1, n_nodes_batch1[0]);
    node_batch2 = arma::randperm(N_NODE2, n_nodes_batch2[0]);
  } else{ // N_TIME  >= 2
    arma::uvec node_batch1_t, node_batch2_t;
    arma::uword ind1 = 0, ind2 = 0;
    for(arma::uword t = 0; t < N_TIME; ++t){
      node_batch1_t = arma::randperm(n_nodes_time1[t], n_nodes_batch1[t]);
      node_batch2_t = arma::randperm(n_nodes_time2[t], n_nodes_batch2[t]);
      for(arma::uword p = 0; p < node_batch1_t[t]; ++p){
        node_batch1[ind1] = node_id_period1[t][node_batch1_t[p]];
        ++ind1;
      }
      for(arma::uword q = 0; q < node_batch2_t[t]; ++q){
        node_batch2[ind2] = node_id_period2[t][node_batch2_t[q]];
        ++ind2;
      }
    }
  }
  
  for(arma::uword p = 0; p < N_NODE1; ++p){
    node_in_batch1[p] = arma::any(node_batch1 == p) ? 1 : 0;
  }
  
  for(arma::uword p = 0; p < N_NODE2; ++p){
    node_in_batch2[p] = arma::any(node_batch2 == p) ? 1 : 0;
  }

  for(arma::uword d = 0; d < N_DYAD; ++d){
    dyad_in_batch[d] = (arma::any(node_batch1 == node_id_dyad(d, 0)) //family 1 node | family 2 node
                          | arma::any(node_batch2 == node_id_dyad(d, 1))) ? 1 : 0;
  }
  reweightFactor = (1. * N_DYAD) / arma::sum(dyad_in_batch);
  step_size = 1.0 / pow(delay + iter, forget_rate);
}

/**
 * CONVERGENCE CHECKER
 */
void MMModelB::convCheck(bool& conv,
                        const arma::cube& beta1_new,
                        const arma::cube& beta1_old,
                        const arma::cube& beta2_new,
                        const arma::cube& beta2_old,
                        const arma::mat& b_new,  
                        const arma::mat& b_old, 
                        const arma::vec& gamma_new, 
                        const arma::vec& gamma_old, 
                        const double& tol)
{
  //double cor_val = arma::as_scalar(arma::cor(ll, arma::regspace(1, ll.n_elem)));
  // Rcpp::Rcout << b_old << std::endl;
  // Rcpp::Rcout << b_new << std::endl;
  // double c1 = arma::as_scalar(arma::cor(arma::vectorise(beta_new), arma::vectorise(beta_old)));
  // double c2 = arma::as_scalar(arma::cor(arma::vectorise(b_new), arma::vectorise(b_old)));
  // double c3 = arma::as_scalar(arma::cor(gamma_new, gamma_old));
  // Rcpp::NumericVector cor_vec = {c1, c2, c3};
  // //Rcpp::Rcout << cor_vec << std::endl;
  // Rcpp::LogicalVector nans = Rcpp::is_nan(cor_vec);
  // cor_vec[nans] = 1.0;
  // conv = Rcpp::is_true(Rcpp::all(cor_vec > (1.0 - tol)));
  
  arma::cube::const_iterator beta1_old_it = beta1_old.begin(),
    beta2_old_it = beta2_old.begin(),
    beta1_new_it = beta1_new.begin(),
    beta2_new_it = beta2_new.begin(),
    beta1_old_end = beta1_old.end(),
    beta2_old_end = beta2_old.end(),
    beta1_new_end = beta1_new.end(),
    beta2_new_end = beta2_new.end();;
  
  arma::mat::const_iterator b_old_it = b_old.begin(),
    b_new_it = b_new.begin(),
    b_old_end = b_old.end(),
    b_new_end = b_new.end();
  
  arma::vec::const_iterator gamma_old_it = gamma_old.begin(),
    gamma_new_it = gamma_new.begin(),
    gamma_old_end = gamma_old.end(),
    gamma_new_end = gamma_new.end();
  
  
  conv = true;
  
  for( ; beta1_new_it != beta1_new_end; ++beta1_new_it,++beta1_old_it){
    if(fabs(*beta1_new_it - *beta1_old_it) > tol){
      conv = false;
      return;
    }  
  }
  
  for( ; beta2_new_it != beta2_new_end; ++beta2_new_it,++beta2_old_it){
    if(fabs(*beta2_new_it - *beta2_old_it) > tol){
      conv = false;
      return;
    }  
  }
  
  for( ; b_new_it != b_new_end; ++b_new_it,++b_old_it){
    if(fabs(*b_new_it - *b_old_it) > tol){
      conv = false;
      return;
    }  
  }
  for( ; gamma_new_it != gamma_new_end; ++gamma_new_it,++gamma_old_it){
    if(fabs(*gamma_new_it - *gamma_old_it) > tol){
      conv = false;
      return;
    }  
  }
  
  //Rprintf("Cor is %f\n", cor_val);
  //Rcpp::Rcout << ll << std::endl;
  //conv = (fabs(arma::as_scalar((ll.tail(1)-ll.head(1))/ll.head(1))) < tol);
  //conv = (cor_val > (1.0-tol));
  // if((arma::stddev(ll) == 0.0) | (fabs(cor_val) <= tol)){
  //   conv = true;
  // }
}

/**
   GETTER FUNCTIONS
*/

arma::mat MMModelB::getB()
 {
   return(b_t);
 }

void MMModelB::getB(arma::mat& res) {
   std::copy(b_t.begin(), b_t.end(), res.begin());
 }

arma::cube MMModelB::getBeta1(){
   return(beta1);
 }
void MMModelB::getBeta1(arma::cube& res) {
   std::copy(beta1.begin(), beta1.end(), res.begin());
 }
arma::cube MMModelB::getBeta2(){
  return(beta2);
}
void MMModelB::getBeta2(arma::cube& res){
  std::copy(beta2.begin(), beta2.end(), res.begin());
}

arma::cube MMModelB::getAlpha1(){
  return(alpha1);
}
arma::cube MMModelB::getAlpha2(){
  return(alpha2);
}

arma::vec MMModelB::getPostMM1(arma::uword p)// one p
{
  arma::vec e_alpha(N_BLK1);
  arma::vec res(N_BLK1);
  if(node_est1[p]){
    e_alpha.zeros();
    double row_total = 0.0;
    for(arma::uword g = 0; g < N_BLK1; ++g){
      for(arma::uword m = 0; m < N_STATE; ++m){
        e_alpha[g] += alpha1(g, p, m) * kappa_t(m, time_id_node1[p]);
      }
      row_total += e_alpha[g] + e_c_t1(g, p);
    }
    for(arma::uword g = 0; g < N_BLK1; ++g){
      res(g) = (e_alpha[g] + e_c_t1(g, p)) / row_total;
    }
  } else {
    res = e_c_t1.col(p)/arma::sum(e_c_t1.col(p));
  }
  return res;
}

arma::mat MMModelB::getPostMM1() //all p
{
  arma::mat res(N_BLK1, N_NODE1);
  arma::vec e_alpha(N_BLK1);
  double row_total;
  for(arma::uword p = 0; p < N_NODE1; ++p){
    if(node_est1[p]){
      e_alpha.zeros();
      row_total = 0.0;
      for(arma::uword g = 0; g < N_BLK1; ++g){
        for(arma::uword m = 0; m < N_STATE; ++m){
          e_alpha[g] += alpha1(g, p, m) * kappa_t(m, time_id_node1[p]);
        }
        row_total += e_alpha[g] + e_c_t1(g, p);
      }
      for(arma::uword g = 0; g < N_BLK1; ++g){
        res(g, p) = (e_alpha[g] + e_c_t1(g, p)) / row_total;
      }
    } else{
      res.col(p) = e_c_t1.col(p)/arma::sum(e_c_t1.col(p));
    }
  }
  return res;
}

arma::vec MMModelB::getPostMM2(arma::uword p)
{
  arma::vec e_alpha(N_BLK2);
  arma::vec res(N_BLK2);
  if(node_est2[p]){
    e_alpha.zeros();
    double row_total = 0.0;
    for(arma::uword g = 0; g < N_BLK2; ++g){
      for(arma::uword m = 0; m < N_STATE; ++m){
        e_alpha[g] += alpha2(g, p, m) * kappa_t(m, time_id_node2[p]);
      }
      row_total += e_alpha[g] + e_c_t2(g, p);
    }
    for(arma::uword g = 0; g < N_BLK2; ++g){
      res(g) = (e_alpha[g] + e_c_t2(g, p)) / row_total;
    }
  } else {
    res = e_c_t2.col(p)/arma::sum(e_c_t2.col(p));
  }
  return res;
}

arma::mat MMModelB::getPostMM2()
{
  arma::mat res(N_BLK2, N_NODE2);
  arma::vec e_alpha(N_BLK2);
  double row_total;
  for(arma::uword p = 0; p < N_NODE2; ++p){
    if(node_est2[p]){
      e_alpha.zeros();
      row_total = 0.0;
      for(arma::uword g = 0; g < N_BLK2; ++g){
        for(arma::uword m = 0; m < N_STATE; ++m){
          e_alpha[g] += alpha2(g, p, m) * kappa_t(m, time_id_node2[p]);
        }
        row_total += e_alpha[g] + e_c_t2(g, p);
      }
      for(arma::uword g = 0; g < N_BLK2; ++g){
        res(g, p) = (e_alpha[g] + e_c_t2(g, p)) / row_total;
      }
    } else{
      res.col(p) = e_c_t2.col(p)/arma::sum(e_c_t2.col(p));
    }
  }
  return res;
}

arma::mat MMModelB::getC(bool mode2
                           //,arma::uword tmpBLK
                           //,arma::uword tmpNODE
                            )
{
  if(mode2){
    return e_c_t1.t();
  }else{
    return e_c_t2.t();
  }
  //arma::mat& tmpE_C_T = mode2 ? e_c_t2 : e_c_t1; 
  // arma::mat res(tmpBLK, tmpNODE);
  // double row_total;
  // for(arma::uword p = 0; p < tmpNODE; ++p){
  //   row_total = 0.0;
  //   for(arma::uword g = 0; g < tmpBLK; ++g){
  //     row_total += tmpE_C_T(g, p);
  //   }
  //   for(arma::uword g = 0; g < tmpBLK; ++g){
  //     res(g, p) = tmpE_C_T(g, p) / row_total;
  //   }
  // }
  // return res;
}

arma::mat MMModelB::getPhi(bool send){
  if(send){
    return(send_phi);
  } else {
    return(rec_phi);
  }
}

arma::uvec MMModelB::getN(bool send)
{
  if(send){
  return(tot_nodes1);
  } else{
    return(tot_nodes2);
  }
}

arma::vec MMModelB::getGamma(){
  return(gamma);
}
void  MMModelB::getGamma(arma::vec& res){
  std::copy(gamma.begin(), gamma.end(), res.begin());
}

arma::mat MMModelB::getKappa(){
  return(kappa_t);
}

arma::mat MMModelB::getWmn()
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
