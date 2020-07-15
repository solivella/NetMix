#include "MMModelClassB.h"

// [[Rcpp::export(mmsbm_fit)]]
Rcpp::List mmsbm_fit(const arma::mat& z_t,
                     const arma::mat& x1_t,//for bipartite 2 X structures
                     const arma::mat& x2_t,
                     const arma::vec& y,
                     const arma::uvec& time_id_dyad,
                     const arma::uvec& time_id_node1,
                     const arma::uvec& time_id_node2,
                     const arma::uvec& nodes_per_period,//stacked, in the case of bipartite length N_TIME+N_TIME
                     const arma::uvec& nodes_per_period1,
                     const arma::uvec& nodes_per_period2,
                     const arma::umat& node_id_dyad,
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
                     arma::cube& beta_init1,
                     arma::cube& beta_init2,
                     arma::vec& gamma_init,
                     arma::cube& theta_init,
                     const arma::vec& phi_order,
                     Rcpp::List control
)
{
  //Create model instance
  MMModelB Model(z_t,
                x1_t,//two sets, X1 and X2
                x2_t,
                y,
                time_id_dyad,
                time_id_node1,
                time_id_node2,
                nodes_per_period,//stacked family 1 nodes per period, then family 2 nodes per period
                nodes_per_period1,
                nodes_per_period2,
                node_id_dyad,
                mu_b,
                var_b,
                mu_beta1,
                var_beta1,
                mu_beta2,
                var_beta2,
                mu_gamma,
                var_gamma,
                phi_init1,
                phi_init2,
                kappa_init_t,
                b_init_t,
                beta_init1,
                beta_init2,
                gamma_init,
                theta_init,
                phi_order,
                control
		);

  // VARIATIONAL EM
  arma::uword iter = 0,
    EM_ITER = control["em_iter"],
    BIPART = control["bipartite"],
    N_NODE1 = sum(nodes_per_period1), 
    N_NODE2 = sum(nodes_per_period2), 
    N_BLK1 = control["blocks1"],
    N_BLK2 = control["blocks2"],
    N_DYAD_PRED = z_t.n_rows,
    N_MONAD_PRED1 = x1_t.n_rows,//two sets
    N_MONAD_PRED2 = control["npred2"],
    N_STATE = control["states"];
  
  bool conv = false,
    verbose = Rcpp::as<bool>(control["verbose"]);

  double tol = Rcpp::as<double>(control["conv_tol"])
    ,newLL, oldLL, change_LL;
  
  if(verbose){
    Rprintf("Estimating model...\n");
  }
  oldLL = Model.cLB();
  newLL = 0.0;
  
  while(iter < EM_ITER && conv == false){
    Rcpp::checkUserInterrupt();
    // Sample batch of dyads for stochastic VI
    Model.sampleDyads(iter);
    // E-STEP
    //Rprintf("Flag 1\n");
    Model.updatePhi();
    
    if(N_STATE > 1){
          Model.updateKappa();
       }
    //
    //M-STEP
    //Rprintf("Flag 2\n");
    Model.getB();
    //Rprintf("Flag 3\n");
    Model.getGamma();
    //Rprintf("Flag 4\n");
    Model.optim_ours(true); //optimize beta1/beta2 (alphaLB)
    //Rprintf("Flag 5\n");
    Model.optim_ours(false);//theta,gamma LB
    //
    //Check convergence
    newLL = Model.cLB();
    change_LL = fabs((newLL-oldLL)/oldLL);
    if(change_LL < tol){
      conv = true;
    } else {
      oldLL = newLL;
    }
    if(verbose){
      if((iter+1) % 50 == 0) {
        Rprintf("Iter %i, change in LB: %f\n", iter + 1, change_LL);
      }
    }
    ++iter;
  }
  if((conv == false) & verbose)
    Rprintf("Warning: model did not converge after %i iterations.\n", iter);
  else if (verbose)
    Rprintf("done after %i iterations.\n", iter);
  
  
  //Form return objects: 
  arma::mat phi_res1 = Model.getC(false, N_BLK1, N_NODE1);
  arma::mat phi_res2 = Model.getC(true, N_BLK2, N_NODE2);
  arma::mat send_phi = Model.getPhi(true);
  arma::mat rec_phi = Model.getPhi(false);
  arma::mat A = Model.getWmn(); 
  arma::mat kappa_res = Model.getKappa();
  arma::mat B = Model.getB();
  arma::vec gamma_res = Model.getGamma();
  arma::cube beta1_res = Model.getBeta1();
  arma::cube beta2_res = Model.getBeta2();
  //arma::cube Alpha1 = Model.getAlpha1();
  //arma::cube Alpha2 = Model.getAlpha2();
  
  Rcpp::List res;
  res["MixedMembership 1"] = phi_res1; 
  res["MixedMembership 2"] = phi_res2;
  res["SenderPhi"] = send_phi;
  res["ReceiverPhi"] = rec_phi;
  
  res["BlockModel"] = B;
  res["DyadCoef"] = gamma_res;
  res["TransitionKernel"] = A;
  res["MonadCoef1"] = beta1_res;
  res["MonadCoef2"] = beta2_res;
  //res["Alpha1"] = Alpha1;
  //res["Alpha2"] = Alpha2;
  
  res["Kappa"] = kappa_res;
  res["n_states"] = N_STATE;
  res["n_blocks1"] = N_BLK1;//for bipartite: 2x
  res["n_blocks2"] = N_BLK2;
  res["LowerBound"] = newLL;
  res["niter"] = iter;
  res["converged"] = conv;
  
  return res;
}
