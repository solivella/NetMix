//' @name mmsbm_fit
//' @title Fitter Function for dynamic MMSBM Model
//' 
//' @description This is the interface to the C++ fitter for the dynamic mixed-membership
//' stochastic blockmodel for network regression.
//' 
//' @param z_t Numeric matrix; transpose of monadic design matrix. Should not include intercept row.
//' @param x_t Numeric matrix; transpose of dyadic design matrix.
//' @param y Numeric vector; vector of edge values. Must have same number of elements as \code{ncol(x_t)}
//' @param time_id_dyad Integer vector; zero-based time-period identifier for each dyad.
//' @param time_id_dyad Integer vector; zero-based time-period identifier for each node.
//' @param nodes_per_period Integer vector; total number of unique nodes observed in each time period.
//' @param node_id_dyad Integer matrix; zero-based sender and receiver identifier per dyad.
//' @param mu_b Numeric matrix; matrix of prior means for elements in blockmodel matrix.
//' @param var_b Numeric matrix; matrix of prior variances for elements in blockmodel matrix.
//' @param pi_init Numeric matrix; matrix of initial mixed-memberships. Nodes along columns.
//' @param kappa_init_t Numeric matrix; matrix of initial marginal HMM state probabilities. Time-periods along columns.
//' @param b_init_t Numeric matrix; square matrix of initial values of blockmodel.
//' @param beta_init Numeric vector; flat array (column-major order) of initial values of monadic coefficients.
//' @param gamma_init Numeric vector; vector of initial values of dyadic coefficients
//' @param control List; see the \code{mmsbm.control} argument of \code{\link{mmsbm}}
//' 
//' @return Unclassed list with named components; see \code{Value} of \code{\link{mmsbm}}
//' @section Warning:
//'          This function is for internal use only. End-users should always resort to \code{\link{mmsbm}}.
//'          In particular, that interface post-processes the return value of this internal in important ways. 
//'          
//' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)



#include "MMModelClass.h"



// [[Rcpp::export(mmsbm_fit)]]
Rcpp::List mmsbm_fit(const arma::mat& z_t,
                     const arma::mat& z_t_ho,
                     const arma::mat& x_t,
                     const arma::vec& y,
                     const arma::vec& y_ho,
                     const arma::uvec& time_id_dyad,
                     const arma::uvec& time_id_node,
                     const arma::uvec& nodes_per_period,
                     const arma::umat& node_id_dyad,
                     const arma::umat& node_id_dyad_ho,
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
                     double sparsity,
                     Rcpp::List& control)
{
  //Create model instance
  MMModel Model(z_t,
                z_t_ho,
                x_t,
                y,
                y_ho,
                time_id_dyad,
                time_id_node,
                nodes_per_period,
                node_id_dyad,
                node_id_dyad_ho,
                node_id_period,
                mu_b,
                var_b,
                mu_beta,
                var_beta,
                mu_gamma,
                var_gamma,
                pi_init,
                kappa_init_t,
                b_init_t,
                beta_init_r,
                gamma_init_r,
                sparsity,
                control
  );

  // VARIATIONAL EM
  arma::uword iter = 0, nworse = 0,
    win_size = control["conv_window"],
    VI_ITER = control["vi_iter"],
    N_BLK = control["blocks"],
    N_STATE = control["states"];
  
  bool conv = false,
    verbose = Rcpp::as<bool>(control["verbose"]);
  
  double tol = Rcpp::as<double>(control["conv_tol"]),
     newLL, oldLL;
  
  oldLL = Model.llho();
  newLL = 0.0;
  
  //arma::mat oldMM = Model.getPostMM(), newMM, corMM;
  arma::vec running_ll(win_size, arma::fill::zeros);
  std::vector<double> ll_vec;
  
  
  while(iter < VI_ITER && conv == false){
    Rcpp::checkUserInterrupt();
    // Sample batch of dyads for stochastic VI

    Model.sampleDyads(iter);
    
    // E-STEP
    Model.updatePhi();
    
    if(N_STATE > 1){
      Model.updateKappa();
    }
    // 
    // 
    // //M-STEP
    Model.optim_ours(true); //optimize alphaLB
    Model.optim_ours(false); //optimize thetaLB
    // 
    //Check convergence
    newLL = Model.llho();
    ll_vec.push_back(newLL);

    std::rotate(running_ll.begin(), running_ll.begin() + 1, running_ll.end());
    running_ll[win_size - 1] = newLL;
    if(iter > win_size) {
      Model.convCheck(nworse, conv,running_ll, tol);
    }
    oldLL = newLL;
    
    if(verbose){
      if((iter+1) % 1 == 0) {
        Rprintf("Iter: %i, held-out ll: %f\r", iter + 1, newLL);
      }
    }
    ++iter;
  }
  if(verbose){
      Rprintf("Final held-out loglik at sparsity: %f.                     \n", iter+1, newLL);
  }
  
  
  //Form return objects
  arma::mat C_res = Model.getC();
  arma::mat postmm_res = Model.getPostMM();
  arma::mat send_phi = Model.getPhi(true);
  arma::mat rec_phi = Model.getPhi(false);
  arma::uvec tot_nodes = Model.getN();
  arma::mat A = Model.getWmn();
  arma::mat kappa_res = Model.getKappa();
  arma::mat B = Model.getB();
  arma::vec gamma_res = Model.getGamma();
  arma::cube beta_res = Model.getBeta();
  
  
  Rcpp::List res;
  res["MixedMembership"] = postmm_res;
  res["CountMatrix"] = C_res;
  res["SenderPhi"] = send_phi;
  res["ReceiverPhi"] = rec_phi;
  res["TotNodes"] = tot_nodes;
  res["BlockModel"] = B;
  res["DyadCoef"] = gamma_res;
  res["TransitionKernel"] = A;
  res["MonadCoef"] = beta_res;
  res["Kappa"] = kappa_res;
  res["n_states"] = N_STATE;
  res["n_blocks"] = N_BLK;
  res["LowerBound"] = newLL;
  res["niter"] = iter + 1;
  res["converged"] = conv;
  res["ll"] = Rcpp::wrap(ll_vec);
  
  
  return res;
}
