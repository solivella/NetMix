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
//' @param phi_init Numeric matrix; matrix of initial mixed-memberships. Nodes along columns.
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
                     const arma::mat& x_t,
                     const arma::vec& y,
                     const arma::uvec& time_id_dyad,
                     const arma::uvec& time_id_node,
                     const arma::uvec& nodes_per_period,
                     const arma::umat& node_id_dyad,
                     const arma::mat& mu_b,
                     const arma::mat& var_b,
                     const arma::cube& mu_beta,
                     const arma::cube& var_beta,
                     const arma::vec& mu_gamma,
                     const arma::vec& var_gamma,
                     const arma::mat& phi_init,
                     arma::mat& kappa_init_t,
                     arma::mat& b_init_t,
                     arma::cube& beta_init,
                     arma::vec& gamma_init,
                     Rcpp::List control
)
{

  //Create model instance
  MMModel Model(z_t,
                x_t,
                y,
                time_id_dyad,
                time_id_node,
                nodes_per_period,
                node_id_dyad,
                mu_b,
                var_b,
                mu_beta,
                var_beta,
                mu_gamma,
                var_gamma,
                phi_init,
                kappa_init_t,
                b_init_t,
                beta_init,
                gamma_init,
                control
  );

  // VARIATIONAL EM
  arma::uword iter = 0,
    EM_ITER = control["em_iter"],
                     N_BLK = control["blocks"],
                                    N_STATE = control["states"];

  bool conv = false,
    verbose = Rcpp::as<bool>(control["verbose"]);
  
  double tol = Rcpp::as<double>(control["conv_tol"])
    ,newLL, oldLL, change_LL;
  
  oldLL = Model.cLL();
  newLL = 0.0;

  while(iter < EM_ITER && conv == false){
    Rcpp::checkUserInterrupt();
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
    newLL = Model.cLL();
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
  res["niter"] = iter;
  res["converged"] = conv;
  
  
  return res;
}
