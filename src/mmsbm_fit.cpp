//' @name .mmsbm_fit
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
//' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)



#include "MMModelClass.h"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::checkUserInterrupt;
using Rcpp::as;


// [[Rcpp::export(.mmsbm_fit)]]
List mmsbm_fit(const NumericMatrix& z_t,
               const NumericMatrix& x_t,
               const IntegerVector& y,
               const IntegerVector& time_id_dyad,
               const IntegerVector& time_id_node,
               const IntegerVector& nodes_per_period,
               const IntegerMatrix& node_id_dyad,
               const NumericMatrix& mu_b,
               const NumericMatrix& var_b,
               const NumericMatrix& phi_init,
               NumericMatrix& kappa_init_t,
               NumericMatrix& b_init_t,
               NumericVector& beta_init,
               NumericVector& gamma_init,
               List control
)
{
  //Obtain nr. of cores
  int N_THREADS = 1;
  int requested = as<int>(control["threads"]);
#ifdef _OPENMP
{
  omp_set_num_threads(requested);
  N_THREADS = requested;
}    
#endif

//Create model instance
MMModel Model(z_t,
              x_t,
              y,
              N_THREADS,
              time_id_dyad,
              time_id_node,
              nodes_per_period,
              node_id_dyad,
              mu_b,
              var_b,
              phi_init,
              kappa_init_t,
              b_init_t,
              beta_init,
              gamma_init,
              control
);

// VARIATIONAL EM
int iter = 0,
  EM_ITER = as<int>(control["em_iter"]),
  N_BLK = as<int>(control["blocks"]),
  N_STATE = as<int>(control["states"]);

bool conv = false,
  verbose = as<bool>(control["verbose"]);

double tol = as<double>(control["conv_tol"]);
  //,newLL, oldLL;
  

NumericVector old_beta(N_BLK * N_STATE * x_t.nrow()),
old_gamma(z_t.nrow()),
//old_c(N_BLK * x_t.ncol()),
old_bm(N_BLK * N_BLK);

// oldLL = Model.cLL();
// newLL = 0.0;
while(iter < EM_ITER && conv == false){
  checkUserInterrupt();
  
  // E-STEP
  //Model.getC(old_c);
  Model.updatePhi();
  
  
  
  if(N_STATE > 1){
    Model.updateKappa();
  }
  
  
  
  //M-STEP
  Model.getBeta(old_beta);
  Model.optim_ours(true); //optimize alphaLB
  
  Model.getGamma(old_gamma);
  Model.getB(old_bm);
  Model.optim_ours(false); //optimize thetaLB
  
  
  
  
  //Check convergence
  if(//Model.checkConv('c', old_c, tol) &&
     Model.checkConv('b', old_bm, tol) &&
     Model.checkConv('e', old_beta, tol) &&
     Model.checkConv('g', old_gamma, tol)){
    conv = true;
  }
  
  
  // newLL = Model.cLL();
  // if(fabs((newLL-oldLL)/oldLL) < tol){
  //   conv = true;
  // } else {
  //   oldLL = newLL;
  // }
  if(verbose){
    if((iter+1) % 50 == 0) {
      Rprintf("Iter %i\n", iter + 1);
    }
  }
    
  ++iter;
}


//Form return objects
NumericMatrix phi_res = Model.getC();
NumericMatrix A = Model.getWmn();
NumericMatrix kappa_res = Model.getKappa();
NumericMatrix B = Model.getB();
NumericVector gamma_res = Model.getGamma();
List beta_res = Model.getBeta();


List res;
res["MixedMembership"] = phi_res;
res["BlockModel"] = B;
res["DyadCoef"] = gamma_res;
res["TransitionKernel"] = A;
res["MonadCoef"] = beta_res;
res["Kappa"] = kappa_res;
res["n_states"] = N_STATE;
res["n_blocks"] = N_BLK;
//res["LowerBound"] = newLL;
res["niter"] = iter;
res["converged"] = conv;


return res;
}
