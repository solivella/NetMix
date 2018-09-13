#include "mmsbm_fit.hpp"

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::checkUserInterrupt;
using Rcpp::as;


// [[Rcpp::export]]
List mmsbm_fit(const NumericMatrix& z_t,
	       const NumericMatrix& x_t,
	       const NumericVector& y,
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
#ifdef _OPENMP
  int requested = as<int>(control["threads"]);
  omp_set_num_threads(requested);
  N_THREADS = omp_get_max_threads();
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
    EM_ITER = as<int>(control["max_em_iter"]),
    N_NODE = sum(nodes_per_period),
    N_BLK = as<int>(control["blocks"]),
    N_DYAD_PRED = Rcpp::sum(z_t(0, Rcpp::_)) == 0 ? 0 : z_t.nrow(),
    N_MONAD_PRED = x_t.nrow(),
    N_STATE = as<int>(control["states"]);

  int TOT_BETA = N_MONAD_PRED * N_BLK * N_STATE,
    TOT_B = N_BLK * N_BLK;
  
  bool conv = false, gamma_conv = true,
    verbose = as<bool>(control["verbose"]);

  double tol = as<double>(control["conv_tol"]);

      
  NumericVector Old_B(TOT_B), Old_Gamma,
    Old_Beta(TOT_BETA);
  NumericMatrix Old_C(N_BLK, N_NODE); 
  if(N_DYAD_PRED > 0){
    Old_Gamma = NumericVector(N_DYAD_PRED);
  }

  if(verbose){
    Rprintf("Estimating model...\n");
  }
  double oldLL = Model.cLL();
  //Rprintf("\tInitial LB: %f\n",oldLL);
  while(iter < EM_ITER && conv == false){
    checkUserInterrupt();

    // E-STEP
    Model.getC(Old_C);
    
     Model.updatePhi();
     if(N_STATE > 1){
        Model.updateKappa();
      }
     // newLL = Model.cLL();
     // if(verbose){
     //   Rprintf("\tLB after E kappa %i: %f\n", iter + 1, newLL);
     // }
     
     // newLL = Model.cLL();
     // if(verbose){
     //   Rprintf("\tLB after E phi %i: %f\n", iter + 1, newLL);
     // }
     
     
    //M-STEP
    Model.getB(Old_B);
    if(N_DYAD_PRED > 0){
      Model.getGamma(Old_Gamma);
    }
    Model.getBeta(Old_Beta);
    Model.optim(true); //optimize alphaLB
    // newLL = Model.cLL();
    // if(verbose){
    //   Rprintf("\tLB after M alpha %i: %f\n", iter + 1, newLL);
    // }
    Model.optim(false); //optimize thetaLB
    
    
    //Check convergence
    double newLL = Model.cLL();
    // if(verbose){
    //   Rprintf("\tLB after M theta %i: %f\n", iter + 1, newLL);
    // }
    if(verbose){
      Rprintf("\tLB %i: %f\n", iter + 1, newLL);
    }
    gamma_conv = N_DYAD_PRED > 0 ?
      Model.checkConvChng(Old_Gamma.begin(), Old_Gamma.end(), 0, tol) :
      true;
    if(//fabs((oldLL - newLL)/oldLL) < tol
        //  &&
        (gamma_conv &&
        Model.checkConvChng(Old_B.begin(), Old_B.end(), 1, tol) &&
        Model.checkConvChng(Old_Beta.begin(), Old_Beta.end(), 2, tol)) &&
        Model.checkConvChng(Old_C.begin(), Old_C.end(), 3, tol)
         ){
      conv = true;
    }
    oldLL = newLL;
    ++iter;
  }
  if(conv == false){
    Rprintf("Warning: model did not converge after %i iterations.\n", iter);
  } else if (verbose) {
    Rprintf("...converged after %i iterations.\n", iter - 1);
  }
  
  //Form return objects
  NumericMatrix pi_res = Model.getC();
  NumericMatrix phi_send_res = Model.getPhi(true);
  NumericMatrix phi_rec_res = Model.getPhi(false);
  NumericMatrix A = Model.getWmn();
  NumericMatrix kappa_res = Model.getKappa();
  NumericMatrix B = Model.getB();
  NumericVector gamma_res = Model.getGamma(); 
  List beta_res = Model.getBeta();
  double conc = Model.getConcentration();


  List res;
  res["MixedMembership"] = pi_res;
  res["PhiSend"] = phi_send_res;
  res["PhiRec"] = phi_rec_res;
  res["BlockModel"] = B;
  res["DyadCoef"] = gamma_res;
  res["MonadCoef"] = beta_res;
  res["MMConcentration"] = conc;
  res["TransitionKernel"] = A;
  res["Kappa"] = kappa_res;
  res["n_states"] = as<int>(control["states"]);
  res["n_blocks"] = as<int>(control["blocks"]);
  res["LowerBound"] = oldLL;
  res["niter"] = iter;


  return res;
}
