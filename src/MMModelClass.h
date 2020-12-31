#ifndef MMMODEL_CLASS
#define MMMODEL_CLASS

#include <vector>

// #ifndef DEBUG_MODE
// #define DEBUG_MODE
// #define ARMA_EXTRA_DEBUG
// #endif 

#include <RcppArmadillo.h>
#include "AuxFuns.h"
#ifdef _OPENMP
#include <omp.h>
#endif



class MMModel
{
public:
  MMModel(const arma::mat& z_t,
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
          Rcpp::List& control
  );
  ~MMModel();
  
  
  void sampleDyads(arma::uword iter);
  void updatePhi();
  void updateKappa();
  void optim_ours(bool);
  double LL();
  double LB();
  //double llho();
  
  
  arma::mat getPostMM();
  arma::vec getPostMM(arma::uword);
  arma::mat getC();
  arma::mat getPhi(bool);
  arma::uvec getN();
  arma::mat getWmn();
  arma::mat getKappa();
  arma::mat getB();
  void getB(arma::mat&);
  arma::vec getGamma();
  void getGamma(arma::vec&);
  arma::cube getBeta();
  void getBeta(arma::cube&);
  void convCheck(bool& conv,
                 const arma::cube& beta_new,
                 const arma::cube& beta_old,
                 const arma::mat& b_new,  
                 const arma::mat& b_old, 
                 const arma::vec& gamma_new, 
                 const arma::vec& gamma_old, 
                 const double& tol);
  
  
  
private:
  
  const arma::uword N_NODE,
  N_DYAD,
  N_BLK,
  N_STATE,
  N_TIME,
  N_MONAD_PRED,
  N_DYAD_PRED,
  N_B_PAR,
  OPT_ITER,
  N_NODE_BATCH,
  N_THREAD;
  //N_DYAD_HO;
  
  const double eta,
  forget_rate,
  delay;//,
  //sparsity;
  
  const arma::vec var_gamma,
  mu_gamma;
  const arma::cube var_beta,
  mu_beta;
  
  double fminAlpha,
  fminTheta,
  reweightFactor,
  step_size;
  
  
  int fncountAlpha,
  fncountTheta,
  grcountAlpha,
  grcountTheta,
  m_failAlpha,
  m_failTheta;
  
  
  bool verbose,
  directed;
  
  const arma::vec y;// y_ho;
  
  const arma::uvec time_id_dyad,
  time_id_node,
  n_nodes_time,
  n_nodes_batch,
  node_est;
  
  arma::uvec tot_nodes,
  node_in_batch,
  dyad_in_batch,
  node_batch;
  
  std::vector<int> maskalpha,
  masktheta;
  
  arma::field<arma::uvec> node_id_period;
  
  arma::vec theta_par, thetaold,
  e_wm,
  alpha_gr, theta_gr,
  gamma,
  gamma_init;
  
  const arma::umat node_id_dyad;// node_id_dyad_ho; //matrix (column major)
  arma::umat par_ind;
  
  const arma::mat x_t, //matrix (column major)
  z_t,
  //z_t_ho,
  mu_b_t,
  var_b_t;
  
  arma::mat kappa_t,
  b_t,
  alpha_term,
  send_phi,
  rec_phi,
  e_wmn_t,
  e_c_t;
  
  arma::cube alpha, //3d array (column major)
  theta,
  beta, betaold,
  beta_init;
  
  //std::vector< Array<double> > new_e_c_t; //For reduce op.
  //arma::mat new_e_c_t;
  
  arma::cube::iterator beta_end;
  arma::vec::iterator theta_par_end;
  
  
  void computeAlpha(bool= false);
  void computeTheta(bool = false);
  double alphaLB(bool = false);
  static double alphaLBW(int, double*, void*);
  void alphaGr(int, double*);
  static void alphaGrW(int, double*, double*, void*);
  double thetaLB(bool = false, bool = false);
  static double thetaLBW(int, double*, void*);
  void thetaGr(int, double*);
  static void thetaGrW(int, double*, double*, void*);
  
  void updatePhiInternal(arma::uword, arma::uword,
                         double*,
                         double* ,
                         double*,
                         arma::uword* );
  
};

#endif //MMMODEL_CLASS