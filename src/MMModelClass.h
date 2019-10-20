#ifndef MMMODEL_CLASS
#define MMMODEL_CLASS

#include <vector>

// #ifndef DEBUG_MODE
// #define DEBUG_MODE
// #define ARMA_EXTRA_DEBUG
// #endif 

#include <RcppArmadillo.h>
#include "AuxFuns.h"



class MMModel
{
public:
  MMModel(const arma::mat& z_t,
          const arma::mat& x_t,
          const arma::vec& y,
          const arma::uvec& time_id_dyad,
          const arma::uvec& time_id_node,
          const arma::uvec& nodes_per_period,
          const arma::umat& node_id_dyad,
          const arma::mat& mu_b,
          const arma::mat& var_b,
          const arma::mat& phi_init,
          arma::mat& kappa_init_t,
          arma::mat& b_init_t,
          arma::cube& beta_init_r,
          arma::vec& gamma_init_r,
          Rcpp::List& control
  );
  ~MMModel();
  
  
  void updatePhi();
  void updateKappa();
  void optim_ours(bool);
  double cLL();
  
  
  arma::mat getC();
  void getC(arma::mat& res);
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
  // bool checkConv(char p,
  //                arma::vec& old,
  //                double tol);
  // bool maxDiffCheck(Array<double>& current,
  //                   arma::vec& old,
  //                   double tol);
  
  
private:
  
  const arma::uword N_NODE,
  N_DYAD,
  N_BLK,
  N_STATE,
  N_TIME,
  N_MONAD_PRED,
  N_DYAD_PRED,
  N_B_PAR,
  OPT_ITER;
  //N_THREAD;
  
  const double eta;
  const arma::vec var_gamma,
  mu_gamma;
  const arma::cube var_beta,
  mu_beta;
  
  double fminAlpha,
  fminTheta;
  
  
  int fncountAlpha,
  fncountTheta,
  grcountAlpha,
  grcountTheta,
  m_failAlpha,
  m_failTheta;
  
  
  bool verbose,
  directed;
  
  const arma::vec y;
  
  const arma::uvec time_id_dyad,
  time_id_node,
  n_nodes_time;
  
  arma::uvec tot_nodes; 
  std::vector<int> maskalpha, 
  masktheta;
  
  arma::vec theta_par,
  e_wm,
  gamma,
  gamma_init;
  
  const arma::umat node_id_dyad; //matrix (column major)
  arma::umat par_ind;
  
  const arma::mat x_t, //matrix (column major)
  z_t,
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
  beta,
  beta_init;
  
  //std::vector< Array<double> > new_e_c_t; //For reduce op.
  arma::mat new_e_c_t;
  
  
  void computeAlpha();
  
  void computeTheta();
  double alphaLB();
  static double alphaLBW(int, double*, void*);
  void alphaGr(int, double*);
  static void alphaGrW(int, double*, double*, void*);
  double thetaLB(bool);
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
