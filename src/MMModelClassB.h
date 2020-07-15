#ifndef MMMODEL_CLASS
#define MMMODEL_CLASS

#include <vector>

// #ifndef DEBUG_MODE
// #define DEBUG_MODE
// #define ARMA_EXTRA_DEBUG
// #endif 

#include <RcppArmadillo.h>
#include "AuxFuns.h"



class MMModelB
{
public:
  MMModelB(const arma::mat& z_t,
            const arma::mat& X1_t, // Add second family data frame X
            const arma::mat& X2_t,
            const arma::vec& y,
            const arma::uvec& time_id_dyad,
            const arma::uvec& time_id_node1_r, //should include for both family 1 and family 2
            const arma::uvec& time_id_node2_r, 
            const arma::uvec& nodes_per_period, //stacked family 1 nodes per period, family 2 nodes per period
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
            const arma::mat& phi_init1, // two sets, family 1, family 2
            const arma::mat& phi_init2,
            arma::mat& kappa_init_t,
            arma::mat& b_init_t,
            arma::cube& beta1_init_r,//two sets, family 1, family 2
            arma::cube& beta2_init_r,
            arma::vec& gamma_init_r,
            arma::cube& theta_init_r,
            const arma::vec& phi_order_r,
            Rcpp::List& control
            
            
  );
  ~MMModelB();
  
  void sampleDyads(arma::uword iter);
  void updatePhi();
  void updateKappa();
  void optim_ours(bool);
  double cLB();
  
  
  arma::mat getC(bool, arma::uword, arma::uword);
  void getC(arma::mat& res);
  arma::mat getPhi(bool);
  //arma::uvec getN();
  arma::mat getWmn();
  arma::mat getKappa();
  arma::mat getB();
  void getB(arma::mat&);
  arma::vec getGamma();
  void getGamma(arma::vec&);
  arma::cube getBeta1();
  void getBeta1(arma::cube&);
  arma::cube getBeta2();
  void getBeta2(arma::cube&);
  arma::cube getAlpha1();
  arma::cube getAlpha2();
  // bool checkConv(char p,
  //                arma::vec& old,
  //                double tol);
  // bool maxDiffCheck(Array<double>& current,
  //                   arma::vec& old,
  //                   double tol);
  
  
private:
  
  const arma::uword BIPART, N_NODE1, N_NODE2,
  N_DYAD,
  N_BLK1, N_BLK2,
  N_STATE,
  N_TIME,
  N_MONAD_PRED1, N_MONAD_PRED2,
  N_DYAD_PRED,
  N_B_PAR,
  OPT_ITER,
  N_NODE_BATCH1,
  N_NODE_BATCH2;
  
  const double eta,
  forget_rate,
  delay;
  
  const arma::vec var_gamma,
  mu_gamma;
  const arma::cube var_beta1,
  mu_beta1, var_beta2, mu_beta2;
  
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
  
  const arma::vec y;
  
  const arma::uvec time_id_dyad,
  time_id_node1, time_id_node2,
  n_nodes_time1, n_nodes_time2;
  
  arma::uvec tot_nodes1, tot_nodes2,
  node_in_batch1,node_in_batch2,
  dyad_in_batch,
  node_batch1,node_batch2; 
  std::vector<int> maskalpha1, maskalpha2, 
  masktheta;
  
  arma::vec phi_order, theta_par, thetaold,
  e_wm,
  gamma,
  gamma_init;
  
  const arma::umat node_id_dyad; //matrix (column major)
  arma::umat par_ind;
  
  const arma::mat x1_t, x2_t, //matrix (column major)
  z_t,
  mu_b_t,
  var_b_t;
  
  arma::mat kappa_t,
  b_t,
  alpha_term1, alpha_term2,
  send_phi,
  rec_phi,
  e_wmn_t,
  e_c_t1, e_c_t2;
  
  arma::cube alpha1,alpha2, //3d array (column major)
  theta,
  beta1, beta2, beta1old, beta2old,
  beta_init1, beta_init2;
  
  //std::vector< Array<double> > new_e_c_t; //For reduce op.
  arma::mat new_e_c_t1, new_e_c_t2;
  
  
  void computeAlpha(const arma::uword,
                    const arma::uword,
                    const arma::uword,
                    bool,
                    const arma::mat&,
                    arma::cube&,
                    arma::cube&,
                    arma::mat&,
                    arma::mat&,
                    const arma::uvec,
                    arma::uvec,
                    arma::uvec,
                    const arma::uword);
  
  void computeTheta(bool);
  double alphaLBInternal(const arma::uword,
                 const arma::uword,
                 const arma::uword,
                 bool,
                 const arma::mat,
                 arma::cube&,
                 const arma::cube,
                 const arma::cube,
                 arma::cube&,
                 arma::mat&,
                 arma::mat&,
                 const arma::uvec,
                 arma::uvec,
                 arma::uvec,
                 const arma::uword);
  double alphaLB(bool, bool);
  static double alphaLBWMode1(int, double*, void*);
  static double alphaLBWMode2(int, double*, void*);
  void alphaGr(bool, int, double*);
  void alphaGrInternal(  int, double*, bool, const arma::uword, 
                         const arma::uword, const arma::uword, 
                         const arma::mat, arma::cube&, const arma::cube, const arma::cube, arma::cube&,
                         arma::mat&, arma::mat&, const arma::uvec, arma::uvec,
                         arma::uvec, const arma::uword);
  static void alphaGrWMode1(int, double*, double*, void*);
  static void alphaGrWMode2(int, double*, double*, void*);
  double thetaLB(bool, bool);
  static double thetaLBW(int, double*, void*);
  void thetaGr(int, double*);
  static void thetaGrW(int, double*, double*, void*);
  
  void updatePhiInternal(arma::uword,
                         arma::uword,
                         double*,
                         double*,
                         double*,
                         arma::uword*,
                         arma::cube&,
                         arma::cube&,
                         arma::mat&,
                         const arma::uvec&,
                         const arma::vec&,
                         arma::mat&,
                         const arma::umat&,
                         arma::uword,
                         arma::uword,
                         arma::uword
                         );
  
  // Add to check if something off with family 2 phi updates: updatePhiInternal1 and updatePhiInternal2
  void updatePhiInternal1(arma::uword,
                         arma::uword,
                         double*,
                         double*,
                         double*,
                         arma::uword*,
                         arma::cube&,
                         arma::cube&,
                         arma::mat&,
                         const arma::uvec&,
                         const arma::vec&,
                         arma::mat&,
                         const arma::umat&,
                         arma::uword,
                         arma::uword,
                         arma::uword
  );
  void updatePhiInternal2(arma::uword,
                          arma::uword,
                          double*,
                          double*,
                          double*,
                          arma::uword*,
                          arma::cube&,
                          arma::cube&,
                          arma::mat&,
                          const arma::uvec&,
                          const arma::vec&,
                          arma::mat&,
                          const arma::umat&,
                          arma::uword,
                          arma::uword,
                          arma::uword
  );
  
};

#endif //MMMODEL_CLASS
