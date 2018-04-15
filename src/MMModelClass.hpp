#ifndef MMMODEL_CLASS
#define MMMODEL_CLASS

#include <vector>
#include <Rcpp.h>
#include "Aux.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


class MMModel
{
public:
  MMModel(const Rcpp::NumericMatrix& z_t,
	  const Rcpp::NumericMatrix& x_t,
	  const Rcpp::IntegerVector& y,
	  const int N_THREAD,
	  const Rcpp::IntegerVector& time_id_dyad,
	  const Rcpp::IntegerVector& time_id_node,
	  const Rcpp::IntegerVector& nodes_per_period,
	  const Rcpp::IntegerMatrix& node_id_dyad,
	  const Rcpp::NumericMatrix& mu_b,
	  const Rcpp::NumericMatrix& var_b,
	  const Rcpp::NumericMatrix& phi_init,
	  Rcpp::NumericMatrix& kappa_init_t,
	  Rcpp::NumericMatrix& b_init_t,
	  Rcpp::NumericVector& beta_init,
	  Rcpp::NumericVector& gamma_init,
	  Rcpp::List& control
	  );


  void updatePhi();
  void updateKappa();
  void optim(bool);
  double cLL();


  Rcpp::NumericMatrix getC();
  Rcpp::NumericMatrix getPhi(bool);
  Rcpp::NumericMatrix getWmn();
  Rcpp::NumericMatrix getKappa();
  Rcpp::NumericMatrix getB();
  void getB(Rcpp::NumericVector&);
  Rcpp::NumericVector getGamma();
  void getGamma(Rcpp::NumericVector&);
  Rcpp::List getBeta();
  void getBeta(Rcpp::NumericVector&);
  int checkConvChng(Rcpp::NumericVector::iterator,
		    Rcpp::NumericVector::iterator, 
		    int,
		    double);


private:

  const int N_NODE,
    N_DYAD,
    N_BLK,
    N_STATE,
    N_TIME,
    N_MONAD_PRED,
    N_DYAD_PRED,
    N_B_PAR,
    OPT_ITER,
    N_THREAD;

  const double eta,
    var_gamma,
    var_beta;

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

  const Array<int> y, //vectors
    time_id_dyad,
    time_id_node,
    n_nodes_time;

  std::vector<int> maskalpha, //vectors
    masktheta;

  Array<double> theta_par,
    e_wm,
    gamma;

  const Array<int> node_id_dyad; //matrix (column major)
  Array<int> par_ind;

  const Array<double> x_t, //matrix (column major)
    z_t,
    mu_b_t,
    var_b_t;

  Array<double> kappa_t,
    b_t,
    alpha_term1,
    alpha_term2,
    send_phi,
    rec_phi,
    e_wmn_t,
    e_c_t;

  Array<double> alpha, //3d array (column major)
    theta,
    beta;

  std::vector< Array<double> > new_e_c_t; //For reduce op.


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

  void updatePhiInternal(int,
			 int,
			 double*,
			 double*,
			 double*);

};

#endif //MMMODEL_CLASS
