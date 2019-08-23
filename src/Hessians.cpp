#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

inline double ftrigamma(double x)	
{	
  double p;	
  int i;	
  
  x+=6;	
  p=1.0/(x*x);	
  p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)	
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;	
  for (i=0; i<6 ;i++)	
  {	
    x-=1;	
    p+=1.0/(x*x);	
  }	
  return p;	
}	

inline double fdigamma(double x)
{
  double p;
  x+=6;
  p=1.0/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
    0.008333333333333)*p-0.083333333333333)*p;
  p+=log(x)-0.5/x-1./(x-1)-1./(x-2)-1./(x-3)-1./(x-4)-1./(x-5)-1./(x-6);
  return p;
}

//[[Rcpp::export(.calcHessBeta)]]
Rcpp::NumericMatrix calcHessBeta(Rcpp::NumericMatrix A, // alpha, nodes by blocks. Assumes subset to nodes in a time sampled to be in state m
                                 Rcpp::NumericMatrix C, // C, nodes by blocks
                                 Rcpp::NumericMatrix X, // X, nodes by preds
                                 Rcpp::NumericVector Xi, // sum_alpha, nodes
                                 Rcpp::NumericVector N, // tot_nodes interacted with, nodes 
                                 Rcpp::NumericMatrix beta, //beta coefs, preds by blocks 
                                 double var_beta)
{
  int NNODE = X.nrow();
  int NPRED = X.ncol();
  int NBLK = A.ncol();
  int NCELL = pow(NPRED * NBLK, 2.0);
  Rcpp::NumericMatrix hess(NPRED * NBLK, NPRED * NBLK);
  int k_row = 0, k_col = 0, j_row = 0, j_col = 0;
  double res = 0.0;
  for(int c = 0; c < NCELL; ++c, ++j_row){
    if(j_row == NPRED){
      j_row = 0;
      k_row++;
    }
    if(k_row == NBLK){
      k_row = 0;
      j_col++;
    }
    if(k_col == NBLK){
      k_col = 0;
    }
    if(j_col == NPRED){
      j_col = 0;
      k_col++;
    }
    if(((k_row == k_col) && (j_col > j_row)) || (k_col > k_row)){
      continue;
    }
    if(k_row == k_col){
      res = 0.0;
      for(int i = 0; i < NNODE; ++i){
        res -= A(i,k_row) * X(i, j_row) * X(i, j_col)
        * (fdigamma(A(i,k_row)) - fdigamma(A(i,k_row) + C(i,k_row))
             - fdigamma(Xi[i]) + fdigamma(Xi[i] + N[i]) +
             A(i, k_row) * (ftrigamma(A(i,k_row)) - ftrigamma(A(i,k_row) + C(i,k_row))
                              - ftrigamma(Xi[i]) + ftrigamma(Xi[i] + N[i])));
      }
      hess[c] = res;
      if(j_row == j_col){
        hess[c] -= beta(j_row, k_row) / var_beta;
      }
    } else {
      res = 0.0;
      for(int i = 0; i < NNODE; ++i) {
        res += A(i, k_row) * A(i, k_col) * X(i, j_row) * X(i, j_col)
        * (ftrigamma(Xi[i]) - ftrigamma(Xi[i] + N[i]));
      }
      hess[c] = res;
    }
  }
  for(int col = 0; col < NPRED * NBLK; ++col){
    for(int row = 0; row < col; ++row){
      hess(row, col) = hess(col, row);
    }
  }
  return(hess);
}

//[[Rcpp::export(.calcHessTheta)]]
Rcpp::NumericMatrix calcHessTheta(Rcpp::IntegerMatrix z,
                                  Rcpp::IntegerMatrix w,
                                  Rcpp::NumericVector theta,
                                  Rcpp::NumericMatrix d,
                                  Rcpp::NumericVector gamma,
                                  Rcpp::NumericMatrix B,
                                  double var_gamma,
                                  Rcpp::NumericVector var_b_vec)
{
  int NBLK2 = B.size();
  int NBLK = B.nrow();
  int NGAMMA = gamma.size();
  int NDYAD = d.nrow();
  int NPAR = NBLK2 + NGAMMA;
  double prob = 0.0, var_b = 1.0, res = 0.0; 
  int blk1_b = 0, blk2_b = 0, blk1_g = 0, blk2_g = 0;
  Rcpp::NumericMatrix hess(NPAR, NPAR);
  
  for(int col = 0; col < NPAR; ++col){
    for(int row = 0; row < NPAR; ++row){
      if(col > row){
        hess(row, col) = hess(col, row);
        continue;
      }
      if((row < NBLK2) && (col < NBLK2)){
        if(row == col){
          res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif          
          for(int i = 0; i < NDYAD; ++i){
            prob = theta[blk1_b + NBLK * (blk2_b + NBLK * i)]; 
            res += (1.0 - prob) * (prob) * z(i, blk1_b) * w(i, blk2_b);
          }
          hess(row, col) -= res;
          var_b = blk1_b == blk2_b ? var_b_vec[0] : var_b_vec[1];
          hess(row, col) -= B(blk1_b, blk2_b) / var_b;
          blk1_b++;
          if(blk1_b == NBLK){
            blk1_b = 0;
            blk2_b++;
          }
        }
      } else if((row >= NBLK2) && (col < NBLK2)){
        res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif          
        for(int i = 0; i < NDYAD; ++i){
          prob = theta[blk1_g + NBLK * (blk2_g + NBLK * i)];
          res += d(i, (row - NBLK2)) * (1.0 - prob) * prob * z(i, blk1_g) * w(i, blk2_g);
        }
        hess(row, col) -= res;
        if((row + 1) == NPAR){
          blk1_g++;
        }
        if(blk1_g == NBLK){
          blk1_g = 0;
          blk2_g++;
        }
      } else {
        res = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:res)
#endif          
        for(int i = 0; i < NDYAD; ++i){
          for(int g = 0; g < NBLK; ++g){
            for(int h = 0; h < NBLK; ++h){
              prob = theta[g + NBLK * (h + NBLK * i)];
              res += d(i, (row - NBLK2)) * d(i, (col - NBLK2)) * (1.0 - prob) * prob * z(i, g) * w(i, h);
            }
          }
        }
        hess(row, col) -= res;
        if(row == col){
          hess(row, col) -= gamma[row - NBLK2] / var_gamma;
        }
      }
      
    }
  }
  return(hess);
} 
