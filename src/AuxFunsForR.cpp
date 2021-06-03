#include <RcppArmadillo.h>

//' @rdname auxfuns
// [[Rcpp::export(approxB)]]
Rcpp::NumericMatrix approxB(Rcpp::NumericVector y,
                            Rcpp::IntegerMatrix d_id,
                            Rcpp::NumericMatrix pi1_mat,
			    Rcpp::Nullable<Rcpp::NumericMatrix> pi2_mat_tmp = R_NilValue, 
                            bool directed = true)
{
  bool bipart;
  Rcpp::NumericMatrix pi2_mat;
  if(pi2_mat_tmp.isNotNull()){
    pi2_mat = Rcpp::as<Rcpp::NumericMatrix>(pi2_mat_tmp);
    bipart = true;
  } else {
    pi2_mat = pi1_mat;
    bipart = false;
  }
   int N_BLK1 = pi1_mat.nrow();
  int N_BLK2 = pi2_mat.nrow();
  int N_DYAD = d_id.nrow();
  Rcpp::NumericMatrix den(N_BLK2, N_BLK1), num(N_BLK2, N_BLK1), B_t(N_BLK2, N_BLK1);
  int s = 0, r = 0;
  double prob_temp;
  for(int d = 0; d < N_DYAD; ++d){
    s = d_id(d, 0);
    r = d_id(d, 1);
    for(int g = 0; g < N_BLK1; ++g){
      for(int h = 0; h < N_BLK2; ++h){
        if((g <= h) | directed | bipart){
          prob_temp = pi1_mat(g, s) * pi2_mat(h, r);
          num(h, g) += y[d] * prob_temp;
          den(h, g) += prob_temp;
        } else{
          num(h, g) = num(g, h);
          den(h, g) = den(g, h);
        }
      }
    }
  }
  std::transform(num.begin(), num.end(),
                 den.begin(),
                 B_t.begin(),
                 std::divides<double>());
  return B_t;
}

//' @rdname auxfuns
//[[Rcpp::export(getZ)]]
Rcpp::IntegerMatrix getZ(Rcpp::NumericMatrix pi_mat)
{
  int NROW = pi_mat.nrow();
  int NCOL = pi_mat.ncol();
  int mflag, bloc; 
  double u, acc;
  Rcpp::NumericVector cprob(NROW); 
  Rcpp::IntegerMatrix res(NROW, NCOL);
  for(int i = 0; i < NCOL; ++i){
    u = R::runif(0, 1);
    acc = 0.0;
    for(int j = 0; j < NROW; ++j){
      acc += pi_mat(j, i);
      cprob[j] = acc;
    }
    bloc = findInterval(&(cprob[0]), NROW, u, FALSE, FALSE, 0, &mflag);
    res(bloc, i) = 1;
  }
  return(res);
}

//' @rdname auxfuns
// [[Rcpp::export]]
arma::mat vcovGamma(const arma::mat& X,
                    const arma::vec& probs, //must be computed using offset given by blockmodel and estimated mm's
                    const arma::vec& pen) {
  arma::uword N_PRED = X.n_cols;
  arma::mat hess(N_PRED, N_PRED, arma::fill::zeros);
  arma::vec d_elem = probs % (1.0 - probs);
  for(arma::uword dc = 0; dc < N_PRED; ++dc){
    for(arma::uword dr = dc; dr < N_PRED; ++dr){
      hess(dr, dc) = -arma::dot(X.col(dr),(d_elem % X.col(dc)));
      hess(dc, dr) = hess(dr, dc);
    }
  }
  hess.diag() -= (1.0) / pen;
  return hess.i();
}

//' @rdname auxfuns
// [[Rcpp::export]]
arma::mat vcovBeta(const arma::mat& X,
                   const arma::mat& pi_mat,
                   const arma::mat& alpha,
                   const arma::vec& alpha_sum,
                   const arma::vec& pen) {
  arma::uword N_PRED = X.n_cols;
  arma::uword N_BLK = alpha.n_cols;
  arma::mat hess(N_PRED*N_BLK, N_PRED*N_BLK, arma::fill::zeros);
  arma::mat di_alpha = alpha, tri_alpha = alpha;
  di_alpha.for_each([](arma::mat::elem_type& val){val=R::digamma(val);});
  tri_alpha.for_each([](arma::mat::elem_type& val){val=R::trigamma(val);});
  arma::vec di_alpha_sum = alpha_sum, tri_alpha_sum = alpha_sum;
  di_alpha_sum.for_each([](arma::vec::elem_type& val){ val=R::digamma(val);});
  tri_alpha_sum.for_each([](arma::vec::elem_type& val){ val=R::trigamma(val);});
  arma::mat l_pi = log(pi_mat + 1e-6);
  arma::uword ind = 0;
  for(arma::uword blkc = 0; blkc < N_BLK; ++blkc){
    for(arma::uword blkr = 0; blkr < N_BLK; ++blkr){
      for(arma::uword dc = 0; dc < N_PRED; ++dc){
        for(arma::uword dr = 0; dr < N_PRED; ++dr){
          ind = (dr + N_PRED * blkr) + (dc + N_PRED * blkc) * (N_PRED * N_BLK);
          if(blkr == blkc){
            hess.at(ind) =  arma::accu((X.col(dr) % X.col(dc) % alpha.col(blkr)) %
              (l_pi.col(blkr) + di_alpha_sum - di_alpha.col(blkr)
                 + alpha.col(blkr) % (tri_alpha_sum - tri_alpha.col(blkr))));
          } else {
            hess.at(ind) = arma::dot((X.col(dr) % X.col(dc)), 
                    alpha.col(blkr) % alpha.col(blkc) % tri_alpha_sum);
          }
        }
      } 
    }
  }
  //hess = arma::symmatl(hess);
  hess.diag() -= (1.0) / pen;
  return (-hess).i();
}

//' @rdname auxfuns
// [[Rcpp::export()]]
double alphaLBound(arma::vec par,
               arma::uvec tot_nodes,
               arma::umat c_t, 
               arma::mat x_t,
               arma::umat s_mat,
               arma::uvec t_id,
               arma::cube var_beta,
               arma::cube mu_beta)
{
  arma::uword N_NODE = x_t.n_cols, N_BLK = c_t.n_rows,
    N_MONAD_PRED = x_t.n_rows,  N_STATE = s_mat.n_rows;
  double linpred = 0.0, row_sum = 0.0, res = 0.0, res_int = 0.0;
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword p = 0; p < N_NODE; ++p){
      row_sum = 0.0;
      res_int = 0.0;
      for(arma::uword g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
          linpred += x_t(x, p) * par[x + N_MONAD_PRED * (g + N_BLK * m)];
        }
        linpred = exp(linpred);
        row_sum += linpred;
        res_int += (lgamma(linpred + c_t(g, p)) - lgamma(linpred));
      }
      res_int += (lgamma(row_sum) - lgamma(row_sum + tot_nodes[p]));
      res += res_int * s_mat(m, t_id[p]);
    }
  }
  
  // Prior for beta
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
        res -= 0.5 * pow(par[x + N_MONAD_PRED * (g + N_BLK * m)] - mu_beta(x, g, m), 2.0) / var_beta(x, g, m);
      }
    } 
  }
  
  
  return -res;
}

//' @rdname auxfuns
// [[Rcpp::export()]]
arma::vec alphaGrad(arma::vec par,
               arma::uvec tot_nodes,
               arma::umat c_t, 
               arma::mat x_t,
               arma::umat s_mat,
               arma::uvec t_id,
               arma::cube var_beta,
               arma::cube mu_beta)
{
  arma::uword N_NODE = x_t.n_cols, N_BLK = c_t.n_rows,
    N_MONAD_PRED = x_t.n_rows,  N_STATE = s_mat.n_rows;
  double res=0.0, prior_gr=0.0, linpred = 0.0;
  arma::uword U_NPAR = par.n_elem;
  
  arma::vec gr(U_NPAR, arma::fill::zeros);
  
  arma::cube alpha (N_BLK, N_NODE, N_STATE, arma::fill::zeros);
  arma::mat alpha_row(N_NODE, N_STATE, arma::fill::zeros);
  
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword p = 0; p < N_NODE; ++p){
      for(arma::uword g = 0; g < N_BLK; ++g){
        linpred = 0.0;
        for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
          linpred += x_t(x, p) * par[x + N_MONAD_PRED * (g + N_BLK * m)];
        }
        alpha(g, p, m) = exp(linpred);
        alpha_row(p, m) += alpha(g, p, m);
      }
    }
  }
  
  
  for(arma::uword m = 0; m < N_STATE; ++m){
    for(arma::uword g = 0; g < N_BLK; ++g){
      for(arma::uword x = 0; x < N_MONAD_PRED; ++x){
        res = 0.0;
#pragma omp parallel for reduction(+:res)
        for(arma::uword p = 0; p < N_NODE; ++p){
            res += (R::digamma(alpha_row(p, m)) - R::digamma(alpha_row(p,m) + tot_nodes[p])
                      + R::digamma(alpha(g, p, m) + c_t(g, p)) - R::digamma(alpha(g, p, m)))
              * s_mat(m, t_id[p]) * alpha(g, p, m) * x_t(x, p);
          }
        prior_gr = (par[x + N_MONAD_PRED * (g + N_BLK * m)] - mu_beta(x, g, m)) / var_beta(x, g, m);
        gr[x + N_MONAD_PRED * (g + N_BLK * m)] = -(res - prior_gr);
      }
    }
  }
  return(gr);
}


//' @rdname auxfuns
// [[Rcpp::export(vertboot_matrix_rcpp2)]]
Rcpp::List vertboot_matrix_rcpp2(Rcpp::IntegerMatrix m1,
				 Rcpp::IntegerVector blist1,
				 Rcpp::IntegerVector blist2){
  int a;
  int b;
  int num1=m1.nrow();
  int num2=m1.ncol();
  Rcpp::IntegerMatrix x1(num1,num2);
  for(int k=0;k<num1;k++){
    for(int q=0;q<num2;q++){
      x1(k,q)=0;
    }
  }
  for(int i=0;i<num1;i++){
    a=blist1[i];
    for(int j=0;j<num2;j++){
      b=blist2[j];
      //if(a!=b){
        x1(i,j)=m1(a,b);
      // }
      // else{
      //   a=round(R::runif(-0.49,num-0.5));
      //   b=round(R::runif(-0.49,num-0.5));
      //   while(a==b){
      //     b=round(R::runif(-0.49,num-0.5));
      //   }
      //   x1(i,j)=m1(a,b);
      //   x1(j,i)=m1(a,b);
      // }
    }
  }
  Rcpp::List res;
  res["x"] = x1;
  res["index1"] = blist1;
  res["index2"] = blist2;
  return (res);
}

