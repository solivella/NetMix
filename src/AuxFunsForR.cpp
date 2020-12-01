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
// [[Rcpp::export(alphaLB)]]
double alphaLB(arma::vec par,
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

