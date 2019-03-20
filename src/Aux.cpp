#include "Aux.hpp"


// [[Rcpp::export]]
Rcpp::NumericMatrix approxB(Rcpp::NumericVector y,
                            Rcpp::IntegerMatrix d_id,
                            Rcpp::NumericMatrix pi_mat)
{
  int N_BLK = pi_mat.nrow();
  int N_DYAD = d_id.nrow();
  Rcpp::NumericMatrix den(N_BLK, N_BLK), num(N_BLK, N_BLK), B_t(N_BLK, N_BLK);
  int s, r;
  double prob_temp;
  for(int d = 0; d < N_DYAD; ++d){
    s = d_id(d, 0);
    r = d_id(d, 1);
    for(int g = 0; g < N_BLK; ++g){
      for(int h = 0; h < N_BLK; ++h){
        prob_temp = pi_mat(g, s) * pi_mat(h, r);
        num(h, g) += y[d] * prob_temp;
        den(h, g) += prob_temp;
      }
    }
  }
  std::transform(num.begin(), num.end(),
                 den.begin(),
                 B_t.begin(),
                 std::divides<double>());
  return B_t;
}

//[[Rcpp::export]]
Rcpp::IntegerMatrix getZ(Rcpp::NumericMatrix pmat)
{
  int NROW = pmat.nrow();
  int NCOL = pmat.ncol();
  int mflag, bloc;
  double u, acc;
  Rcpp::NumericVector cprob(NROW); 
  Rcpp::IntegerMatrix res(NROW, NCOL);
  for(int i = 0; i < NCOL; ++i){
    u = R::runif(0, 1);
    acc = 0.0;
    for(int j = 0; j < NROW; ++j){
      acc += pmat(j, i);
      cprob[j] = acc;
    }
    bloc = findInterval(&(cprob[0]), NROW, u, FALSE, FALSE, 0, &mflag);
    res(bloc, i) = 1;
  }
  return(res);
}

// Coming from Abramowitz and Stegun 6.4.13 and 6.4.6
inline double tetragamma(double x)
{
  x += 6;
  double p = 1.0/(x*x);
  p = p*(-1-p*(0.5+p*(0.16666666666667
		      - p*(0.16666666666667
			   + p*(0.3-p*(0.83333333333333
				       + p*3.2904761904762))))))-1.0/(x*x*x);
  for(int i = 0; i<6; ++i)
    {
      x -= 1;
      p -= 2.0/(x*x*x);
    }
  return p;
}

// Coming from Blei's lda implementation in C (github.com/blei-lab/lda-c)
inline double trigamma(double x)
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

// inline double digamma(double x)
// {
//   double p;
//   x+=6;
//   p=1.0/(x*x);
//   p=(((0.004166666666667*p-0.003968253986254)*p+
//       0.008333333333333)*p-0.083333333333333)*p;
//   p+=log(x)-0.5/x-1./(x-1)-1./(x-2)-1./(x-3)-1./(x-4)-1./(x-5)-1./(x-6);
//   return p;
// }

double logSumExp(const std::vector<double>& invec)
{
  double offset = *std::max_element(invec.begin(), invec.end());
  double res = std::accumulate(invec.begin(), invec.end(),
			       0.0,
			       [&offset](double a, const double b){
				 return a + exp(b - offset);
			       });
  return offset + log(res);
}

/*
// Adaptation of vmmin in optim.c to
// enable use in threaded application
 */
void vmmin_ours(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
      int maxit, int trace, int *mask,
      double abstol, double reltol, int nREPORT, void *ex,
      int *fncount, int *grcount, int *fail)
{
    bool accpoint, enough;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n;

    if (maxit <= 0) {
	*fail = 0;
	*Fmin = fminfn(n0, b, ex);
	*fncount = *grcount = 0;
	return;
    }

    std::vector<int> l(n0, 1);
    n = 0;
    for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
    std::vector<double> g(n0, 0.0);
    std::vector<double> t(n, 0.0);
    std::vector<double> X(n, 0.0);
    std::vector<double> c(n, 0.0);
    std::vector< std::vector<double> > B(n, std::vector<double>(n));
    f = fminfn(n0, b, ex);
    if (!R_FINITE(f))
      Rcpp::stop("initial value in 'vmmin' is not finite");
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(n0, b, &g[0], ex);
    iter++;
    ilast = gradcount;
    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = FALSE;
	    do {
		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (10.0 + X[i] == 10.0 + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
		    f = fminfn(n0, b, ex);
		    funcount++;
		    accpoint = R_FINITE(f) &&
			(f <= *Fmin + gradproj * steplength * 0.0001);
		    if (!accpoint) {
			steplength *= 0.2;
		    }
		}
	    } while (!(count == n || accpoint));
	    enough = (f > abstol) &&
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }
	    if (count < n) {/* making progress */
		*Fmin = f;
		fmingr(n0, b, &g[0], ex);
		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
}
