#include "Aux.hpp"
#include "MMModelClass.hpp"


double digamma_approx(double x)
{
  double p;
  x=x+6;
  p=1/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
    0.008333333333333)*p-0.083333333333333)*p;
  p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
  return p;
}

double lgamma_approx(double x)
{
  double z=1/(x*x);
  
  x=x+6;
  z=(((-0.000595238095238*z+0.000793650793651)
        *z-0.002777777777778)*z+0.083333333333333)/x;
  z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
  log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
  return z;
}

double lgammaDiff(double alpha, double C) {
  return lgamma(alpha + C) - lgamma(alpha);
}
double digammaDiff(double alpha, double C) {
  return R::digamma(alpha + C) - R::digamma(alpha);
}


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
    
    if (gradproj < 0.0) {	// search direction is downhill 
      steplength = 1.0;
      accpoint = FALSE;
      do {
        count = 0;
        for (i = 0; i < n; i++) {
          b[l[i]] = X[i] + steplength * t[i];
          if (10.0 + X[i] == 10.0 + b[l[i]]) // no change
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
      // stop if value if small or if relative change is low 
      if (!enough) {
        count = n;
        *Fmin = f;
      }
      if (count < n) {
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
        } else {	// D1 < 0 
          ilast = gradcount;
        }
      } else {	// no progress
        if (ilast < gradcount) {
          count = 0;
          ilast = gradcount;
        }
      }
    } else {		// uphill search 
      count = 0;
      if (ilast == gradcount) count = n;
      else ilast = gradcount;
      // Resets unless has just been reset 
    }
    if (trace && (iter % nREPORT == 0))
      Rprintf("iter%4d value %f\n", iter, f);
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n)
      ilast = gradcount;	// periodic restart 
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
