#include "AuxFuns.h"

double logSumExp(const arma::vec& invec)
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
 // enable use in threaded call
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
  if (!R_FINITE(f)){
    *fail = 1;
    return;
  }
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
    if (iter >= maxit) break;
    if (gradcount - ilast > 2 * n)
      ilast = gradcount;	/* periodic restart */
  } while (count != n || ilast != gradcount);
  *fail = (iter < maxit) ? 0 : 1;
  *fncount = funcount;
  *grcount = gradcount;
}
