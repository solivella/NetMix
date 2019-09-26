#ifndef AUX_HPP
#define AUX_HPP

#include <vector>
#include <initializer_list>
#include <functional>
#include <numeric>
#include <RcppArmadillo.h>



double logSumExp(const arma::vec& invec);

typedef double optimfn(int, double*, void*);
typedef void optimgr(int, double*, double*, void*);

void vmmin_ours(int,
		double*,
		double*,
		optimfn,
		optimgr,
		int,
		int,
		int*,
		double,
		double,
		int,
		void*,
		int*,
		int*,
		int*);



#endif // AUX_HPP
