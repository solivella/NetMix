#ifndef AUX_HPP
#define AUX_HPP

#include <vector>
#include <initializer_list>
#include <functional>
#include <numeric>
#include <Rcpp.h>
#include <R_ext/Utils.h>


template<typename T>
class Array //traverse in order of indeces (i.e. i fastest)
{
public:
  template <typename Source>
  Array(std::initializer_list<int> dim, const Source& source)
  : dims(dim),
    data(source.begin(), source.end())
  {
  }
  Array(std::initializer_list<int> dim, T val)
    : dims(dim),
      data(std::accumulate(dims.begin(),dims.end(), 1, std::multiplies<int>()), val)
  {
  }
  //typedef T* iterator;
typename std::vector<T>::iterator begin(){
    return data.begin();
  }
typename std::vector<T>::iterator end(){
    return data.end();
  }
  //1d
  T& operator[](int i){
    return (data[i]);
  }
const T& operator[](int i) const {
    return (data[i]);
  }
  //2d
  T& operator()(int i, int j){//j varies most slowly.
    return (data[i + dims[0] * j]);
  }
  const T& operator()(int i, int j) const { 
    return (data[i + dims[0] * j]);
  }
  //3d
  T& operator()(int i, int j, int k){
    return (data[i + dims[0] * (j + dims[1] * k)]);
  }
  const T& operator()(int i, int j, int k) const {
    return (data[i + dims[0] * (j + dims[1] * k)]);
  }

private:
  std::vector<int> dims;
  std::vector<T> data;
};

Rcpp::NumericMatrix approxBdyad(Rcpp::NumericVector y,
                                Rcpp::IntegerMatrix node_id,
                                Rcpp::NumericMatrix pi_mat,
                                bool directed);
  
Rcpp::NumericMatrix approxB(Rcpp::NumericMatrix y,
                            Rcpp::NumericMatrix pi_mat,
                            bool directed);


double digamma_approx(double x);
double lgamma_approx(double x);
double lgammaDiff(double alpha, double C);
double digammaDiff(double alpha, double C);

double logSumExp(const std::vector<double>& invec);

template <typename T> int sgn(T val);

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
