/* --------------------------------------------------------------------------------
#
#   Sample from inhomogeneous Poisson Process (piecewise constant)
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "inhomPP.hpp"

/* sample piecewise constant PP */

//' Sample an Inhomogeneous Poisson Process with Piecewise Constant Intensity
//'
//' Sample an inhomogeneous Poisson process withe piecewise constant intensity using an algorithm presented in:
//'   * Leemis, Lawrence M., and Stephen Keith Park. Discrete-event simulation: A first course. Upper Saddle River, NJ: Pearson Prentice Hall, 2006.
//'
//' @param tvec times at which the intensity changes (in j=0,1,...,k)
//' @param lambdavec values at the left endpoints of intensity (in j=0,1,...,k)
//' @param tmax maximum simulation time (usually equal to value of \code{tvec} at k)
//' @param first if \code{TRUE}, stop early and return the first event in the process
//'
//' @return a vector of event times in the counting process
//'
//'
//' @export
// [[Rcpp::export]]
std::vector<double> inhomPP_piecewiseconst(const Rcpp::NumericVector& tvec, const Rcpp::NumericVector& lambdavec, const double tmax, const bool first){

  // precalculate the cumulative intensity and inverse fn
  int k = tvec.size() - 1;
  if(k+1 != lambdavec.size()){
    Rcpp::stop("input vector 'tvec' and 'lambdavec' must be the same length");
  }

  if(tmax < tvec[k]){
    Rcpp::Rcout << " --- warning, tmax is less than the end of 'tvec'; please make sure this is what you want --- \n";
  }

  std::vector<double> Lambda(k+1,0.);
  for(int j=1; j<k+1; j++){
    Lambda[j] = Lambda[j-1] + 0.5*(lambdavec[j] + lambdavec[j-1])*(tvec[j] - tvec[j-1]);
  }

  // output and storage
  std::vector<double> a;
  a.push_back(0.);
  double u{0.};
  int n{0};
  int j{0};

  while(a[n] < tmax){

    u += R::rexp(1.);
    while( (Lambda[j+1] < u) && (j < k) ){
      j++;
    }

    double a_n1 = tvec[j] + ((u - Lambda[j]) / lambdavec[j]);
    a.push_back(a_n1); /* Lambda[j] <= u < Lambda[j+1] */
    n++;
    // only want the first event of the process
    if(first && n==1){
      a.erase(a.begin());
      return a;
    }

  }

  a.erase(a.begin());
  if(a.size()>1){
    a.erase(a.end()-1);
  }
  return a;
};



/* sample first event time via rejection algorithm */

//' Sample First Event Time Inhomogeneous Poisson Process with Piecewise Constant Intensity
//'
//' Sample an inhomogeneous Poisson process withe piecewise constant intensity using an accept-reject algorithm.
//' If the first event occurs after the time that the intensity is defined via \code{lambdavec}, return a \code{NaN}
//'
//' @param tvec times at which the intensity changes (in j=0,1,...,k)
//' @param lambdavec values at the left endpoints of intensity (in j=0,1,...,k)
//'
//' @return a vector of number of rejections and the time of first event
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector inhomPP_piecewiseconst_reject(const Rcpp::NumericVector& tvec, const Rcpp::NumericVector& lambdavec){

  int k = tvec.size();

  double tnow = tvec[0];
  double tnext{0.};
  double putative{0.};
  int nreject{0};

  for(int j=1; j<k; j++){
    tnext = tvec[j];
    putative = tnow + R::rexp(1./lambdavec[j-1]);
    if(putative < tnext){
      return Rcpp::NumericVector::create(
        Rcpp::Named("nreject") = nreject,
        Rcpp::Named("time") = putative
      );
    }
    nreject++;
    tnow = tnext;
  }

  return Rcpp::NumericVector::create(
    Rcpp::Named("nreject") = nreject,
    Rcpp::Named("time") = R_NaN
  );
};
