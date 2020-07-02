/* --------------------------------------------------------------------------------
#
#   Discretize continuous time output
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "discretise.hpp"

//' Discretize Output
//'
//' Modified from \code{smfsb} package
//'
//' @param out a matrix
//' @param out dt the discretization lattice
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix discretise(const Rcpp::NumericMatrix& out, const double dt){

  int ncol = out.ncol();
  int events = out.nrow();

  double end = out.at(events-1,0);
  double start = out.at(0,0);

  // int len = ((int)(end - start)) / dt + 1;
  int len = (int)floor((end - start) / dt) + 1;
  Rcpp::NumericMatrix x(len,ncol);

  double target{0.};
  int j{0};

  for(int i=0; i<events; i++){
    while(out.at(i,0) >= target){
      x.at(j,0) = target;
      for(int k=1; k<ncol; k++){
        x.at(j,k) = out.at(i,k);
      }
      j += 1;
      target += dt;
    }
  }

  Rcpp::CharacterVector outnames = Rcpp::colnames(out);
  Rcpp::colnames(x) = outnames;
  return x;
};
