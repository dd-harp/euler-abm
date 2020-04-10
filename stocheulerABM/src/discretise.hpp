/* --------------------------------------------------------------------------------
#
#   Discretize continuous time output
#   Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef DISCRETISE_HPP
#define DISCRETISE_HPP

#include <Rcpp.h>

Rcpp::NumericMatrix discretise(const Rcpp::NumericMatrix& out, const int dt);


#endif
