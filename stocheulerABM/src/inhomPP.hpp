/* --------------------------------------------------------------------------------
#
#   Sample from inhomogeneous Poisson Process (piecewise constant)
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef INHOMPP_HPP
#define INHOMPP_HPP

#include <vector>

#include <Rcpp.h>

/* sample piecewise constant PP */
std::vector<double> inhomPP_piecewiseconst(const Rcpp::NumericVector& tj, const Rcpp::NumericVector& lambdaj, const double tmax, const bool first);

#endif
