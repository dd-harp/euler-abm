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
std::vector<double> inhomPP_piecewiseconst(const Rcpp::NumericVector& tvec, const Rcpp::NumericVector& lambdavec, const double tmax, const bool first);

/* sample first event time via rejection algorithm */
Rcpp::NumericVector inhomPP_piecewiseconst_reject(const Rcpp::NumericVector& tvec, const Rcpp::NumericVector& lambdavec);

#endif
