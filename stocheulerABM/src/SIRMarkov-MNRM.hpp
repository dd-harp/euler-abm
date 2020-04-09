/* --------------------------------------------------------------------------------
#
#   Markovian SIR model via Anderson's MNRM
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef SIRMARKOV_MNRM
#define SIRMARKOV_MNRM

#include <vector>
#include <array>

#include <Rcpp.h>

Rcpp::NumericMatrix SIRMarkov_MNRM(const double tmax, const int S, const int I, const int R, const double beta, const double gamma, const bool verbose);


#endif
