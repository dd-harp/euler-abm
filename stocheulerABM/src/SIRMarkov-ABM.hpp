/* --------------------------------------------------------------------------------
#
#   Markovian SIR model via ABM
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef SIRMARKOV_ABM
#define SIRMARKOV_ABM

#include <vector>
#include <memory>
#include <string>

#include <Rcpp.h>


/* --------------------------------------------------------------------------------
#   human
-------------------------------------------------------------------------------- */

typedef struct human_markov {

  double        tnow;       /* time at which I entered my current state */
  std::string   state;      /* my current state */

  double        tnext;      /* time of my next state transition */
  std::string   nextstate;  /* my current state */

  // constructor & destructor
  human_markov(const std::string state_);
  ~human_markov();
} human_markov;


/* function to sample a trajectory */
Rcpp::NumericMatrix SIRMarkov_ABM(const double dt, const double tmax, const int S, const int I, const int R, const double beta, const double gamma, const bool verbose);


#endif
