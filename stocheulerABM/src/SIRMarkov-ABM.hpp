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
#include <limits>
#include <string>
#include <sstream>
#include <unordered_map>

#include <Rcpp.h>


/* --------------------------------------------------------------------------------
#   human struct
-------------------------------------------------------------------------------- */

typedef struct human_markov {

  double        tnow;       /* time at which I entered my current state */
  std::string   state;      /* my current state */

  double        tnext;      /* time of my next state transition */
  std::string   nextstate;  /* my current state */

  double        FOI;        /* my personal force of infection */
  double        gamma;      /* recovery rate */

  // constructor & destructor
  human_markov(const std::string state_, const double gamma_);
  ~human_markov();
} human_markov;

using human_markov_ptr = std::unique_ptr<human_markov>;
using humans_markov_vec = std::vector<human_markov_ptr>;


/* --------------------------------------------------------------------------------
#   human simulation
-------------------------------------------------------------------------------- */

// susceptible
void S_compartment(human_markov_ptr& human, const double tmax);

// infected
void I_compartment(human_markov_ptr& human, const double tmax);

// recovered
void R_compartment(human_markov_ptr& human, const double tmax);

/* updates for a human over a time step */
void simulate_human(human_markov_ptr& human, const double tmax);

/* calculates the FOI */
void calc_FOI(humans_markov_vec& humans, const double beta);


/* --------------------------------------------------------------------------------
#   function to sample a trajectory
-------------------------------------------------------------------------------- */

Rcpp::List SIRMarkov_ABM(const double dt, const double tmax, const int S, const int I, const int R, const double beta, const double gamma, const bool verbose);


#endif
