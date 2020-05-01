/* --------------------------------------------------------------------------------
#
#   non-Markovian SIR model via ABM
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
-------------------------------------------------------------------------------- */

#ifndef SIRNONMARKOV_ABM
#define SIRNONMARKOV_ABM

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

typedef struct human_nonmarkov {

  double        tnow;       /* time at which I entered my current state */
  std::string   state;      /* my current state */

  double        tnext;      /* time of my next state transition */
  std::string   nextstate;  /* my current state */

  double        FOI;        /* my personal force of infection */
  double        gamma_shape;
  double        gamma_scale;

  // constructor & destructor
  human_nonmarkov(const std::string state_, const double gamma_shape_, const double gamma_scale_);
  ~human_nonmarkov();
} human_nonmarkov;

using human_nonmarkov_ptr = std::unique_ptr<human_nonmarkov>;
using humans_nonmarkov_vec = std::vector<human_nonmarkov_ptr>;


/* --------------------------------------------------------------------------------
#   human simulation
-------------------------------------------------------------------------------- */

// susceptible
void S_compartment(human_nonmarkov_ptr& human, const double tmax);

// infected
void I_compartment(human_nonmarkov_ptr& human, const double tmax);

// recovered
void R_compartment(human_nonmarkov_ptr& human, const double tmax);

/* updates for a human over a time step */
void simulate_human(human_nonmarkov_ptr& human, const double tmax);

/* calculates the FOI */
void calc_FOI(humans_nonmarkov_vec& humans, const double beta);


/* --------------------------------------------------------------------------------
#   function to sample a trajectory
-------------------------------------------------------------------------------- */

Rcpp::List SIRnonMarkov_ABM(const double dt, const double tmax, const int S, const int I, const int R, const double beta, const double gamma, const bool verbose);


#endif
