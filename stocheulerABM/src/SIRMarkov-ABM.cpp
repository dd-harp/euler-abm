/* --------------------------------------------------------------------------------
#
#   Markovian SIR model via ABM
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "SIRMarkov-ABM.hpp"

// a very small quantity
static const double eps = std::numeric_limits<double>::epsilon();


/* --------------------------------------------------------------------------------
#   human struct
-------------------------------------------------------------------------------- */

// ctor
human_markov::human_markov(const std::string state_, const double gamma_) :
  tnow{0.}, state(state_),
  tnext{0.}, nextstate(state_),
  FOI{0.}, gamma(gamma_)
{};

// dtor
human_markov::~human_markov(){};


/* --------------------------------------------------------------------------------
#   human simulation
-------------------------------------------------------------------------------- */

// put the functions in a hash table
static const std::unordered_map<std::string, std::function<void(human_markov_ptr&, const double)> > compartment_functions  = {
  {"S",S_compartment},
  {"I",I_compartment},
  {"R",R_compartment}
};

// S: susceptible
void S_compartment(human_markov_ptr& human, const double tmax){

  // update my current state and time
  human->tnow = human->tnext;
  human->state = human->nextstate;

  // putative infection time
  double tau = human->tnow + R::rexp(1./human->FOI);

  // if it occurs on this time step, go ahead and accept, otherwise reject
  if(tau < tmax){
    human->tnext = tau;
    human->nextstate = "I";
  } else {
    human->tnext = tmax + eps;
    human->nextstate = "S";
  }

};

// I: infected
void I_compartment(human_markov_ptr& human, const double tmax){

  // update my current state and time
  human->tnow = human->tnext;
  human->state = human->nextstate;

  /* sample transition to R */
  human->tnext = human->tnow + R::rexp(1./human->gamma);
  human->nextstate = "R";

};

// R: recovered
void R_compartment(human_markov_ptr& human, const double tmax){

  // update my current state and time
  human->tnow = human->tnext;
  human->state = human->nextstate;

  /* stay in R forever */
  human->tnext = std::numeric_limits<double>::infinity();
  human->nextstate = "R";

};

/* updates for a human over a time step */
void simulate_human(human_markov_ptr& human, const double tmax){
  while(human->tnext < tmax){
    compartment_functions.at(human->nextstate)(human,tmax);
  }
};

/* calculates the FOI */
void calc_FOI(humans_markov_vec& humans, const double beta){

  // compute the FOI
  int I{0};
  for(auto& human : humans){
    if(human->state.compare("I") == 0){
      I += 1;
    }
  }
  double FOI = beta * (double)I;

  // update it for everyone
  for(auto& human : humans){
    if(human->state.compare("S") == 0){
      human->FOI = FOI;
    }
  }
};


/* --------------------------------------------------------------------------------
#   function to sample a trajectory
-------------------------------------------------------------------------------- */

//' Simulate Markovian SIR Model via Agent-based Model (ABM)
//'
//' Sample a trajectory from the Markovian SIR model using the approximate ABM.
//'
//' @param dt the time step
//' @param tmax maximum time of simulation
//' @param S initial number of susceptible individuals
//' @param I initial number of infected & infectious individuals
//' @param R initial number of recovered individuals
//' @param beta the product of transmission probability and contact rate
//' @param gamma inverse of the duration of infectiousness
//' @param verbose print extra information?
//'
//' @return a list (use \code{do.call(cbind,out)} to convert to \code{matrix})
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::List SIRMarkov_ABM(
  const double dt,
  const double tmax,
  const int S,
  const int I,
  const int R,
  const double beta,
  const double gamma,
  const bool verbose
){

  int N = S+I+R;

  int nstep = tmax / dt + 1;
  int nout{0};
  double tnow{0.};

  // trajectory of the process
  std::vector<double> t_hist(nstep,0.);
  std::vector<int> S_hist(nstep,0);
  std::vector<int> I_hist(nstep,0);
  std::vector<int> R_hist(nstep,0);

  // t=0 output
  t_hist[nout] = tnow;
  S_hist[nout] = S;
  I_hist[nout] = I;
  R_hist[nout] = R;
  nout += 1;

  // initialize humans
  if(verbose){Rcpp::Rcout << " --- initializing humans --- \n";}

  humans_markov_vec humans;
  humans.reserve(N);
  for(int i=0; i<N; i++){
    std::string state("");
    /* initial susceptible */
    if(i < S){
      state = "S";
    /* initial recovered */
    } else if(i >= S && i < (S+R)){
      state = "R";
    /* initial infected */
    } else if(i >= (S+R)){
      state = "I";
    }
    humans.emplace_back(std::make_unique<human_markov>(state,gamma));
  }


  if(verbose){Rcpp::Rcout << " --- beginning simulation --- \n";}

  /* initial FOI */
  calc_FOI(humans,beta);

  /* main simulation loop */
  while(nout < nstep){

    if(nout % 10 == 0){
      Rcpp::checkUserInterrupt();
    }

    tnow += dt;

    /* simulate humans */
    for(auto& human : humans){
      simulate_human(human,tnow);
    }

    /* update FOI */
    calc_FOI(humans,beta);

    /* push output */
    t_hist[nout] = tnow;
    for(auto& human : humans){
      if(human->state.compare("S") == 0){
        S_hist[nout] += 1;
      } else if(human->state.compare("I") == 0){
        I_hist[nout] += 1;
      } else if(human->state.compare("R") == 0){
        R_hist[nout] += 1;
      } else {
        std::stringstream msg;
        msg << "human state is illegal value: " << human->state;
        Rcpp::stop(msg.str());
      }
    }

    /* if epidemic is over return early */
    if(R_hist[nout] == N){
      if(verbose){Rcpp::Rcout << " --- all individuals in 'R', returning early --- \n";}
      break;
    }

    nout += 1;
  }

  /* resize output and return */
  t_hist.resize(nout-1);
  S_hist.resize(nout-1);
  I_hist.resize(nout-1);
  R_hist.resize(nout-1);

  return Rcpp::List::create(
    Rcpp::Named("time") = t_hist,
    Rcpp::Named("S") = S_hist,
    Rcpp::Named("I") = I_hist,
    Rcpp::Named("R") = R_hist
  );
};
