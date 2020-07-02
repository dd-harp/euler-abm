/* --------------------------------------------------------------------------------
#
#   Markovian SIR model via Anderson's MNRM
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "SIRMarkov-MNRM.hpp"

const static double eps = 1.E-9;


//' Simulate Markovian SIR Model via Modified Next Reaction Method (MNRM)
//'
//' Sample a trajectory from the Markovian SIR model using the MNRM algorithm presented in:
//'   * Anderson, D. F. (2007). A modified next reaction method for simulating chemical systems with time dependent propensities and delays. Journal of Chemical Physics, 127(21). \url{https://doi.org/10.1063/1.2799998}
//'
//' @param tmax maximum time of simulation
//' @param S initial number of susceptible individuals
//' @param I initial number of infected & infectious individuals
//' @param R initial number of recovered individuals
//' @param beta the product of transmission probability and contact rate
//' @param gamma inverse of the duration of infectiousness
//' @param verbose print extra information?
//'
//' @return a matrix
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix SIRMarkov_MNRM(
  const double tmax,
  const int S,
  const int I,
  const int R,
  const double beta,
  const double gamma,
  const bool verbose
){

  // system state
  std::array<int,3> X{S,I,R};

  double tnow{0.};
  const int outsize = 1E5;
  int i{1};

  // trajectory of the process
  std::vector<double> t_hist;
  std::vector<double> S_hist;
  std::vector<double> I_hist;
  std::vector<double> R_hist;
  t_hist.reserve(outsize);
  S_hist.reserve(outsize);
  I_hist.reserve(outsize);
  R_hist.reserve(outsize);
  t_hist.push_back(tnow);
  S_hist.push_back(X[0]);
  I_hist.push_back(X[1]);
  R_hist.push_back(X[2]);

  // 1. initialize
  std::array<double,2> Pk{0.,0.};   // next internal firing time of Poisson process (Pk > Tk)
  std::array<double,2> Tk{0.,0.};   // internal time of Poisson process (integrated propensity)
  std::array<double,2> delta_t{0.,0.};   // absolute time to fire
  std::array<double,2> ak{0.,0.};        // propensity functions

  int mu{0};
  double Delta{0.};

  // 2. calculate propensities
  ak[0] = beta * static_cast<double>(X[0]) * static_cast<double>(X[1]);
  ak[1] = gamma * static_cast<double>(X[1]);

  // 3-4: draw internal jump times
  Pk[0] = log(1. / R::runif(0.,1.));
  Pk[1] = log(1. / R::runif(0.,1.));

  if(verbose){
    Rcpp::Rcout << " --- beginning simulation --- \n";
  }

  while(tnow < tmax){

    if(verbose){
      if(i % 100 == 0){
        Rcpp::Rcout << " --- simulated " << i << " reactions at time " << tnow << " --- \n";
      }
    }

    // 5. set absolute times to fire
    delta_t[0] = (Pk[0] - Tk[0]) / ak[0];
    delta_t[1] = (Pk[1] - Tk[1]) / ak[1];

    // 6. find minimum
    mu = 0;
    if(delta_t[1] < delta_t[0]){
      mu = 1;
    }
    Delta = delta_t[mu];

    // 7. set t += delta and update state
    tnow += Delta;
    // check if we can break early
    if(tnow > tmax){
      break;
    }
    if(mu==0){
      // S -> I
      X[0] -= 1;
      X[1] += 1;
    } else if(mu==1){
      // I -> R
      X[1] -= 1;
      X[2] += 1;
    } else {
      Rcpp::stop("illegal value of minimum reaction 'mu'");
    }

    // 8. update Tk
    Tk[0] += (ak[0]*Delta);
    Tk[1] += (ak[1]*Delta);

    // 9. update P_mu
    Pk[mu] += log(1./R::runif(0.,1.));

    // 10. recalculate propensities
    ak[0] = beta * static_cast<double>(X[0]) * static_cast<double>(X[1]);
    ak[1] = gamma * static_cast<double>(X[1]);

    // store output
    t_hist.push_back(tnow);
    S_hist.push_back(X[0]);
    I_hist.push_back(X[1]);
    R_hist.push_back(X[2]);
    i++;

    // check to see if we can end early
    if((ak[0] < eps) && (ak[1] < eps)){
      Rcpp::Rcout << " --- all propensities approximately zero; returning output early --- \n";
      int k = t_hist.size();
      Rcpp::NumericMatrix out(k,4);
      Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I","R");
      for(int j=0; j<k; j++){
        out(j,0) = t_hist[j];
        out(j,1) = S_hist[j];
        out(j,2) = I_hist[j];
        out(j,3) = R_hist[j];
      }
      return out;
    }
  }

  int k = t_hist.size();
  Rcpp::NumericMatrix out(k,4);
  Rcpp::colnames(out) = Rcpp::CharacterVector::create("time","S","I","R");
  for(int j=0; j<k; j++){
    out(j,0) = t_hist[j];
    out(j,1) = S_hist[j];
    out(j,2) = I_hist[j];
    out(j,3) = R_hist[j];
  }
  return out;
};
