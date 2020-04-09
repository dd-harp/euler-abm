# --------------------------------------------------------------------------------
#
#   Master equations for simple Markov SIR model
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
# --------------------------------------------------------------------------------

SIR_getC <- function(N){
  ((N+1)*(N+2))/2
}

SIR_getix <- function(S,I,N){
  (N*S - 1/2 * (S*(S-3) + I + 1))
}
