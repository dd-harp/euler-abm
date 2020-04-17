# --------------------------------------------------------------------------------
#
#   Closed form solutions for final epidemic size distributions for SIR models
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
# --------------------------------------------------------------------------------

#' Compute Epidemic Final Size Distribution for SIR Model
#'
#' Compute the final epidemic size of a stochastic SIR epidemic with arbitrary distribution of infectious period using the closed form method of:
#'    * Ball, Frank. "A unified approach to the distribution of total size and total area under the trajectory of infectives in epidemic models." Advances in Applied Probability 18.2 (1986): 289-310.
#'
#'  This method becomes numerically unstable for N larger than around 50-60.
#'
#'
#' @param N initial numer of susceptibles
#' @param m initial number of infectives
#' @param lambda the effective contact rate for a force of infection term following the mass action law
#' @param phi the MGF of the infectious period distribution (see \url{https://en.wikipedia.org/wiki/Moment-generating_function#Examples} for MGFs of various named distributions)
#'
#' @examples
#' \dontrun{
#' # this reproduces figure 2.2 from the book: Andersson, Hakan, and Tom Britton. Stochastic epidemic models and their statistical analysis. Vol. 151. Springer Science & Business Media, 2012.
#' # N=50 susceptibles and m=1 infected at time 0
#'  N <- 50
#'  lambda <- 1.5
#'  m <- 1
#' # duration of infectiousness is deterministic and lasts 1 day
#'  phi <- function(b){
#'    exp(b)
#'  }
#'
#'  probs <- SIR_finalsize(N = N,m = m,lambda = lambda,phi = phi)
#'  plot(0:N,probs,pch=16)
#'  abline(h=0)
#' }
#'
#' @export
SIR_finalsize <- function(N,m,lambda,phi){
  probs <- rep(NaN,N+1)
  k <- 0
  probs[k+1] <- phi(-lambda)^m
  repeat{
    k <- k + 1
    if(k+1>N+1){
      break()
    }
    b <- ((N-k)*lambda) / N
    term1 <- choose(N,k)*(phi(-b)^(k+m))
    term2 <- 0
    for(i in 0:(k-1)){
      term2 <- term2 + (choose(N-i,k-i) * (phi(-b)^(k-i)) * probs[i+1])
    }
    probs[k+1] <- term1 - term2
  }
  return(probs)
}
