% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{SIRnonMarkov_ABM}
\alias{SIRnonMarkov_ABM}
\title{Simulate non-Markovian SIR Model via Agent-based Model (ABM)}
\usage{
SIRnonMarkov_ABM(dt, tmax, S, I, R, beta, gamma_shape, gamma_scale, verbose)
}
\arguments{
\item{dt}{the time step}

\item{tmax}{maximum time of simulation}

\item{S}{initial number of susceptible individuals}

\item{I}{initial number of infected & infectious individuals}

\item{R}{initial number of recovered individuals}

\item{beta}{the product of transmission probability and contact rate}

\item{gamma_shape}{shape parameter of Gamma distributed infectious period}

\item{gamma_scale}{scale parameter of Gamma distributed infectious period}

\item{verbose}{print extra information?}
}
\value{
a list (use \code{do.call(cbind,out)} to convert to \code{matrix})
}
\description{
In this non-Markovian variant of the SIR model, the infectious period has a Gamma distribution. Sample a trajectory
from it using the approximate ABM.
}
