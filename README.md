# euler-abm
scripts and code to accompany the manuscript "Principled Simulation of Agent Based Models in Epidemiology"

## Directory structure
  * **stocheulerABM**: contains the R package `stocheulerABM` that has the simulation algorithms and analytic epidemic final size algorithm needed to produce the results of the manuscript. It requires `Rcpp` with a C++11 compiler to build.
  * **figs**: contains R scripts to produce all figures used in the manuscript as well as a Keynote document (*KeynoteFigs.key*) for graphics not produced in R. Subdirectory **scripts** contains R scripts to produce all figures, listed below.
    * *fig_deltalatticeMarkov.R*: produces Figure 5
    * *fig_finalSizeMarkov.R*: produces Figure 6
    * *fig_finalSizeNonMarkov.R*: produces Figure 8
    * *fig_inhomintensity.R*: produces Figure S2
    * *fig_sampleinhom.R*: produces Figure 2
    * *fig_trajectory.R*: produces Figure S3
    * *fig_transientMarkov.R*: produces Figure 4
    * *fig_transientNonMarkov.R*: produces Figure 7
  * **notebooks**: contains a Mathematica notebook used to check a few calculations.
