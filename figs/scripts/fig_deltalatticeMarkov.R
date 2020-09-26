# --------------------------------------------------------------------------------
#
#   Figure: evaluate Markov SIR ABM over lattice of delta t
#   Sean L. Wu (slwu89@berkeley.edu)
#   September 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(MultiBD)
# library(entropy)
library(here)
# library(data.table)
# library(ggplot2)
# library(viridis)
# library(gridExtra)
library(parallel)
library(foreach)
library(doSNOW)

# cores for parallel
ncores <- 16

# parameters
S0 = 60 # S_0
I0 = 10 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/3.5
beta <- R0 * (gamma/N)
tmax <- 5

# mc reps
nrep <- 5e5

# lattice of time points
abm_dt <- c(0.001,0.005,0.01,0.025,0.05,0.075,0.1,0.5,1)


# --------------------------------------------------------------------------------
#   Draw (S,I) pair from ABM
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(ncores)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 186337428L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

mc_reps <- expand.grid(rep = 1:nrep,dt = abm_dt)

deltaMarkov_ABM <- foreach(i = mc_reps$rep, dt = mc_reps$dt,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

   out <- stocheulerABM::SIRMarkov_ABM(
     tmax = tmax,
     dt = dt,
     S = S0,
     I = I0,
     R = 0,
     beta = beta,
     gamma = gamma,
     verbose = FALSE
   )

    c("dt" = dt,"S"=tail(out$S,1),"I"=tail(out$I,1))
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

saveRDS(object = transraw_ABM,file = here::here("/figs/deltaMarkov_ABM.rds"),compress = TRUE)

