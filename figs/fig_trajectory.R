# --------------------------------------------------------------------------------
#
#   Figure: compare trajectories
#   Sean L. Wu (slwu89@berkeley.edu)
#   June 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(parallel)
library(foreach)
library(doSNOW)

# parameters
S0 = 100 # S_0
I0 = 1 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/5
beta <- R0 * (gamma/N)

nrep <- 1e3


# --------------------------------------------------------------------------------
#   MNRM trajectories
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 9009543L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

traj_MNRM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRMarkov_MNRM(
    tmax = 200,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma = gamma,
    verbose = FALSE
  )

  out <- stocheulerABM::discretise(out,dt=0.25)
  out <- cbind(out,mc=i)
  return(out)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()
