# --------------------------------------------------------------------------------
#
#   Figure: evaluate non-Markov SIR ABM over lattice of delta t
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
library(doParallel)

# cores for parallel
ncores <- 16

# parameters
S0 = 60 # S_0
I0 = 10 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/3.5
beta <- R0 * (gamma/N)

# non-Markov parameters
gamma_mean <- 1/gamma
gamma_var <- 0.5^2

gamma_shape <- (gamma_mean^2)/gamma_var
gamma_scale <- gamma_var/gamma_mean

# stop time
tmax <- 5

# mc reps
nrep <- 2e5

# lattice of time points
abm_dt <- c(0.001,0.005,0.01,0.025,0.05,0.075,0.1,0.5,1)


# --------------------------------------------------------------------------------
#   Draw (S,I) pair from ABM
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 257152385L)

# parallel options
opts <- list(preschedule=FALSE)

mc_reps <- expand.grid(rep = 1:nrep,dt = abm_dt)

deltaNonMarkov_ABM <- foreach(i = mc_reps$rep, dt = mc_reps$dt,.combine = "rbind",.options.multicore=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_ABM(
    tmax = tmax,
    dt = dt_ABM,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
  )

    c("dt" = dt,"S"=tail(out$S,1),"I"=tail(out$I,1))
}

# clean up the parallel cluster and remove it
parallel::stopCluster(cl);rm(cl);gc()

saveRDS(object = deltaNonMarkov_ABM,file = here::here("/figs/deltaNonMarkov_ABM.rds"),compress = TRUE)
