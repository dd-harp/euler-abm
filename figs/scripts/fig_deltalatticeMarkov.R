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
library(doParallel)

# cores for parallel
ncores <- 16

fig_directory <- here::here("/figs")
if (!dir.exists(fig_directory)) {
  dir.create(fig_directory)
}

# parameters
S0 = 60 # S_0
I0 = 10 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/3.5
beta <- R0 * (gamma/N)
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
parallel::clusterSetRNGStream(cl = cl,iseed = 186337428L)

# parallel options
opts <- list(preschedule=FALSE)

mc_reps <- expand.grid(rep = 1:nrep,dt = abm_dt)

deltaMarkov_ABM <- foreach(i = mc_reps$rep, dt = mc_reps$dt,.combine = "rbind",.options.multicore=opts,.packages = c("stocheulerABM")) %dopar% {

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
parallel::stopCluster(cl);rm(cl);gc()
result_file <- paste(fig_directory, "/deltaMarkov_ABM.rds", sep = "", collapse = "")
cat(paste("writing to", result_file, "\n"))
saveRDS(object = deltaMarkov_ABM,file = result_file, compress = TRUE)
