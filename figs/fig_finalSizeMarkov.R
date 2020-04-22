# --------------------------------------------------------------------------------
#
#   Figure: compare final epidemic sizes in Markov SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(parallel)
library(foreach)
library(doSNOW)


# parameters
S0 = 50 # S_0
I0 = 1 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/5
beta <- R0 * (gamma/N)
nrep <- 1e4


# --------------------------------------------------------------------------------
#   MNRM final size
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(3)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 6978149L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

final_size_markovMNRM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRMarkov_MNRM(
    tmax = 1e4,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma = gamma,
    verbose = TRUE
  )

   tail(out,1)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()


# --------------------------------------------------------------------------------
#   ABM final size
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(3)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 9782101L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

final_size_markovABM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRMarkov_ABM(
    dt = 0.01,
    tmax = 1e4,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma = gamma,
    verbose = TRUE
  )

   n <- length(out$time)
   cbind(time=out$time[n],S=out$S[n],I=out$I[n],R=out$R[n])
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()


# --------------------------------------------------------------------------------
#   analytic final size
# --------------------------------------------------------------------------------

# MGF of exponential distribution
phi_exp <- function(b){
  1 / (1 - b*(1/gamma))
}

final_size_markovExact <- stocheulerABM::SIR_finalsize(N = S0, m = I0, lambda = beta*N, phi = phi_exp)


# --------------------------------------------------------------------------------
#   figure plot (8 x 14 PDF)
# --------------------------------------------------------------------------------

par(mfrow=c(1,2))

plot(0:S0,final_size_markovExact,
     xlab = "Final Number Infected",ylab = "Probability",
     main = "Analytic Solution vs. Exact Simulation",
     col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25,
     cex=1.25,
     type="b",pch=16,ylim = c(0,.3)
)
lines(
  x = 0:S0,
  y = table(as.vector(final_size_markovMNRM[,"R"])-1)/nrep,
  col = adjustcolor("steelblue",alpha.f = 0.95),
  type="b",pch=16,cex=1.25
)


plot(0:S0,final_size_markovExact,
     xlab = "Final Number Infected",ylab = "Probability",
     main = "Analytic Solution vs. Approximate ABM",
     col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25,
     cex=1.25,
     type="b",pch=16,ylim = c(0,.3)
)
lines(
  x = 0:S0,
  y = table(as.vector(final_size_markovABM[,"R"])-1)/nrep,
  col = adjustcolor("steelblue",alpha.f = 0.95),
  type="b",pch=16,cex=1.25
)

par(mfrow=c(1,1))
