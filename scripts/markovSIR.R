# --------------------------------------------------------------------------------
#
#   Markovian SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(MultiBD)
library(stocheulerABM)

tList <- c(.1, .2, .25, .3 ,.35, .4, .5, .6, .7, .8, .9, 1)
gridLength = 128
S0 = 110 # S_0
I0 = 15 # I_0

A=0
B = gridLength - 1
gamma = 3.2 #3.2 #this is death rate
beta = .025 #.019 #this is transition or infection rates
brates1=function(S,I){0}
drates1=function(S,I){0}
brates2=function(S,I){0}
drates2=function(S,I){gamma*I}
trans=function(S,I){beta*S*I}

tmax <- 100

system.time(SIRexact <- MultiBD::dbd_prob(
  t=tmax,
  a0=S0,
  b0=I0,
  mu1=drates1,
  lambda2=brates2,
  mu2=drates2,
  gamma=trans,
  a=0,
  B=B,
  computeMode=1
))


# packages to run in parallel
library(parallel)
library(foreach)
library(doSNOW)

tmax <- 100
S0 = 100 # S_0
I0 = 5 # I_0
N <- S0 + I0
R0 <- 1.85
gamma <- 1/5
beta <- R0 * (gamma/N)

out <- stocheulerABM::SIRMarkov_MNRM(
  tmax = tmax,
  S = S0,
  I = I0,
  R = 0,
  beta = beta,
  gamma = gamma,
  verbose = TRUE
)

out <- stocheulerABM::discretise(out = out,dt = 1)

library(reshape2)
library(ggplot2)

outdf <- melt(as.data.frame(out),id.vars="time")

ggplot(outdf) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()


nrep <- 1e3

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(2)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 32423423L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

mc_sims <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM","reshape2")) %dopar% {

  out <- stocheulerABM::SIRMarkov_MNRM(
    tmax = tmax,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma = gamma,
    verbose = TRUE
  )

  out <- stocheulerABM::discretise(out = out,dt = 1)
  outdf <- reshape2::melt(as.data.frame(out),id.vars="time")
  outdf$mc <- i
  storage.mode(outdf$time) <- "integer"
  storage.mode(outdf$value) <- "integer"
  outdf$variable <- as.character(outdf$variable)
  colnames(outdf) <- c("time","state","count","mc")
  return(outdf)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

ggplot(mc_sims) +
  geom_path(aes(x=time,y=count,color=state,group=interaction(mc,state)),alpha=0.15) +
  theme_bw()

ggplot(mc_sims[which(mc_sims$state=="I"),]) +
  geom_path(aes(x=time,y=count,group=mc),color="darkorchid3",alpha=0.15) +
  theme_bw()

# compute means and 95% quantiles
mc_sims %>%
  as_tibble %>%
  group_by(state,time) %>%
  summarise(mean = mean(count),
            q_lo = quantile(count, 0.025),
            q_hi = quantile(count, 0.975)) -> mc_simsQuant

ggplot(mc_simsQuant) +
  geom_path(aes(x=time,y=mean,color=state),alpha=0.85) +
  geom_ribbon(aes(x=time,ymin=q_lo,ymax=q_hi,fill=state),alpha=0.45) +
  theme_bw()


# # exact
# lambda <- 1.5
# I <- 1
# N <- 50
# pk <- rep(0,50)
# pk[1] <- exp(-lambda)
# pk[2] <- N*exp(-((N-1)*lambda)/N) * (exp(-((N-1)*lambda)/N) - pk[1])
#
# # assumes I(0)=1
# exact_probs <- function(N,lambda,phi){
#   probs <- rep(NaN,N+1)
#   k <- 0
#   probs[k+1] <- phi(-lambda)
#   repeat{
#     k <- k + 1
#     if(k+1>N+1){
#       break()
#     }
#     b <- ((N-k)*lambda) / N
#     term1 <- choose(N,k)*(phi(-b)^(k+1))
#     term2 <- 0
#     for(i in 0:(k-1)){
#       term2 <- term2 + (choose(N-i,k-i) * (phi(-b)^(k-i)) * probs[i+1])
#     }
#     probs[k+1] <- term1 - term2
#   }
#   return(probs)
# }
#
# phi <- function(b){
#   exp(b)
# }
#
# N <- 50
# lambda <- 1.5
#
# probs <- exact_probs(N = N,lambda = lambda,phi = phi)
# plot(0:N,probs,pch=16)
# abline(h=0)



out <- stocheulerABM::SIRMarkov_ABM(
  dt = 0.01,
  tmax = tmax,
  S = S0,
  I = I0,
  R = 0,
  beta = beta,
  gamma = gamma,
  verbose = TRUE
)
out <- as.data.frame(do.call(cbind,out))
out <- melt(out,id.vars="time")

ggplot(out) +
  geom_line(aes(x=time,y=value,color=variable)) +
  theme_bw()



# FINAL SIZES

S0 = 50 # S_0
I0 = 1 # I_0
N <- S0 + I0
R0 <- 2.5
gamma <- 1/5
beta <- R0 * (gamma/N)

# --------------------------------------------------------------------------------
#   MNRM final size
# --------------------------------------------------------------------------------

nrep <- 1e4

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

nrep <- 1e4

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

phi_exp <- function(b){
  1 / (1 - b*(1/gamma))
}

final_size_markovExact <- stocheulerABM::SIR_finalsize(N = S0, m = I0, lambda = beta*N, phi = phi_exp)


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
