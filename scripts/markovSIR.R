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


# exact
lambda <- 1.5
I <- 1
N <- 50
pk <- rep(0,50)
pk[1] <- exp(-lambda)
pk[2] <- N*exp(-((N-1)*lambda)/N) * (exp(-((N-1)*lambda)/N) - pk[1])

# assumes I(0)=1
exact_probs <- function(N,lambda,phi){
  probs <- rep(NaN,N+1)
  k <- 0
  probs[k+1] <- phi(-lambda)
  repeat{
    k <- k + 1
    if(k+1==N+1){
      break()
    }
    b <- ((N-k)*lambda) / N
    term1 <- choose(N,k)*(phi(-b)^(k+1))
    term2 <- 0
    for(i in 0:(k-1)){
      term2 <- term2 + (choose(N-i,k-i) * (phi(-b)^(k-i)) * probs[i+1])
    }
    probs[k+1] <- term1 - term2
  }
  
  # while(k <= N+1){
  #   k <- k + 1
  #   b <- ((N-k)*lambda) / N
  #   term1 <- choose(N,k)*(phi(-b)^(k+1))
  #   term2 <- 0
  #   for(i in 0:(k-1)){
  #     term2 <- term2 + (choose(N-i,k-i) * (phi(-b)^(k-i)) * probs[i+1])
  #   }
  #   probs[k+1] <- term1 - term2
  # }
  return(probs)
}

phi <- function(b){
  exp(b)
}

probs <- exact_probs(N = N,lambda = lambda,phi = phi)

debug(exact_probs)
probs <- exact_probs(N = N,lambda = lambda,phi = phi)
