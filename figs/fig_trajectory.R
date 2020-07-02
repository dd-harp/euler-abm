# --------------------------------------------------------------------------------
#
#   Figure: compare trajectories
#   Sean L. Wu (slwu89@berkeley.edu)
#   June 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(data.table)
library(ggplot2)
library(gridExtra)
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

nrep <- 1e4

# non-Markov parameters
gamma_mean <- 1/gamma
gamma_var <- 0.5^2

gamma_shape <- (gamma_mean^2)/gamma_var
gamma_scale <- gamma_var/gamma_mean


# --------------------------------------------------------------------------------
#   Markov SIR: MNRM trajectories
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 458946L)

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

traj_MNRM <- as.data.table(traj_MNRM)
traj_MNRM <- melt(traj_MNRM,id.vars = c("time","mc"),measure.vars = c("S","I","R"),variable.name = "state",value.name = "count")

traj_MNRM <- traj_MNRM[,
          list(
            mean = mean(count),
            lo = quantile(count,probs = 0.025),
            hi = quantile(count,probs = 0.975)
          ),
          ,by = c("time","state")
          ]


# --------------------------------------------------------------------------------
#   Markov SIR: ABM trajectories
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 6859423L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

traj_ABM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRMarkov_ABM(
    dt = 0.005,
    tmax = 200,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma = gamma,
    verbose = FALSE
  )
  out <- do.call(cbind,out)

  out <- stocheulerABM::discretise(out,dt=0.25)
  out <- cbind(out,mc=i)
  return(out)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

traj_ABM <- as.data.table(traj_ABM)
traj_ABM <- melt(traj_ABM,id.vars = c("time","mc"),measure.vars = c("S","I","R"),variable.name = "state",value.name = "count")

traj_ABM <- traj_ABM[,
          list(
            mean = mean(count),
            lo = quantile(count,probs = 0.025),
            hi = quantile(count,probs = 0.975)
          ),
          ,by = c("time","state")
          ]


traj_ABM$model <- "ABM"
traj_MNRM$model <- "SSA"

traj_markov <- rbind(traj_ABM,traj_MNRM)


# --------------------------------------------------------------------------------
#   non-Markovian SIR: MNRM trajectories
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 695941L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

traj_MNRM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_MNRM(
    tmax = 200,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
  )

  out <- stocheulerABM::discretise(out,dt=0.25)
  out <- cbind(out,mc=i)
  return(out)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

traj_MNRM <- as.data.table(traj_MNRM)
traj_MNRM <- melt(traj_MNRM,id.vars = c("time","mc"),measure.vars = c("S","I","R"),variable.name = "state",value.name = "count")

traj_MNRM <- traj_MNRM[,
          list(
            mean = mean(count),
            lo = quantile(count,probs = 0.025),
            hi = quantile(count,probs = 0.975)
          ),
          ,by = c("time","state")
          ]


# --------------------------------------------------------------------------------
#   non-Markovian SIR: ABM trajectories
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 977173L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

traj_ABM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_ABM(
    dt = 0.005,
    tmax = 200,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
  )
  out <- do.call(cbind,out)

  out <- stocheulerABM::discretise(out,dt=0.25)
  out <- cbind(out,mc=i)
  return(out)
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

traj_ABM <- as.data.table(traj_ABM)
traj_ABM <- melt(traj_ABM,id.vars = c("time","mc"),measure.vars = c("S","I","R"),variable.name = "state",value.name = "count")

traj_ABM <- traj_ABM[,
          list(
            mean = mean(count),
            lo = quantile(count,probs = 0.025),
            hi = quantile(count,probs = 0.975)
          ),
          ,by = c("time","state")
          ]


traj_ABM$model <- "ABM"
traj_MNRM$model <- "SSA"

traj_nonmarkov <- rbind(traj_ABM,traj_MNRM)


# --------------------------------------------------------------------------------
#   plot
# --------------------------------------------------------------------------------

p_markov <- ggplot(data = traj_markov[time<=80,]) +
  geom_line(aes(x=time,y=mean,color=state,linetype=model)) +
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=state,color=state,linetype=model),alpha=0.2) +
  guides(linetype=FALSE,color=FALSE,fill=FALSE) +
  xlab("Time") + ylab("Count") + labs(title="A. Markovian SIR") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(2.5)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

p_nonmarkov <- ggplot(data = traj_nonmarkov[time<=40,]) +
  geom_line(aes(x=time,y=mean,color=state,linetype=model)) +
  geom_ribbon(aes(x=time,ymin=lo,ymax=hi,fill=state,color=state,linetype=model),alpha=0.2) +
  guides(linetype=FALSE,color=FALSE,fill=FALSE) +
  xlab("Time") + ylab("Count") + labs(title="B. non-Markovian SIR") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(2.5)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

pp <- grid.arrange(p_markov,p_nonmarkov,nrow=1)
# save as 8 x 14 landscape PDF