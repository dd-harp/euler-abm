# --------------------------------------------------------------------------------
#
#   Figure: compare final epidemic sizes in non-Markovian SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(PropCIs)
library(parallel)
library(foreach)
library(doSNOW)


# parameters
S0 = 50 # S_0
I0 = 1 # I_0
N <- S0 + I0
R0 <- 1.85

gamma_mean <- 5
gamma_var <- 0.5^2

gamma_shape <- (gamma_mean^2)/gamma_var
gamma_scale <- gamma_var/gamma_mean

beta <- R0 * (1/gamma_mean/N)
nrep <- 1e4


# --------------------------------------------------------------------------------
#   MNRM final size
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 191459L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

final_size_nonmarkovMNRM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_MNRM(
    tmax = 1e4,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
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
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 83595489L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

final_size_nonmarkovABM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_ABM(
    dt = 0.01,
    tmax = 1e4,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
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

# MGF of gamma distribution
phi_gamma <- function(b){
  stopifnot(b < 1/gamma_scale)
  1 / ( (1 - (gamma_scale*b))^gamma_shape )
}

final_size_nonmarkovExact <- stocheulerABM::SIR_finalsize(N = S0, m = I0, lambda = beta*N, phi = phi_gamma)


# --------------------------------------------------------------------------------
#   figure plot (10 x 16 PDF)
# --------------------------------------------------------------------------------

# probs MNRM
freq_MNRM <- rep(0,N)
tab_MNRM <- as.data.frame(table(as.vector(final_size_nonmarkovMNRM[,"R"])-I0))
freq_MNRM[as.integer(tab_MNRM$Var1)] <- tab_MNRM$Freq
probs_MNRM <- freq_MNRM/nrep

# probs_MNRM <- table(as.vector(final_size_nonmarkovMNRM[,"R"])-I0)/nrep

# probs ABM
freq_ABM <- rep(0,N)
tab_ABM <- as.data.frame(table(as.vector(final_size_nonmarkovABM[,"R"])-I0))
freq_ABM[as.integer(tab_ABM$Var1)] <- tab_ABM$Freq
probs_ABM <- freq_ABM/nrep

# probs_ABM <- table(as.vector(final_size_nonmarkovABM[,"R"])-I0)/nrep

# compute pointwise CIs
cis_MNRM <- lapply(X = freq_MNRM,FUN = function(x){
  PropCIs::scoreci(x = x,n = nrep,conf.level = 0.99)
})
cis_lo_MNRM <- sapply(X = cis_MNRM,FUN = function(x){
  x$conf.int[1]
})
cis_hi_MNRM <- sapply(X = cis_MNRM,FUN = function(x){
  x$conf.int[2]
})

cis_ABM <- lapply(X = freq_ABM,FUN = function(x){
  PropCIs::scoreci(x = x,n = nrep,conf.level = 0.99)
})
cis_lo_ABM <- sapply(X = cis_ABM,FUN = function(x){
  x$conf.int[1]
})
cis_hi_ABM <- sapply(X = cis_ABM,FUN = function(x){
  x$conf.int[2]
})


# par(mfrow=c(1,2))
#
# xadj <- 0.5
#
# # MNRM
# plot(0:S0,final_size_nonmarkovExact,
#      xlab = "Final Number Infected",ylab = "Probability",
#      main = "Analytic Solution vs. Exact Simulation (non-Markovian)",
#      col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25,
#      cex=1.1,
#      type="p",pch=16,ylim = c(0,.2)
# )
# invisible(mapply(FUN = function(x,y){
#   segments(
#     x0 = x-xadj,
#     y0 = y,
#     x1 = x+xadj,
#     y1 = y,
#     col = "darkorchid3",
#     lwd = 1.85,lend = 2
#   )
# },x=0:S0,y=probs_MNRM))
# rect(
#   xleft = (0:S0)-xadj,
#   ybottom = cis_lo_MNRM,
#   xright = (0:S0)+xadj,
#   ytop = cis_hi_MNRM,
#   border = NA,
#   col = adjustcolor("steelblue",alpha.f = 0.45)
# )
#
# # ABM
# plot(0:S0,final_size_nonmarkovExact,
#      xlab = "Final Number Infected",ylab = "Probability",
#      main = "Analytic Solution vs. Approximate ABM (non-Markovian)",
#      col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25,
#      cex=1.1,
#      type="p",pch=16,ylim = c(0,.2)
# )
# invisible(mapply(FUN = function(x,y){
#   segments(
#     x0 = x-xadj,
#     y0 = y,
#     x1 = x+xadj,
#     y1 = y,
#     col = "darkorchid3",
#     lwd = 1.85,lend = 2
#   )
# },x=0:S0,y=probs_ABM))
# rect(
#   xleft = (0:S0)-xadj,
#   ybottom = cis_lo_ABM,
#   xright = (0:S0)+xadj,
#   ytop = cis_hi_ABM,
#   border = NA,
#   col = adjustcolor("steelblue",alpha.f = 0.45)
# )
#
# par(mfrow=c(1,1))

library(ggplot2)
library(gridExtra)

dat_exact <- data.frame(x=0:S0,y=final_size_nonmarkovExact)
dat_mnrm_ci <- data.frame(x=0:S0,ylo=cis_lo_MNRM,yhi=cis_hi_MNRM)
dat_mnrm <- data.frame(x=0:S0,y=probs_MNRM)

plot_mnrm <- ggplot(data = dat_exact) +
  geom_point(aes(x=x,y=y),color="firebrick3",cex=2) +
  geom_errorbar(data=dat_mnrm_ci,aes(x=x,ymin=ylo,ymax=yhi),color="darkorchid3",alpha=0.85) +
  geom_point(data = dat_mnrm,aes(x=x,y=y),color="darkorchid3",alpha=0.5) +
  xlab("Final Number Infected") + ylab("Probability") + labs(title="A. Analytic Solution vs. Exact Simulation (non-Markovian)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.75)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

dat_abm_ci <- data.frame(x=0:S0,ylo=cis_lo_ABM,yhi=cis_hi_ABM)
dat_abm <- data.frame(x=0:S0,y=probs_ABM)

plot_abm <- ggplot(data = dat_exact) +
  geom_point(aes(x=x,y=y),color="firebrick3",cex=2) +
  geom_errorbar(data=dat_abm_ci,aes(x=x,ymin=ylo,ymax=yhi),color="darkorchid3",alpha=0.85) +
  geom_point(data = dat_abm,aes(x=x,y=y),color="darkorchid3",alpha=0.5) +
  xlab("Final Number Infected") + ylab("Probability") + labs(title="B. Analytic Solution vs. Approximate ABM (non-Markovian)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.75)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

grid.arrange(plot_mnrm,plot_abm,nrow=1)
# save as 8x16 pdf: nonmarkovSIR_finalsize.pdf
