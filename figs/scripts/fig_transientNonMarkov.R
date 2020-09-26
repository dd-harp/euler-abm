# --------------------------------------------------------------------------------
#
#   Figure: compare transient transition distributions of non-Markovian SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   June 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(stocheulerABM)
library(data.table)
library(ggplot2)
library(viridis)
library(gridExtra)
library(parallel)
library(foreach)
library(doSNOW)

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

nrep <- 1e5
dt_ABM <- 0.01

tmax <- 5
abm_dt <- c(0.01)

# --------------------------------------------------------------------------------
#   Draw (S,I) pair from MNRM
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 158229L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

transraw_MNRM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

  out <- stocheulerABM::SIRnonMarkov_MNRM(
    tmax = tmax,
    S = S0,
    I = I0,
    R = 0,
    beta = beta,
    gamma_shape = gamma_shape,
    gamma_scale = gamma_scale,
    verbose = FALSE
  )

   tail(out,1)[,c("S","I")]
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()


trans_MNRM <- matrix(data = 0, nrow = S0+I0+1, ncol = S0+I0+1)
for(i in 1:nrep){
  row_i <- transraw_MNRM[i,] + 1
  trans_MNRM[row_i["S"],row_i["I"]] <- trans_MNRM[row_i["S"],row_i["I"]] + 1
}
trans_MNRM <- trans_MNRM / sum(trans_MNRM)

dimnames(trans_MNRM) <- list(as.character(0:(nrow(trans_MNRM)-1)),as.character(0:(ncol(trans_MNRM)-1)))


# --------------------------------------------------------------------------------
#   Draw (S,I) pair from ABM
# --------------------------------------------------------------------------------

# set up cluster and source the file on each core
cl <- snow::makeSOCKcluster(4)
doSNOW::registerDoSNOW(cl)

# use parallel RNG for multiple streams of random numbers
parallel::clusterSetRNGStream(cl = cl,iseed = 82459129L)

# progress bar in parallel foreach (see: https://blog.revolutionanalytics.com/2015/10/updates-to-the-foreach-package-and-its-friends.html)
pb <- txtProgressBar(max=nrep, style=3)
progress <- function(n){setTxtProgressBar(pb, n)}
opts <- list(progress=progress)

transraw_ABM <- foreach(i = 1:nrep,.combine = "rbind",.options.snow=opts,.packages = c("stocheulerABM")) %dopar% {

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

    c("S"=tail(out$S,1),"I"=tail(out$I,1))
}

# clean up the parallel cluster and remove it
close(pb)
snow::stopCluster(cl);rm(cl);gc()

trans_ABM <- matrix(data = 0, nrow = S0+I0+1, ncol = S0+I0+1)
for(i in 1:nrep){
  row_i <- transraw_ABM[i,] + 1
  trans_ABM[row_i["S"],row_i["I"]] <- trans_ABM[row_i["S"],row_i["I"]] + 1
}
trans_ABM <- trans_ABM / sum(trans_ABM)

dimnames(trans_ABM) <- list(as.character(0:(nrow(trans_ABM)-1)),as.character(0:(ncol(trans_ABM)-1)))


# --------------------------------------------------------------------------------
#   format output and plot
# --------------------------------------------------------------------------------

trans_MNRM_dt <- as.data.table(as.table(trans_MNRM))
trans_MNRM_dt$V1 <- as.integer(trans_MNRM_dt$V1)
trans_MNRM_dt$V2 <- as.integer(trans_MNRM_dt$V2)
trans_MNRM_dt$N <- as.numeric(trans_MNRM_dt$N)
colnames(trans_MNRM_dt) <- c("S","I","density")
trans_MNRM_dt$model <- "MNRM"

trans_ABM_dt <- as.data.table(as.table(trans_ABM))
trans_ABM_dt$V1 <- as.integer(trans_ABM_dt$V1)
trans_ABM_dt$V2 <- as.integer(trans_ABM_dt$V2)
trans_ABM_dt$N <- as.numeric(trans_ABM_dt$N)
colnames(trans_ABM_dt) <- c("S","I","density")
trans_ABM_dt$model <- "ABM"

trans_nonmarkov <- rbind(trans_MNRM_dt,trans_ABM_dt)

ggplot(data = trans_nonmarkov) +
  geom_contour(aes(x=S,y=I,z=density,colour=after_stat(level),group=model,linetype=model),size=0.85,alpha=1) +
  scale_color_viridis(option = "D") +
  scale_linetype_manual(values = c("ABM"=2,"MNRM"=1)) +
  guides(linetype=FALSE,colour=FALSE) +
  xlab("S") + ylab("I") + labs(title="Non-Markovian SIR (MNRM vs. ABM)") +
  theme_bw() +
  theme(
    panel.background=element_rect(fill="grey90"),
    panel.grid=element_blank(),
    plot.title = element_text(size = rel(2.5)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )


# save as 8 x 10 landscape PDF: contour_nonMarkov.pdf
