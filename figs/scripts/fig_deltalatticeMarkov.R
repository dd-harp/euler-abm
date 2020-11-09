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
library(data.table)
library(ggplot2)
library(viridis)
# library(gridExtra)
library(parallel)
library(foreach)
library(doParallel)

# cores for parallel
ncores <- 16

fig_directory <- here::here("figs")
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

# save file
result_file <- paste(fig_directory, "/deltaMarkov_ABM.rds", sep = "", collapse = "")
cat(paste("writing to", result_file, "\n"))
saveRDS(object = deltaMarkov_ABM,file = result_file, compress = TRUE)


# --------------------------------------------------------------------------------
#   Compute the master equation directly
# --------------------------------------------------------------------------------

brates1 <- function(a,b){0}
drates1 <- function(a,b){0}
brates2 <- function(a,b){0}
drates2 <- function(a,b){gamma*b}
trans <- function(a,b){beta*a*b}

# dimensions give ending states
# rows S: a:a0 (0:S0)
# cols S: 0:B (0:S0+I0)
trans_dbd <- MultiBD::dbd_prob(
   t = tmax,a0 = S0,b0 =  I0,
   mu1 = drates1,lambda2 = brates2,mu2 = drates2,gamma = trans,
   a = 0, B = S0+I0,
   computeMode = 4,tol = 1e-16,nblocks = 512
)

dbd_rows <- as.integer(rownames(trans_dbd)) + 1
dbd_cols <- as.integer(colnames(trans_dbd)) + 1

trans_KFE <- matrix(data = 0, nrow = S0+I0+1, ncol = S0+I0+1)
trans_KFE[dbd_rows,dbd_cols] <- trans_dbd
trans_KFE <- trans_KFE / sum(trans_KFE)

dimnames(trans_KFE) <- list(as.character(0:(nrow(trans_KFE)-1)),as.character(0:(ncol(trans_KFE)-1)))


# --------------------------------------------------------------------------------
#   plot of error
# --------------------------------------------------------------------------------

# load MC data
deltaMarkov_ABM <- readRDS(file = file.path(fig_directory,"deltaMarkov_ABM.rds"))

trans_ABM <- replicate(
   n = length(abm_dt),
   expr = {matrix(data = 0, nrow = S0+I0+1, ncol = S0+I0+1,dimnames = list(0:(S0+I0), 0:(S0+I0)))},
   simplify = FALSE
)

for(dt in seq_along(abm_dt)){
   dt_ix <- which(deltaMarkov_ABM[,"dt"] == abm_dt[dt])
   for(j in seq_along(dt_ix)){
      this_row <- deltaMarkov_ABM[dt_ix[j],] + 1
      trans_ABM[[dt]][this_row["S"],this_row["I"]] <- trans_ABM[[dt]][this_row["S"],this_row["I"]] + 1
   }
   trans_ABM[[dt]] <- trans_ABM[[dt]] / sum(trans_ABM[[dt]])
}

# calculate differences
differences <- lapply(X = seq_along(trans_ABM),FUN = function(x){
   dif <- abs(trans_KFE - trans_ABM[[x]])
   dif <- as.data.table(as.table(dif))
   dif$V1 <- as.integer(dif$V1)
   dif$V2 <- as.integer(dif$V2)
   dif$N <- as.numeric(dif$N)
   dif$dt <- abm_dt[x]
   colnames(dif) <- c("S","I","density","dt")
   return(dif)
})
differences <- do.call(rbind,differences)

# absolute error plot
ggplot(data = differences[I<=40 & S <= 60,]) +
   geom_raster(aes(x=S,y=I,fill=density)) +
   scale_fill_viridis(option = "A") +
   facet_wrap(. ~ dt,labeller = labeller(dt = function(x){paste0("dt: ",x)})) +
   xlab("S") + ylab("I") + 
   labs(title="Absolute Error (ABM vs. KFE)") +
   guides(fill = guide_legend(title = "Error")) + 
   theme_bw() +
   theme(
      panel.background=element_rect(fill="grey90"),
      panel.grid=element_blank(),
      plot.title = element_text(size = rel(2.5)),
      axis.title = element_text(size = rel(1.5)),
      axis.text = element_text(size = rel(1.25)),
      strip.text.x = element_text(size = rel(1.2))
   )




# absolute diff
abs_dif <- sapply(X = seq_along(trans_ABM),FUN = function(x){
   sum(abs(trans_KFE - trans_ABM[[x]]))
})

# # means and variances
# n_mc <- 1e1
# blocks <- lapply(X = (1:n_mc),FUN = function(x){
#    if(x == 1){
#       return(c(1,x*(nrep/n_mc)))
#    } else {
#       return(c(
#          ((x-1)*(nrep/n_mc))+1,
#          x*(nrep/n_mc)
#       ))
#    }
# })
# 
# 
# dt_idx <- lapply(X = abm_dt,FUN = function(x){
#    unname(which(deltaMarkov_ABM[,"dt"] == x))
# })
# 
# 
# system.time(abs_dif_block <- parallel::mclapply(X = dt_idx,FUN = function(ix){
#    out <- rep(0,length(blocks))
#    mat_dt <- deltaMarkov_ABM[ix,]
#    mc_probs <- matrix(data = 0, nrow = S0+I0+1, ncol = S0+I0+1,dimnames = list(0:(S0+I0), 0:(S0+I0)))
#    for(i in seq_along(blocks)){
#       mc_probs <- mc_probs * 0
#       for(j in blocks[[i]][1]:blocks[[i]][2]){
#          this_row <- mat_dt[j,] + 1
#          mc_probs[this_row["S"],this_row["I"]] <- mc_probs[this_row["S"],this_row["I"]] + 1
#       }
#       mc_probs <- mc_probs/sum(mc_probs)
#       out[i] <- sum(abs(trans_KFE - mc_probs))
#    }
#    return(out)
# }))
# 
# abs_dif_block_mu <- sapply(abs_dif_block,mean)
# 
# abs_dif_block_ci_lo <- abs_dif_block_mu - (1.96 * sapply(abs_dif_block,sd)/sqrt(n_mc))
# abs_dif_block_ci_hi <- abs_dif_block_mu + (1.96 * sapply(abs_dif_block,sd)/sqrt(n_mc))
# 
# rbind(abs_dif_block_ci_lo,abs_dif_block_mu,abs_dif_block_ci_hi)
# 
# # sqrt( (trans_ABM[[1]])*(1 - trans_ABM[[1]])/N )
# # 
# # # 9 largest transition probabilities
# # largest <- order(trans_KFE,decreasing = TRUE)[1:9]
# # pos <- arrayInd(ind = largest,.dim = dim(trans_KFE),useNames = TRUE)
# # 
# # 
# # largest_ABM <- lapply(X = seq_along(abm_dt),FUN = function(x){
# #    probs <- trans_ABM[[x]][pos]
# #    probs_sd <- sqrt( (probs)*(1 - probs)/nrep )
# #    probs_lo <- probs - (probs_sd * 1.96)
# #    probs_hi <- probs + (probs_sd * 1.96)
# #    cbind(probs,probs_sd,probs_lo,probs_hi,dt=abm_dt[x])
# # })
# # # largest_ABM <- do.call(rbind,largest_ABM)
# # 
# # 
# # library(plotrix)
# # i = 7
# # plotrix::plotCI(x = 1:9,y = largest_ABM[[i]][,"probs"],ui = largest_ABM[[i]][,"probs_hi"],li = largest_ABM[[i]][,"probs_lo"],ylim=c(0,0.02))
# # points(x = 1:9,y = sort(trans_KFE,decreasing = T)[1:9],pch=16,col="steelblue")
