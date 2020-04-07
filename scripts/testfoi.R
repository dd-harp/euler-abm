rm(list=ls());gc()
library(parallel)

simmer <- function(tend,rate,dt){
  tnow <- tnext <- 0
  nreject <- 0
  while(tnow < tend){
    tnext <- tnow + dt
    putative <- tnow + rexp(1,rate)
    if(putative < tnext){
      return(setNames(object = c(nreject,putative),nm = c("nreject","time")))
    }
    nreject <- nreject + 1
    tnow <- tnext
  }
  return(NaN)
}

# debug(simmer)
# simmer(tend = 1e4,rate = 1/20,dt = 1)

n <- 2e6
lambda <- 1/20
dt <- 0.5

out_05 <- mclapply(X = 1:n,FUN = function(i){simmer(tend = 1e4,rate = lambda,dt = dt)})
out_05 <- do.call(rbind,out_05)

out_1 <- mclapply(X = 1:n,FUN = function(i){simmer(tend = 1e4,rate = lambda,dt = 1)})
out_1 <- do.call(rbind,out_1)
# out <- replicate(1e5,simmer(tend = n,rate = 1/20,dt = 1))

# colMeans(out_05)
# 
# p_dt <- 1-exp(-lambda*dt)
# 1/p_dt
# mean(rgeom(n = n,p_dt)+1)
# 
# (1/p_dt)*0.5
# 
# (1/lambda) + (dt / (1-exp(lambda*dt)))
# 
# 
# b = 1.5
# lambda <- 1/2
# mean(rtrunc(n = 4e6,spec = "exp",b = b,rate = lambda))
# 
# (1/lambda) + (b / (1-exp(lambda*b)))


# NHPP: piecewise constant intensity

findInterval()

x <- seq(0,48,by=0.01)
cyclicIntensity <- function(t){
  # t <- t %% 24
  (1/8)*sin((t - 6)*2*pi/24)+((1/8)+.Machine$double.eps)
}

cyclicApprox <- approx(x = 0:23,y = cyclicIntensity(0:23),method = "constant",n=length(0:23))

plot(x,cyclicIntensity(x),type="l",ylim = c(0,0.25),
     main = "Inhomogeneous Intensity Function",xaxt="n",
     xlab = "Time (hours)",ylab = "Intensity",col="firebrick3",lwd=1.5)

for(i in 0:48){
  
  j <- (i %% 24) + 1
  
  segments(x0 = i,x1 = i+1,
           y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j],
           col = "darkorchid3",lwd=1.5)
  points(x=i,y=cyclicApprox$y[j],pch=16,col="darkorchid3",cex=0.5)
  points(x=i+1,y=cyclicApprox$y[j],pch=1,col="darkorchid3",cex=0.5)
  if(i != 48){
    segments(x0 = i+1,x1 = i+1,y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j+1],
             lty = 2,col="darkorchid3") 
  }
  
}
axis(side = 1,at = c(0,cumsum(rep(6,48/6))),
     labels = c(0,as.character(rep(cumsum(rep(6,24/6)),times=2))))


piecewise_approx <- approx(x = 0:(24*14),y = cyclicIntensity(0:(24*14)),method = "constant",n=length(0:(24*14)))

NHPP_inversion_1st <- function(approx){
  t <- approx$x
  lambda <- approx$y
  j <- 1
  u <- rexp(n=1)
  Lambda <- rep(0,2)
  Lambda[1] <- 0 # lambda_{j}
  Lambda[2] <- Lambda[1] + 0.5*(lambda[j+1] + lambda[j])*(t[j+1] - t[j]) # lambda_{j+1}
  while(Lambda[2] < u){
    j <- j + 1
    if(j+1 > length(t)){
      return(NaN)
    }
    Lambda[1] <- Lambda[2]
    Lambda[2] <- Lambda[1] + 0.5*(lambda[j+1] + lambda[j])*(t[j+1] - t[j]) # lambda_{j+1}
  }
  tout <- t[j-1] + ((u - Lambda[1])/lambda[j-1]) #  Lambda_{j} <= u < Lambda_{j+1}
  return(tout)
}

times <- unlist(parallel::mclapply(X = 1:1e5,FUN = function(x){NHPP_inversion_1st(approx = piecewise_approx)}))

hist_out <- hist(times[times<48],breaks = 48)

# NHPP_rejection <- function(){
#   
# }

library(stocheulerABM)

out <- stocheulerABM::inhomPP_piecewiseconst(
  tvec = piecewise_approx$x,
  lambdavec = piecewise_approx$y,
  tmax = tail(piecewise_approx$x,1),
  first = F
)

times2 <- unlist(parallel::mclapply(X = 1:1e5,FUN = function(x){stocheulerABM::inhomPP_piecewiseconst(
  tvec = piecewise_approx$x,
  lambdavec = piecewise_approx$y,
  tmax = tail(piecewise_approx$y,1),
  first = T
)}))

hist(times2[times2<48],breaks = 48)

# NHPP_inversion_all <- function(approx,tmax){
#   a <- 0
#   u <- 0
#   j <- 1
#   Lambda <- c(0,0)
#   
#   t <- approx$x
#   lambda <- approx$y
#   
#   while(tail(a,1) < tmax){
#     
#     u <- u + rexp(n=1)
#     Lambda[j] <- 0 # lambda_{j}
#     Lambda[j+1] <- Lambda[1] + 0.5*(lambda[j+1] + lambda[j])*(t[j+1] - t[j]) # lambda_{j+1}
#     
#     
#   }
#   
#   
#   t <- approx$x
#   lambda <- approx$y
#   j <- 1
#   u <- rexp(n=1)
#   Lambda <- rep(0,2)
#   Lambda[1] <- 0 # lambda_{j}
#   Lambda[2] <- Lambda[1] + 0.5*(lambda[j+1] + lambda[j])*(t[j+1] - t[j]) # lambda_{j+1}
#   while(Lambda[2] < u){
#     j <- j + 1
#     if(j+1 > length(t)){
#       return(NaN)
#     }
#     Lambda[1] <- Lambda[2]
#     Lambda[2] <- Lambda[1] + 0.5*(lambda[j+1] + lambda[j])*(t[j+1] - t[j]) # lambda_{j+1}
#   }
#   tout <- t[j-1] + ((u - Lambda[1])/lambda[j-1]) #  Lambda_{j} <= u < Lambda_{j+1}
#   return(tout)
# }