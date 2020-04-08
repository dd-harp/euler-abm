rm(list=ls());gc()
library(parallel)

# --------------------------------------------------------------------------------
#   test rejection for constant FOI
# --------------------------------------------------------------------------------

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

# n <- 2e6
# lambda <- 1/20
# dt <- 0.5
#
# out_05 <- mclapply(X = 1:n,FUN = function(i){simmer(tend = 1e4,rate = lambda,dt = dt)})
# out_05 <- do.call(rbind,out_05)
#
# out_1 <- mclapply(X = 1:n,FUN = function(i){simmer(tend = 1e4,rate = lambda,dt = 1)})
# out_1 <- do.call(rbind,out_1)
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



# --------------------------------------------------------------------------------
# NHPP: piecewise constant intensity
# --------------------------------------------------------------------------------

# findInterval()

x <- seq(0,48,by=0.01)
cyclicIntensity <- function(t,fudge=T){
  # t <- t %% 24
  if(fudge){
    (1/8)*sin((t - 6)*2*pi/24)+((1/8)+.Machine$double.eps)
  } else {
    (1/8)*sin((t - 6)*2*pi/24)+(1/8)
  }
}

cyclicApprox <- approx(x = 0:23,y = cyclicIntensity(0:23),method = "constant",n=length(0:23))

# Intensity
plot(x,cyclicIntensity(x),type="l",ylim = c(0,0.25),
     main = "Inhomogeneous Intensity Function",xaxt="n",
     xlab = "Time (hours)",ylab = "Intensity",col="firebrick3",lwd=2.15)

for(i in 0:48){

  j <- (i %% 24) + 1

  segments(x0 = i,x1 = i+1,
           y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j],
           col = "darkorchid3",lwd=1.85)
  points(x=i,y=cyclicApprox$y[j],pch=16,col="darkorchid3",cex=0.85)
  points(x=i+1,y=cyclicApprox$y[j],pch=1,col="darkorchid3",cex=0.85)
  if(i != 48){
    segments(x0 = i+1,x1 = i+1,y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j+1],
             lty = 2,col="darkorchid3")
  }

}
axis(side = 1,at = c(0,cumsum(rep(6,48/6))),
     labels = c(0,as.character(rep(cumsum(rep(6,24/6)),times=2))))


# time to 1st event
lambda <- cyclicIntensity(t = x,fudge = FALSE)
lambda_discrete <- approx(x = x,y = cyclicIntensity(x,F),method = "constant",n = length(x))

IntegratedCyclicIntensity <- function(t){
  (1/8) * (t - ( (12 * sin( (pi*t) / 12 ) ) / pi ))
}

Lambda <- IntegratedCyclicIntensity(t = x)
Lambda_discrete <- cumsum(lambda_discrete$y*0.01)

pdf_1st <- lambda * exp(-Lambda)
pmf_1st <- lambda_discrete$y * exp(-Lambda_discrete)

# plot(x,pdf_1st,type="l",
#      xlab = "Time (hours)",ylab = "Density",
#      main = "Time to First Event",col = "firebrick3",lwd=1.5)
# lines(x,pmf_1st,col="darkorchid3",lwd=1.5,lty=1)




# --------------------------------------------------------------------------------
# NHPP: 1st event time in R by inversion
# --------------------------------------------------------------------------------

# for 14 days
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
  tmax = tail(piecewise_approx$x,1),
  first = T
)}))

hist(times2[times2<48],breaks = 48)

times_lambda <- unlist(parallel::mclapply(X = 1:1e5,FUN = function(x){stocheulerABM::inhomPP_piecewiseconst(
  tvec = lambda_discrete$x,
  lambdavec = lambda_discrete$y,
  tmax = tail(lambda_discrete$x,1),
  first = T
)}))

times_lambda_rej <- parallel::mclapply(X = 1:1e5,FUN = function(x){stocheulerABM::inhomPP_piecewiseconst_reject(
  tvec = lambda_discrete$x,
  lambdavec = lambda_discrete$y
)})

par(mfrow=c(1,2))
plot(x,pdf_1st,type="l",
     xlab = "Time (hours)",ylab = "Density",
     main = "Time to First Event (Sampling Integrated Hazard)",col = "firebrick3",lwd=2.15)
hist(times_lambda,breaks=48,probability = T,add=T,col = adjustcolor("steelblue",alpha.f = 0.5),border=grey(0.5,0.5))

plot(x,pdf_1st,type="l",
     xlab = "Time (hours)",ylab = "Density",
     main = "Time to First Event (Rejection Sampling)",col = "firebrick3",lwd=2.15)
hist(sapply(times_lambda_rej,function(x){x[["time"]]}),breaks=48,probability = T,add=T,col = adjustcolor("darkorchid3",alpha.f = 0.5),border=grey(0.5,0.5))
par(mfrow=c(1,1))
