rm(list=ls());gc()
library(parallel)

# --------------------------------------------------------------------------------
# NHPP: piecewise constant intensity
# --------------------------------------------------------------------------------


x <- seq(0,48,by=0.01)
cyclicIntensity <- function(t,fudge=T){
  # t <- t %% 24
  if(fudge){
    (1/8)*sin((t - 6)*2*pi/24)+((1/8)+.Machine$double.eps)
  } else {
    (1/8)*sin((t - 6)*2*pi/24)+(1/8)
  }
}

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



# --------------------------------------------------------------------------------
# NHPP: 1st event time in R by inversion
# --------------------------------------------------------------------------------

n <- 1e6

times_lambda <- unlist(parallel::mclapply(X = 1:n,FUN = function(x){stocheulerABM::inhomPP_piecewiseconst(
  tvec = lambda_discrete$x,
  lambdavec = lambda_discrete$y,
  tmax = tail(lambda_discrete$x,1),
  first = T
)}))

times_lambda_rej <- parallel::mclapply(X = 1:n,FUN = function(x){stocheulerABM::inhomPP_piecewiseconst_reject(
  tvec = lambda_discrete$x,
  lambdavec = lambda_discrete$y
)})

par(mfrow=c(1,2))
plot(x,pdf_1st,type="l",
     xlab = "Time (hours)",ylab = "Density",
     main = "Time to First Event (Sampling Integrated Hazard)",
     col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25)
hist(times_lambda,breaks=48,probability = T,add=T,col = adjustcolor("steelblue",alpha.f = 0.5),border=grey(0.5,0.5))

plot(x,pdf_1st,type="l",
     xlab = "Time (hours)",ylab = "Density",
     main = "Time to First Event (Rejection Sampling)",
     col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25)
hist(sapply(times_lambda_rej,function(x){x[["time"]]}),breaks=48,probability = T,add=T,col = adjustcolor("darkorchid3",alpha.f = 0.5),border=grey(0.5,0.5))
par(mfrow=c(1,1))
