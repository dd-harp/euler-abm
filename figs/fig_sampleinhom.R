# --------------------------------------------------------------------------------
#
#   Figure: compare sampling of NHPP first event times via inversion vs. accept-reject
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
# --------------------------------------------------------------------------------

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

# par(mfrow=c(1,2))
# plot(x,pdf_1st,type="l",
#      xlab = "Time (hours)",ylab = "Density",
#      main = "Time to First Event (Sampling Integrated Hazard)",
#      col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25)
# hist(times_lambda,breaks=48,probability = T,add=T,col = adjustcolor("steelblue",alpha.f = 0.5),border=grey(0.5,0.5))
# 
# plot(x,pdf_1st,type="l",
#      xlab = "Time (hours)",ylab = "Density",
#      main = "Time to First Event (Rejection Sampling)",
#      col = "firebrick3",lwd=2.15,cex.lab=1.45,cex.main=1.25)
# hist(sapply(times_lambda_rej,function(x){x[["time"]]}),breaks=48,probability = T,add=T,col = adjustcolor("darkorchid3",alpha.f = 0.5),border=grey(0.5,0.5))
# par(mfrow=c(1,1))

library(ggplot2)
library(gridExtra)
library(data.table)

dat_int <- data.table(time=x,density=pdf_1st)

plot_int <- ggplot(data = dat_int) +
  geom_line(aes(x=time,y=density),color="firebrick3",lwd=1.5) +
  geom_histogram(data=data.frame(x=times_lambda),aes(x=x,y=after_stat(density)),color=grey(0.5,0.5),fill="steelblue",bins=48,alpha=0.5) +
  xlab("Time (Hours)") + ylab("Density") + labs(title="A. Time to First Event (Sampling Integrated Hazard)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(2)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

plot_rej <- ggplot(data = dat_int) +
  geom_line(aes(x=time,y=density),color="firebrick3",lwd=1.5) +
  geom_histogram(data=data.frame(x=sapply(times_lambda_rej,function(x){x[["time"]]})),aes(x=x,y=after_stat(density)),color=grey(0.5,0.5),fill="darkorchid3",bins=48,alpha=0.5) +
  xlab("Time (Hours)") + ylab("Density") + labs(title="B. Time to First Event (Rejection Sampling)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(2)),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5))
  )

plot_both <- grid.arrange(plot_int,plot_rej,nrow=1)
# save as 8 x 16 landscape PDF: sampleinhom.pdf