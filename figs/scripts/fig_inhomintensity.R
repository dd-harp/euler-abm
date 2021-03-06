# --------------------------------------------------------------------------------
#
#   Figure: show inhomogeneous intensity function and piecewise approx
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()

mar_default <- c(5.1, 4.1, 4.1, 2.1)
par(mar=c(4.5,4.5,2,2))

dt <- 0.01
x <- seq(0,48,by=dt)

# the inhomogeneous intensity function
cyclicIntensity <- function(t,fudge=T){
  if(fudge){
    (1/8)*sin((t - 6)*2*pi/24)+((1/8)+.Machine$double.eps)
  } else {
    (1/8)*sin((t - 6)*2*pi/24)+(1/8)
  }
}

cyclicApprox <- approx(x = 0:23,y = cyclicIntensity(0:23),method = "constant",n=length(0:23))


# # Intensity
# plot(x,cyclicIntensity(x),type="l",ylim = c(0,0.25),xaxt="n",
#      xlab = "Time (hours)",ylab = expression(paste(lambda,"(t)")),
#      col="firebrick3",lwd=2.15,cex.lab=1.45)
# 
# for(i in 0:48){
# 
#   j <- (i %% 24) + 1
# 
#   segments(x0 = i,x1 = i+1,
#            y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j],
#            col = "darkorchid3",lwd=1.85)
#   points(x=i,y=cyclicApprox$y[j],pch=16,col="darkorchid3",cex=0.85)
#   points(x=i+1,y=cyclicApprox$y[j],pch=1,col="darkorchid3",cex=0.85)
#   if(i != 48){
#     segments(x0 = i+1,x1 = i+1,y0 = cyclicApprox$y[j],y1 = cyclicApprox$y[j+1],
#              lty = 2,col="darkorchid3")
#   }
# 
# }
# axis(side = 1,at = c(0,cumsum(rep(6,48/6))),
#      labels = c(0,as.character(rep(cumsum(rep(6,24/6)),times=2))))

library(ggplot2)
library(data.table)

cyclicApprox <- approx(x = 0:48,y = cyclicIntensity(0:48),method = "constant",n=length(0:48))

plotdat_c <- data.table(time=x,intensity=cyclicIntensity(x),type="cont")
plotdat_s <- data.table(time=cyclicApprox$x,intensity=cyclicApprox$y,type="step")


plot_int <- ggplot(data = plotdat_c) +
  geom_line(aes(x=time,y=intensity),color="firebrick3") +
  geom_segment(data = plotdat_s,aes(x=time,y=intensity,xend=shift(time,type = "lead"),yend=intensity),color="darkorchid3") +
  geom_segment(data = plotdat_s,aes(x=shift(time,type="lead"),y=intensity,xend=shift(time,type = "lead"),yend=shift(intensity,type = "lead")),color="darkorchid3",linetype=2) +
  geom_point(data = plotdat_s,aes(x=time,y=intensity),color="darkorchid3",shape=16,cex=2) +
  geom_point(data = plotdat_s,aes(x=shift(time,type = "lead"),y=intensity),color="darkorchid3",shape=1,cex=2) +
  xlab("Time (Hours)") + ylab(expression(paste(lambda,"(t)"))) + labs(title="Cyclic Intensity Function") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(2.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(2)),
    axis.text = element_text(size = rel(1.5))
  )

# save as 8 x 10 landscape PDF: inhomintensity.pdf


# IntegratedCyclicIntensity <- function(t){
#   (1/8) * (t - ( (12 * sin( (pi*t) / 12 ) ) / pi ))
# }
#
# lambda <- cyclicIntensity(t = x,fudge = FALSE)
# lambda_discrete <- approx(x = x,y = cyclicIntensity(x,F),method = "constant",n = length(x)/50)
#
# Lambda <- IntegratedCyclicIntensity(t = x)
# Lambda_discrete <- cumsum(lambda_discrete$y*dt*50)
#
# plot(x,Lambda,type="l",col="firebrick3",lwd=2.5,xaxt="n")
#
#
# for(i in 1:(length(lambda_discrete$x)-1)){
#   segments(x0 = lambda_discrete$x[i],x1 = lambda_discrete$x[i+1],
#            y0 = Lambda_discrete[i],y1 = Lambda_discrete[i],
#            col = "darkorchid3",lwd=1.85)
#   points(x=lambda_discrete$x[i],y=Lambda_discrete[i],pch=16,col="darkorchid3",cex=0.85)
#   points(x=lambda_discrete$x[i+1],y=Lambda_discrete[i],pch=1,col="darkorchid3",cex=0.85)
#   if(i != length(lambda_discrete$x)-1){
#     segments(x0 = lambda_discrete$x[i+1],x1 = lambda_discrete$x[i+1],y0 = Lambda_discrete[i],y1 = Lambda_discrete[i+1],
#              lty = 2,col="darkorchid3")
#   }
# }
#
# axis(side = 1,at = c(0,cumsum(rep(6,48/6))),
#      labels = c(0,as.character(rep(cumsum(rep(6,24/6)),times=2))))
#
# par(mfrow=c(1,1))
