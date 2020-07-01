# transient probability distributions
rm(list=ls());gc()
library(MultiBD)
library(matrixStats)

# --------------------------------------------------------------------------------
#   parameters
# --------------------------------------------------------------------------------

#Choose initial S and I population here
S <- 140
I <- 10
beta = .5/(S+I)
gamma = .1
N <- 1000 #number of MC realizations: increase N in practice
t.end <- .5 #time interval length


# --------------------------------------------------------------------------------
#   Branching process approximation
# --------------------------------------------------------------------------------

# expression for the phi1 cross-partials
A <- function(t, m, n, k, j, beta, gam){
  if(j > m-k){return(0)} else{
    return( factorial(m) * exp(-k*beta*n*t) *
              ( 1 - beta*n*exp(-gam*t)/(beta*n-gam) - exp(-beta*n*t) *
                  (1 - beta*n/(beta*n-gam)  ) )^(m-k-j) *
              ( beta*n*(exp(-gam*t) - exp(-beta*n*t)) /
                  (beta*n - gam) )^(j) / factorial(m-k-j) )
  }
}

# expression for the phi2 partials
B <- function(t, n, j, gam){
  if(j>n){return(0)} else{
    return( factorial(n)* (1-exp(-gam*t))^(n-j) *
              exp(-gam*j*t)/factorial(n-j) )
  }
}

# remove if statements from A(), B(), so that vector arguments work
# input j must be greater than m-k
A.vectorized <- function(t, m, n, k, j, beta, gam){
  return( factorial(m) * exp(-k*beta*n*t) *
            ( 1 - beta*n*exp(-gam*t)/(beta*n-gam) - exp(-beta*n*t) *
                (1 - beta*n/(beta*n-gam)  ) )^(m-k-j) *
            ( beta*n*(exp(-gam*t) - exp(-beta*n*t)) /
                (beta*n - gam) )^(j) / factorial(m-k-j) )
}

#log version of previous, to handle larger numbers
A.vec.log <- function(t, m, n, k, j, beta, gam){
  return( lgamma(m+1) - lgamma(m-k-j+1) - (k*beta*n*t) +
            (m-k-j)* log( ( 1 - beta*n*exp(-gam*t)/(beta*n-gam) -
                              exp(-beta*n*t) *
                              (1 - beta*n/(beta*n-gam)  ) ) )
          + j * log( ( beta*n*(exp(-gam*t) - exp(-beta*n*t)) /
                         (beta*n - gam) ) )
  )
}

B.vectorized <- function(t, n, j, gam){
  return( factorial(n)* (1-exp(-gam*t))^(n-j) * exp(-gam*j*t) /
            factorial(n-j) )
}

B.vec.log <- function(t, n, j, gam){
  return( lgamma(n+1) - lgamma(n-j+1) + (n-j)*log(1 - exp(-gam*t))
          - gam*j*t  )
}

# note: of course k can never be greater than m
# uses log versions to handle large populations
# put k and l first for use with outer()
TransProb_mnkl <- function(k,l, t, m, n, beta, gam){
  if(k>m){return(0)}
  logAA <- logBB <- rep(0,l+1)
  aj <- rev(seq(0,min(m-k,l)))
  bj <- seq(0, min(n,l) )
  c <- lgamma(l+1) - lgamma(seq(0,l) + 1) - lgamma( l+1 - seq(0,l))
  logAA[( l+1-min(m-k,l) ):(l+1)] <- A.vec.log(t, m, n, k, aj, beta, gam)
  logBB[ 1: (min(n,l)+1) ] <- B.vec.log(t, n, bj, gam)
  #  print(AA)
  #  print(BB)
  #  print(c)
  term <- c + logAA + logBB
  return( exp( matrixStats::logSumExp(term) - lgamma(k+1) - lgamma(l+1) )  )
}

#vectorize the previous for use with outer
vectorized <- Vectorize(TransProb_mnkl, vectorize.args = c('k','l'))

# S_0, I_0 are m,n: fixes those and returns a grid of transitions
# over different end states for comparison with getTransProbs()
#vectorized version
getTransProbsClosed <- function(t, gridLength, beta, gam, S0, I0){
  return( outer(0:gridLength, 0:gridLength, vectorized,
                t = t, beta = beta, gam = gam, m = S0, n = I0))
}


# --------------------------------------------------------------------------------
#   Gillespie SIR
# --------------------------------------------------------------------------------

# simple simulation of true SIR model over a time interval
simulateSIR <- function(t.end, S, I, beta, gam, maxEvents = 99999999){
  t.cur <- 0
  for(i in 1:maxEvents){
    if( S<0 || I <0 ){
      print("Negative population? Error")
      return(-99) #error code
    }
    if( S == 0 || I ==0){
      #print("S or I is zero, end epidemic")
      return( c(S,I) )   #end epidemic
    }
    infectRate <- S*I*beta
    recovRate <- I*gam
    rates <- c(infectRate, recovRate)

    t.next <- rexp(1, sum(rates)) #time until next event
    t.cur <- t.cur+t.next
    if(t.cur > t.end){            #end of simulation period
      return( c(S,I) )
    }
    #sample the type of next event
    decision <- rbinom(1, 1, infectRate/sum(rates))
    if(decision == 1) {   #infection
      S <- S - 1; I <- I + 1
    } else {             #recovery
      I <- I - 1
    }
  }
  return(-99) #error code for testing
}

#run simulation once with error catch
sim.once <- function(t.end, S, I, beta, gam, maxEvents = 99999999){
  res = -99 # error catch
  while(res[1] == -99){
    res <- simulateSIR(t.end, S, I, beta, gam, maxEvents)  }
  return(res)
}

getTrans.MC <- function(N, t.end, S, I, beta, gam){
  result <- replicate(N, sim.once(t.end, S, I, beta, gam))
  #make big enough to account for all events: count end states
  trans.count <- matrix(0, S+I,S+I)
  for(i in 1:N){
    id <- result[,i]+1
    #indices in the resulting transition count: ie if you end at (1,1),
    # you add a count to the (2,2) entry of the count matrix, etc
    trans.count[id[1], id[2]] = trans.count[id[1], id[2]] + 1
  }
  tpm <- trans.count/sum(trans.count)
  return(tpm)
}


## ----multiBDsetup-----------------------------------------------------------------------
# tList  <- c(.1, .2, .25, .3 ,.35, .4, .5, .6, .7, .8, .9, 1)
tList <- c(0.1, 0.5, 1)
gridLength = 128
a0 = 110 # S_0
b0 = 15 # I_0
A = 0
B = gridLength - 1
alpha = 3.2 #3.2 #this is death rate
beta = .025 #.019 #this is transition or infection rates
nSim = 1e3 #number of MC simulations

brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){alpha*b}
trans=function(a,b){beta*a*b}

#indexed by time, type of computation, and dimensions of the tpm
tpmArray <- array(NA, dim = c(length(tList), 3, 52, 25 )) #store a subset of transition probabilities

for(i in 1:length(tList)){
  t.end <- tList[i]

  tpm.Closed <- getTransProbsClosed(t.end, gridLength, beta, alpha, a0, b0)
  tpm1 = tpm.Closed[1:(a0+1),] #using 2-type branching approximation

  #using continued fractions via MultiBD
  tpm2 <- dbd_prob(t.end, a0, b0, drates1, brates2, drates2, trans, a=A, B)
  #MC simulation "ground truth"
  tpm.MC <- getTrans.MC(nSim, t.end, a0, b0, beta, alpha)
  tpm3 <- tpm.MC[1:(a0+1), ]

  #store subset of matrices containing about 99 percent of the mass:
  tpmArray[i,1,,] <- tpm1[60:(a0+1),1:25]
  tpmArray[i,2,,] <- tpm2[60:(a0+1),1:25]
  tpmArray[i,3,,] <- tpm3[60:(a0+1),1:25]
}


## ---------------------------------------------------------------------------------------
#for example, look at the ones with t.end = .5
ix <- which(tList == 0.5)
small1 <- tpmArray[ix,1,,]
small2 <- tpmArray[ix,2,,]
small3 <- tpmArray[ix,3,,]

#they comprise most of transition probability mass:
sum(small1); sum(small2); sum(small3)

# mean errors
mean(abs(small1- small3 ) ) #2-type vs MC
mean(abs(small2 - small3) ) #Continued Frac vs MC

# scaled heatmap images to compare tpm visually
par(mfrow=(c(3,1)))
image(small1, main = "Two-type branching approximation")
image(small2, main = "Continued Fraction expansion")
image(small3, main = "Monte Carlo estimates")

## ---------------------------------------------------------------------------------------
par(mfrow=(c(3,1)))
image(tpmArray[3,1,,], main = "Two-type branching approximation")
image(tpmArray[3,2,,], main = "Continued Fraction expansion")
image(tpmArray[3,3,,], main = "Monte Carlo estimates")

## ----transitionCompare------------------------------------------------------------------
library(plotrix)
inds <- t(which(tpmArray[3,2,,] >= sort(tpmArray[3,2,,], decreasing=T)[16],
                arr.ind=TRUE))
#ind1 <- sample(52,25, replace=T); ind2 <- sample(25,25,replace=T)
par(mfrow = c(4,4), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
for(i in 1:16){
  plot(tList, tpmArray[,2,inds[1,i], inds[2,i] ], pch = 17, col = 'red',
       ylim = c(0,max(tpmArray[,,inds[1,i], inds[2,i]])),
       yaxt = 'n', xlab = "dt")
  MCp <- tpmArray[,3,inds[1,i], inds[2,i] ] #MC prob
  plotCI(tList, MCp, pch = 4, col = 'green',
         ui = MCp+1.96*sqrt(MCp*(1-MCp)/nSim),
         li <- MCp-1.96*sqrt(MCp*(1-MCp)/nSim), add = TRUE)
  points(tList, tpmArray[,1,inds[1,i], inds[2,i] ],
         col='purple', pch = 16)
}
