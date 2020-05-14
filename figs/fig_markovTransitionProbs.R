# --------------------------------------------------------------------------------
#
#   Figure: compare final epidemic sizes in Markov SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   May 2020
#
# --------------------------------------------------------------------------------


library(MultiBD)
tList  <- c(.1, .2, .25, .3 ,.35, .4, .5, .6, .7, .8, .9, 1)
gridLength = 128
a0 = 110 # S_0
b0 = 15 # I_0
A = 0
B = gridLength - 1

# gamma
alpha = 3.2 #3.2 #this is death rate

# beta
beta = .025 #.019 #this is transition or infection rates
nSim = 4000 #number of MC simulations

brates1=function(a,b){0}
drates1=function(a,b){0}
brates2=function(a,b){0}
drates2=function(a,b){alpha*b}
trans=function(a,b){beta*a*b}





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













#indexed by time, type of computation, and dimensions of the tpm
tpmArray <- array(NA, dim = c(length(tList), 3, 52, 25 )
) #store a subset of transition probabilities
time1 <- rep(0, length(tList)); time2 <- rep(0, length(tList))

for(i in 1:length(tList)){
  t.end <- tList[i]
  time1[i] <- system.time(
    tpm.Closed <- getTransProbsClosed(t.end, gridLength,
                                      beta, alpha, a0, b0) )
  tpm1 = tpm.Closed[1:(a0+1),] #using 2-type branching approximation

  #using continued fractions via MultiBD
  time2[i] <- system.time(
    tpm2 <- dbd_prob(t.end, a0, b0, drates1, brates2, drates2, trans,
                     a=A, B))#, computeMode=2))
  #MC simulation "ground truth"
  tpm.MC <- getTrans.MC(nSim, t.end, a0, b0, beta, alpha)
  tpm3 <- tpm.MC[1:(a0+1), ]

  #store subset of matrices containing about 99 percent of the mass:
  tpmArray[i,1,,] <- tpm1[60:(a0+1),1:25]
  tpmArray[i,2,,] <- tpm2[60:(a0+1),1:25]
  tpmArray[i,3,,] <- tpm3[60:(a0+1),1:25]
}
save(tpmArray, time1, time2, file = "../inst/vignetteCache/tpm_array.RData")
