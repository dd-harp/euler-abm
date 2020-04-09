# --------------------------------------------------------------------------------
#
#   Markovian SIR
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
# --------------------------------------------------------------------------------

rm(list=ls());gc()
library(MultiBD)

tList <- c(.1, .2, .25, .3 ,.35, .4, .5, .6, .7, .8, .9, 1)
gridLength = 128
S0 = 110 # S_0
I0 = 15 # I_0

A=0
B = gridLength - 1
gamma = 3.2 #3.2 #this is death rate
beta = .025 #.019 #this is transition or infection rates nSim = 4000 #number of MC simulations
brates1=function(S,I){0}
drates1=function(S,I){0}
brates2=function(S,I){0}
drates2=function(S,I){gamma*I}
trans=function(S,I){beta*S*I}

tmax <- 100

system.time(SIRexact <- MultiBD::dbd_prob(
  t=tmax,
  a0=S0,
  b0=I0,
  mu1=drates1,
  lambda2=brates2,
  mu2=drates2,
  gamma=trans,
  a=0,
  B=B,
  computeMode=1
))
