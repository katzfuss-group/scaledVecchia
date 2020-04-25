
####################################################################
##############  collection of test functions   #####################
####################################################################


##  all functions from Virtual Library of Simulation Experiments
#
# Authors: Sonja Surjanovic, Simon Fraser University
#          Derek Bingham, Simon Fraser University
# Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
#
# Copyright 2013. Derek Bingham, Simon Fraser University.
#
# THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
# FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
# derivative works, such modified software should be clearly marked.
# Additionally, this program is free software; you can redistribute it 
# and/or modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation; version 2.0 of the License. 
# Accordingly, this program is distributed in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# For function details and reference information, see:
# http://www.sfu.ca/~ssurjano/


testfuns=list()


###  function that scales functions to unit input

unit.scale=function(fun,ranges){
  scaled.fun=function(x){
    xx = x*(ranges[,2]-ranges[,1])+ranges[,1]
    fun(xx)
  }
  return(scaled.fun)
}



######################   borehole arm   ############################

###   borehole function  ###

borehole <- function(xx){
  # OUTPUT AND INPUT:
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}


bore.ranges=matrix(c(.05,.15,
                100,50000,
                63070,115600,
                990,1110,
                63.1,116,
                700,820,
                1120,1680,
                9855,12045),ncol=2,byrow=TRUE)

bore=list(fun=unit.scale(borehole,bore.ranges),d=8,name='borehole')
testfuns$bore=bore





######################   robot arm   ############################

robotfun <- function(xx)
{
  # OUTPUT AND INPUTS:
  # y = distance from the end of the arm to the origin
  # xx = c(theta1, theta2, theta3, theta4, L1, L2, L3, L4)
  
  theta <- xx[1:4]
  L     <- xx[5:8]
  
  thetamat <- matrix(rep(theta,times=4), 4, 4, byrow=TRUE)
  thetamatlow <- thetamat
  thetamatlow[upper.tri(thetamatlow)] <- 0
  sumtheta <- rowSums(thetamatlow)
  
  u <- sum(L*cos(sumtheta))
  v <- sum(L*sin(sumtheta))
  
  y <- (u^2 + v^2)^(0.5)
  return(y)
}


robot.ranges=matrix(c(rep(c(0,2*pi),4),rep(c(0,1),4)),ncol=2,byrow=TRUE)

robot=list(fun=unit.scale(robotfun,robot.ranges),d=8,name='robot arm')
testfuns$robot=robot






######################   Piston    ############################

pistonfun <- function(xx)
{
  # OUTPUT AND INPUT:
  # C = cycle time
  # xx = c(M, S, V0, k, P0, Ta, T0)

  M  <- xx[1]
  S  <- xx[2]
  V0 <- xx[3]
  k  <- xx[4]
  P0 <- xx[5]
  Ta <- xx[6]
  T0 <- xx[7]
  
  Aterm1 <- P0 * S
  Aterm2 <- 19.62 * M
  Aterm3 <- -k*V0 / S
  A <- Aterm1 + Aterm2 + Aterm3
  
  Vfact1 <- S / (2*k)
  Vfact2 <- sqrt(A^2 + 4*k*(P0*V0/T0)*Ta)
  V <- Vfact1 * (Vfact2 - A)
  
  fact1 <- M
  fact2 <- k + (S^2)*(P0*V0/T0)*(Ta/(V^2))
  
  C <- 2 * pi * sqrt(fact1/fact2)
  return(C)
}

piston.ranges=matrix(c(30, 60,
                       0.005, 0.020,
                       0.002, 0.010,
                       1000, 5000,
                       90000, 110000,
                       290, 296,
                       340, 360),ncol=2,byrow=TRUE)

piston=list(fun=unit.scale(pistonfun,piston.ranges),d=7,name='piston')
testfuns$piston=piston



