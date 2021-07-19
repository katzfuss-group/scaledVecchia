
#######   methods for comparison: low-rank and hybrid laGP


### prediction fct for lowrank

predictions_LR <- function(fit,locs_pred,X_pred=matrix(1,nrow(locs_pred),1),m){
  
  y_obs = fit$y
  locs_obs = fit$locs
  X_obs = fit$X
  beta = fit$betahat
  covparms = fit$covparms
  covfun_name = fit$covfun_name
  n_obs <- nrow(locs_obs)
  n_pred <- nrow(locs_pred)
  n.all=n_obs+n_pred
  
  # get orderings
  ord1=1:n_obs
  ord2=1:n_pred
  
  # reorder stuff
  X_obs <- as.matrix(X_obs)
  X_pred <- as.matrix(X_pred)
  Xord_obs  <- X_obs[ord1,,drop=FALSE]
  Xord_pred <- X_pred[ord2,,drop=FALSE]
  yord_obs  <- y_obs[ord1]
  
  # put all coordinates together
  locs_all <- rbind( locs_obs[ord1,,drop=FALSE], locs_pred[ord2,,drop=FALSE] )
  inds1 <- 1:n_obs
  inds2 <- (n_obs+1):(n_obs+n_pred)
  
  # get nearest neighbor array
  NNarray_all <- find_ordered_nn_pred(locs_all,m,fix.first=n_obs)
  fm = NNarray_all[m+1,2:(m+1)]
  NNarray_all[(m+2):n.all,2:(m+1)]=matrix(rep(fm,n.all-m-1),byrow=TRUE,ncol=m)
  
  # get entries of Linv for obs locations and pred locations
  Linv_all <- GpGp::vecchia_Linv(covparms,covfun_name,locs_all,NNarray_all,n_obs+1)
  Linv_all[1:n_obs,1] <- 1.0
  
  ## prediction mean
  y_withzeros <- c(yord_obs - Xord_obs %*% beta, rep(0,n_pred) )
  v1 <- GpGp::Linv_mult(Linv_all, y_withzeros, NNarray_all )
  v1[inds1] <- 0
  v2 <- -GpGp::L_mult(Linv_all,v1,NNarray_all)
  
  condexp <- c(v2[inds2] + Xord_pred %*% beta)
  condexp[ord2] <- condexp
  return(condexp)
  
}




## implementation of hybrid laGP -- adapted from Bobby Gramacy
library(laGP)
library(tictoc)
global_local_laGP <- function(xtrain, ytrain, xtest, m, nth=4,ma=100)
{
  ## fixing a tiny nugget is very helpful on this problem
  g <- 1/10000000
  
  ## macro-scale analysis on a random subset of the data
  tic()
  N=nrow(xtrain)
  n <- min(1000,N)
  subs <- sample(1:N, n, replace = FALSE)
  d2 <- darg(list(mle = TRUE, max = ma), xtrain)
  gpsepi <- newGPsep(xtrain[subs, ], ytrain[subs], 
                     rep(d2$start, ncol(xtrain)), g=g, dK=TRUE)
  that <- mleGPsep(gpsepi, param="d", tmin=d2$min, tmax=d2$max, 
                   ab=d2$ab, maxit=200)
  # p <- predGPsep(gpsepi, xtest, lite=TRUE)
  deleteGPsep(gpsepi)
  temp=toc()
  time.train=temp$toc-temp$tic
  
  ## scale the inputs according to the macro-analysis lengthscales
  tic()
  scale <- sqrt(that$d)
  xs <- xtrain
  xpreds <- xtest
  for(j in 1:ncol(xs)) {
    xs[,j] <- xs[,j] / scale[j]
    xpreds[,j] <- xpreds[,j] / scale[j]
  }
  
  ## laGP analysis on scaled inputs
  out <- aGP(xs, ytrain, xpreds, d=list(start=1, max=20), start=min(5,m-1),
             end=m, g=g, omp.threads=nth, verb=0)
  temp=toc()
  time.pred=temp$toc-temp$tic  
  
  out$time.train=time.train
  out$time.pred=time.pred
  return(out)
}



## compute several scores assessing UQ
compute_scores=function(dat,mu,sig2,samples,alpha=.05){
  
  scores=rep(NA,times=6)
  
  # interval coverage
  scores[1]=int.cov(dat,mu,sig2,alpha=alpha)
  
  # interval width
  scores[2]=mean(qnorm(1-alpha/2)*sqrt(sig2))
  
  # interval score
  scores[3]=IS.normal(dat,mu,sig2,alpha=alpha)
  
  # log score
  scores[4]=-mean(dnorm(dat,mu,sqrt(sig2),log=TRUE))
  
  # CRPS
  scores[5]=mean(scoringRules::crps_norm(dat,mu,sqrt(sig2)))
  
  # energy score
  if(!missing(samples)) scores[6]=scoringRules::es_sample(dat,samples)
  
  return(scores)
  
}

## interval score for normal distribution
IS.normal=function(dat,mu,sig2,alpha=.05){
  q=qnorm(1-alpha/2)
  l=mu-q*sqrt(sig2)
  r=mu+q*sqrt(sig2)
  width=r-l
  outside.penalty=(l-dat)*(dat<l)+(dat-r)*(dat>r)
  mean(width+(2/alpha)*outside.penalty)
}

## interval coverage for normal distribution
int.cov=function(dat,mu,sig2,alpha=.05){
  q=qnorm(1-alpha/2)
  mean(abs(dat-mu)<=q*sqrt(sig2))
}