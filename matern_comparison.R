#####   comparison on simulated matern data  #####

# setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R')


range.imp=.05
range.minor=5


####### simulation for increasing m  #######


###  dimensions
d=10
n=5000
n.test=2000

### scale covariates according to importance
d.imp=2
d.minor=d-d.imp
ranges=c(rep(range.imp,d.imp),rep(range.minor,d.minor),rep(Inf,d-d.imp-d.minor))

### true parameter values
nu=3.5
variance=1
nugget=0
covfun=paste0("matern",nu*10,"_scaledim")
covparms=c(variance,ranges,nugget)


### simulation settings
reps=10
ms=c(5,10,15,20,30,50)


##### simulation for increasing m  ###
els.exact=rep(NA,reps)
els=array(dim=c(2,3,length(ms),reps))
mse=array(dim=c(3,length(ms),reps))
par.ests=array(dim=c(3,length(ms),reps,d))
for(rep in 1:reps){
  
  ### true cov and responses
  inputs=lhs::randomLHS(n,d)
  inputs.test=matrix(runif(n.test*d),n.test)
  Sigma.true=get(covfun)(covparms,rbind(inputs,inputs.test))
  Sigma.chol=t(chol(Sigma.true))
  y.all=as.numeric(Sigma.chol%*%rnorm(n+n.test))
  y=y.all[1:n]
  y.test=y.all[n+(1:n.test)]
  els.exact[rep]=-sum(log(diag(Sigma.chol[1:n,1:n])))-.5*n*log(2*pi)
  
  for(i.m in 1:length(ms)){
    
    m=ms[i.m]
    print(paste0('rep=',rep,', m=',m))
    
    
    ##### scaled vecchia
    scales=1/covparms[1+(1:d)]
    ord=order_maxmin_exact(t(t(inputs)*scales))
    inputs.ord=inputs[ord,]
    NNarray=find_ordered_nn(t(t(inputs.ord)*scales),m)
    els[1,1,i.m,rep]=vecchia_meanzero_loglik(covparms,covfun,
                                             rep(0,n),inputs.ord,NNarray)[[1]]
    
    ## fit parameters
    fit=fit_scaled(y,inputs,ms=c(4,m),nu=nu)
    par.ests[1,i.m,rep,]=fit$covparms[1+(1:d)]
    NNarray=find_ordered_nn(t(t(fit$locs)*1/fit$covparms[1+(1:d)]),m)
    els[2,1,i.m,rep]=vecchia_meanzero_loglik(fit$covparms,covfun,
                                             rep(0,n),fit$locs,NNarray)[[1]]
    
    # predictions
    preds=predictions_scaled(fit=fit,locs_pred=inputs.test,m=2*m)
    mse[1,i.m,reps]=mean((preds-y.test)^2)
    
    
    ##### vecchia w/o scaling
    ord=order_maxmin_exact(inputs)
    inputs.ord=inputs[ord,]
    NNarray=find_ordered_nn(inputs.ord,m)
    els[1,2,i.m,rep]=vecchia_meanzero_loglik(covparms,covfun,
                                             rep(0,n),inputs.ord,NNarray)[[1]]
    
    ## fit parameters
    y.ord=y[ord]
    start_parms=c(var(y),rep(.2,d),0)
    fit=fit_model(y.ord,inputs.ord,NNarray=NNarray,m_seq=m,
                  start_parms=start_parms,covfun_name=covfun,silent=TRUE,
                  reorder=FALSE,fixed_parms=2+d)
    par.ests[2,i.m,rep,]=fit$covparms[1+(1:d)]
    els[2,2,i.m,rep]=vecchia_meanzero_loglik(fit$covparms,covfun,
                                             rep(0,n),inputs.ord,NNarray)[[1]]
    
    # prediction
    preds=predictions(fit=fit,locs_pred=inputs.test,
                      X_pred=matrix(0,n.test,1),m=2*m)
    mse[2,i.m,reps]=mean((preds-y.test)^2)
    
    
    ##### low rank
    fm = NNarray[m+1,2:(m+1)]
    NNarray[(m+2):n,2:(m+1)]=matrix(rep(fm,n-m-1),byrow=TRUE,ncol = m)
    els[1,3,i.m,rep]=vecchia_meanzero_loglik(covparms,covfun,
                                             rep(0,n),inputs.ord,NNarray)[[1]]
    
    ## fit parameters    
    fit=fit_model(y.ord,inputs.ord,NNarray=NNarray,m_seq=m,
                  start_parms=start_parms,covfun_name=covfun,silent=TRUE,
                  reorder=FALSE,fixed_parms=2+d)
    par.ests[3,i.m,rep,]=fit$covparms[1+(1:d)]
    els[2,3,i.m,rep]=vecchia_meanzero_loglik(fit$covparms,covfun,
                                             rep(0,n),inputs.ord,NNarray)[[1]]
    
    # prediction
    preds=predictions_LR(fit,inputs.test,m=2*m,X_pred=matrix(0,n.test,1))
    mse[3,i.m,reps]=mean((preds-y.test)^2)
    
    
  }
}

# save(ms,d,els,els.exact,par.ests,mse,n,file='results/maternSimPred.RData')




################   plot results   ################


# load(file='results/maternSimPred.RData')


##### log score

els.avg=apply(els,1:3,mean,na.rm=TRUE)
kls=mean(els.exact,na.rm=TRUE)-els.avg
tp=c('fixed','est')

y.range=range(kls,na.rm=TRUE)
y.ticks=c(1000,2000,5000,10000,15000,20000)

for(k in 1:2){
  pdf(file=paste0('plots/matern_',tp[k],'_ylog.pdf'),width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
  matplot(ms,log10(t(kls[k,,])),type='l',ylab='dLS /1000',xlab='m',lwd=2,
          axes=FALSE,ylim=log10(y.range))
  axis(1,at=seq(10,50,by=10),labels=seq(10,50,by=10),lwd=0,lwd.ticks=1)
  axis(2,at=log10(y.ticks),labels=y.ticks/1000,lwd=0,lwd.ticks=1)
  box()
  legend('right',c('SVecchia','Vecchia','LowRank'),col=1:3,lty=1:3,lwd=2,bg='white')
  dev.off()
}


##### RMSE
rmse=sqrt(apply(mse,1:2,mean,na.rm=TRUE))
y.ticks=c(.01,.02,.05,.1,.2,.5,1,2)
pdf(file=paste0('plots/matern_rmse.pdf'),width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
matplot(ms,log10(t(rmse)),type='l',ylab='RMSE',xlab='m',lwd=2,axes=FALSE)
axis(1,at=seq(10,50,by=10),labels=seq(10,50,by=10),lwd=0,lwd.ticks=1)
axis(2,at=log10(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
box()
dev.off()
