#####   extensive comparison on the borehole function   ####

# setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R')


### load borehole function
source('code/test_functions.R')
bore=testfuns[[1]]$fun
d=testfuns[[1]]$d




##########   comparison for different m and n   #########

ns=c(100,200,400,1000,4000,10000)
ms=c(5,10,15,20,30,50)
reps=10
nu.fixed=3.5
covfun=paste0("matern",nu.fixed*10,"_scaledim")

### initial parameter values
fit.sv=fit.v=fit.lr=list()
var.ini=2000
ra.ini=rep(.2,d)
fit.sv$covparms=fit.v$covparms=fit.lr$covparms=c(var.ini,ra.ini,nu.fixed,0)

### start comparison
mse=array(dim=c(reps,4,length(ms),length(ns)))
par.ests=array(dim=c(reps,3,length(ms),length(ns),d+1))

for(i.n in 1:length(ns)){
  
  n=ns[i.n]

  for(rep in 1:reps){
    
    ### training data
    inputs=lhs::randomLHS(n,d)
    y=apply(inputs,1,bore)
    
    ### test data
    n.test=2000
    inputs.test=matrix(runif(n.test*d),n.test)
    y.test=apply(inputs.test,1,bore)
    
    for(i.m in 1:length(ms)){
      
      m=ms[i.m]
      print(paste0('n=',n,', rep=',rep,', m=',m))
      m.seq=m #c(min(round(m/2),10),m)
      
      ### fit using iterative scaling
      fit.sv=fit_scaled(y,inputs,ms=m.seq,nu=nu.fixed,select=1e3)
      par.ests[rep,1,i.m,i.n,]=fit.sv$covparms[1:(d+1)]
      preds=predictions_scaled(fit=fit.sv,locs_pred=inputs.test,m=2*m)
      mse[rep,1,i.m,i.n]=mean((preds-y.test)^2)
      
      ### vecchia w/o scaling
      start.parms=c(var(y),rep(.2,d),0)
      ord=order_maxmin_exact(inputs)
      inputs.ord=inputs[ord,]
      y.ord=y[ord]
      NNarray=find_ordered_nn(inputs.ord,m)
      fit.v=fit_model(y.ord,inputs.ord,NNarray=NNarray,m_seq=m.seq,
                     X=as.matrix(sample(c(-1,1),n,replace=TRUE)),
                     start_parms=start.parms,max_iter=50,
                     covfun_name=covfun,silent=TRUE,
                     reorder=FALSE,fixed_parms=2+d)
      par.ests[rep,2,i.m,i.n,]=fit.v$covparms[1:(d+1)]
      preds=predictions2(fit=fit.v,locs_pred=inputs.test,
                        X_pred=matrix(0,n.test,1),m=2*m)
      mse[rep,2,i.m,i.n]=mean((preds-y.test)^2)
      
      ### low rank
      fm = NNarray[m+1,2:(m+1)]
      NNarray[(m+2):n,2:(m+1)]=matrix(rep(fm,n-m-1),byrow=TRUE,ncol = m)
      fit.lr=fit_model(y.ord,inputs.ord,NNarray=NNarray,m_seq=m.seq,
                     X=as.matrix(sample(c(-1,1),n,replace=TRUE)),
                     start_parms=start.parms,max_iter=50,
                     covfun_name=covfun,silent=TRUE,
                     reorder=FALSE,fixed_parms=2+d)
      par.ests[rep,3,i.m,i.n,]=fit.lr$covparms[1:(d+1)]
      preds=predictions_LR(fit.lr,inputs.test,m=2*m,X_pred=matrix(0,n.test,1))
      mse[rep,3,i.m,i.n]=mean((preds-y.test)^2)
      
      
      ### laGP
      la.pred=laGP::aGP(inputs,y,inputs.test,start=4,end=min(2*m,n-1),
                          method="alcray",verb=0)
      mse[rep,4,i.m,i.n]=mean((la.pred$mean-y.test)^2)
  
      ### save  results
      save(ms,ns,mse,par.ests,file='results/bore_mse_comp.RData')
    
    }
  }
}




###########   plot results   ########


load(file='results/bore_mse_comp.RData')
rmse=sqrt(apply(mse,2:4,mean,na.rm=TRUE))

ms.ind=c(1,3,6)
y.range=range(rmse[,ms.ind,],na.rm=TRUE)
# y.ticks=c(round(y.range[1],2),10^seq(ceiling(log10(y.range[1])),
#                                      floor(log10(y.range[2]))),round(y.range[2]))
y.ticks=c(.02,.05,.1,.5,1,5,10,15)


#### comparison for increasing n
pdf(file='plots/bore_mse_comp_logn.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.5)) # bltr
plot(1,1,col='white',xlim=range(log10(ns)),ylim=log10(y.range),
     ylab='RMSE',xlab='n',axes=FALSE)
for(i.m in 1:length(ms.ind)) matplot(log10(ns),log10(t(rmse[,ms.ind[i.m],])),
                                     type='l',add=TRUE,lwd=2,lty=i.m)
axis(1,at=log10(ns),labels=ns,lwd=0,lwd.ticks=1)
axis(2,at=log10(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
box()
legend('bottomleft',legend=ms[ms.ind],col=1,lty=1:length(ms.ind),lwd=2,
       title=expression(bold('m')),bg='white')
dev.off()



#### comparison for increasing m
ns.ind=c(1,4,6)
pdf(file='plots/bore_mse_comp.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.5)) # bltr
plot(1,1,col='white',xlim=range(ms),ylim=log10(y.range),
     ylab='RMSE',xlab='m',axes=FALSE)
for(i.n in 1:length(ns.ind)) matplot(ms,log10(t(rmse[,,ns.ind[i.n]])),type='l',
                            add=TRUE,lwd=2,lty=i.n)
axis(1,at=seq(10,50,by=10),labels=seq(10,50,by=10),lwd=0,lwd.ticks=1)
axis(2,at=log10(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
box()
legend('topright',legend=ns[ns.ind],col=1,lty=1:length(ns.ind),lwd=2,
       title=expression(bold('n')),bg='white')
dev.off()


### methods legend
pdf(file='plots/bore_mse_comp_legend.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center',c('SVecchia','Vecchia','LowRank','laGP'),col=1:4,lty=1,lwd=2,
       title=expression(bold('Method')),bg='white')
dev.off()


#### parameter estimates
mean.scales=array(dim=c(length(ms),length(ns),d))
for(i.n in 1:length(ns)){
  for(i.m in 1:length(ms)){
    mean.scales[i.m,i.n,]=colMeans(par.ests[,1,i.m,i.n,2]/par.ests[,1,i.m,i.n,-1])
  }
}

i.n=6; i.m=6
boxplot(par.ests[,1,i.m,i.n,2]/par.ests[,1,i.m,i.n,-1])
boxplot(log10(par.ests[,1,i.m,i.n,2]/par.ests[,1,i.m,i.n,-1]))

i.m=5
matplot(log10(ns),mean.scales[i.m,,],type='l')

i.n=6
matplot(ms,mean.scales[,i.n,],type='l')

