######   comparison on several test functions   ####

# setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R') # SVecchia
source('code/other_methods.R') # laGP and lowRank
source('code/test_functions.R')

library(tictoc)
library(laGP)

### settings
n.train=1e5
n.test=2e4
n.est=3000
ms=c(30,50)
m.pred=140
reps=5
methods=c('SVecchia','Vecchia','laGP','H-laGP',
          'laGP4c','H-laGP4c','trivial')



### loop over different settings and test functions
for(i.m in 1:length(ms)){
  
  m.est=m.lagp=m.hybrid=ms[i.m]
  
  for(i.f in 1:length(testfuns)){
    
    testfun=testfuns[[i.f]]
    d=testfun$d
    fun=testfun$fun
    print(testfun$name)
    
    rmse=array(dim=c(length(methods),reps))
    time=array(dim=c(length(methods),2,reps))
    relevance=array(dim=c(d,reps))
    
    ### repeat data simulation reps times for each function
    for(rep in 1:reps){
      
      ### simulate data
      inputs.train=lhs::randomLHS(n.train,d)
      y.train=apply(inputs.train,1,fun)
      inputs.test=matrix(runif(n.test*d),n.test)
      y.test=apply(inputs.test,1,fun)
      
      #### Vecchia methods
      for(meth in 1:2){
        
        print(methods[meth])
        
        ## how to scale inputs depends on method
        if(methods[meth]=='SVecchia'){ scale='parms'
        } else if(methods[meth]=='Vecchia'){ scale='ranges'
        }
        
        ## estimate parameters
        tic()
        fit=fit_scaled(y.train,inputs.train,ms=m.est,scale=scale,n.est=n.est)
        temp=toc()
        time[meth,1,rep]=temp$toc-temp$tic
        # summary.GpGp_fit(fit)
        
        ## make prediction
        tic()
        preds=predictions_scaled(fit=fit,locs_pred=inputs.test,m=m.pred,scale=scale)
        temp=toc()
        time[meth,2,rep]=temp$toc-temp$tic
        
        ## compute rmse
        rmse[meth,rep]=sqrt(mean((preds-y.test)^2))
        
        ## save estimated ranges/relevances for SVecchia
        if(methods[meth]=='SVecchia'){
          input.ranges=apply(inputs.train,2,function(x) diff(range(x)))
          relevance[,rep]=input.ranges/fit$covparms[1+(1:d)]
        }
        
      }
      
      #### laGP methods
      
      g <- 1/10000000  # tiny nugget as suggested by bobby
      
      # ## laGP
      # tic()
      # la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
      #                   verb=0,omp.threads=1,g=g)
      # temp=toc()
      # time[3,1,rep]=0
      # time[3,2,rep]=temp$toc-temp$tic
      # rmse[3,rep]=sqrt(mean((la.pred$mean-y.test)^2))
      # 
      # ## hybrid laGP
      # h.pred=global_local_laGP(inputs.train,y.train,inputs.test,m.hybrid,nth=1)
      # time[4,1,rep]=h.pred$time.train
      # time[4,2,rep]=h.pred$time.pred
      # rmse[4,rep]=sqrt(mean((h.pred$mean-y.test)^2))
      
      ## laGP - 4 cores
      tic()
      la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
                        verb=0,omp.threads=4,g=g)
      temp=toc()
      time[5,1,rep]=0
      time[5,2,rep]=temp$toc-temp$tic
      rmse[5,rep]=sqrt(mean((la.pred$mean-y.test)^2))
      
      ## hybrid laGP - 4 cores
      h.pred=global_local_laGP(inputs.train,y.train,inputs.test,m.hybrid,nth=4)
      time[6,1,rep]=h.pred$time.train
      time[6,2,rep]=h.pred$time.pred
      rmse[6,rep]=sqrt(mean((h.pred$mean-y.test)^2))
      
      
      ##### naive
      tic()
      rmse[7,rep]=sqrt(mean((mean(y.train)-y.test)^2))
      temp=toc()
      time[7,2,rep]=temp$toc-temp$tic
      
      
      #### print and save
      print(rmse)
      print(time)
      save(time,rmse,relevance,file=paste0('results/testfun_m',ms[i.m],'_',
                                           names(testfuns)[i.f],'.RData'))
      
    }
  }
}



######  table of prediction results
tabs=array(dim=c(length(ms),length(methods),length(testfuns)*2))
tab.data=array(dim=c(length(testfuns),length(methods),2))
fun.names=numeric(length=length(testfuns))
for(i.m in 1:length(ms)){
  for(i.f in 1:length(testfuns)){
    load(file=paste0('results/testfun_m',ms[i.m],'_',names(testfuns)[i.f],'.RData'))
    if(i.f<=2) mult=1e2 else mult=1e5
    tab.data[i.f,,1]=round(sqrt(rowMeans(rmse^2,na.rm=TRUE))*mult,1)
    # tab.data[i.f,5:6,1]=NA
    time.train=round(rowMeans(time[,1,],na.rm=TRUE)/60,1)
    time.pred=round(rowMeans(time[,2,],na.rm=TRUE)/60,1)
    tab.data[i.f,,2]=paste0(time.train,'+',time.pred,'=',time.train+time.pred)
    tab.data[i.f,,2]=paste0(time.train,'$+$',time.pred,'$=$',
                            time.train+time.pred)
    fun.names[i.f]=testfuns[[i.f]]$name
  }
  tabs[i.m,,]=t(rbind(tab.data[1,,1],tab.data[1,,2],tab.data[2,,1],tab.data[2,,2],
              tab.data[3,,1],tab.data[3,,2]))
}
tab=rbind(tabs[1,c(1:3,5,4,6),],tabs[2,c(1:3,5,4,6,7),])
tab=tab[-c(3,5,9,11),]
rownames(tab)=c(rep(c('SVecchia','Vecchia','laGP','H-laGP'),2),'trivial')
colnames(tab)=c('RMSE x 1e2','time (sec)','RMSE x 1e2',
                'time (sec)','RMSE x 1e5','time (sec)')

## prepare table for latex
library(xtable)
print(xtable(tab),sanitize.text.function=identity)
fun.names



######  plots of prediction results
results=array(dim=c(length(testfuns),length(methods)*2,3))
fun.names=numeric(length=length(testfuns))
for(i.f in 1:length(testfuns)){
  for(i.m in 1:length(ms)){
    load(file=paste0('results/testfun_m',ms[i.m],'_',names(testfuns)[i.f],'.RData'))
    ind=(1:length(methods))+length(methods)*(i.m-1)
    results[i.f,ind,1]=rowMeans(time[,1,],na.rm=TRUE)/60
    results[i.f,ind,2]=results[i.f,ind,1]+rowMeans(time[,2,],na.rm=TRUE)/60
    results[i.f,ind,3]=sqrt(rowMeans(rmse^2,na.rm=TRUE))
  }
  fun.names[i.f]=testfuns[[i.f]]$name
}
results=results[,-c(7,14),] # remove trivial predictor
results[,c(3,5,9,11),1]=NA # remove 0 training times for laGP
for(i.f in 1:length(testfuns)) # average H-laGP training times
  results[i.f,c(4,6,10,12),1]=mean(results[i.f,c(4,6,10,12),1],na.rm=TRUE)

## plot settings
cols=c(1,2,4,5,4,5)
symbs=matrix(c(1,16,2,17,0,15,5,18,0,15,5,18),ncol=2,byrow=TRUE)

## actual plotting
for(i.f in 1:length(testfuns)){
  pdf(file=paste0('plots/testfun_earl_',names(testfuns)[i.f],'.pdf'),
                  width=4.0,height=4.0)
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.5)) # bltr
  if(i.f<=2) mult=2 else mult=5
  dat=log2(cbind(results[i.f,,1:2],results[i.f,,3]*10^mult))
  xl=range(c(dat[,1:2]),na.rm=TRUE)
  yl=range(dat[,3],na.rm=TRUE)
  x.ticks=2^seq(floor(xl[1]),ceiling(xl[2]))
  y.ticks=2^seq(floor(yl[1]),ceiling(yl[2]))
  plot(NULL,xaxt='n',yaxt='n',ylim=yl,xlim=xl,
       ylab=paste0('RMSE x 10^',mult),xlab='time (min)')
  axis(1,at=log2(x.ticks),labels=x.ticks,lwd=0,lwd.ticks=1)
  axis(2,at=log2(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
  for(i.t in 1:2) points(dat[,i.t],dat[,3],pch=symbs[,i.t],col=cols,cex=i.t+.5)
  segments(dat[,1],dat[,3],dat[,2],dat[,3],col=cols,lwd=1)
  dev.off()
}

### methods legend
meths.leg=c('SVecchia  ','Vecchia ','laGP','H-laGP')
pdf(file='plots/testfun_legend.pdf',width=5.0,height=1.0)
par(mgp = c(0,0,0), mar=c(0,0,0,0)) # bltr
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)
legend('center',meths.leg,col=cols,pch=symbs[1:length(meths.leg),2],
       bg='white',horiz=TRUE,pt.cex=1.5)
dev.off()

pdf(file='plots/testfun_legend2.pdf',width=2.5,height=1.0)
par(mgp = c(0,0,0), mar=c(0,0,0,0)) # bltr
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=c(.44,.53))
points(.25,.5,pch=symbs[1,1],cex=1.5)
points(.75,.5,pch=symbs[1,2],cex=2.5)
segments(.25,.5,.75,.5)
text(c(.25,.75),rep(.475,2),c('train','total'),cex=1.5)
dev.off()



########  plot relevance
for(i.f in 1:length(testfuns)){
  load(file=paste0('results/testfun_reps_',names(testfuns)[i.f],'.RData'))
  matplot(log10(relevance),xlab='input variable',main=names(testfuns)[i.f])
}





################   assess UQ   ###############

methods=c('SVecchia','Vecchia','laGP','H-laGP')
m.est=m.lagp=m.hybrid=30

i.f=3
testfun=testfuns[[i.f]]
d=testfun$d
fun=testfun$fun
print(testfun$name)

## scores: int.cov, int.width, int.score, log.score, CRPS, energy.score
scores=array(dim=c(length(methods),6,reps))

### repeat data simulation reps times for each function
for(rep in 1:reps){ # rep=1
  
  ### simulate data
  inputs.train=lhs::randomLHS(n.train,d)
  y.train=apply(inputs.train,1,fun)
  inputs.test=matrix(runif(n.test*d),n.test)
  y.test=apply(inputs.test,1,fun)
  
  #### Vecchia methods
  for(meth in 1:2){ # meth=1
    
    print(methods[meth])
    
    ## how to scale inputs depends on method
    if(methods[meth]=='SVecchia'){ scale='parms'
    } else if(methods[meth]=='Vecchia'){ scale='ranges'
    }
    
    ## estimate parameters
    fit=fit_scaled(y.train,inputs.train,ms=m.est,scale=scale,n.est=n.est,
                   find.vcf=TRUE)

    ## make prediction
    preds=predictions_scaled(fit=fit,locs_pred=inputs.test,m=m.pred,nsims=50,
                             scale=scale)

    ## compute scores
    vars=apply(preds$samples,1,var)
    scores[meth,,rep]=compute_scores(y.test,preds$means,vars,preds$samples)
    
  }
  
  
  #### laGP methods
  
  g <- 1/10000000  # tiny nugget as suggested by bobby
  
  ## laGP - 4 cores
  la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
                    verb=0,omp.threads=4,g=g)
  scores[3,,rep]=compute_scores(y.test,la.pred$mean,la.pred$var)

  
  ## hybrid laGP - 4 cores
  h.pred=global_local_laGP(inputs.train,y.train,inputs.test,m.hybrid,nth=4)
  scores[4,,rep]=compute_scores(y.test,h.pred$mean,h.pred$var)
  
  
  #### print and save
  print(scores[,,rep])
  save(scores,file=paste0('results/testUQ_m',m.est,'_',
                                       names(testfuns)[i.f],'.RData'))
  
}


### print UQ table
i.f=3
load(file=paste0('results/testUQ_m',m.est,'_',names(testfuns)[i.f],'.RData'))
mults=c(1e2,1e5,1e5,1e1,1e5,1e5)
tab=round(t(t(apply(scores,1:2,mean,na.rm=TRUE))*mults),1)
rownames(tab)=methods
colnames(tab)=c('ICov','IWidth','IScore','LogScore','CRPS','Energy')
print(xtable::xtable(tab,digits=1))




################   estimation of the nugget variance  ###############

i.f=3
testfun=testfuns[[i.f]]
fun=testfun$fun

### simulate data with artifical noise
tau2=.02^2
inputs.train=lhs::randomLHS(n.train,d)
y.train=apply(inputs.train,1,fun)+rnorm(n.train,0,sqrt(tau2))


### fit SVecchia
m.est=30
fit=fit_scaled(y.train,inputs.train,ms=m.est,nug=NULL)
sqrt(fit$covparms[length(fit$covparms)]*fit$covparms[1])
