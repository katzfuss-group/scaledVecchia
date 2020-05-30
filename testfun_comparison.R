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
m.est=30
m.pred=140
m.lagp=30
m.hybrid=30
reps=5
methods=c('SVecchia','Vecchia','laGP','H-laGP',
          'laGP4c','H-laGP4c','trivial')



### loop over test functions
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
    
    ## laGP
    tic()
    la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
                      verb=0,omp.threads=1)
    temp=toc()
    time[3,1,rep]=0
    time[3,2,rep]=temp$toc-temp$tic
    rmse[3,rep]=sqrt(mean((la.pred$mean-y.test)^2))
    
    ## hybrid laGP
    h.pred=global_local_laGP(inputs.train,y.train,inputs.test,m.hybrid,nth=1)
    time[4,1,rep]=h.pred$time.train
    time[4,2,rep]=h.pred$time.pred
    rmse[4,rep]=sqrt(mean((h.pred$mean-y.test)^2))
    
    ## laGP - 4 cores
    tic()
    la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
                      verb=0,omp.threads=4)
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
    save(time,rmse,relevance,file=paste0('results/testfun_m30_',
                                         names(testfuns)[i.f],'.RData'))
    
  }
}



######  look at prediction results
tab.data=array(dim=c(length(testfuns),length(methods),2))
fun.names=numeric(length=length(testfuns))
for(i.f in 1:length(testfuns)){
  load(file=paste0('results/testfun_m30_',names(testfuns)[i.f],'.RData'))
  if(i.f<=2) mult=1e2 else mult=1e5
  tab.data[i.f,,1]=round(sqrt(rowMeans(rmse^2,na.rm=TRUE))*mult,1)
  tab.data[i.f,5:6,1]=NA
  time.train=round(rowMeans(time[,1,],na.rm=TRUE)/60,1)
  time.pred=round(rowMeans(time[,2,],na.rm=TRUE)/60,1)
  tab.data[i.f,,2]=paste0(time.train,'+',time.pred,'=',time.train+time.pred)
  tab.data[i.f,,2]=paste0(time.train,'$+$',time.pred,'$=$',
                          time.train+time.pred)
  fun.names[i.f]=testfuns[[i.f]]$name
}

tab=t(rbind(tab.data[1,,1],tab.data[1,,2],tab.data[2,,1],tab.data[2,,2],
            tab.data[3,,1],tab.data[3,,2]))
rownames(tab)=methods
colnames(tab)=c('RMSE x 1e2','time (sec)','RMSE x 1e2',
                'time (sec)','RMSE x 1e5','time (sec)')

### prepare table for latex
library(xtable)
print(xtable(tab[c(1:3,5,4,6,7),]),sanitize.text.function=identity)
fun.names


## plot relevance
for(i.f in 1:length(testfuns)){
  load(file=paste0('results/testfun_reps_',names(testfuns)[i.f],'.RData'))
  matplot(log10(relevance),xlab='input variable',main=names(testfuns)[i.f])
}
