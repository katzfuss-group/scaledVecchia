######   comparison on several test functions   ####

# setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R')
source('code/test_functions.R')

library(tictoc)


### settings
n.train=1e5
n.test=2e4
n.est=3000
m.est=30
m.pred=140
m.lagp=100
reps=5
methods=c('SVecchia','Vecchia','laGP','trivial')


### loop over test functions
for(i.f in 1:length(testfuns)){
  
  testfun=testfuns[[i.f]]
  d=testfun$d
  fun=testfun$fun
  print(testfun$name)
  
  rmse=time=array(dim=c(length(methods),reps))
  relevance=array(dim=c(d,reps))
  
  ### repeat data simulation reps times for each function
  for(rep in 1:reps){
    
    ### simulate data
    inputs.train=lhs::randomLHS(n.train,d)
    y.train=apply(inputs.train,1,fun)
    inputs.test=matrix(runif(n.test*d),n.test)
    y.test=apply(inputs.test,1,fun)
    
    ### estimate parameters on subset
    ind.est=sample(1:n.train,n.est)
    y.est=y.train[ind.est]
    inputs.est=inputs.train[ind.est,]
    
    for(meth in 1:2){
      
      print(methods[meth])
      
      ## how to scale inputs depends on method
      if(methods[meth]=='SVecchia'){ scale='parms'
      } else if(methods[meth]=='Vecchia'){ scale='ranges'
      }
      
      ## estimate parameters on subset
      tic()
      fit=fit_scaled(y.est,inputs.est,ms=c(m.est),scale=scale)
      # summary.GpGp_fit(fit)
    
      if(methods[meth]=='SVecchia'){
        input.ranges=apply(inputs.train,2,function(x) diff(range(x)))
        relevance[,rep]=input.ranges/fit$covparms[1+(1:d)]
      }
      
      ## overwrite training data for prediction on full training set
      fit$y=y.train
      fit$locs=inputs.train
      if(fit$trend=='zero'){ fit$X=as.matrix(rep(0,n.train)) 
      } else if(fit$trend=='intercept'){ fit$X=as.matrix(rep(1,n.train)) 
      } else stop()
    
      ## make prediction
      preds=predictions_scaled(fit=fit,locs_pred=inputs.test,m=m.pred,scale=scale)
      temp=toc()
      time[meth,rep]=temp$toc-temp$tic
      
      ## compute rmse
      rmse[meth,rep]=sqrt(mean((preds-y.test)^2))
      
    }
    
    ## laGP
    tic()
    la.pred=laGP::aGP(inputs.train,y.train,inputs.test,end=m.lagp,
                      method="alcray",verb=0)
    temp=toc()
    time[3,rep]=temp$toc-temp$tic
    rmse[3,rep]=sqrt(mean((la.pred$mean-y.test)^2))
    
    ## naive
    tic()
    rmse[4,rep]=sqrt(mean((mean(y.train)-y.test)^2))
    temp=toc()
    time[4,rep]=temp$toc-temp$tic
    
    ## print and save
    print(rmse)
    print(time)
    save(time,rmse,relevance,file=paste0('results/testfun_reps_',
                                         names(testfuns)[i.f],'.RData'))
    
  }
}

  


######  look at prediction results
tab.data=array(dim=c(length(testfuns),length(methods),2))
fun.names=numeric(length=length(testfuns))
for(i.f in 1:length(testfuns)){
  load(file=paste0('results/testfun_reps_',names(testfuns)[i.f],'.RData'))
  if(i.f<=2) mult=1e2 else mult=1e5
  tab.data[i.f,,1]=rowMeans(rmse)*mult
  tab.data[i.f,,2]=rowMeans(time)
  fun.names[i.f]=testfuns[[i.f]]$name
}

tab=t(rbind(tab.data[1,,1],tab.data[1,,2],tab.data[2,,1],
            tab.data[2,,2],tab.data[3,,1],tab.data[3,,2]))
rownames(tab)=methods
colnames(tab)=c('RMSE x 1e2','time (sec)','RMSE x 1e2',
                'time (sec)','RMSE x 1e5','time (sec)')

### prepare table for latex
library(xtable)
xtable(tab,digits=c(0,rep(c(1,0),3)))
fun.names


## plot relevance
for(i.f in 1:length(testfuns)){
  load(file=paste0('results/testfun_reps_',names(testfuns)[i.f],'.RData'))
  matplot(log10(relevance),xlab='input variable',main=names(testfuns)[i.f])
}
