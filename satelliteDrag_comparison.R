
#########  comparison of GP methods on satellite drag data

##  data and laGP results are from
## https://bitbucket.org/gramacylab/tpm/src/master/
# To reproduce the results, save the following files
#    from tpm/data/HST to a local subfolder 'HSTdata/':
#   - hstXX.dat for XX in (O,O2,N,N2,He,H)
#     (combine hstHe and hstHe2 into a single file hstHe)
#   - hstQ_05.dat, hstQ_joint_05.csv, hstA_05.dat, hstA_joint_05.csv
#   - for laGP CV scores, request sat_hstXX_laGP_cv.RData from Furong Sun


###  load functions and packages

# setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R')
library(tictoc)




##########    cross-validation experiment   ########

species=c('O','O2','N','N2','He','H')
methods=c('SVecchia','Vecchia')
m.est=30
m.pred=140
n.est=10000  # size of random subset for estimation
n.all=2e6
d=8

### split into 10 subsets for crossvalidation
set.seed(999)
folds=10
cv.inds=matrix(sample(1:n.all,n.all),nrow=folds)
n.test=n.all/folds
n.train=n.all-n.test

### initialize output
time=mse=array(dim=c(length(methods),length(species),2,folds))
par.ests=array(dim=c(length(methods),length(species),d+1,folds))

### loop over the 6 chemical species
for(i.s in 1:length(species)){
  
  print(species[i.s])
  
  ### load data
  sat.data=read.table(paste0('HSTdata/hst',species[i.s],'.dat'),header=TRUE)

  ### loop over CV folds
  for(fold in 1:folds){
    
    print(paste0('fold=',fold))
    
    ## training and test data
    ind.test=cv.inds[fold,]
    y.train=as.numeric(sat.data[-ind.test,9])
    inputs.train=as.matrix(sat.data[-ind.test,1:8],rownames=FALSE)
    y.test=as.numeric(sat.data[ind.test,9])
    inputs.test=as.matrix(sat.data[ind.test,1:8],rownames=FALSE)
    
    ## estimation data set (smaller, for speed)
    ind.est=sample(1:n.train,n.est)
    y.est=y.train[ind.est]
    inputs.est=inputs.train[ind.est,]

    ## loop over GP methods
    for(meth in 1:length(methods)){
    
      print(methods[meth])
      
      ## how to scale inputs depends on method
      if(methods[meth]=='SVecchia'){ scale='parms'
      } else if(methods[meth]=='Vecchia'){ scale='ranges'
      }
      
      ## estimate parameters on subset
      tic()
      fit=fit_scaled(y.est,inputs.est,ms=c(m.est),scale=scale)
      temp=toc()
      time[meth,i.s,1,fold]=temp$toc-temp$tic
      par.ests[meth,i.s,,fold]=fit$covparms[1:(d+1)]
      
      ## overwrite training data for prediction on full training set
      fit$y=y.train
      fit$locs=inputs.train
      if(fit$trend=='zero'){ fit$X=as.matrix(rep(0,n.train)) 
      } else if(fit$trend=='intercept'){ fit$X=as.matrix(rep(1,n.train)) 
      } else stop()
    
      ## make prediction
      tic()
      preds=predictions_scaled(fit=fit,locs_pred=inputs.test,m=m.pred,scale=scale)
      temp=toc()
      time[meth,i.s,2,fold]=temp$toc-temp$tic
      
      ## compute rmse
      mse[meth,i.s,1,fold]=mean((preds-y.test)^2)
      mse[meth,i.s,2,fold]=mean((100*(preds - y.test)/y.test)^2)
      
      ## print new results and save
      print(sqrt(mse[meth,i.s,,fold])); print(time[meth,i.s,,fold])
      save(mse,time,par.ests,file=paste0('results/satelliteCV.RData'))
      
    }
  }
}


## CV results from laGP (files from Furong Sun)
# rmspe.laGP=array(dim=c(folds,19,length(species)))
# for(i.s in 1:length(species)){
#   load(paste0('HSTdata/sat_hst',species[i.s],'_laGP_cv.RData'))
#   rmspe.laGP[,,i.s]=as.matrix(rmspe)
# }
# meth.names.laGP=colnames(rmspe)
# save(rmspe.laGP,meth.names.laGP,file=paste0('results/satelliteCV_laGP.RData'))
load(file=paste0('results/satelliteCV_laGP.RData'))
rmse.laGP=sqrt(apply(rmspe.laGP^2,2:3,mean,na.rm=TRUE))



####  prediction results

load(paste0('results/satelliteCV.RData'))

## prediction rmse
rmse=sqrt(apply(mse,1:3,mean,na.rm=TRUE))
rmse.all=rbind(rmse[,,2],rmse.laGP)

## plotting details
cols=c(1,2,rep(4,7),rep(5,12))
pchs=c(1,2,rep(0,7),rep(4,12))
y.ticks=c(.5,1,2,4,8,16)

## save the plot
pdf(file='plots/satellite_CV.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.3)) # bltr
matplot(log10(t(rmse.all)),col=cols,pch=pchs,cex=1.5,xaxt='n',yaxt='n',
        ylab='RMSPE')
axis(1,at=1:length(species),labels=species,lwd=0,lwd.ticks=1)
axis(2,at=log10(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
abline(h=0,col=8)
legend('topright',c(methods,'laGP','H-laGP'),col=c(1,2,4,5),pch=c(1,2,0,4),
       pt.cex=1.3,bg='white')
dev.off()


## RMSPE of best laGP is XX% higher than of SVecchia
round((apply(rmse.all[-(1:2),],2,min)/rmse.all[1,]-1)*100,1)

# naive
pred.naive=mean(y.train)
sqrt(mean((100*(pred.naive - y.test)/y.test)^2)) # around 50%


#### time
par(mfrow=c(1,3))
boxplot(apply(time[,,1,],1,cbind)/60,ylab='min',main='estimation')
boxplot(apply(time[,,2,],1,cbind)/60,ylab='min',main='prediction')
boxplot(apply(apply(time,c(1,2,4),sum),1,cbind)/60,ylab='min',main='total')
par(mfrow=c(1,1))


#### estimated parameters
matplot(log10(t(par.ests[1,i.s,,])))
# variability over folds is very small




##########    trajectories   ########


#####   estimate parameters based on subset of full training set

species=c('O','O2','N','N2','He','H')
methods=c('SVecchia','Vecchia')
m.est=30
n.est=10000  # size of random subset for estimation

for(i.s in 1:length(species)){
  
  print(species[i.s])
  
  ### load data
  sat.data=read.table(paste0('HSTdata/hst',species[i.s],'.dat'),header=TRUE)
  y.train=as.numeric(sat.data[,9])
  inputs.train=as.matrix(sat.data[,1:8],rownames=FALSE)
  n.train=nrow(inputs.train)
  
  ## estimation data set (smaller, for speed)
  ind.est=sample(1:n.train,n.est)
  y.est=y.train[ind.est]
  inputs.est=inputs.train[ind.est,]
  
  ## loop over GP methods
  for(meth in 1:length(methods)){
    
    print(methods[meth])
    
    ## how to scale inputs depends on method
    if(methods[meth]=='SVecchia'){ scale='parms'
    } else if(methods[meth]=='Vecchia'){ scale='ranges'
    }
    
    ## estimate parameters on subset
    fit=fit_scaled(y.est,inputs.est,ms=c(m.est),scale=scale)[c(1,2,10,14)]
      # covparms betahat covfun_name trend

    save(fit,file=paste0('results/traj_fit_',species[i.s],'_',methods[meth],'.RData'))
    
  }
}




#### plot estimated parameters

## compute relevance
d=8
parms=array(dim=c(length(species),d+1))
relevance=array(dim=c(length(species),d))
for(i.s in 1:length(species)){
  load(file=paste0('results/traj_fit_',species[i.s],'_',methods[1],'.RData'))
  parms[i.s,]=fit$covparms[1:(d+1)]
  sat.data=read.table(paste0('HSTdata/hst',species[i.s],'.dat'),header=TRUE)
  input.ranges=apply(sat.data[,1:8],2,function(x) diff(range(x)))
  relevance[i.s,]=input.ranges/fit$covparms[1+(1:d)]
}


## save the plot
pdf(file='plots/satellite_relevance.pdf',width=4.0,height=4.0)
par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.3)) # bltr
y.ticks=10^((-2):1) #2^c(-4,-2,0,2,4)
matplot(log10(t(relevance)),pch=length(species):1,cex=1.5,yaxt='n',
        ylab='relevance',xlab='input variable',
        ylim=log10(range(c(y.ticks,relevance))))
axis(2,at=log10(y.ticks),labels=y.ticks,lwd=0,lwd.ticks=1)
legend('bottomright',species,col=1:length(species),pch=length(species):1,
       pt.cex=1.3,bg='white')
dev.off()




#####   make predictions for trajectory

regimes=c('Q','A') # quiet vs active
m.pred=140
n.test=8600
species=c('O','O2','N','N2','He','H')
methods=c('SVecchia','Vecchia')

### initialize output
time=array(dim=c(length(methods),length(species),length(regimes)))
pred.s=array(dim=c(length(methods),length(regimes),n.test,length(species)))
pred.w=array(dim=c(length(methods),length(regimes),n.test))
rmse.percent=array(dim=c(length(methods),length(regimes)))

### loop over the two regimes (quiet vs active)
for(i.r in 1:length(regimes)){

  ### load trajectory data
  ds=read.table(paste0('HSTdata/hst',regimes[i.r],'_05.dat'),header=TRUE,nrows=n.test)
  moles=data.frame(O=ds$X_O,O2=ds$X_O2,N=ds$X_N,N2=ds$X_N2,He=ds$X_He, H=ds$X_H)
  X=data.frame(Umag=ds$Umag.m.s.,Ts=ds$T_s.K.,Ta=ds$T_a.K.,alphan=ds$alpha_n.n.a., 
                  sigmat=ds$sigma_t.n.a.,theta=ds$Yaw.rad.,phi=ds$Pitch.rad.,panang=50)
  inputs.test=as.matrix(X,rownames=FALSE)
  y.test=as.numeric(ds$Cd)

  ### weights for each species
  w=matrix(nrow=n.test, ncol=ncol(moles))
  for(i in 1:n.test){
    mf <- as.numeric(moles[i,])
    pm <- c(2.65676, 5.31352, 2.32586, 4.65173, 0.664648, 0.167356)
    w[i,] <- mf * pm/sum(mf * pm)   ## weight each location
  }
  
  
  ### loop over the 6 chemical species
  for(i.s in 1:length(species)){
    
    print(species[i.s])
    
    ### load training data
    sat.data=read.table(paste0('HSTdata/hst',species[i.s],'.dat'),header=TRUE)
    y.train=as.numeric(sat.data[,9])
    inputs.train=as.matrix(sat.data[,1:8],rownames=FALSE)
    n.train=nrow(inputs.train)
    
    ## loop over GP methods
    for(meth in 1:length(methods)){
      
      print(methods[meth])
      
      ## load parameter estimates
      load(paste0('results/traj_fit_',species[i.s],'_',methods[meth],'.RData'))
      
      ## specify training data
      fit$y=y.train
      fit$locs=inputs.train
      if(fit$trend=='zero'){ fit$X=as.matrix(rep(0,n.train)) 
      } else if(fit$trend=='intercept'){ fit$X=as.matrix(rep(1,n.train)) 
      } else stop()
      
      ## add nugget for numerical stability
      fit$covparms[length(fit$covparms)]=fit$covparms[1]*1e-12
      
      ## how to scale inputs depends on method
      if(methods[meth]=='SVecchia'){ scale='parms'
      } else if(methods[meth]=='Vecchia'){ scale='ranges'
      }
      
      ## make prediction
      tic()
      pred.s[meth,i.r,,i.s]=predictions_scaled(fit=fit,locs_pred=inputs.test,
                                      m=m.pred,scale=scale)
      temp=toc()
      time[meth,i.s,i.r]=temp$toc-temp$tic
      
    }
  }
  
  
  ### compute weighted prediction over species
  for(meth in 1:length(methods)){
    for(j in 1:n.test){
      pred.w[meth,i.r,j]=sum(pred.s[meth,i.r,j,]*w[j,])
    }
    rmse.percent[meth,i.r]=sqrt(mean((100*(pred.w[meth,i.r,]-y.test)/y.test)^2))
  }
  
  save(time,pred.s,pred.w,rmse.percent,file=paste0('results/traj_pred.RData'))
  
}



####  plot and print trajectory results

# combine ours with gramacy results
load(paste0('results/traj_pred.RData'))
rmse.percent.all=array(dim=c(length(methods)+15,length(regimes)))
rmse.percent.all[1:length(methods),]=rmse.percent
for(i.r in 1:length(regimes)){
  gram=read.csv(paste0('HSTdata/hst',regimes[i.r],'_joint_05.csv'),header=TRUE)
  rmspes.gram.all=gram[,6:10]
  rmse.percent.all[-(1:length(methods)),i.r]=sqrt(colMeans(rmspes.gram.all^2))
}


# plot prediction accuracy
matplot(t(rmse.percent.all))

## RMSPE of best laGP is XX% higher than of SVecchia
round((apply(rmse.percent.all[-(1:2),],2,min)/rmse.percent.all[1,]-1)*100,1)


### time
boxplot(apply(time,1,cbind),ylab='sec')


