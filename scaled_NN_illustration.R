####  create plot to illustrate ordering and NN in scaled space

setwd("G:/My Drive/projects/computerExperiments")
source('code/vecchia_scaled.R')

## settings
n=500
d=2
range.mult=2
ranges=c(1/range.mult,range.mult)
i=28
m=4
set.seed(9)

## generate inputs
inputs=lhs::randomLHS(n,2)

inputs.ord=NN=list()

## euclidean ordering, original space
ord=order_maxmin_exact(inputs)
inputs.ord[[1]]=inputs[ord,]
NN[[1]]=find_ordered_nn(inputs.ord[[1]],m)[,-1]

## euclidean ordering, scaled space
inputs.ord[[2]]=t(t(inputs.ord[[1]])/ranges)
NN[[2]]=NN[[1]]

## scaled ordering, scaled space
inputs.scaled=t(t(inputs)/ranges)
ord=order_maxmin_exact(inputs.scaled)
inputs.ord[[3]]=inputs.scaled[ord,]
NN[[3]]=find_ordered_nn(inputs.ord[[3]],m)[,-1]

## scaled ordering, original space
inputs.ord[[4]]=t(t(inputs.ord[[3]])*ranges)
NN[[4]]=NN[[3]]



### illustration plots

cols=c(2,2,1,1)
lab.scale=c('','scaled ','scaled ','')
mult=c(1,range.mult,range.mult,1)
y.adj=c(0,1,1,0)*.08/range.mult
for(j in 1:4){
  pdf(file=paste0('plots/orderillus_',i,'_',j,'.pdf'),
      width=3.7*mult[j],height=3.6/mult[j])
  par(mgp = c(1.6,.5,0), mar=c(2.6,2.6,.3,.1)) # bltr
  plot(inputs.ord[[j]][-(1:i),],pch='.',col='grey90',cex=2.5,
       xlab=paste0(lab.scale[j],'input 1'),xaxt='n',
       ylab=paste0(lab.scale[j],'input 2'),yaxt='n',
       xlim=c(0,mult[j]),ylim=c(0-y.adj[j],1/mult[j]+y.adj[j]))
  axis(1,at=seq(0,mult[j],by=.5))
  axis(2,at=seq(0,1/mult[j],by=.5))
  text(inputs.ord[[j]][1:(i-1),],col=cols[j],cex=1,labels=as.character(1:(i-1)))
  points(inputs.ord[[j]][NN[[j]][i,],],col=cols[j],pch=1,cex=3)
  text(inputs.ord[[j]][i,,drop=FALSE],col=cols[j],cex=1,labels=as.character(i))
  points(inputs.ord[[j]][i,,drop=FALSE],col=4,pch=0,cex=3)
  dev.off()
}


