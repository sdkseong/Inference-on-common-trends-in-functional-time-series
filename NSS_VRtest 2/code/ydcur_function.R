
if (trend==0){
  
  load("cvalue/criticalVR.RData")
  load("cvalue/criticalVR12.RData")
  fdata = (ydcurve) - rowMeans(ydcurve)  
  
  xcoef = t(fdata)
} else {
  
  load("cvalue/criticalVRTrend.RData")
  load("cvalue/criticalVR12Trend.RData")
  
  
  expl = seq(1,tobs,1)
  expl = cbind(1, expl)
  fdata = ydcurve  -  t(expl%*%solve(crossprod(expl))%*%t(expl)%*%t(ydcurve) )
  xcoef = t(fdata)
}
  

lrx0 = lrvar(xcoef, aband = 0, bband = 0 ,kernel = ker  )
eig.lrx0 = eigen(lrx0$omega/(tobs^2))$vectors 

fdtmp0 = t(xcoef%*%eig.lrx0)  


##Projection for VR(2,0), VR(1,0), VR(0,2)##
lrx1 = lrvar(xcoef, aband = 1, bband = hr1 ,kernel = ker )
eig.lrx1 = eigen(lrx1$omega/(tobs^2))$vectors 

fdtmp1 = t(xcoef%*%eig.lrx1) 


##Projection for  VR(1,2)##
lrx2 = lrvar(xcoef, aband = 1, bband = hr1 ,kernel = ker )
eig.lrx2 = eigen(lrx2$omega/(tobs^2))$vectors 

fdtmp2 = t(xcoef%*%eig.lrx2) 

dimest.vr12.t <-  NULL   


dimest.vr12 =  t(sapply(c(0:max.xdim + adddim), function(x){
  fdtmp.t = fdtmp2[1:x,]
  lftmp = t(apply(fdtmp.t,1,cumsum)); rftmp =  fdtmp.t
  vrtmp = inv.vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x-adddim, al = al_inv, ar = 1, hl=hl_inv, hr = hr_inv, kn = ker)
  dimest.vr12.tmp =  c(vrtmp$vr.t >=CMatVR12t[(x-adddim+1),adddim,2] , vrtmp$vr.m >=CMatVR12m[(x-adddim+1),adddim,2])    
  return( c(vrtmp$vr.t,vrtmp$vr.m,dimest.vr12.tmp )   )}))

est.dim.inv = apply(rbind(dimest.vr12[,3:4],0),2,which.min) -1 + aux.add -1


dimest.vr12.t = cbind(dimest.vr12.t,dimest.vr12[,1],rowSums(matrix(rep(dimest.vr12[,1],3),ncol = 3) > CMatVR12t[1:(max.xdim+1),adddim,]))
 
if (is.null(smax_est) == T){smax_est = est.dim.inv[1]}

 
dimest.vr20.t<- dimest.vr20.m <-dimest.vr21.t<-dimest.vr21.m <-dimest.vr10.t<-dimest.vr10.m <- NULL       
for (tdim in 1:smax_est){
  pdim   = tdim + adddim 
  fdptmp0 = fdtmp0[1:pdim,]
  fdptmp1 = fdtmp1[1:pdim,]
  fdptmp2 = fdtmp2[1:pdim,]
  
  
  
  lftmp = t(apply(fdptmp0,1,cumsum)); rftmp =  fdptmp0
  vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker1) 
  dimest.vr21.t <- rbind(dimest.vr21.t,vrtmp$vr.t); dimest.vr21.m <- rbind(dimest.vr21.m,vrtmp$vr.m);
  
  if (tdim == smax_est){
    eigestvr21 = eigest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker) 
  }
  lftmp = t(apply(fdptmp1,1,cumsum))[,2:tobs]; rftmp = fdptmp1[,2:tobs] - fdptmp1[,1:(tobs-1)]  
  vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 0, al = al_20, ar = 1, hl=hlspecified, hr = hr2, kn = ker, optimal = optm)
  dimest.vr20.t <- rbind(dimest.vr20.t,vrtmp$vr.t); dimest.vr20.m <- rbind(dimest.vr20.m,vrtmp$vr.m);
  
  lftmp = fdptmp1[,2:tobs]
  vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 1, rfo = 0, al = al_10, ar = 1, hl=hlspecified, hr = hr2, kn = ker, optimal = optm)
  dimest.vr10.t <- rbind(dimest.vr10.t,vrtmp$vr.t); dimest.vr10.m <- rbind(dimest.vr10.m,vrtmp$vr.m);
}

dimest.vr21.t = cbind(dimest.vr21.t,rowSums(matrix(rep(dimest.vr21.t,3),ncol = 3) > CMatVR21t[1:smax_est,]))
dimest.vr21.m = cbind(dimest.vr21.m,rowSums(matrix(rep(dimest.vr21.m,3),ncol = 3) > CMatVR21m[1:smax_est,])) 
dimest.vr20.t = cbind(dimest.vr20.t,rowSums(matrix(rep(dimest.vr20.t,3),ncol = 3) > CMatVR20t[1:smax_est,]))
dimest.vr20.m = cbind(dimest.vr20.m,rowSums(matrix(rep(dimest.vr20.m,3),ncol = 3) > CMatVR20m[1:smax_est,]))
dimest.vr10.t = cbind(dimest.vr10.t,rowSums(matrix(rep(dimest.vr10.t,3),ncol = 3) > CMatVR10t[1:smax_est,]))
dimest.vr10.m = cbind(dimest.vr10.m,rowSums(matrix(rep(dimest.vr10.m,3),ncol = 3) > CMatVR10m[1:smax_est,]))


 
#########LRS(2020)#########
pcacov = eigen(crossprod(xcoef), only.values = TRUE)
eiglist = tobs*pcacov$values 
del.value = 1e-4
eiglist[which(eiglist/max(eiglist) < del.value)] = 0
eiglist = eiglist[2:length(eiglist)]/eiglist[1:(length(eiglist)-1)]
eiglist[which(is.nan(eiglist)==TRUE)] = 1


dimest.lrs = eiglist[1:smax_est]