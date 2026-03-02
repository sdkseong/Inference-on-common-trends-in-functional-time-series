
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



##Projection for  VR(1,2)##
lrx2 = lrvar(xcoef, aband = 1, bband = hr1 ,kernel = ker )
eig.lrx2 = eigen(lrx2$omega/(tobs^2))$vectors 

fdtmp2 = t(xcoef%*%eig.lrx2) 

dimest.vr12.t <-  NULL   


dimest.vr12 =  t(sapply(c(0:max.xdim.inv + adddim.inv), function(x){
  fdtmp.t = fdtmp2[1:x,]
  if (x==1){
    lftmp = cumsum(fdtmp.t); rftmp =  fdtmp.t
    vrtmp = inv.vrtest(lfdata = t((lftmp)), rfdata = t(rftmp),  testdim = x-adddim.inv, al = al_inv, ar = 1, hl=hl_inv, hr = hr_inv, kn = ker)
    
  }else {
    lftmp = t(apply(fdtmp.t,1,cumsum)); rftmp =  fdtmp.t
    vrtmp = inv.vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x-adddim.inv, al = al_inv, ar = 1, hl=hl_inv, hr = hr_inv, kn = ker)
    
  
  }
  dimest.vr12.tmp =  c(vrtmp$vr.t >=CMatVR12t[(x-adddim.inv+1),adddim.inv,2] , vrtmp$vr.m >=CMatVR12m[(x-adddim.inv+1),adddim.inv,2])    
  return( c(vrtmp$vr.t,vrtmp$vr.m,dimest.vr12.tmp )   )}))

est.dim.inv = apply(rbind(dimest.vr12[,3:4],0),2,which.min) -1 + aux.add -1


dimest.vr12.t = cbind(dimest.vr12.t,dimest.vr12[,1],rowSums(matrix(rep(dimest.vr12[,1],3),ncol = 3) > CMatVR12t[1:(max.xdim.inv+1),adddim.inv,]))
