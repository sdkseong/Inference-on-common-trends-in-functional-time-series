
################################################################
############ Simulation code to replicate Table 2 ##############
#### "Inference on common trends in functional time series" #### 
################################################################

rm(list = ls())
library(fda)
library(tseries)
library(sandwich)
library(sde)
library(variables)
library(basefun)
library(polynom)
library(geigen) 
library(ftsa)

today = "result20260217/" #Result output folder


setwd("/");  
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

if(file.exists(paste0("result/",today))==FALSE){
  dir.create(paste0("result/",today))
}



## Critical values ## 
load("cvalue/criticalVR.RData")  #Load Critical values for VR tests
load("cvalue/criticalVR3.RData")  #Load Critical values for VR Higher order
load("cvalue/criticalVR12.RData") #Local Critical values for Inverse VR


## Simulation-specific parameters ##
nnbasis=31  ## the number of basis functions
tsim = 1e4  ## the number of iteration
ARDfac=0.9  ## parameter in a simulation setup (decaying rate of coefficient)
aux.add = 5 ## the number of additional hypotheses for UD

iniseed=123456789 ## random see for simulation
smax=5                          ## s_{max} = s_0 + smax   
adddim = 2                      ## dimension allowance for the slack extractor
bmin=-0.8; bmax=0.8             ## simulation AR(1) operator coefficients ; a_j ~unif[-0.8, 0.8] for each iteration
iterset =c(seq(from = 7, 1, by = -2),0) ## dimension in the null of interest


discard=200

optm = TRUE; 
hr1<-1/3;     # bandwidth for slack extractor
hr2<-1/5;     # bandwidth for h_R If optm =TRUE, this line will be ignored.
hr_inv <- 1/5; # h_R for inverse VR test
hl_inv <- hr1; # h_L for inverse VR test
al_inv <- 1 ; # 0 if you would prefer to use the mere covariance for the inverse test.
hlspecified<-hr1 # bandwidth for h_L


nobslist = c(250,500)    ## The list of sample size
tobs =max(nobslist)


ker=5;

nt = 200;t = (0:(nt-1))/(nt-1); # Support of the function
neigen = 20   #number of basis used to generate process
adddim = 2    #slack extractor; additional dimension
lbnumber=40   # LBB generation
maxalterdim = 2         # When you study power; how many null hypotheses would you like to test?
 


lbb<- Legendre_basis(numeric_var("x", support = c(0, 1)),order = lbnumber)


source("code/vrtest_fun.R")  #Obtain the functions for computing test statistics


tobs =max(nobslist)
nt = 200;t = (0:(nt-1))/(nt-1); # Support of the function
neigen = 20   #number of basis used to generate process
adddim = 2    #slack extractor; additional dimension
lbnumber=40   # LBB generation

lbb<- Legendre_basis(numeric_var("x", support = c(0, 1)),order = lbnumber)
lb=lbb(t)
LB=matrix(0,nrow=length(t),ncol=lbnumber)
for(i in 2:lbnumber){  
  for(j in 1:i)  { 
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber){
  LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
} 

  
# Fourier Basis Function
lbnumber2=40
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}
LBF=cbind(rep(1,length(t)),LBF)



for (x_dim in iterset){
  
  set.seed(iniseed)
  
  dimestmat.10 <-dimestmat.20 <-dimestmat.10A <-dimestmat.20A <-dimestmat.21 <- dimestmat.31 <-dimestmat.30 <-dimestmat.41 <-dimestmat.40 <-array(NA, dim = c(length(nobslist), tsim, 4))
  dimestmat.ud.10 <-dimestmat.ud.20 <-dimestmat.ud.10A <-dimestmat.ud.20A <-dimestmat.ud.21 <-dimestmat.ud.31 <- dimestmat.ud.30 <-dimestmat.ud.41 <- dimestmat.ud.40 <-array(NA, dim = c(length(nobslist), tsim, 4))
  dimestmat.12 <- dimestmat.12A <-array(NA, dim = c(length(nobslist), tsim , 2))
  dimestmat.lrs <- array(NA, dim = c(length(nobslist), tsim , 3))
  
  for (iter in 1:tsim){
    eigenxx=runif(12,bmin,bmax)*0.9^(0:11)
    eigen_x  = append(runif(x_dim,bmin,bmax),eigenxx) 
    cutt = x_dim+3
    mulnumber=30
    armu= rnorm(mulnumber,0,1)*0.9^(0:(mulnumber-1))
    
    mu=LBF[,1:mulnumber]%*%as.matrix(armu) 
    rpert=append(sample(1:cutt,size=x_dim),sample((15):30,size=12))
    
    basis_x = LBF[,rpert] 
    
    ttobs = tobs+discard
    
    eta = matrix(NA, nrow = nt, ncol = tobs+discard)
    eta[,1] = LBF%*%cbind(rnorm(ncol(LBF),0,1)*ARDfac^(0:(ncol(LBF)-1)))
    
    x_mat <- pnx_mat <- psx_mat <- eta 
    x_mat[,1] =   eta[,1]
    
    dim_null = x_dim
    for (jiter in 2:ttobs){
      
      eta[,jiter] = LBF%*%cbind(rnorm(ncol(LBF),0,1)*ARDfac^(0:(ncol(LBF)-1))) 
      x_mat[,jiter] = operator(basis_x,eigen_x, x_mat[,jiter-1] ,t) + eta[,jiter] 
      aaaa=0
      bbbb=0      
      if (x_dim>0){
        for ( jjij in 1:x_dim)
        { 
          aaaa=aaaa+inner(x_mat[,jiter],basis_x[,jjij],t)*basis_x[,jjij]
        }}else{aaaa=0}
      pnx_mat[,jiter] = aaaa
      psx_mat[,jiter] = x_mat[,jiter] - aaaa    + mu  
    }
    nx_mat=t(apply(pnx_mat, 1, cumsum))
    x_mat= nx_mat+psx_mat
    
    X_mat=x_mat[,201:ttobs]
    dimest.20 <-dimest.21 <-dimest.10  <-dimest.ud.20 <-
      dimest.10A  <-dimest.ud.10A <-
      dimest.20A  <-dimest.ud.20A <-
      dimest.ud.21 <-dimest.ud.10  <- dimestlrs <- dimest.12 <-dimest.12A <-
      dimest.31  <-dimest.ud.31 <- dimest.30  <-dimest.ud.30 <-
      dimest.41  <-dimest.ud.41 <- dimest.40  <-dimest.ud.40 <-
      dimest.02 <-dimest.lrs<-NULL
    for(nobs in nobslist){
      x_mat=X_mat[,1:nobs] 
      xx_mat=x_mat-rowMeans(x_mat)       
      hh2=t(LB[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1])
      xcoef=t(hh2)
      
      
      ##Slack Extractor; VR test d_R=1
      lrx0 = lrvar(xcoef, aband = 0, bband = 0 ,kernel = ker  )
      eig.lrx0 = eigen(lrx0$omega/(nobs^2))$vectors 
      
      pdim = x_dim + adddim
      fdtmp0 = t(xcoef%*%eig.lrx0)
      fdtmp0 = fdtmp0 - rowMeans(fdtmp0)
      
      
      ##Slack extractor for VR test; d_R=0
      lrx1 = lrvar(xcoef, aband = 1, bband = hr1 ,kernel = ker )
      eig.lrx1 = eigen(lrx1$omega/(nobs^2))$vectors 
      
      fdtmp1 = t(xcoef%*%eig.lrx1)
      fdtmp1 = fdtmp1 - rowMeans(fdtmp1)
      
      
      ##Slack extractor for inverse test
      lrx2 = lrvar(xcoef, aband = 1, bband = hr1 ,kernel = ker )
      eig.lrx2 = eigen(lrx2$omega/(nobs^2))$vectors 
      
      fdtmp2 = t(xcoef%*%eig.lrx2)
      fdtmp2 = fdtmp2 - rowMeans(fdtmp2)
      
      a.fac = smax
      max.xdim =dim_null + a.fac
      a.dim =dim_null + a.fac - 1  ##a.fac - a.dim > 0 ##
      
      dimest.vr21<- dimest.vr30<- dimest.vr31 <-
        dimest.vr40<- dimest.vr41 <-   dimest.vr20 <-     dimest.vr10 <-dimest.ud.vr20<-dimest.ud.vr10 <-
        dimest.vr20A <-     dimest.vr10A <-dimest.ud.vr20A<-dimest.ud.vr10A <-
        dimest.ud.vr21     <-dimest.vr12<-
        dimest.ud.vr21A     <-dimest.vr12A<-
        dimestlrs <- NULL       
      
      ## Inverse test for the initial dimension conjecture
      dimest.vr12 =  t(sapply(c(max.xdim:0 + adddim), function(x){
        fdtmp.t = fdtmp2[1:x,]
        lftmp = t(apply(fdtmp.t,1,cumsum)); rftmp =  fdtmp.t
        vrtmp = inv.vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x-adddim, al = al_inv, ar=1,  hl=hl_inv, hr = hr_inv, kn = ker)
        dimest.vr12.t =  c(vrtmp$vr.t >CMatVR12t[(x-adddim+1),adddim,2] , vrtmp$vr.m >CMatVR12m[(x-adddim+1),adddim,2])   
        
        return( c(dimest.vr12.t )   )}))
      
      est.dim.inv = apply(rbind(dimest.vr12[(max.xdim+1):1, ],0),2,which.min) -1
      
      
     dimest.12  = cbind(dimest.12,est.dim.inv)
      
     
     dimest.vr12A =  t(sapply(c(max.xdim:0 + adddim), function(x){
       fdtmp.t = fdtmp2[1:x,]
       lftmp = t(apply(fdtmp.t,1,cumsum)); rftmp =  fdtmp.t
       vrtmp = inv.vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x-adddim, al = 0, ar=1,  hl=hl_inv, hr = hr_inv, kn = ker)
       dimest.vr12.t =  c(vrtmp$vr.t >CMatVR12t[(x-adddim+1),adddim,2] , vrtmp$vr.m >CMatVR12m[(x-adddim+1),adddim,2])   
       
       return( c(dimest.vr12.t )   )}))
     
     
     est.dim.inv = apply(rbind(dimest.vr12A[(max.xdim+1):1, ],0),2,which.min) -1
     
     dimest.12A  = cbind(dimest.12A,est.dim.inv)
     
     
      for (tdim in max.xdim:0){
        
        pdim   = tdim + adddim # slack extractor dimension for the alternative hypothesis
        
        fdptmp0 = fdtmp0[1:pdim,];
        fdptmp1 = fdtmp1[1:pdim,];      
        
        if (tdim >0){
          
          # VR test; d_R=1
          
          lftmp = t(apply(fdptmp0,1,cumsum)); rftmp =  fdptmp0
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker) 
          dimest.vr21 = rbind(dimest.vr21, c(vrtmp$vr.t >CMatVR21t[tdim,2] , vrtmp$vr.m >CMatVR21m[tdim,2], vrtmp$vr.est )   )
          
          
          lftmp = t(apply(lftmp,1,cumsum)); rftmp = fdptmp0
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 3, rfo = 1, al = 0, ar = 0, kn = ker)
          dimest.vr31 = rbind(dimest.vr31, c(vrtmp$vr.t >CMatVR31t[tdim,2] , vrtmp$vr.m >CMatVR31m[tdim,2], vrtmp$vr.est )   )
          
          lftmp = t(apply(lftmp,1,cumsum)); 
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 4, rfo = 1, al = 0, ar = 0, kn = ker)
          dimest.vr41 = rbind(dimest.vr41, c(vrtmp$vr.t >CMatVR41t[tdim,2] , vrtmp$vr.m >CMatVR41m[tdim,2], vrtmp$vr.est )   )
          
          # VR test; d_R = 0
          
          lftmp = t(apply(fdptmp1,1,cumsum))[,2:nobs]; rftmp = fdptmp1[,2:nobs] - fdptmp1[,1:(nobs-1)]  
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 0, al = 1,  kn = ker, ar=1,hl = hlspecified, hr = hr2, optimal = optm)
          dimest.vr20 = rbind(dimest.vr20, c(vrtmp$vr.t >CMatVR20t[tdim,2] , vrtmp$vr.m >CMatVR20m[tdim,2] , vrtmp$vr.est )   )
          
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 0, al = 0,  kn = ker, ar=1,hl = hlspecified, hr = hr2, optimal = optm)
          dimest.vr20A = rbind(dimest.vr20A, c(vrtmp$vr.t >CMatVR20t[tdim,2] , vrtmp$vr.m >CMatVR20m[tdim,2] , vrtmp$vr.est )   )
          
          
          lftmp = t(apply(lftmp,1,cumsum)); 
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 3, rfo = 0, al = 1,  kn = ker, ar=1,hl = hlspecified, hr = hr2,optimal = optm)
          dimest.vr30 = rbind(dimest.vr30, c(vrtmp$vr.t >CMatVR30t[tdim,2] , vrtmp$vr.m >CMatVR30m[tdim,2] , vrtmp$vr.est )   )
          
          lftmp = t(apply(lftmp,1,cumsum)); 
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 4, rfo = 0, al = 1,  kn = ker, ar=1,hl = hlspecified, hr = hr2 ,optimal = optm)
          dimest.vr40 = rbind(dimest.vr40, c(vrtmp$vr.t >CMatVR40t[tdim,2] , vrtmp$vr.m >CMatVR40m[tdim,2] , vrtmp$vr.est )   )
          
          
          lftmp = fdptmp1[,2:nobs]
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 1, rfo = 0, al = 1,  kn = ker, ar=1, hl = hlspecified,hr = hr2, optimal = optm)
          dimest.vr10 = rbind(dimest.vr10, c(vrtmp$vr.t >CMatVR10t[tdim,2] , vrtmp$vr.m >CMatVR10m[tdim,2] , vrtmp$vr.est )   )
          
          vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 1, rfo = 0, al = 0,  kn = ker, ar=1, hl = hlspecified,hr = hr2, optimal = optm)
          dimest.vr10A = rbind(dimest.vr10A, c(vrtmp$vr.t >CMatVR10t[tdim,2] , vrtmp$vr.m >CMatVR10m[tdim,2] , vrtmp$vr.est )   )
          
        }
      }
      
      dimest.21 = cbind(dimest.21,c((max.xdim+1)-apply(rbind(dimest.vr21[,1:2],0),2,which.min), dimest.vr21[c(1,1),3]))
      dimest.31 = cbind(dimest.31,c((max.xdim+1)-apply(rbind(dimest.vr31[,1:2],0),2,which.min), dimest.vr31[c(1,1),3]))
      dimest.30 = cbind(dimest.30,c((max.xdim+1)-apply(rbind(dimest.vr30[,1:2],0),2,which.min), dimest.vr30[c(1,1),3]))
      dimest.41 = cbind(dimest.41,c((max.xdim+1)-apply(rbind(dimest.vr41[,1:2],0),2,which.min), dimest.vr41[c(1,1),3]))
      dimest.40 = cbind(dimest.40,c((max.xdim+1)-apply(rbind(dimest.vr40[,1:2],0),2,which.min), dimest.vr40[c(1,1),3]))
      dimest.20 = cbind(dimest.20,c((max.xdim+1)-apply(rbind(dimest.vr20[,1:2],0),2,which.min), dimest.vr20[c(1,1),3]))
      dimest.10 = cbind(dimest.10,c((max.xdim+1)-apply(rbind(dimest.vr10[,1:2],0),2,which.min), dimest.vr10[c(1,1),3]))
      dimest.20A = cbind(dimest.20A,c((max.xdim+1)-apply(rbind(dimest.vr20A[,1:2],0),2,which.min), dimest.vr20A[c(1,1),3]))
      dimest.10A = cbind(dimest.10A,c((max.xdim+1)-apply(rbind(dimest.vr10A[,1:2],0),2,which.min), dimest.vr10A[c(1,1),3]))
      
      
      aux.max = (est.dim.inv[1]+aux.add)
      if (max.xdim >= aux.max  ){
        
        ## if the upper bound implied by the inverse test is lower than the initial hypothesis considered in the top down procedure, up-down = top-down
        dimest.ud.21 = cbind(dimest.ud.21,c((aux.max+1)-apply(rbind(dimest.vr21[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr21[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.31 = cbind(dimest.ud.31,c((aux.max+1)-apply(rbind(dimest.vr31[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr31[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.30 = cbind(dimest.ud.30,c((aux.max+1)-apply(rbind(dimest.vr30[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr30[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.41 = cbind(dimest.ud.41,c((aux.max+1)-apply(rbind(dimest.vr41[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr41[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.40 = cbind(dimest.ud.40,c((aux.max+1)-apply(rbind(dimest.vr40[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr40[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.20 = cbind(dimest.ud.20,c((aux.max+1)-apply(rbind(dimest.vr20[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr20[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.10 = cbind(dimest.ud.10,c((aux.max+1)-apply(rbind(dimest.vr10[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr10[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.20A = cbind(dimest.ud.20A,c((aux.max+1)-apply(rbind(dimest.vr20A[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr20A[c(1,(max.xdim - aux.max+1)),3]))
        dimest.ud.10A = cbind(dimest.ud.10A,c((aux.max+1)-apply(rbind(dimest.vr10A[(max.xdim - aux.max+1):max.xdim,1:2],0),2,which.min), dimest.vr10A[c(1,(max.xdim - aux.max+1)),3]))
        
      } else {
        
        for (tdim in (max.xdim+1):aux.max ){
          pdim   = tdim + adddim 
          fdptmp0 = fdtmp0[1:pdim,]
          fdptmp1 = fdtmp1[1:pdim,];     
          
          ##VR12##
          
          if (tdim >0){
            lftmp = t(apply(fdptmp0,1,cumsum)); rftmp =  fdptmp0
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker) 
            dimest.vr21 = rbind( c(vrtmp$vr.t >CMatVR21t[tdim,2] , vrtmp$vr.m >CMatVR21m[tdim,2], vrtmp$vr.est )  ,dimest.vr21 )
            
            lftmp = t(apply(lftmp,1,cumsum)); rftmp = fdptmp0
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 3, rfo = 1, al = 0, ar = 0, kn = ker)
            dimest.vr31 = rbind( c(vrtmp$vr.t >CMatVR31t[tdim,2] , vrtmp$vr.m >CMatVR31m[tdim,2], vrtmp$vr.est )  ,dimest.vr31 )
            
            lftmp = t(apply(lftmp,1,cumsum)); rftmp = fdptmp0
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 4, rfo = 1, al = 0, ar = 0, kn = ker)
            dimest.vr41 = rbind( c(vrtmp$vr.t >CMatVR41t[tdim,2] , vrtmp$vr.m >CMatVR41m[tdim,2], vrtmp$vr.est )  ,dimest.vr41 )
            
            
            lftmp = t(apply(fdptmp1,1,cumsum))[,2:nobs]; rftmp = fdptmp1[,2:nobs] - fdptmp1[,1:(nobs-1)]  
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 0, al = 1,   kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr20 = rbind(c(vrtmp$vr.t >CMatVR20t[tdim,2] , vrtmp$vr.m >CMatVR20m[tdim,2] , vrtmp$vr.est ),dimest.vr20   )
           
             vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 2, rfo = 0, al = 0,   kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr20A = rbind(c(vrtmp$vr.t >CMatVR20t[tdim,2] , vrtmp$vr.m >CMatVR20m[tdim,2] , vrtmp$vr.est ),dimest.vr20A   )
            
            lftmp =  t(apply(lftmp,1,cumsum))
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 3, rfo = 0, al = 1,   kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr30 = rbind(c(vrtmp$vr.t >CMatVR30t[tdim,2] , vrtmp$vr.m >CMatVR30m[tdim,2] , vrtmp$vr.est ),dimest.vr30   )
            
            lftmp =  t(apply(lftmp,1,cumsum))
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 4, rfo = 0, al = 1,   kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr40 = rbind(c(vrtmp$vr.t >CMatVR40t[tdim,2] , vrtmp$vr.m >CMatVR40m[tdim,2] , vrtmp$vr.est ),dimest.vr40   )
            
            lftmp = fdptmp1[,2:nobs]
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 1, rfo = 0, al = 1,  kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr10 = rbind( c(vrtmp$vr.t >CMatVR10t[tdim,2] , vrtmp$vr.m >CMatVR10m[tdim,2] , vrtmp$vr.est ) ,dimest.vr10  )
            
            vrtmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = tdim, lfo = 1, rfo = 0, al = 0,  kn = ker,  ar=1,hl=hlspecified, hr = hr2, optimal = optm)
            dimest.vr10A = rbind( c(vrtmp$vr.t >CMatVR10t[tdim,2] , vrtmp$vr.m >CMatVR10m[tdim,2] , vrtmp$vr.est ) ,dimest.vr10A  )
            
          }
        }
        
        
        dimest.ud.21 = cbind(dimest.ud.21,c((aux.max+1)-apply(rbind(dimest.vr21[,1:2],0),2,which.min), dimest.vr21[c(1,1),3]))
        dimest.ud.31 = cbind(dimest.ud.31,c((aux.max+1)-apply(rbind(dimest.vr31[,1:2],0),2,which.min), dimest.vr31[c(1,1),3]))
        dimest.ud.30 = cbind(dimest.ud.30,c((aux.max+1)-apply(rbind(dimest.vr30[,1:2],0),2,which.min), dimest.vr30[c(1,1),3]))
        dimest.ud.41 = cbind(dimest.ud.41,c((aux.max+1)-apply(rbind(dimest.vr41[,1:2],0),2,which.min), dimest.vr41[c(1,1),3]))
        dimest.ud.40 = cbind(dimest.ud.40,c((aux.max+1)-apply(rbind(dimest.vr40[,1:2],0),2,which.min), dimest.vr40[c(1,1),3]))
        dimest.ud.20 = cbind(dimest.ud.20,c((aux.max+1)-apply(rbind(dimest.vr20[,1:2],0),2,which.min), dimest.vr20[c(1,1),3]))
        dimest.ud.10 = cbind(dimest.ud.10,c((aux.max+1)-apply(rbind(dimest.vr10[,1:2],0),2,which.min), dimest.vr10[c(1,1),3]))
        dimest.ud.20A = cbind(dimest.ud.20A,c((aux.max+1)-apply(rbind(dimest.vr20A[,1:2],0),2,which.min), dimest.vr20A[c(1,1),3]))
        dimest.ud.10A = cbind(dimest.ud.10A,c((aux.max+1)-apply(rbind(dimest.vr10A[,1:2],0),2,which.min), dimest.vr10A[c(1,1),3]))
        
        
      }
      
      
      
      
      
      
      
      if (x_dim >0){
        del.value = 1e-4
        #########LRS(2020)#########
        pcacov = eigen(crossprod(xcoef), only.values = TRUE)
        eiglist = nobs*pcacov$values
        eiglist[which(eiglist/max(eiglist) < del.value)] = 0
        eiglist = eiglist[2:length(eiglist)]/eiglist[1:(length(eiglist)-1)]
        eiglist[which(is.nan(eiglist)==TRUE)] = 1
        dimest.lrs = cbind(dimest.lrs, c(which.min(eiglist[1:max.xdim]), which.min(eiglist[1:(x_dim+5)]), which.min(eiglist[1:(x_dim+3)]) ) )
      }
      
    } 
    dimestmat.12[,iter,] = t(dimest.12)
    dimestmat.12A[,iter,] = t(dimest.12A) 
    ;dimestmat.21[,iter,] = t(dimest.21);
    dimestmat.31[,iter,] = t(dimest.31);dimestmat.30[,iter,] = t(dimest.30);
    dimestmat.41[,iter,] = t(dimest.41);dimestmat.40[,iter,] = t(dimest.40);
    dimestmat.10[,iter,] = t(dimest.10); dimestmat.20[,iter,] =t(dimest.20)
    dimestmat.10A[,iter,] = t(dimest.10A); dimestmat.20A[,iter,] =t(dimest.20A)
    
    dimestmat.ud.21[,iter,] = t(dimest.ud.21); 
    dimestmat.ud.31[,iter,] = t(dimest.ud.31); 
    dimestmat.ud.30[,iter,] = t(dimest.ud.30); 
    dimestmat.ud.41[,iter,] = t(dimest.ud.41); 
    dimestmat.ud.40[,iter,] = t(dimest.ud.40); 
    dimestmat.ud.10[,iter,] = t(dimest.ud.10); 
    dimestmat.ud.20[,iter,] =t(dimest.ud.20)
    dimestmat.ud.10A[,iter,] = t(dimest.ud.10A); 
    dimestmat.ud.20A[,iter,] =t(dimest.ud.20A)
    
    
    if (x_dim >0) {
      dimestmat.lrs[,iter,] = t(dimest.lrs)}
    
    if(iter%%15==0){ 
      printarg = c(iter,x_dim,bmin,bmax,
                   mean(dimestmat.21[length(nobslist),1:iter,1]==x_dim),
                   mean(dimestmat.31[length(nobslist),1:iter,1]==x_dim),
                   mean(dimestmat.30[length(nobslist),1:iter,1]==x_dim),
                   mean(dimestmat.41[length(nobslist),1:iter,1]==x_dim),
                   mean(dimestmat.40[length(nobslist),1:iter,1]==x_dim)
                   ,mean(dimestmat.20[length(nobslist),1:iter,1]==x_dim)
                   ,mean(dimestmat.12[length(nobslist),1:iter,1]==x_dim)
                   ,mean(dimestmat.ud.21[length(nobslist),1:iter,1]==x_dim)
                   ,mean(dimestmat.ud.20[length(nobslist),1:iter,1]==x_dim))
      printarg1 = c(iter,x_dim,bmin,bmax
                    ,mean(dimestmat.21[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.31[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.30[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.41[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.40[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.20[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.12[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.ud.21[length(nobslist)-1,1:iter,1]==x_dim)
                    ,mean(dimestmat.ud.20[length(nobslist)-1,1:iter,1]==x_dim))			 
      printarg2 = c(iter,x_dim,bmin,bmax
                    ,mean(dimestmat.21[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.31[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.30[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.41[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.40[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.20[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.12[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.ud.21[1,1:iter,1]==x_dim)
                    ,mean(dimestmat.ud.20[1,1:iter,1]==x_dim))          
      
      
      print( round(printarg2,digits=3))  
      print( round(printarg1,digits=3))  
      print(  round(printarg,digits=3)) 
    }
  }
  save(file = paste0("result/",today,"Nobs",nobs,"x_dim",x_dim,"adddim",adddim,"Dimest","HR_Slack",round(hr1, digits=2),"HR2",round(hr2, digits=2),
                     "OPTM", optm,"AL_INV",al_inv, 
                     ".RData"),
       dimestmat.12,dimestmat.10,dimestmat.20,dimestmat.12A,
       dimestmat.21,dimestmat.ud.10,dimestmat.ud.20,dimestmat.10A,dimestmat.20A,
       dimestmat.ud.10A,dimestmat.ud.20A,
       dimestmat.31,dimestmat.30,      dimestmat.41,dimestmat.40,
       dimestmat.ud.21, dimestmat.ud.31, dimestmat.ud.30,dimestmat.ud.41, dimestmat.ud.40,  dimestmat.lrs)
  
  
}





