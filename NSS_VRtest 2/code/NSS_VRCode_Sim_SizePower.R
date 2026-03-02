
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
aux.add = 2 ## the number of additional hypotheses for UD

iniseed=123456789 ## random see for simulation
smax=5                          ## s_{max} = s_0 + smax   
adddim = 2                      ## dimension allowance for the slack extractor
bmin=-0.8; bmax=0.8             ## simulation AR(1) operator coefficients ; a_j ~unif[-0.8, 0.8] for each iteration
iterset =c(seq(from = 7, 1, by = -2),0) ## dimension in the null of interest


discard=200

hr1<- 1/3   # bandwidth for slack extractor for d_R = 0
hr2<- 1/5   # right bandwidth for VR test with d_R = 0  # If obtm = TRUE, it will not be used.
hlspecified <- hr1 # left bandwidth for VR test with d_R = 0 
optm <- TRUE # optimal or no for the right bandwidth

 
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
   
  
  vr21result.t <-  vr20result.t <- vr10result.t <-   vr20Aresult.t <- vr10Aresult.t <- array(NA, dim = c(length(nobslist),tsim, 1+maxalterdim) )   
  vr21result.m <-  vr20result.m <- vr10result.m <-  vr20Aresult.m <- vr10Aresult.m <-   array(NA, dim = c(length(nobslist),tsim, 1+maxalterdim) )   
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
    
    x_mat = eta 
    pnx_mat = eta
    psx_mat = eta
    x_mat[,1] = eta[,1]
    
    dim_null = x_dim
    for (jiter in 2:(tobs+discard)){
      eta[,jiter] = LBF%*%cbind(rnorm(ncol(LBF),0,1)*ARDfac^(0:(ncol(LBF)-1))) 
      
      x_mat[,jiter] = operator(basis_x,eigen_x, x_mat[,jiter-1] ,t) + eta[,jiter] 
      aaaa=0;      bbbb=0   	
      
      if (x_dim>0){
        for ( jjij in 1:x_dim) { 
          aaaa=aaaa+inner(x_mat[,jiter],basis_x[,jjij],t)*basis_x[,jjij]
        }		   } else{aaaa=0}
      
      pnx_mat[,jiter] = aaaa
      psx_mat[,jiter] = x_mat[,jiter] - aaaa    + mu
    }
    
    nx_mat=t(apply(pnx_mat, 1, cumsum))
    x_mat= nx_mat+psx_mat
    
    X_mat=x_mat[,(discard+1):ncol(x_mat)] 
    
    vr21 <-vr20 <- vr10 <- vr20A <- vr10A <-  array(NA, dim = c(2, length(nobslist), maxalterdim+1) )  
    for(nobs in nobslist){
      x_mat=X_mat[,1:nobs] 
      xx_mat=x_mat-rowMeans(x_mat)        #Demean the time series
      
      hh2=t(LB[2:(nt),])%*%xx_mat[2:(nt),]*(t[2]-t[1]) # Inner product representation of X_t with respect to basis functions on the regular interval. It's taken to simplify computation.
      xcoef=t(hh2)
      
      pdim = x_dim + adddim  # Slack extractor dimension
      
      
      ## Slack Extractor for VR test d_R = 1
      lrx0 = lrvar(xcoef, aband = 0, bband = 0 ,kernel = ker  )
      eig.lrx0 = eigen(lrx0$omega/(nobs^2))$vectors 
      
      
      fdtmp0 = t(xcoef%*%eig.lrx0[,1:pdim])
      fdtmp0 = fdtmp0 - rowMeans(fdtmp0)
      
       
      ## Slack Extractor for VR test d_R = 0
      lrx1 = lrvar(xcoef, aband = 1, bband =  hr1 ,kernel = ker )
      eig.lrx1 = eigen(lrx1$omega/(nobs^2))$vectors 
      
      fdtmp = t(xcoef%*%eig.lrx1[,1:pdim])
      fdtmp = fdtmp - rowMeans(fdtmp)
      
# Make the below active if you would like to use the optimal bandwidth for the slack extractor (d_R = 0)      
#      op_tmp = opbw(fdtmp[,2:nobs] - fdtmp[,1:(nobs-1)]  )
#      lrx1 = lrvar(xcoef, aband = op_tmp, bband =  1/5 ,kernel = ker )
#      eig.lrx1 = eigen(lrx1$omega/(nobs^2))$vectors 
      
      fdtmp = t(xcoef%*%eig.lrx1[,1:pdim])
      fdtmp = fdtmp - rowMeans(fdtmp)
      
      lftmp = t(apply(fdtmp,1,cumsum)); rftmp = fdtmp 

      ##Size of the top-down procedure is obtainable only when x_dim >0##
      
      if (x_dim>0){ 
        #VR test with d_R =1; No bandwidth on both sides
        
        lftmp0 = t(apply(fdtmp0,1,cumsum)); rftmp0 = fdtmp0
        vr21tmp = vrtest(lfdata = lftmp0, rfdata = rftmp0, testdim = x_dim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker)
        vr21[,which(nobs==nobslist),1] = c(vr21tmp$vr.t, vr21tmp$vr.m)
        
        # VR test with d_R = 0; optimal bandwidth on the right
        
        lftmp = lftmp[,2:nobs]; rftmp = fdtmp[,2:nobs] - fdtmp[,1:(nobs-1)]  
        vr20tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x_dim, lfo = 2, rfo = 0, al = 1, hl=hlspecified,    kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr20[,which(nobs==nobslist),1] = c(vr20tmp$vr.t, vr20tmp$vr.m)
         
        vr20tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x_dim, lfo = 2, rfo = 0, al = 0,    kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr20A[,which(nobs==nobslist),1] = c(vr20tmp$vr.t, vr20tmp$vr.m)
        
        
 
        lftmp = fdtmp[,2:nobs]
        vr10tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x_dim, lfo = 1,rfo = 0, al = 1,  hl=hlspecified,  kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr10[,which(nobs==nobslist),1] = c(vr10tmp$vr.t, vr10tmp$vr.m)
        
        
        vr10tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = x_dim, lfo = 1,rfo = 0, al = 0,   kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr10A[,which(nobs==nobslist),1] = c(vr10tmp$vr.t, vr10tmp$vr.m)
      }
      
      
      
      
      
      for (p_dim in (x_dim+1):(x_dim+maxalterdim )){
        p_loc = p_dim-x_dim+1
        dim_y = p_dim+adddim 
        
        fdtmp0 = t(xcoef%*%eig.lrx0[,1:dim_y])
        fdtmp0 = fdtmp0 - rowMeans(fdtmp0)
        
        fdtmp = t(xcoef%*%eig.lrx1[,1:dim_y])
        fdtmp = fdtmp - rowMeans(fdtmp)
        
        
        ##VR test with d_R = 1;no bandwidth on both sides
        lftmp0 = t(apply(fdtmp0,1,cumsum)); rftmp0 = fdtmp0
        vr21tmp = vrtest(lfdata = lftmp0, rfdata = rftmp0, testdim = p_dim, lfo = 2, rfo = 1, al = 0, ar = 0, kn = ker)
        vr21[,which(nobs==nobslist),p_loc] = c(vr21tmp$vr.t, vr21tmp$vr.m)
        
      
        ## VR test with d_R=0; hl = hlspecified; h_R = h_opt
        
        lftmp = t(apply(fdtmp,1,cumsum)); rftmp = fdtmp 
        lftmp = lftmp[,2:nobs]; rftmp = fdtmp[,2:nobs] - fdtmp[,1:(nobs-1)]  
        vr20tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = p_dim, lfo = 2, rfo = 0,  al = 1, hl=hlspecified,  kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr20[,which(nobs==nobslist),p_loc] = c(vr20tmp$vr.t, vr20tmp$vr.m)
        
  
        vr20tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = p_dim, lfo = 2, rfo = 0,  al = 0,  kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr20A[,which(nobs==nobslist),p_loc] = c(vr20tmp$vr.t, vr20tmp$vr.m)
        
        lftmp = fdtmp[,2:nobs]
        vr10tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = p_dim, lfo = 1, rfo = 0,  al = 1, hl=hlspecified, kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr10[,which(nobs==nobslist),p_loc] = c(vr10tmp$vr.t, vr10tmp$vr.m)
        
        vr10tmp = vrtest(lfdata = lftmp, rfdata = rftmp, testdim = p_dim, lfo = 1, rfo = 0,  al = 0, kn = ker,  ar=1, hr = hr2, optimal = optm)
        vr10A[,which(nobs==nobslist),p_loc] = c(vr10tmp$vr.t, vr10tmp$vr.m)
      }
      
      
      
      
      
      
    } 
    
    vr21result.t[,iter,] =  vr21[1,,] ; vr21result.m[,iter,] =  vr21[2,,]
    vr20result.t[,iter,] =  vr20[1,,] ; vr20result.m[,iter,] =  vr20[2,,]
    vr10result.t[,iter,] =  vr10[1,,] ; vr10result.m[,iter,] =  vr10[2,,] 
    vr20Aresult.t[,iter,] =  vr20A[1,,] ; vr20Aresult.m[,iter,] =  vr20A[2,,]
    vr10Aresult.t[,iter,] =  vr10A[1,,] ; vr10Aresult.m[,iter,] =  vr10A[2,,] 
   
    
    lnobs=length(nobslist)
    if(iter%%60==0){
      if (x_dim==0){
        print( cbind(iter,bmin,bmax,hr1,x_dim,rbind( 
          c(mean(vr21result.t[1,1:iter,2] >CMatVR21t[1,2]),
            mean(vr20result.t[1,1:iter,2] >CMatVR20t[1,2]) ,
            mean(vr20Aresult.t[1,1:iter,2] >CMatVR20t[1,2]),
            mean(vr10result.t[1,1:iter,2] >CMatVR10t[1,2])
            ),
          c(mean(vr21result.t[2,1:iter,2] >CMatVR21t[1,2]),
            mean(vr20result.t[2,1:iter,2] >CMatVR20t[1,2]) ,
            mean(vr20Aresult.t[2,1:iter,2] >CMatVR20t[1,2]),
            mean(vr10result.t[2,1:iter,2] >CMatVR10t[1,2]) ),
          c(mean(vr21result.t[(lnobs),1:iter,2] >CMatVR21t[1,2]),
            mean(vr20result.t[(lnobs),1:iter,2] >CMatVR20t[1,2]) ,
            mean(vr20Aresult.t[(lnobs),1:iter,2] >CMatVR20t[1,2]),
            mean(vr10result.t[(lnobs),1:iter,2] >CMatVR10t[1,2]) )	
        ))
        ,digits=3)
      } else{
        print( cbind(iter,bmin,bmax,hr1,x_dim,rbind( 
          c(rowMeans(t(vr21result.t[1,1:iter,1:2]) >CMatVR21t[x_dim:(x_dim+1),2]),
            rowMeans(t(vr20result.t[1,1:iter,1:2]) >CMatVR20t[x_dim:(x_dim+1),2]) ,
            rowMeans(t(vr20Aresult.t[1,1:iter,1:2]) >CMatVR20t[x_dim:(x_dim+1),2]),
            rowMeans(t(vr10result.t[1,1:iter,1:2]) >CMatVR10t[x_dim:(x_dim+1),2])  ),
          c(rowMeans(t(vr21result.t[2,1:iter,1:2]) >CMatVR21t[x_dim:(x_dim+1),2]),
            rowMeans(t(vr20result.t[2,1:iter,1:2]) >CMatVR20t[x_dim:(x_dim+1),2]) ,
            rowMeans(t(vr20Aresult.t[2,1:iter,1:2]) >CMatVR20t[x_dim:(x_dim+1),2]) ,
            rowMeans(t(vr10result.t[2,1:iter,1:2]) >CMatVR10t[x_dim:(x_dim+1),2])  ) 
        ))
        ,digits=3)
      }
    }
  }
  
  save(file=paste0("result/",today,"Nobs",nobs,"x_dim",x_dim,"adddim",adddim,"Size_Power_","HR_Slack",round(hr1, digits = 2),"HR2",round(hr2, digits=2),"OPTM", optm,  ".RData"),
       vr21result.t,vr21result.m,
       vr20result.t,vr20result.m,
       vr20Aresult.t,vr20Aresult.m,
       vr10result.t,vr10result.m,
       vr10Aresult.t,vr10Aresult.m )
  
  
}





