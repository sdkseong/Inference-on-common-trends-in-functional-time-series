
vrtest = function(lfdata, rfdata, lfo, rfo, testdim, al=0, ar=0, hl=1/3, hr=1/3, kn, optimal = FALSE){
  nt = ncol(lfdata)
  hc = 2*lfo - rfo - 1
  if (optimal == FALSE){
    Lmat = lrvar(t(lfdata), aband = al, bband = hl, kernel = kn)
    Rmat = lrvar(t(rfdata), aband = ar, bband = hr, kernel = kn)
  } else{ 
    Lmat = lrvar(t(lfdata), aband = al, bband = hl, kernel = kn)
    Rmat = lrvar(t(rfdata), aband =  opbw(rfdata), bband = 1/5, kernel = kn)
  }
  lfac = (Lmat$intker*(Lmat$bandw)*(al>0) +1*(al==0))
  rfac = 1*(rfo==0) + (rfo>0)*(Rmat$intker*(Rmat$bandw)*(ar>0) +1*(ar==0)  )
  nt = (nt^hc)*lfac/rfac
  vr21tmp = nt*sort(eigen(solve(Lmat$omega)%*%Rmat$omega,only.values = TRUE)$values)
  est.rank = which.max(vr21tmp[2:length(vr21tmp)]/vr21tmp[1:(length(vr21tmp)-1)])
  return(list(vr.t = sum(vr21tmp[1:testdim]), vr.m =  max(vr21tmp[1:testdim])  , vr.est = est.rank ))
}

inv.vrtest = function(lfdata, rfdata, testdim, al=0, ar=1, hl=1/3, hr=1/3,kn){
  nt = ncol(lfdata) ; nk = nrow(lfdata)
  Lmat = lrvar(t(lfdata), aband = al, bband = hl, kernel = kn)
  Rmat = lrvar(t(rfdata), aband = ar, bband = hr, kernel = kn)
  
  lfac = (Lmat$intker*(Lmat$bandw)*(al>0) +1*(al==0)) 
  nt = (nt)*lfac 
  vr21tmp = 1/(nt*sort(eigen(solve(Lmat$omega)%*%Rmat$omega,only.values = TRUE)$values))
  if (testdim==0){
    vr.t =  sum(vr21tmp)
    vr.m =  max(vr21tmp)
  } else {
  vr.t =  sum(vr21tmp[-(1:testdim)])  
  vr.m =  max(vr21tmp[-(1:testdim)])  
  }
  return(list(vr.t = vr.t, vr.m = vr.m, nt = nt ))
}



## Identical to VRtest function but it only reports the estimated eigenvalue ratio in detail.
eigest = function(lfdata, rfdata, lfo, rfo, testdim, al=0, ar=0, hl=1/3, hr=1/3, kn){
  nt = ncol(lfdata)
  hc = 2*lfo - rfo - 1
  Lmat = lrvar(t(lfdata), aband = al, bband = hl, kernel = kn)
  Rmat = lrvar(t(rfdata), aband = ar, bband = hr, kernel = kn)
  lfac = (Lmat$intker*(Lmat$bandw)*(al>0) +1*(al==0))
  rfac = 1*(rfo==0) + (rfo>0)*(Rmat$intker*(Rmat$bandw)*(ar>0) +1*(ar==0)  )
  nt = (nt^hc)*lfac/rfac
  vr21tmp = nt*sort(eigen(solve(Lmat$omega)%*%Rmat$omega,only.values = TRUE)$values)
  est.rank =  (vr21tmp[2:length(vr21tmp)]/vr21tmp[1:(length(vr21tmp)-1)])
  return( est.rank )
}



inner = function(f,g,grid){
  h = f*g
  return(sum((0.5*h[1:(length(grid)-1)] + 0.5*h[2:(length(grid))])*(grid[2] - grid[1])))
}


opbw = function(X){
  p = dim(X)[1]
  n = dim(X)[2]
  
  rho = rep(0, p)
  s2 = rep(0, p)
  for (i in c(1:p)){
    atemp = X[i, 1:(n-1)]
    rho[i] = t(X[i, 2:n]) %*% (atemp)/(t(atemp) %*% (atemp))
    btemp = X[i, 2:n] - rho[i] * X[i, 1:(n-1)]
    s2[i] = t(btemp) %*% btemp/(n-1)
  }
  num = sum(4*rho^2 * s2^2/((1-rho)^8))
  den = sum(s2^2/(1 - rho)^4)
  op_bw_fac = 1.7462*(num/den)^(1/5) 
  
  return(op_bw_fac)
  
}


operator = function(basis_temp,eigen,f,grid){
  temp = 0 
  for (k in 1:length(eigen)){
    temp = temp + eigen[k]*inner(basis_temp[,k],f,grid)*basis_temp[,k]
  }
  return(temp)
}

lrvar = function(u, kernel = 2, aband = 1, bband = 1/3, urprint = 1, white = 0 ){
    tu <- nrow(u)
    p <- ncol(u)
    bandw = round(aband*(tu^(bband)))
    if (white == 1){
      te <- tu-1
      au <- qr.solve(as.matrix(u[1:te,]),as.matrix(u[2:tu,]))
      e <- as.matrix(u[2:tu,]) - as.matrix(u[1:te,])%*%au
    }else{
      e <- u
      te <- tu
    }
    
    # Estimate Covariances #
    if (bandw >0){
      if (kernel == 1){                             # Parzen kernel #
        tm <- floor(bandw)
        if (tm > 0){
          jb <- as.matrix(seq(1,tm,1)/bandw)
          kern <- (1 - (jb^2)*6 + (jb^3)*6)*(jb <= .5)
          kern <- kern + ((1-jb)^3)*(jb > .5)*2
        }
        intker = 0.75
      } else if (kernel == 2){                       ####Biweight
        tm <- floor(bandw)
        if (tm >0){
          jb = as.matrix(seq(1,tm,1)/bandw)
          #   kern <- .75*(1-jb^2)
          kern <- (1-jb^2)^2
        }
        intker = 16/15
      }else if (kernel == 3){                 ####Triweight
        tm <- floor(bandw)
        if (tm >0){
          jb = as.matrix(seq(1,tm,1)/bandw) 
          kern <- (1-jb^2)^3
        }
        intker = 32/35
      }else if (kernel == 4){                 ####Cosine
        tm <- floor(bandw)
        if (tm >0){
          jb = as.matrix(seq(1,tm,1)/bandw) 
          kern <-  cos(jb*pi/2)  
        }
        intker = 4/pi
      }else if (kernel == 5){                ####Tukey Hanning 
        tm <- floor(bandw)
        if (tm >0){
          jb = as.matrix(seq(1,tm,1)/bandw) 
          kern <-  (1+cos(jb*pi))/2  
        }
        intker = 1
      }else if (kernel == 6){               ###quadweight   
        tm <- floor(bandw)
        if (tm >0){
          jb = as.matrix(seq(1,tm,1)/bandw) 
          kern <- (1-jb^2)^4
        }
        intker = 256/315
      }
      
      
      lam <- matrix(0,p,p)
      for (j in 1:tm){
        kj <- kern[j]
        lam <- lam + (t(as.matrix(e[1:(te-j),]))%*%as.matrix(e[(1+j):te,]))*kj
      }
      omega <- (t(e)%*%e + lam + t(lam)) 
      
      if (white == 1){
        eau <- solve(diag(p) - au)
        omega <- t(eau)%*%omega%*%eau
      }
    } else {
      omega = crossprod(e)
      intker = 0
    }
    list(omega=omega,bandw=bandw,intker = intker)
  }




##Orthonormalization of the Legendre polynomial basis functions
lb=lbb(t)
LB=matrix(0,nrow=length(t),ncol=lbnumber)
for(i in 2:lbnumber){  
  for(j in 1:i)  { 
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],t)/inner(lb[,j],lb[,j],t))*lb[,j]  }}}

for(i in 1:lbnumber){
  LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],t)))
} 

rm(lbb)
rm(lb)


sigfunc = function(x){
  if(x[2]==3){s_level = "$^{**}$"} else if(x[2]==2){s_level = "$^{*}$ "}  else if (x[2]==1){s_level = "$^{\\dagger}$"} else {s_level = ""}
  paste0(formatC(x[1], digits =2 , format = "f"),s_level)
}


