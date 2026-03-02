rm(list=ls())
library(fda);library(AER);library(sandwich);library(geigen)



setwd("/");  
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))



optm = TRUE; 
hr1<-1/3;     # bandwidth for slack extractor
hr2<-1/5;     # bandwidth for h_R (Top-Down and Up-Down Tests) If optm =TRUE, this line will be ignored.
hr_inv <- 1/5; # h_R for inverse VR test
hl_inv <- hr1; # h_L for inverse VR test
al_inv <- 1 ; # 0 if you would prefer to use the mere covariance for the inverse test.
al_20 <- 1;  # 1 for longrun variance in VR(2,0)
al_10 <- 1;  # 1 for longrun variance in VR(1,0)
hlspecified<-hr1 # bandwidth for h_L
ker <- 5;  # Choice of kernel function

#######################################################
##The data is available at https://home.treasury.gov ## 
#######################################################
month.list = seq.Date(from = as.Date("1984/1/1"), by = "month", to = as.Date("2022/6/1"))

datmat = NULL
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_84_88_0.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 
datmat = as.matrix(dat[,3:ncol(dat)])

dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_89_93.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_94_98.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_99_03.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_04_08.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_09_13.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_14_18.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
dat = as.data.frame(readxl::read_xls(path = "data/hqmeom_19_23.xls",sheet = 1, 
                                     skip = 6, col_names = F)) 

datmat = cbind(datmat,as.matrix(dat[,3:ncol(dat)]))
colnames(datmat) = month.list
t.grid = dat$...1


rm(dat);dat = datmat;rm(datmat)

 
tobs = dim(dat)[2]
datt = dat 
t.grid = (t.grid )/(max(t.grid) )


## Generate the basis functions
library(polynom); library(basefun);
lbnumber= 40
lbb<- Legendre_basis(numeric_var("x", support = c(min(t.grid), max(t.grid))),order = lbnumber)
t = t.grid     ## for the consistency of the names in simulations


source("code/vrtest_fun.R") 
ydcurve = apply(dat,2,function(x){return(apply(LB  , 2, inner,x,t.grid)  )})  

 
trend = 1        ## trend = 0 if the model has no linear trend but a constant trend
max.xdim = 10    ## Maximum dimension to begin with
adddim = 2       ## allowance for the slack extractor
aux.add = 2      ## the number of additional hypotheses for UD procedure 
smax_est  = NULL ## chosen by the inverse procedure


## Load Critical Values and Remove Deterministic Trends
if (trend==0){   load("cvalue/criticalVR.RData")
  load("cvalue/criticalVR12.RData")
  fdata = (ydcurve) - rowMeans(ydcurve)  
  xcoef = t(fdata)
} else {   load("cvalue/criticalVRTrend.RData")
  load("cvalue/criticalVR12Trend.RData")
  expl = cbind(1,seq(1,tobs,1))   
  fdata = ydcurve  -  t(expl%*%solve(crossprod(expl))%*%t(expl)%*%t(ydcurve) )
  xcoef = t(fdata)
}


## VR Testing Procedure 
 
source("code/vrtest_procedure.R")

##################################################################################################
###############################REPORT THE ESTIMATED RESULTS#######################################
##################################################################################################
 
library(xtable)

ind.max = est.dim.inv[1] 


Rmat = apply(dimest.vr12.t,1,sigfunc) ; TDmat =  c(apply(dimest.vr21.t,1,sigfunc)) ; TDmat = rbind(TDmat,  c(apply(dimest.vr20.t,1,sigfunc) ))  ; TDmat = rbind(TDmat,  c(apply(dimest.vr10.t,1,sigfunc) )  ) 
rown.list = NULL
for (sss in 0:ind.max){    rown.list = append(rown.list,paste0("$\\smalls_0=",sss,"$"))  }

Rmat =  Rmat[1:(ind.max+1)] 
res.table = rbind(rown.list,Rmat,cbind(NA,TDmat))
rownames(res.table)<- c("", "Inv.VR ","VR(2,1)", "VR(2,0)","VR(1,0)") 




rown.list = NULL
for (sss in 1:(ind.max) ){    rown.list = append(rown.list,paste0("$j=",sss,"$")) }
RRmat = rbind(rown.list,  paste0(formatC(eigestvr21[1:ind.max], digits=2, format="f"), "") 
              ,paste0(formatC(dimest.lrs[1:ind.max], digits=2, format="f"), "") )
RRmat = cbind(NA, RRmat)
rownames(RRmat) <- c("Statistic", "NSS$_j$", "LRS$_j$") 

print(xtable(rbind(RRmat,res.table)) ,sanitize.text.function = identity
      ,include.rownames = T, file=paste0("result/main_est_op.txt")) 








#######################################################################################################################
###############################        Replication of Section 7.2      ###############################################
###############################  Table 4 corresponds to the result in Model 4##########################################
#######################################################################################################################
lamb.value.seq = c(1.37); ## Smoothing parameter for Nielsen-Siegel model

dat = datt; m.point = 30; tt.grid = t.grid[1:which(100*t.grid==m.point)] ## Truncation of maturities upto 30 years


## Generate Legendre Basis Functions on the restricted support

lbnumber= 40
lbb<- Legendre_basis(numeric_var("x", support = c(min(tt.grid), max(tt.grid))),order = lbnumber)
lb=lbb(tt.grid) 
LB=matrix(0,nrow=length(tt.grid),ncol=lbnumber)
for(i in 2:lbnumber)  {
  for(j in 1:i)     { 
    if (j != i) {lb[,i] = lb[,i]-(inner(lb[,i],lb[,j],tt.grid)/inner(lb[,j],lb[,j],tt.grid))*lb[,j]
    }}} 
for(i in 1:lbnumber)   {
  LB[,i] = lb[,i]/(sqrt(inner(lb[,i],lb[,i],tt.grid)))
} 



trend = 1 ;max.xdim = 5 ; smax_est=5 ; max.xdim.inv = 0

est.mat<-est.mat.1<- NULL
modeltyp = 1:7

## Model type 1:
## Model type 2: sp\{\varsigma_0}
## Model type 3: sp\{\varsigma_0. \varsigma_1}
## Model type 4: sp\{\varsigma_0. \varsigma_1, \varsigma_2}
## Model type 5: sp\{\varsigma_1}
## Model type 6: sp\{\varsigma_2}
## Model type 7: sp\{\varsigma_1, \varsigma_2}

dat2=dat[1:which(100*t.grid==m.point),]  ;res.vr12<- array(NA, dim = c(length(lamb.value.seq), 4,length(modeltyp)))
res.vr12.pv <- matrix(NA,nrow = length(lamb.value.seq), ncol = length(modeltyp) )
res.vr21<-res.vr20<-res.vr10<-res.eig.lrs <- res.eig.vr<- array(NA, dim = c(length(lamb.value.seq), 5,length(modeltyp)))
for (model in modeltyp){
  if (model==1){
    ydcurve =apply(dat2,2,function(x){return(apply(LB, 2, inner,x,tt.grid)  )}) 
    source("code/ydcur_function.R")
    res.vr12[,,model] = t(matrix(rep(dimest.vr12.t[1:4,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr21[,,model] = t(matrix(rep(dimest.vr21.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr20[,,model] = t(matrix(rep(dimest.vr20.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr10[,,model] = t(matrix(rep(dimest.vr10.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.eig.lrs[,,model] = t(matrix(rep(dimest.lrs, length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.eig.vr[,,model] = t(matrix(rep(eigestvr21[1:5], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    
    adddim.inv = 1
    ydcurve =t(matrix(apply(dat2,2,inner, LB, tt.grid))) 
    source("code/ydcur_function_invonly.R")
    res.vr12.pv[,model] = dimest.vr12.t[1,1]
    
  }else  if (model==2){
    ## Because the basis functions are orthogonalized, it is equivalent to exclude the first score
    
    ydcurve =apply(dat2,2,function(x){return(apply(LB[,-1] , 2, inner,x,tt.grid)  )}) 
    source("code/ydcur_function.R")
    res.vr12[,,model] = t(matrix(rep(dimest.vr12.t[1:4,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr21[,,model] = t(matrix(rep(dimest.vr21.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr20[,,model] = t(matrix(rep(dimest.vr20.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.vr10[,,model] = t(matrix(rep(dimest.vr10.t[,1], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.eig.lrs[,,model] = t(matrix(rep(dimest.lrs, length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    res.eig.vr[,,model] = t(matrix(rep(eigestvr21[1:5], length(lamb.value.seq)), ncol = length(lamb.value.seq)))
    
    adddim.inv = 1
    ydcurve =t(matrix(apply(dat2,2,inner, LB[,1], tt.grid))) 
    source("code/ydcur_function_invonly.R")
    res.vr12.pv[,model] = dimest.vr12.t[1,1]
    
  } else if (model ==3 | model==4) {
    
    ## test.mat includes basis functions spanning H under the null
    for (lamb.value in lamb.value.seq  ){
      lam.order= which(lamb.value==lamb.value.seq)
      slop.com =  lamb.value*(1-exp(-(100*tt.grid)/lamb.value))/(100*tt.grid) # Slope component in NS model
      curv.com =  lamb.value*(1-exp(-(100*tt.grid)/lamb.value))/(100*tt.grid)  -exp(-(100*tt.grid)/lamb.value)   #Curvature component in NS model
      
      if (model==3){  test.mat = c(cbind(slop.com))  } else if (model==4){   test.mat = cbind(slop.com, curv.com)  }

      # representation of the testing component into Legendre basis functions      
      if (is.vector(test.mat)){
        Bmat =apply(LB[,-1],2,inner,test.mat, tt.grid)
      }else{
        Bmat =  t(apply(LB[,-1],2,function(x){return(apply((test.mat),2,inner,x,tt.grid))} )) 
      }

      proj.mat =   (Bmat)%*%solve((crossprod((Bmat)) ))%*%t(Bmat) # Projection matrix
      ydcurve =apply(dat2,2,function(x){return(apply(LB[,-1], 2, inner,x,tt.grid)  )})
      ydcurve = (diag(nrow(proj.mat)) - proj.mat)%*%ydcurve # (I-P)X; orthogonalized curve represented by basis functions
      
      
      ## Apply the test to the orthogonalized process
      smax_est=5
      source("code/ydcur_function.R")
      
      res.vr12[lam.order,,model] = dimest.vr12.t[1:4,1] ; res.vr21[lam.order,,model] =  dimest.vr21.t[,1] 
      res.vr20[lam.order,,model] = dimest.vr20.t[,1] ;res.vr10[lam.order,,model] =  dimest.vr10.t[,1] 
      res.eig.lrs[lam.order,,model]=  dimest.lrs ; res.eig.vr[lam.order,,model] =  eigestvr21[1:5] 
      
      
      # Inverse Test applied to P*X
      
      adddim.inv = 1 ;      ydcurve =apply(dat2,2,function(x){return(apply(LB[,-1], 2, inner,x,tt.grid)  )}) ;      ydcurve = (proj.mat)%*%ydcurve
      
      source("code/ydcur_function_invonly.R")
      res.vr12.pv[which(lamb.value==lamb.value.seq),model] = dimest.vr12.t[1,1]
      
    }    }else  {
    for (lamb.value in lamb.value.seq  ){
      lam.order= which(lamb.value==lamb.value.seq)
      slop.com =  lamb.value*(1-exp(-(100*tt.grid)/lamb.value))/(100*tt.grid) 
      curv.com =  lamb.value*(1-exp(-(100*tt.grid)/lamb.value))/(100*tt.grid)  -exp(-(100*tt.grid)/lamb.value)   
      
      if (model==5){
        test.mat = slop.com
        Bmat =  as.matrix(apply(LB,2,function(x){return(apply(as.matrix(test.mat),2,inner,x,tt.grid))} )) 
      } else if (model==6){
        test.mat = curv.com
        Bmat =  as.matrix(apply(LB,2,function(x){return(apply(as.matrix(test.mat),2,inner,x,tt.grid))} )) 
      } else if (model==7){ 
        test.mat = cbind(slop.com, curv.com)
        Bmat =  t(apply(LB,2,function(x){return(apply((test.mat),2,inner,x,tt.grid))} )) 
      }
      Cmat =  (crossprod((Bmat)) ) 
      
      proj.mat =   (Bmat)%*%solve(Cmat)%*%t(Bmat)
      ydcurve =apply(dat2,2,function(x){return(apply(LB, 2, inner,x,tt.grid)  )})
      ydcurve = (diag(nrow(proj.mat)) - proj.mat)%*%ydcurve
      
      smax_est=5
      source("code/ydcur_function.R")
      
      res.vr12[lam.order,,model] = dimest.vr12.t[1:4,1] ;   res.vr21[lam.order,,model] =  dimest.vr21.t[,1] 
      res.vr20[lam.order,,model] = dimest.vr20.t[,1] ;      res.vr10[lam.order,,model] =  dimest.vr10.t[,1] 
      res.eig.lrs[lam.order,,model]=  dimest.lrs ;    res.eig.vr[lam.order,,model] =  eigestvr21[1:5] 

      # Inverse Test applied to P*X
      adddim.inv = 1;      ydcurve =apply(dat2,2,function(x){return(apply(LB, 2, inner,x,tt.grid)  )}) ;      ydcurve = (proj.mat)%*%ydcurve
      
      source("code/ydcur_function_invonly.R")
      res.vr12.pv[which(lamb.value==lamb.value.seq),model] = dimest.vr12.t[1,1]
      
      
    }    }}




# Significance level function (** 1% significance * 5% significance level)

sigfunc = function(x, dg=2){
  if(x[2]==3){s_level = "$^{**}$"} else if(x[2]==2){s_level = "$^{*}$\\phantom{$^{*}$}"}  else if(x[2]==1){s_level = "$^{\\dagger}$\\phantom{$^{*}$}"} else  {s_level = "\\phantom{$^{**}$}"}
  paste0(formatC(x[1], digits =dg , format = "f"),s_level)
}
 

for (model in modeltyp){
  if (model<=2){
    lamlist = 1
  } else {
    lamlist = 1:length(lamb.value.seq)
  }
  
  for (lambo in lamlist){
    sigl = rowSums(res.vr12[lambo,,model]>(CMatVR12t[1:4,adddim,]))
    Rmat = apply( cbind(res.vr12[lambo,,model], sigl),1,sigfunc )
    
    sigl = rowSums(res.vr21[lambo,,model]>(CMatVR21t[1:5,]))
    TDmat = apply( cbind(res.vr21[lambo,,model], sigl),1,sigfunc ) 
    sigl = rowSums(res.vr20[lambo,,model]>(CMatVR20t[1:5,]))
    TDmat = rbind(TDmat,apply( cbind(res.vr20[lambo,,model], sigl),1,sigfunc ) )
    sigl = rowSums(res.vr10[lambo,,model]>(CMatVR10t[1:5,]))
    TDmat = rbind(TDmat,apply( cbind(res.vr10[lambo,,model], sigl),1,sigfunc ) )
    
    rown.list = NULL
    for (sss in 0:5){
      rown.list = append(rown.list,paste0("$s_0=",sss,"$\\phantom{$^{**}$}"))
    }
    
    res.table = rbind(rown.list,c(Rmat, rep(NA,2)),cbind(NA,TDmat))
    rownames(res.table)<- c("\\midrule
\\multicolumn{7}{c}{Variance Ratio Test Statistics}\\\\
Statistic",
                            "\\midrule Inv.VR ",
                            "VR(2,1)",
                            "VR(2,0)",
                            "VR(1,0)") 
    
    
    
    rown.list = NULL
    for (sss in 1:5){
      rown.list = append(rown.list,paste0("$j=",sss,"$\\phantom{$^{**}$}"))
    }
    RRmat = rbind(rown.list,  paste0(formatC(res.eig.vr[lambo,,model], digits=2, format="f"), "\\phantom{$^{**}$}") 
                  ,paste0(formatC(res.eig.lrs[lambo,,model], digits=2, format="f"), "\\phantom{$^{**}$}") )
    RRmat = cbind(NA, RRmat)
    rownames(RRmat) <- c("Statistic", "\\midrule $\\mu_{j+1}/\\mu_j$", "$\\hat{\\kappa}_j/\\max_j\\{\\hat{\\kappa}_j  \\}$") 
    
    print(xtable(rbind(RRmat,res.table)) ,sanitize.text.function = identity
          ,include.rownames = T, file=paste0("result/eigen_est_model_",model,"lamb",
                                             100* lamb.value.seq[lambo],"00.txt")) 
    
  }}







####################################################################################################################
###################################################RESULTS  for Testing (5.1)#######################################
###############################################Testing Stationarity of the Projected Time Series#####################
#####################################################################################################################

modelseq = c(2,5,6,3,7,4)

for (llseq in 1:length(lamb.value.seq)){
  
  sigl = sapply(res.vr12.pv[llseq,modelseq],function(x){return(sum(x>CMatVR12t[1,adddim.inv,]))})
  Rmat = apply( cbind(res.vr12.pv[llseq,modelseq], sigl),1,sigfunc, dg = 4 )
  
  rown.list = c("$\\spn\\{\\varsigma_0\\}$","$\\spn\\{\\varsigma_1\\}$","$\\spn\\{\\varsigma_2\\}$",
                "$\\spn\\{\\varsigma_0,\\varsigma_1\\}$",
                "$\\spn\\{\\varsigma_1,\\varsigma_2\\}$",
                "$\\spn\\{\\varsigma_0,\\varsigma_1,\\varsigma_2\\}$")
  
  res.table = rbind(rown.list,c(Rmat))
  rownames(res.table)<- c("\\midrule
Statistic",
                          "\\midrule Inv.VR ") 
  
  
  print(xtable(rbind(res.table)) ,sanitize.text.function = identity
        ,include.rownames = T, file=paste0("result/eigen_est_","lamb",
                                           100* lamb.value.seq[llseq],"pv.txt")) 
  
}


