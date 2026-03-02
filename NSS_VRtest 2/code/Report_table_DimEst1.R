rm(list=ls())
library(xtable) 
library(MASS) 
library(ggplot2)
library(scales)



setwd("/");  
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
 
adddim = 2
ARD = 0.5
 


## Reported Ver. ##
today <- "result20260217/"
iniseed=123456789


setup2 = "cosines" # "sines" or "fourier"  
nobslist=c(250,500)
tobs = max(nobslist)

xdimlist = c(0,1,3,5,7)
bmin=-0.8 ; bmax = 0.8
tsim = 1e4
optimal = FALSE


optm = TRUE; 
hr1<-1/3;     # bandwidth for slack extractor
hr2<-1/5;     # bandwidth for h_R If optm =TRUE, this line will be ignored.
hr_inv <- 1/5; # h_R for inverse VR test
hl_inv <- hr1; # h_L for inverse VR test
al_inv <- 1 ; # 0 if you would prefer to use the mere covariance for the inverse test.
hlspecified<-hr1 # bandwidth for h_L


#Reported result : 1/3 and 1/4  
aux.add = 5




TRmat = NULL  
RRmat = NULL
for (nobs in nobslist ){
  
  Rmat = NULL
  for (x_dim in xdimlist){
    load(file = paste0("result/",today,"Nobs",max(nobslist),"x_dim",x_dim,"adddim",adddim,"Dimest","HR_Slack",round(hr1, digits=2),"HR2",round(hr2, digits=2),
                       "OPTM", optm,"AL_INV",al_inv, 
                       ".RData"))
    trace.test.mat = cbind(dimestmat.21[which(nobs==nobslist),,1],
                           dimestmat.ud.21[which(nobs==nobslist),,1],
                           dimestmat.20A[which(nobs==nobslist),,1],
                           dimestmat.10A[which(nobs==nobslist),,1],
                           dimestmat.12A[which(nobs==nobslist),,1],
                           dimestmat.ud.20A[which(nobs==nobslist),,1],
                           dimestmat.ud.10A[which(nobs==nobslist),,1], 
                           dimestmat.20[which(nobs==nobslist),,1],
                           dimestmat.10[which(nobs==nobslist),,1],  
                           dimestmat.12[which(nobs==nobslist),,1], 
                           dimestmat.ud.20[which(nobs==nobslist),,1],
                           dimestmat.ud.10[which(nobs==nobslist),,1])
    tmat = matrix(colMeans(trace.test.mat==x_dim),ncol=1)
 
    emat = matrix(rep(NA, 2), ncol=1)
    if (x_dim>0){
      est.dim.mat = cbind(dimestmat.21[which(nobs==nobslist),,c(3)], dimestmat.lrs[which(nobs==nobslist),,c(1)])
      emat = cbind(colMeans(est.dim.mat<x_dim),colMeans(est.dim.mat==x_dim),colMeans(est.dim.mat>x_dim))
      emat =  matrix(colMeans(est.dim.mat==x_dim), ncol = 1) 
    } 
    
    tmat = rbind(emat,tmat)    
    Rmat = cbind(Rmat, tmat)
  }
  RRmat = cbind(RRmat,NA, Rmat)  
  
}


rown.list = rep(c("$\\ddot{\\mathbbm{d}}$",
                  "$\\widehat{\\mathbbm{d}}_{\\text{\\tiny{LRS}}}$",
                  "$\\widehat{\\mathbbm{d}}_{\\rm TD}$: VR(2,1) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm UD}$: VR(2,1) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm TD}$: VR$_a$(2,0) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm TD}$: VR$_a$(1,0) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm BU}$: Inv.VR$_a$",  
                  "$\\widehat{\\mathbbm{d}}_{\\rm UD}$: VR$_a$(2,0) ", 
                  "$\\widehat{\\mathbbm{d}}_{\\rm UD}$: VR$_a$(1,0) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm TD}$: VR$_b$(2,0) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm TD}$: VR$_b$(1,0) ", 
                  "$\\widehat{\\mathbbm{d}}_{\\rm BU}$: Inv.VR$_b$ ", 
                  "$\\widehat{\\mathbbm{d}}_{\\rm UD}$: VR$_b$(2,0) ",
                  "$\\widehat{\\mathbbm{d}}_{\\rm UD}$: VR$_b$(1,0) "),1)


Rmat = matrix(paste0(" ", formatC(RRmat, digits =3 , format = "f")," "), ncol = ncol(RRmat))

res.table = cbind(rown.list, RRmat)

res.table = cbind(rown.list,matrix(paste0(" ", formatC(RRmat, digits =3 , format = "f")," "), ncol = ncol(RRmat)))

print(xtable(res.table,digits = c(0, 0,rep(3, ncol(res.table)-1))) ,sanitize.text.function = identity
      ,include.rownames = FALSE, file=paste0("result/Dimestimation",".tex"), na.string = "")

