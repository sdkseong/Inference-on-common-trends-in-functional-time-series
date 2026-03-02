rm(list=ls())
library(xtable) 
library(MASS) 
library(ggplot2)
library(scales)


setwd("/");  
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
 
adddim = 2 

load("cvalue/criticalVR.RData")
load("cvalue/criticalVR12.RData")

## Reported Ver. ##
today <- "result20260217/"
iniseed=123456789

nobslist=c(250,500)
tobs = max(nobslist)
 
optm = TRUE; 
hr1<-1/3;     # bandwidth for slack extractor
hr2<-1/5;     # bandwidth for h_R If optm =TRUE, this line will be ignored.
hr_inv <- 1/5; # h_R for inverse VR test
hl_inv <- hr1; # h_L for inverse VR test
al_inv <- 1 ; # 0 if you would prefer to use the mere covariance for the inverse test.
hlspecified<-hr1 # bandwidth for h_L


tsim <-1e4
TRmat = NULL

xdimlist = c(0,1,3,5,7)
for (tp in 3){
  
  RRmat = NULL
  
  
  for (nobs in nobslist[c(1,2)]){
    Rmat <- NULL
    SAP_VR21 <- SAP_VR20<-SAP_VR10<-SAP_VR20B<-SAP_VR10B<-  NULL
    for (x_dim  in xdimlist){
      
      load(file=paste0("result/",today,"Nobs",max(nobslist),"x_dim",x_dim,"adddim",adddim,"Size_Power_","HR_Slack",round(hr1, digits = 2),"HR2",round(hr2, digits=2),"OPTM", optm,  ".RData"))
      if (x_dim==0){
        SAP_VR21 = rbind(SAP_VR21, c(NA,rowMeans(t(vr21result.t[which(nobs==nobslist),,1+1:2])>CMatVR21t[(1+x_dim:(x_dim+1)),2])))
        SAP_VR20 = rbind(SAP_VR20, c(NA,rowMeans(t(vr20result.t[which(nobs==nobslist),,1+1:2])>CMatVR20t[(1+x_dim:(x_dim+1)),2])))
        SAP_VR10 = rbind(SAP_VR10, c(NA,rowMeans(t(vr10result.t[which(nobs==nobslist),,1+1:2])>CMatVR10t[(1+x_dim:(x_dim+1)),2])))
        SAP_VR20B = rbind(SAP_VR20B, c(NA,rowMeans(t(vr20Aresult.t[which(nobs==nobslist),,1+1:2])>CMatVR20t[(1+x_dim:(x_dim+1)),2])))
        SAP_VR10B = rbind(SAP_VR10B, c(NA,rowMeans(t(vr10Aresult.t[which(nobs==nobslist),,1+1:2])>CMatVR10t[(1+x_dim:(x_dim+1)),2])))
        
      }else{
        SAP_VR21 = rbind(SAP_VR21, rowMeans(t(vr21result.t[which(nobs==nobslist),,1:3])>CMatVR21t[x_dim:(x_dim+2),2]))
        SAP_VR20 = rbind(SAP_VR20, rowMeans(t(vr20result.t[which(nobs==nobslist),,1:3])>CMatVR20t[x_dim:(x_dim+2),2]))
        SAP_VR10 = rbind(SAP_VR10, rowMeans(t(vr10result.t[which(nobs==nobslist),,1:3])>CMatVR10t[x_dim:(x_dim+2),2]))
        SAP_VR20B = rbind(SAP_VR20B, rowMeans(t(vr20Aresult.t[which(nobs==nobslist),,1:3])>CMatVR20t[x_dim:(x_dim+2),2]))
        SAP_VR10B = rbind(SAP_VR10B, rowMeans(t(vr10Aresult.t[which(nobs==nobslist),,1:3])>CMatVR10t[x_dim:(x_dim+2),2]))
        
      }
      
    } 
    
    SAP_VR21 = t(SAP_VR21);SAP_VR20 = t(SAP_VR20);SAP_VR10 = t(SAP_VR10);SAP_VR20B = t(SAP_VR20B);SAP_VR10B = t(SAP_VR10B)
    
    
    SAP_VR21tmp= paste0("{", formatC(SAP_VR21[1,2:ncol(SAP_VR21)], digits = 3, format = "f"),"}")
    SAP_VR21tmp= rbind(c(NA,SAP_VR21tmp), formatC(SAP_VR21[2:nrow(SAP_VR21),], digits = 3, format = "f"))
    SAP_VR20tmp= paste0("{", formatC(SAP_VR20[1,2:ncol(SAP_VR20)], digits = 3, format = "f"),"}")
    SAP_VR20tmp= rbind(c(NA,SAP_VR20tmp), formatC(SAP_VR20[2:nrow(SAP_VR20),], digits = 3, format = "f"))
    SAP_VR10tmp= paste0("{", formatC(SAP_VR10[1,2:ncol(SAP_VR10)], digits = 3, format = "f"),"}")
    SAP_VR10tmp= rbind(c(NA,SAP_VR10tmp), formatC(SAP_VR10[2:nrow(SAP_VR10),], digits = 3, format = "f"))
    SAP_VR20Btmp= paste0("{", formatC(SAP_VR20B[1,2:ncol(SAP_VR20B)], digits = 3, format = "f"),"}")
    SAP_VR20Btmp= rbind(c(NA,SAP_VR20Btmp), formatC(SAP_VR20B[2:nrow(SAP_VR20B),], digits = 3, format = "f"))
    SAP_VR10Btmp= paste0("{", formatC(SAP_VR10B[1,2:ncol(SAP_VR10B)], digits = 3, format = "f"),"}")
    SAP_VR10Btmp= rbind(c(NA,SAP_VR10Btmp), formatC(SAP_VR10B[2:nrow(SAP_VR10B),], digits = 3, format = "f"))
    
    rnn =  paste0("$s =",xdimlist,"$") 
    Rmat = cbind(Rmat,rbind(rnn,SAP_VR21tmp,SAP_VR20Btmp,SAP_VR10Btmp,SAP_VR20tmp,SAP_VR10tmp)) 
    if (nobs==nobslist[1]){
      RRmat = cbind(RRmat,  Rmat)
    }else{RRmat = cbind(RRmat,NA, Rmat)}
  }
  
  pdim = 2
  rown = c("${\\smalls}_N$",rep(c("${\\mathbbm{d}_{N}}$",paste0("${\\mathbbm{d}_{N}}+",1:(2),"$")),5))
  rownn = c(NA,"VR(2,1)",rep(NA,pdim),"\\midrule VR$_a$(2,0)",rep(NA,pdim),"\\midrule VR$_a$(1,0)",rep(NA,pdim),"\\midrule VR$_b$(2,0)",rep(NA,pdim),
            "\\midrule VR$_b$(1,0)",rep(NA,pdim) )
  
  RRmat = cbind( rownn,rown,RRmat) 
  
  TRmat = rbind(TRmat, RRmat)
}


print(xtable(TRmat,digits = c(rep(0,3), rep(3, ncol(TRmat)-2))) ,sanitize.text.function = identity
      ,include.rownames = FALSE, file=paste0("result/SizePower",".tex"), na.string = "")
 
