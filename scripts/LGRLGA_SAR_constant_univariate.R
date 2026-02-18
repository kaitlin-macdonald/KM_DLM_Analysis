rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)


SAR_raw<- readRDS("data/df_cbr_LGRLGA_apr.rds")

SAR<- SAR_raw[-28,]
################################
##Salmon Survive test model
year_SAR<- SAR$year
logit.s_SAR<- SAR$logit.s
CUI_z<- scale(SAR$CUI)
CUTI_z<- scale(SAR$CUTI)
BEUTI_z<- scale(SAR$BEUTI)
mldno3_z<- scale(SAR$Apr.mldno3)
CUI_mean_z<- scale(SAR$mean_CUI)
CUTI_mean_z<- scale(SAR$mean_CUTI)
BEUTI_mean_z<- scale(SAR$mean_BEUTI)
mldno3_mean_z<- scale(SAR$mean_mldno3)
outflow_z<- scale(SAR$mean.out)
temp_z<- scale(SAR$mean.temp)
PNI_z<- scale(SAR$Annual.PNI)

##########################
#Run Reference model
buildFun_SAR_ref<- function(parm){
  mod_SAR_ref<- dlmModPoly(1) 
  V(mod_SAR_ref)<- exp(parm[1])
  diag(W(mod_SAR_ref))[1]<- exp(parm[2])
  return(mod_SAR_ref)
  
}

fit_SAR_ref<- dlmMLE(logit.s_SAR, parm = rep(0,2), build=buildFun_SAR_ref, hessian=TRUE)
conv_SAR_ref<- fit_SAR_ref$convergence

dlmLogits_SAR_ref<- buildFun_SAR_ref(fit_SAR_ref$par)

V(dlmLogits_SAR_ref)

W(dlmLogits_SAR_ref)

mod_SAR_ref<-buildFun_SAR_ref(fit_SAR_ref$par)

#Apply Kalman filter to the model
SAR_reffilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_ref)

#Mean absolute percentage error
MAPE_SAR_ref<- mean(abs((SAR_reffilter$f-logit.s_SAR)/logit.s_SAR))
rmse_SAR_ref<- rmse(logit.s_SAR, SAR_reffilter$f)

#dlmLL caculates the neg. LL (lower is better)
loglik_SAR_ref <- dlmLL(logit.s_SAR, dlmModPoly(1))

n.coef <- 2
SAR.ref.aic <- (2 * (loglik_SAR_ref)) + 2 * (sum(n.coef))

resids_ref <- residuals(SAR_reffilter, sd = FALSE)
plot.ts(resids_ref, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_ref)
qqline(resids_ref)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_ref, mu = 0)$p.value

## plot ACF of innovations
acf(resids_ref, lag.max = 10)

hist(resids_ref)

tsdiag(SAR_reffilter)

##########################
#Run CUI model
buildFun_SAR_CUI<- function(parm){
  mod_SAR_CUI<- dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE)
  V(mod_SAR_CUI)<- exp(parm[1])
  diag(W(mod_SAR_CUI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_CUI)
  
}

fit_SAR_CUI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_CUI, hessian=TRUE)
conv_SAR_CUI<- fit_SAR_CUI$convergence

dlmLogits_SAR_CUI<- buildFun_SAR_CUI(fit_SAR_CUI$par)

V(dlmLogits_SAR_CUI)

W(dlmLogits_SAR_CUI)

mod_SAR_CUI<- buildFun_SAR_CUI(fit_SAR_CUI$par)

#Apply Kalman filter to the model
SAR_CUIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUI)

#Mean absolute percentage error
MAPE_SAR_CUI<- mean(abs((SAR_CUIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_CUI<- rmse(logit.s_SAR, SAR_CUIfilter$f)

#dlmLL caculates the neg. LL (lower is better)
loglik_SAR_CUI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE))

Href_SAR_CUI<- 2*(loglik_SAR_ref - loglik_SAR_CUI)

n.coef <- 3
SAR.CUI.aic <- (2 * (loglik_SAR_CUI)) + 2 * (sum(n.coef))  
#r.bic <- (2 * (loglik_SAR_CUI)) + (log(length(logit.s))) * (n.coef)

cui_mat<- as.matrix(cbind(rep(1, 42), CUI))
SAR_CUI_state<- CUIfilter$m

###########################################
buildFun_SAR_CUTI<- function(parm){
  mod_SAR_CUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE)
  V(mod_SAR_CUTI)<- exp(parm[1])
  diag(W(mod_SAR_CUTI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_CUTI)
  
}

fit_SAR_CUTI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_CUTI, hessian=TRUE)
conv_SAR_CUTI<- fit_SAR_CUTI$convergence

dlmLogits_SAR_CUTI<- buildFun_SAR_CUTI(fit_SAR_CUTI$par)

V(dlmLogits_SAR_CUTI)

W(dlmLogits_SAR_CUTI)

mod_SAR_CUTI<- buildFun_SAR_CUTI(fit_SAR_CUTI$par)

SAR_CUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUTI)

#What to do about first obs??
MAPE_SAR_CUTI<- mean(abs((SAR_CUTIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_CUTI<- rmse(logit.s_SAR, SAR_CUTIfilter$f)

loglik_SAR_CUTI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE))

n.coef<- 3
SAR.CUTI.aic <- (2 * (loglik_SAR_CUTI)) + 2 * (sum(n.coef))  

###########################################
buildFun_SAR_BEUTI<- function(parm){
  mod_SAR_BEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE)
  V(mod_SAR_BEUTI)<- exp(parm[1])
  diag(W(mod_SAR_BEUTI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_BEUTI)
  
}

fit_SAR_BEUTI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_BEUTI, hessian=TRUE)
conv_SAR_BEUTI<- fit_SAR_BEUTI$convergence

dlmLogits_SAR_BEUTI<- buildFun_SAR_BEUTI(fit_SAR_BEUTI$par)

V(dlmLogits_SAR_BEUTI)

W(dlmLogits_SAR_BEUTI)

mod_SAR_BEUTI<- buildFun_SAR_BEUTI(fit_SAR_BEUTI$par)

SAR_BEUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_BEUTI)

#What to do about first obs??
MAPE_SAR_BEUTI<- mean(abs((SAR_BEUTIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_BEUTI<- rmse(logit.s_SAR, SAR_BEUTIfilter$f)

loglik_SAR_BEUTI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE))

SAR.BEUTI.aic <- (2 * (loglik_SAR_BEUTI)) + 2 * (sum(n.coef))  

###########################################
buildFun_SAR_mldno3<- function(parm){
  mod_SAR_mldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE)
  V(mod_SAR_mldno3)<- exp(parm[1])
  diag(W(mod_SAR_mldno3))[1:2]<- exp(parm[2:3])
  return(mod_SAR_mldno3)
  
}

fit_SAR_mldno3<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_mldno3, hessian=TRUE)
conv_SAR_mldno3<- fit_SAR_mldno3$convergence

dlmLogits_SAR_mldno3<- buildFun_SAR_mldno3(fit_SAR_mldno3$par)

V(dlmLogits_SAR_mldno3)

W(dlmLogits_SAR_mldno3)

mod_SAR_mldno3<- buildFun_SAR_mldno3(fit_SAR_mldno3$par)

SAR_mldno3filter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mldno3)

#What to do about first obs??
MAPE_SAR_mldno3<- mean(abs((SAR_mldno3filter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_mldno3<- rmse(logit.s_SAR, SAR_mldno3filter$f)

loglik_SAR_mldno3 <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE))

SAR.mldno3.aic <- (2 * (loglik_SAR_mldno3)) + 2 * (sum(n.coef))  

###########################################
buildFun_SAR_mnCUI<- function(parm){
  mod_SAR_mnCUI<- dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE)
  V(mod_SAR_mnCUI)<- exp(parm[1])
  diag(W(mod_SAR_mnCUI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_mnCUI)
  
}

fit_SAR_mnCUI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_mnCUI, hessian=TRUE)
conv_SAR_mnCUI<- fit_SAR_mnCUI$convergence

dlmLogits_SAR_mnCUI<- buildFun_SAR_mnCUI(fit_SAR_mnCUI$par)

V(dlmLogits_SAR_mnCUI)

W(dlmLogits_SAR_mnCUI)

mod_SAR_mnCUI<- buildFun_SAR_mnCUI(fit_SAR_mnCUI$par)

SAR_mnCUIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mnCUI)

#What to do about first obs??
MAPE_SAR_mnCUI<- mean(abs((SAR_mnCUIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_mnCUI<- rmse(logit.s_SAR, SAR_mnCUIfilter$f)

loglik_SAR_mnCUI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE))

SAR.mnCUI.aic <- (2 * (loglik_SAR_mnCUI)) + 2 * (sum(n.coef))  

###########################################
buildFun_SAR_mnCUTI<- function(parm){
  mod_SAR_mnCUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_mean_z, addInt=FALSE)
  V(mod_SAR_mnCUTI)<- exp(parm[1])
  diag(W(mod_SAR_mnCUTI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_mnCUTI)
  
}

fit_SAR_mnCUTI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_mnCUTI, hessian=TRUE)
conv_SAR_mnCUTI<- fit_SAR_mnCUTI$convergence

dlmLogits_SAR_mnCUTI<- buildFun_SAR_mnCUTI(fit_SAR_mnCUTI$par)

V(dlmLogits_SAR_mnCUTI)

W(dlmLogits_SAR_mnCUTI)

mod_SAR_mnCUTI<- buildFun_SAR_mnCUTI(fit_SAR_mnCUTI$par)

SAR_mnCUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mnCUTI)

#What to do about first obs??
MAPE_SAR_mnCUTI<- mean(abs((SAR_mnCUTIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_mnCUTI<- rmse(logit.s_SAR, SAR_mnCUTIfilter$f)

loglik_SAR_mnCUTI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=CUTI_mean_z, addInt=FALSE))

SAR.mnCUTI.aic <- (2 * (loglik_SAR_mnCUTI)) + 2 * (sum(n.coef)) 

###########################################
buildFun_SAR_mnBEUTI<- function(parm){
  mod_SAR_mnBEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE)
  V(mod_SAR_mnBEUTI)<- exp(parm[1])
  diag(W(mod_SAR_mnBEUTI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_mnBEUTI)
  
}

fit_SAR_mnBEUTI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_mnBEUTI, hessian=TRUE)
conv_SAR_mnBEUTI<- fit_SAR_mnBEUTI$convergence

dlmLogits_SAR_mnBEUTI<- buildFun_SAR_mnBEUTI(fit_SAR_mnBEUTI$par)

V(dlmLogits_SAR_mnBEUTI)

W(dlmLogits_SAR_mnBEUTI)

mod_SAR_mnBEUTI<- buildFun_SAR_mnBEUTI(fit_SAR_mnBEUTI$par)

SAR_mnBEUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mnBEUTI)

#What to do about first obs??
MAPE_SAR_mnBEUTI<- mean(abs((SAR_mnBEUTIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_mnBEUTI<- rmse(logit.s_SAR, SAR_mnBEUTIfilter$f)

loglik_SAR_mnBEUTI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE))

SAR.mnBEUTI.aic <- (2 * (loglik_SAR_mnBEUTI)) + 2 * (sum(n.coef)) 

###########################################
buildFun_SAR_mnmldno3<- function(parm){
  mod_SAR_mnmldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE)
  V(mod_SAR_mnmldno3)<- exp(parm[1])
  diag(W(mod_SAR_mnmldno3))[1:2]<- exp(parm[2:3])
  return(mod_SAR_mnmldno3)
  
}

fit_SAR_mnmldno3<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_mnmldno3, hessian=TRUE)
conv_SAR_mnmldno3<- fit_SAR_mnmldno3$convergence

dlmLogits_SAR_mnmldno3<- buildFun_SAR_mnmldno3(fit_SAR_mnmldno3$par)

V(dlmLogits_SAR_mnmldno3)

W(dlmLogits_SAR_mnmldno3)

mod_SAR_mnmldno3<- buildFun_SAR_mnmldno3(fit_SAR_mnmldno3$par)

SAR_mnmldno3filter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mnmldno3)

#What to do about first obs??
MAPE_SAR_mnmldno3<- mean(abs((SAR_mnmldno3filter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_mnmldno3<- rmse(logit.s_SAR, SAR_mnmldno3filter$f)

loglik_SAR_mnmldno3 <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE))

SAR.mnmldno3.aic <- (2 * (loglik_SAR_mnmldno3)) + 2 * (sum(n.coef)) 

###########################################
buildFun_SAR_outflow<- function(parm){
  mod_SAR_outflow<- dlmModPoly(1) + dlmModReg(X=outflow_z, addInt=FALSE)
  V(mod_SAR_outflow)<- exp(parm[1])
  diag(W(mod_SAR_outflow))[1:2]<- exp(parm[2:3])
  return(mod_SAR_outflow)
  
}

fit_SAR_outflow<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_outflow, hessian=TRUE)
conv_SAR_outflow<- fit_SAR_outflow$convergence

dlmLogits_SAR_outflow<- buildFun_SAR_outflow(fit_SAR_outflow$par)

V(dlmLogits_SAR_outflow)

W(dlmLogits_SAR_outflow)

mod_SAR_outflow<- buildFun_SAR_outflow(fit_SAR_outflow$par)

SAR_outflowfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_outflow)

#What to do about first obs??
MAPE_SAR_outflow<- mean(abs((SAR_outflowfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_outflow<- rmse(logit.s_SAR, SAR_outflowfilter$f)

loglik_SAR_outflow <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=outflow_z, addInt=FALSE))

SAR.outlfow.aic <- (2 * (loglik_SAR_outflow)) + 2 * (sum(n.coef))  

###########################################
buildFun_SAR_temp<- function(parm){
  mod_SAR_temp<- dlmModPoly(1) + dlmModReg(X=temp_z, addInt=FALSE)
  V(mod_SAR_temp)<- exp(parm[1])
  diag(W(mod_SAR_temp))[1:2]<- exp(parm[2:3])
  return(mod_SAR_temp)
  
}

fit_SAR_temp<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_temp, hessian=TRUE)
conv_SAR_temp<- fit_SAR_temp$convergence

dlmLogits_SAR_temp<- buildFun_SAR_temp(fit_SAR_temp$par)

V(dlmLogits_SAR_temp)

W(dlmLogits_SAR_temp)

mod_SAR_temp<- buildFun_SAR_temp(fit_SAR_temp$par)

SAR_tempfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_temp)

#What to do about first obs??
MAPE_SAR_temp<- mean(abs((SAR_tempfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_temp<- rmse(logit.s_SAR, SAR_tempfilter$f)

loglik_SAR_temp <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=temp_z, addInt=FALSE))

SAR.temp.aic <- (2 * (loglik_SAR_temp)) + 2 * (sum(n.coef))  

#################################

buildFun_SAR_PNI<- function(parm){
  mod_SAR_PNI<- dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE)
  V(mod_SAR_PNI)<- exp(parm[1])
  diag(W(mod_SAR_PNI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_PNI)
  
}

fit_SAR_PNI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_PNI, hessian=TRUE)
conv_SAR_PNI<- fit_SAR_PNI$convergence

dlmLogits_SAR_PNI<- buildFun_SAR_PNI(fit_SAR_PNI$par)

V(dlmLogits_SAR_PNI)

W(dlmLogits_SAR_PNI)

mod_SAR_PNI<- buildFun_SAR_PNI(fit_SAR_PNI$par)

SAR_PNIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_PNI)

#What to do about first obs??
MAPE_SAR_PNI<- mean(abs((SAR_PNIfilter$f-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SAR_PNI<- rmse(logit.s_SAR, SAR_PNIfilter$f)

loglik_SAR_PNI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE))

SAR.PNI.aic <- (2 * (loglik_SAR_PNI)) + 2 * (sum(n.coef))  

###################################################################################
#linear trend
buildFun<- function(parm){
  mod_SAR_CUI<- dlmModPoly(2) + dlmModReg(X=CUI_z, addInt=FALSE)
  V(mod_SAR_CUI)<- exp(parm[1])
  diag(W(mod_SAR_CUI))[1:3]<- exp(parm[2:4])
  return(mod_SAR_CUI)
  
}

fit_SAR_CUI<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun, hessian=TRUE)
conv_SAR_CUI<- fit_SAR_CUI$convergence

dlmLogits_SAR_CUI<- buildFun(fit_SAR_CUI$par)

V(dlmLogits_SAR_CUI)

W(dlmLogits_SAR_CUI)

mod_SAR_CUI<- buildFun(fit_SAR_CUI$par)

CUIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUI)

loglik_SAR_CUI <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=CUI_z, addInt=FALSE))

n.coef <- 2
r.aic <- (2 * (loglik_SAR_CUI)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik_SAR_CUI)) + (log(length(logit.s))) * (n.coef)

CUIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUI)

cui_mat<- as.matrix(cbind(rep(1, 42), CUI))
SAR_CUI_state<- CUIfilter$m

###########################################
buildFun_SAR_CUTI_lin<- function(parm){
  mod_SAR_CUTI<- dlmModPoly(2) + dlmModReg(X=CUTI_z, addInt=FALSE)
  V(mod_SAR_CUTI)<- exp(parm[1])
  diag(W(mod_SAR_CUTI))[1:3]<- exp(parm[2:4])
  return(mod_SAR_CUTI)
  
}

fit_SAR_CUTI<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun_SAR_CUTI_lin, hessian=TRUE)
conv_SAR_CUTI<- fit_SAR_CUTI$convergence

dlmLogits_SAR_CUTI<- buildFun(fit_SAR_CUTI$par)

V(dlmLogits_SAR_CUTI)

W(dlmLogits_SAR_CUTI)

mod_SAR_CUTI<- buildFun(fit_SAR_CUTI$par)

SAR_CUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUTI)

loglik_SAR_CUTI <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=CUTI_z, addInt=FALSE))

###########################################
buildFun<- function(parm){
  mod_SAR_BEUTI<- dlmModPoly(2) + dlmModReg(X=BEUTI_z, addInt=FALSE)
  V(mod_SAR_BEUTI)<- exp(parm[1])
  diag(W(mod_SAR_BEUTI))[1:3]<- exp(parm[2:4])
  return(mod_SAR_BEUTI)
  
}

fit_SAR_BEUTI<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun, hessian=TRUE)
conv_SAR_BEUTI<- fit_SAR_BEUTI$convergence

dlmLogits_SAR_BEUTI<- buildFun(fit_SAR_BEUTI$par)

V(dlmLogits_SAR_BEUTI)

W(dlmLogits_SAR_BEUTI)

mod_SAR_BEUTI<- buildFun(fit_SAR_BEUTI$par)

SAR_BEUTIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_BEUTI)

loglik_SAR_BEUTI <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=BEUTI_z, addInt=FALSE))

###########################################
buildFun_SAR_mldno3<- function(parm){
  mod_SAR_mldno3<- dlmModPoly(2) + dlmModReg(X=mldno3_z, addInt=FALSE)
  V(mod_SAR_mldno3)<- exp(parm[1])
  diag(W(mod_SAR_mldno3))[1:3]<- exp(parm[2:4])
  return(mod_SAR_mldno3)
  
}

fit_SAR_mldno3<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun_SAR_mldno3, hessian=TRUE)
conv_SAR_mldno3<- fit_SAR_mldno3$convergence

dlmLogits_SAR_mldno3<- buildFun(fit_SAR_mldno3$par)

V(dlmLogits_SAR_mldno3)

W(dlmLogits_SAR_mldno3)

mod_SAR_mldno3<- buildFun(fit_SAR_mldno3$par)

SAR_mldno3filter<- dlmFilter(logit.s_SAR, mod=mod_SAR_mldno3)

loglik_SAR_mldno3 <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=mldno3_z, addInt=FALSE))

###########################################
buildFun_SAR_outflow<- function(parm){
  mod_SAR_outflow<- dlmModPoly(2) + dlmModReg(X=outflow_z, addInt=FALSE)
  V(mod_SAR_outflow)<- exp(parm[1])
  diag(W(mod_SAR_outflow))[1:3]<- exp(parm[2:4])
  return(mod_SAR_outflow)
  
}

fit_SAR_outflow<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun_SAR_outflow, hessian=TRUE)
conv_SAR_outflow<- fit_SAR_outflow$convergence

dlmLogits_SAR_outflow<- buildFun(fit_SAR_outflow$par)

V(dlmLogits_SAR_outflow)

W(dlmLogits_SAR_outflow)

mod_SAR_outflow<- buildFun(fit_SAR_outflow$par)

SAR_outflowfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_outflow)

loglik_SAR_outflow <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=outflow_z, addInt=FALSE))

###########################################
buildFun_SAR_temp<- function(parm){
  mod_SAR_temp<- dlmModPoly(2) + dlmModReg(X=temp_z, addInt=FALSE)
  V(mod_SAR_temp)<- exp(parm[1])
  diag(W(mod_SAR_temp))[1:3]<- exp(parm[2:4])
  return(mod_SAR_temp)
  
}

fit_SAR_temp<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun_SAR_temp, hessian=TRUE)
conv_SAR_temp<- fit_SAR_temp$convergence

dlmLogits_SAR_temp<- buildFun(fit_SAR_temp$par)

V(dlmLogits_SAR_temp)

W(dlmLogits_SAR_temp)

mod_SAR_temp<- buildFun(fit_SAR_temp$par)

SAR_tempfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_temp)

loglik_SAR_temp <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=temp_z, addInt=FALSE))

#################################

buildFun_SAR_PNI<- function(parm){
  mod_SAR_PNI<- dlmModPoly(2) + dlmModReg(X=PNI_z, addInt=FALSE)
  V(mod_SAR_PNI)<- exp(parm[1])
  diag(W(mod_SAR_PNI))[1:3]<- exp(parm[2:4])
  return(mod_SAR_PNI)
  
}

fit_SAR_PNI<- dlmMLE(logit.s_SAR, parm = rep(0,4), build=buildFun_SAR_PNI, hessian=TRUE)
conv_SAR_PNI<- fit_SAR_PNI$convergence

dlmLogits_SAR_PNI<- buildFun(fit_SAR_PNI$par)

V(dlmLogits_SAR_PNI)

W(dlmLogits_SAR_PNI)

mod_SAR_PNI<- buildFun(fit_SAR_PNI$par)

SAR_PNIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_PNI)

loglik_SAR_PNI <- dlmLL(logit.s_SAR, dlmModPoly(2) + dlmModReg(X=PNI_z, addInt=FALSE))