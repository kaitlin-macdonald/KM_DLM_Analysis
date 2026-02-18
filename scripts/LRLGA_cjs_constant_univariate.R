rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)

cjs_raw<- readRDS("data/df_cjs_LGRLGA_apr.rds")
#Remove NA at end of time series
cjs<- cjs_raw[-28,]
################################
##Salmon Survive test model
year_cjs<- cjs$year
logit.s_cjs<- cjs$logit.s
CUI_z<- scale(cjs$CUI)
CUTI_z<- scale(cjs$CUTI)
BEUTI_z<- scale(cjs$BEUTI)
mldno3_z<- scale(cjs$Apr.mldno3)
CUI_mean_z<- scale(cjs$mean_CUI)
CUTI_mean_z<- scale(cjs$mean_CUTI)
BEUTI_mean_z<- scale(cjs$mean_BEUTI)
mldno3_mean_z<- scale(cjs$mean_mldno3)
outflow_z<- scale(cjs$mean.out)
temp_z<- scale(cjs$mean.temp)
PNI_z<- scale(cjs$Annual.PNI)

##########################
#Run Reference model
buildFun_cjs_ref<- function(parm){
  mod_cjs_ref<- dlmModPoly(1) 
  V(mod_cjs_ref)<- exp(parm[1])
  diag(W(mod_cjs_ref))[1]<- exp(parm[2])
  return(mod_cjs_ref)
  
}

fit_cjs_ref<- dlmMLE(logit.s_cjs, parm = rep(0,2), build=buildFun_cjs_ref, hessian=TRUE)
conv_cjs_ref<- fit_cjs_ref$convergence

dlmLogits_cjs_ref<- buildFun_cjs_ref(fit_cjs_ref$par)

V(dlmLogits_cjs_ref)

W(dlmLogits_cjs_ref)

mod_cjs_ref<-buildFun_cjs_ref(fit_cjs_ref$par)

#Apply Kalman filter to the model
cjs_reffilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_ref)

#Mean absolute percentage error
MAPE_cjs_ref<- mean(abs((cjs_reffilter$f-logit.s_cjs)/logit.s_cjs))
rmse_cjs_ref<- rmse(logit.s_cjs, cjs_reffilter$f)

#dlmLL caculates the neg. LL (lower is better)
loglik_cjs_ref <- dlmLL(logit.s_cjs, dlmModPoly(1))

n.coef <- 2
cjs.ref.aic <- (2 * (loglik_cjs_ref)) + 2 * (sum(n.coef))

resids_ref <- residuals(cjs_reffilter, sd = FALSE)
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

tsdiag(cjs_reffilter)

##########################
#Run CUI model
buildFun_cjs_CUI<- function(parm){
  mod_cjs_CUI<- dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE)
  V(mod_cjs_CUI)<- exp(parm[1])
  diag(W(mod_cjs_CUI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_CUI)
  
}

fit_cjs_CUI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_CUI, hessian=TRUE)
conv_cjs_CUI<- fit_cjs_CUI$convergence

dlmLogits_cjs_CUI<- buildFun_cjs_CUI(fit_cjs_CUI$par)

V(dlmLogits_cjs_CUI)

W(dlmLogits_cjs_CUI)

mod_cjs_CUI<- buildFun_cjs_CUI(fit_cjs_CUI$par)

#Apply Kalman filter to the model
cjs_CUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUI)

#Mean absolute percentage error
MAPE_cjs_CUI<- mean(abs((cjs_CUIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_CUI<- rmse(logit.s_cjs, cjs_CUIfilter$f)

#dlmLL caculates the neg. LL (lower is better)
loglik_cjs_CUI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE))

Href_cjs_CUI<- 2*(loglik_cjs_ref - loglik_cjs_CUI)

n.coef <- 3
cjs.CUI.aic <- (2 * (loglik_cjs_CUI)) + 2 * (sum(n.coef))  
#r.bic <- (2 * (loglik_cjs_CUI)) + (log(length(logit.s))) * (n.coef)

cui_mat<- as.matrix(cbind(rep(1, 42), CUI))
cjs_CUI_state<- CUIfilter$m

###########################################
buildFun_cjs_CUTI<- function(parm){
  mod_cjs_CUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE)
  V(mod_cjs_CUTI)<- exp(parm[1])
  diag(W(mod_cjs_CUTI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_CUTI)
  
}

fit_cjs_CUTI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_CUTI, hessian=TRUE)
conv_cjs_CUTI<- fit_cjs_CUTI$convergence

dlmLogits_cjs_CUTI<- buildFun_cjs_CUTI(fit_cjs_CUTI$par)

V(dlmLogits_cjs_CUTI)

W(dlmLogits_cjs_CUTI)

mod_cjs_CUTI<- buildFun_cjs_CUTI(fit_cjs_CUTI$par)

cjs_CUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUTI)

#What to do about first obs??
MAPE_cjs_CUTI<- mean(abs((cjs_CUTIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_CUTI<- rmse(logit.s_cjs, cjs_CUTIfilter$f)

loglik_cjs_CUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE))

n.coef<- 3
cjs.CUTI.aic <- (2 * (loglik_cjs_CUTI)) + 2 * (sum(n.coef))  

###########################################
buildFun_cjs_BEUTI<- function(parm){
  mod_cjs_BEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE)
  V(mod_cjs_BEUTI)<- exp(parm[1])
  diag(W(mod_cjs_BEUTI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_BEUTI)
  
}

fit_cjs_BEUTI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_BEUTI, hessian=TRUE)
conv_cjs_BEUTI<- fit_cjs_BEUTI$convergence

dlmLogits_cjs_BEUTI<- buildFun_cjs_BEUTI(fit_cjs_BEUTI$par)

V(dlmLogits_cjs_BEUTI)

W(dlmLogits_cjs_BEUTI)

mod_cjs_BEUTI<- buildFun_cjs_BEUTI(fit_cjs_BEUTI$par)

cjs_BEUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_BEUTI)

#What to do about first obs??
MAPE_cjs_BEUTI<- mean(abs((cjs_BEUTIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_BEUTI<- rmse(logit.s_cjs, cjs_BEUTIfilter$f)

loglik_cjs_BEUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE))

cjs.BEUTI.aic <- (2 * (loglik_cjs_BEUTI)) + 2 * (sum(n.coef))  

###########################################
buildFun_cjs_mldno3<- function(parm){
  mod_cjs_mldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE)
  V(mod_cjs_mldno3)<- exp(parm[1])
  diag(W(mod_cjs_mldno3))[1:2]<- exp(parm[2:3])
  return(mod_cjs_mldno3)
  
}

fit_cjs_mldno3<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_mldno3, hessian=TRUE)
conv_cjs_mldno3<- fit_cjs_mldno3$convergence

dlmLogits_cjs_mldno3<- buildFun_cjs_mldno3(fit_cjs_mldno3$par)

V(dlmLogits_cjs_mldno3)

W(dlmLogits_cjs_mldno3)

mod_cjs_mldno3<- buildFun_cjs_mldno3(fit_cjs_mldno3$par)

cjs_mldno3filter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mldno3)

#What to do about first obs??
MAPE_cjs_mldno3<- mean(abs((cjs_mldno3filter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_mldno3<- rmse(logit.s_cjs, cjs_mldno3filter$f)

loglik_cjs_mldno3 <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE))

cjs.mldno3.aic <- (2 * (loglik_cjs_mldno3)) + 2 * (sum(n.coef))  

###########################################
buildFun_cjs_mnCUI<- function(parm){
  mod_cjs_mnCUI<- dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE)
  V(mod_cjs_mnCUI)<- exp(parm[1])
  diag(W(mod_cjs_mnCUI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_mnCUI)
  
}

fit_cjs_mnCUI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_mnCUI, hessian=TRUE)
conv_cjs_mnCUI<- fit_cjs_mnCUI$convergence

dlmLogits_cjs_mnCUI<- buildFun_cjs_mnCUI(fit_cjs_mnCUI$par)

V(dlmLogits_cjs_mnCUI)

W(dlmLogits_cjs_mnCUI)

mod_cjs_mnCUI<- buildFun_cjs_mnCUI(fit_cjs_mnCUI$par)

cjs_mnCUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnCUI)

#What to do about first obs??
MAPE_cjs_mnCUI<- mean(abs((cjs_mnCUIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_mnCUI<- rmse(logit.s_cjs, cjs_mnCUIfilter$f)

loglik_cjs_mnCUI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE))

cjs.mnCUI.aic <- (2 * (loglik_cjs_mnCUI)) + 2 * (sum(n.coef))  

###########################################
buildFun_cjs_mnCUTI<- function(parm){
  mod_cjs_mnCUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_mean_z, addInt=FALSE)
  V(mod_cjs_mnCUTI)<- exp(parm[1])
  diag(W(mod_cjs_mnCUTI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_mnCUTI)
  
}

fit_cjs_mnCUTI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_mnCUTI, hessian=TRUE)
conv_cjs_mnCUTI<- fit_cjs_mnCUTI$convergence

dlmLogits_cjs_mnCUTI<- buildFun_cjs_mnCUTI(fit_cjs_mnCUTI$par)

V(dlmLogits_cjs_mnCUTI)

W(dlmLogits_cjs_mnCUTI)

mod_cjs_mnCUTI<- buildFun_cjs_mnCUTI(fit_cjs_mnCUTI$par)

cjs_mnCUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnCUTI)

#What to do about first obs??
MAPE_cjs_mnCUTI<- mean(abs((cjs_mnCUTIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_mnCUTI<- rmse(logit.s_cjs, cjs_mnCUTIfilter$f)

loglik_cjs_mnCUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUTI_mean_z, addInt=FALSE))

cjs.mnCUTI.aic <- (2 * (loglik_cjs_mnCUTI)) + 2 * (sum(n.coef)) 

###########################################
buildFun_cjs_mnBEUTI<- function(parm){
  mod_cjs_mnBEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE)
  V(mod_cjs_mnBEUTI)<- exp(parm[1])
  diag(W(mod_cjs_mnBEUTI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_mnBEUTI)
  
}

fit_cjs_mnBEUTI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_mnBEUTI, hessian=TRUE)
conv_cjs_mnBEUTI<- fit_cjs_mnBEUTI$convergence

dlmLogits_cjs_mnBEUTI<- buildFun_cjs_mnBEUTI(fit_cjs_mnBEUTI$par)

V(dlmLogits_cjs_mnBEUTI)

W(dlmLogits_cjs_mnBEUTI)

mod_cjs_mnBEUTI<- buildFun_cjs_mnBEUTI(fit_cjs_mnBEUTI$par)

cjs_mnBEUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnBEUTI)

#What to do about first obs??
MAPE_cjs_mnBEUTI<- mean(abs((cjs_mnBEUTIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_mnBEUTI<- rmse(logit.s_cjs, cjs_mnBEUTIfilter$f)

loglik_cjs_mnBEUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE))

cjs.mnBEUTI.aic <- (2 * (loglik_cjs_mnBEUTI)) + 2 * (sum(n.coef)) 

###########################################
buildFun_cjs_mnmldno3<- function(parm){
  mod_cjs_mnmldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE)
  V(mod_cjs_mnmldno3)<- exp(parm[1])
  diag(W(mod_cjs_mnmldno3))[1:2]<- exp(parm[2:3])
  return(mod_cjs_mnmldno3)
  
}

fit_cjs_mnmldno3<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_mnmldno3, hessian=TRUE)
conv_cjs_mnmldno3<- fit_cjs_mnmldno3$convergence

dlmLogits_cjs_mnmldno3<- buildFun_cjs_mnmldno3(fit_cjs_mnmldno3$par)

V(dlmLogits_cjs_mnmldno3)

W(dlmLogits_cjs_mnmldno3)

mod_cjs_mnmldno3<- buildFun_cjs_mnmldno3(fit_cjs_mnmldno3$par)

cjs_mnmldno3filter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnmldno3)

#What to do about first obs??
MAPE_cjs_mnmldno3<- mean(abs((cjs_mnmldno3filter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_mnmldno3<- rmse(logit.s_cjs, cjs_mnmldno3filter$f)

loglik_cjs_mnmldno3 <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE))

cjs.mnmldno3.aic <- (2 * (loglik_cjs_mnmldno3)) + 2 * (sum(n.coef)) 

###########################################
buildFun_cjs_outflow<- function(parm){
  mod_cjs_outflow<- dlmModPoly(1) + dlmModReg(X=outflow_z, addInt=FALSE)
  V(mod_cjs_outflow)<- exp(parm[1])
  diag(W(mod_cjs_outflow))[1:2]<- exp(parm[2:3])
  return(mod_cjs_outflow)
  
}

fit_cjs_outflow<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_outflow, hessian=TRUE)
conv_cjs_outflow<- fit_cjs_outflow$convergence

dlmLogits_cjs_outflow<- buildFun_cjs_outflow(fit_cjs_outflow$par)

V(dlmLogits_cjs_outflow)

W(dlmLogits_cjs_outflow)

mod_cjs_outflow<- buildFun_cjs_outflow(fit_cjs_outflow$par)

cjs_outflowfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_outflow)

#What to do about first obs??
MAPE_cjs_outflow<- mean(abs((cjs_outflowfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_outflow<- rmse(logit.s_cjs, cjs_outflowfilter$f)

loglik_cjs_outflow <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=outflow_z, addInt=FALSE))

cjs.outlfow.aic <- (2 * (loglik_cjs_outflow)) + 2 * (sum(n.coef))  

###########################################
buildFun_cjs_temp<- function(parm){
  mod_cjs_temp<- dlmModPoly(1) + dlmModReg(X=temp_z, addInt=FALSE)
  V(mod_cjs_temp)<- exp(parm[1])
  diag(W(mod_cjs_temp))[1:2]<- exp(parm[2:3])
  return(mod_cjs_temp)
  
}

fit_cjs_temp<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_temp, hessian=TRUE)
conv_cjs_temp<- fit_cjs_temp$convergence

dlmLogits_cjs_temp<- buildFun_cjs_temp(fit_cjs_temp$par)

V(dlmLogits_cjs_temp)

W(dlmLogits_cjs_temp)

mod_cjs_temp<- buildFun_cjs_temp(fit_cjs_temp$par)

cjs_tempfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_temp)

#What to do about first obs??
MAPE_cjs_temp<- mean(abs((cjs_tempfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_temp<- rmse(logit.s_cjs, cjs_tempfilter$f)

loglik_cjs_temp <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=temp_z, addInt=FALSE))

cjs.temp.aic <- (2 * (loglik_cjs_temp)) + 2 * (sum(n.coef))  

#################################

buildFun_cjs_PNI<- function(parm){
  mod_cjs_PNI<- dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE)
  V(mod_cjs_PNI)<- exp(parm[1])
  diag(W(mod_cjs_PNI))[1:2]<- exp(parm[2:3])
  return(mod_cjs_PNI)
  
}

fit_cjs_PNI<- dlmMLE(logit.s_cjs, parm = rep(0,3), build=buildFun_cjs_PNI, hessian=TRUE)
conv_cjs_PNI<- fit_cjs_PNI$convergence

dlmLogits_cjs_PNI<- buildFun_cjs_PNI(fit_cjs_PNI$par)

V(dlmLogits_cjs_PNI)

W(dlmLogits_cjs_PNI)

mod_cjs_PNI<- buildFun_cjs_PNI(fit_cjs_PNI$par)

cjs_PNIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_PNI)

#What to do about first obs??
MAPE_cjs_PNI<- mean(abs((cjs_PNIfilter$f-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_cjs_PNI<- rmse(logit.s_cjs, cjs_PNIfilter$f)

loglik_cjs_PNI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE))

cjs.PNI.aic <- (2 * (loglik_cjs_PNI)) + 2 * (sum(n.coef))  

###################################################################################
#linear trend
buildFun<- function(parm){
  mod_cjs_CUI<- dlmModPoly(2) + dlmModReg(X=CUI_z, addInt=FALSE)
  V(mod_cjs_CUI)<- exp(parm[1])
  diag(W(mod_cjs_CUI))[1:3]<- exp(parm[2:4])
  return(mod_cjs_CUI)
  
}

fit_cjs_CUI<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun, hessian=TRUE)
conv_cjs_CUI<- fit_cjs_CUI$convergence

dlmLogits_cjs_CUI<- buildFun(fit_cjs_CUI$par)

V(dlmLogits_cjs_CUI)

W(dlmLogits_cjs_CUI)

mod_cjs_CUI<- buildFun(fit_cjs_CUI$par)

CUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUI)

loglik_cjs_CUI <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=CUI_z, addInt=FALSE))

n.coef <- 2
r.aic <- (2 * (loglik_cjs_CUI)) + 2 * (sum(n.coef))  #dlmLL caculates the neg. LL
r.bic <- (2 * (loglik_cjs_CUI)) + (log(length(logit.s))) * (n.coef)

CUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUI)

cui_mat<- as.matrix(cbind(rep(1, 42), CUI))
cjs_CUI_state<- CUIfilter$m

###########################################
buildFun_cjs_CUTI_lin<- function(parm){
  mod_cjs_CUTI<- dlmModPoly(2) + dlmModReg(X=CUTI_z, addInt=FALSE)
  V(mod_cjs_CUTI)<- exp(parm[1])
  diag(W(mod_cjs_CUTI))[1:3]<- exp(parm[2:4])
  return(mod_cjs_CUTI)
  
}

fit_cjs_CUTI<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun_cjs_CUTI_lin, hessian=TRUE)
conv_cjs_CUTI<- fit_cjs_CUTI$convergence

dlmLogits_cjs_CUTI<- buildFun(fit_cjs_CUTI$par)

V(dlmLogits_cjs_CUTI)

W(dlmLogits_cjs_CUTI)

mod_cjs_CUTI<- buildFun(fit_cjs_CUTI$par)

cjs_CUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUTI)

loglik_cjs_CUTI <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=CUTI_z, addInt=FALSE))

###########################################
buildFun<- function(parm){
  mod_cjs_BEUTI<- dlmModPoly(2) + dlmModReg(X=BEUTI_z, addInt=FALSE)
  V(mod_cjs_BEUTI)<- exp(parm[1])
  diag(W(mod_cjs_BEUTI))[1:3]<- exp(parm[2:4])
  return(mod_cjs_BEUTI)
  
}

fit_cjs_BEUTI<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun, hessian=TRUE)
conv_cjs_BEUTI<- fit_cjs_BEUTI$convergence

dlmLogits_cjs_BEUTI<- buildFun(fit_cjs_BEUTI$par)

V(dlmLogits_cjs_BEUTI)

W(dlmLogits_cjs_BEUTI)

mod_cjs_BEUTI<- buildFun(fit_cjs_BEUTI$par)

cjs_BEUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_BEUTI)

loglik_cjs_BEUTI <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=BEUTI_z, addInt=FALSE))

###########################################
buildFun_cjs_mldno3<- function(parm){
  mod_cjs_mldno3<- dlmModPoly(2) + dlmModReg(X=mldno3_z, addInt=FALSE)
  V(mod_cjs_mldno3)<- exp(parm[1])
  diag(W(mod_cjs_mldno3))[1:3]<- exp(parm[2:4])
  return(mod_cjs_mldno3)
  
}

fit_cjs_mldno3<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun_cjs_mldno3, hessian=TRUE)
conv_cjs_mldno3<- fit_cjs_mldno3$convergence

dlmLogits_cjs_mldno3<- buildFun(fit_cjs_mldno3$par)

V(dlmLogits_cjs_mldno3)

W(dlmLogits_cjs_mldno3)

mod_cjs_mldno3<- buildFun(fit_cjs_mldno3$par)

cjs_mldno3filter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mldno3)

loglik_cjs_mldno3 <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=mldno3_z, addInt=FALSE))

###########################################
buildFun_cjs_outflow<- function(parm){
  mod_cjs_outflow<- dlmModPoly(2) + dlmModReg(X=outflow_z, addInt=FALSE)
  V(mod_cjs_outflow)<- exp(parm[1])
  diag(W(mod_cjs_outflow))[1:3]<- exp(parm[2:4])
  return(mod_cjs_outflow)
  
}

fit_cjs_outflow<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun_cjs_outflow, hessian=TRUE)
conv_cjs_outflow<- fit_cjs_outflow$convergence

dlmLogits_cjs_outflow<- buildFun(fit_cjs_outflow$par)

V(dlmLogits_cjs_outflow)

W(dlmLogits_cjs_outflow)

mod_cjs_outflow<- buildFun(fit_cjs_outflow$par)

cjs_outflowfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_outflow)

loglik_cjs_outflow <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=outflow_z, addInt=FALSE))

###########################################
buildFun_cjs_temp<- function(parm){
  mod_cjs_temp<- dlmModPoly(2) + dlmModReg(X=temp_z, addInt=FALSE)
  V(mod_cjs_temp)<- exp(parm[1])
  diag(W(mod_cjs_temp))[1:3]<- exp(parm[2:4])
  return(mod_cjs_temp)
  
}

fit_cjs_temp<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun_cjs_temp, hessian=TRUE)
conv_cjs_temp<- fit_cjs_temp$convergence

dlmLogits_cjs_temp<- buildFun(fit_cjs_temp$par)

V(dlmLogits_cjs_temp)

W(dlmLogits_cjs_temp)

mod_cjs_temp<- buildFun(fit_cjs_temp$par)

cjs_tempfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_temp)

loglik_cjs_temp <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=temp_z, addInt=FALSE))

#################################

buildFun_cjs_PNI<- function(parm){
  mod_cjs_PNI<- dlmModPoly(2) + dlmModReg(X=PNI_z, addInt=FALSE)
  V(mod_cjs_PNI)<- exp(parm[1])
  diag(W(mod_cjs_PNI))[1:3]<- exp(parm[2:4])
  return(mod_cjs_PNI)
  
}

fit_cjs_PNI<- dlmMLE(logit.s_cjs, parm = rep(0,4), build=buildFun_cjs_PNI, hessian=TRUE)
conv_cjs_PNI<- fit_cjs_PNI$convergence

dlmLogits_cjs_PNI<- buildFun(fit_cjs_PNI$par)

V(dlmLogits_cjs_PNI)

W(dlmLogits_cjs_PNI)

mod_cjs_PNI<- buildFun(fit_cjs_PNI$par)

cjs_PNIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_PNI)

loglik_cjs_PNI <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=PNI_z, addInt=FALSE))