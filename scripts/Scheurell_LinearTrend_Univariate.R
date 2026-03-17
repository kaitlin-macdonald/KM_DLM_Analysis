rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)
library(tidyverse)

cjs_raw<- readRDS("data/df_SAR_LGRLGA_MV.rds") 

cjs<- cjs_raw

SAR<- as.matrix(cjs_raw %>% select(SW_logit_phi, SAR_logit_phi, cjs_logit_phi))
SAR

################################
##Salmon Survive test model
year_cjs<- cjs$year
logit.s_cjs<- SAR
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
SST_pre_z<- scale(cjs$SST_pre1212)
SST_entry_z<- scale(cjs$SST_entry567)
SST_6mo_z<- scale(cjs$SST_6mo)

##########################
#Run Reference model
buildFun_cjs_ref<- function(parm){
  mod_ref<- dlmModPoly(1) 
  mod_ref$FF <- mod_ref$FF %x% diag(3)
  mod_ref$GG <- mod_ref$GG %x% diag(3)
  mod_ref$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_ref$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6])))
  return(mod_ref)
  
}

fit_cjs_ref<- dlmMLE(logit.s_cjs, parm = rep(0,6), build=buildFun_cjs_ref, hessian=TRUE)
conv_cjs_ref<- fit_cjs_ref$convergence

dlmLogits_cjs_ref<- buildFun_cjs_ref(fit_cjs_ref$par)

V(dlmLogits_cjs_ref)

W(dlmLogits_cjs_ref)

mod_cjs_ref<-buildFun_cjs_ref(fit_cjs_ref$par)

#Apply Kalman filter to the model
cjs_reffilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_ref)

#Mean absolute percentage error
MAPE_cjs_ref<- mean(na.omit(abs((cjs_reffilter$f-logit.s_cjs)/logit.s_cjs)))
rmse_cjs_ref<- sqrt(mean(na.omit((cjs_reffilter$f - logit.s_cjs)^2)))

#dlmLL caculates the neg. LL (lower is better)
loglik_cjs_ref <- dlmLL(logit.s_cjs, dlmModPoly(1))

n.coef <- 6
cjs.ref.aic <- (2 * (loglik_cjs_ref)) + 2 * (sum(n.coef))
cjs.ref.aicc <- (2 * (loglik_cjs_ref)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)

mod_sel_df<- data.frame(Covariate=NA, Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_ref, RMSE=rmse_cjs_ref, Href=0, AICc=cjs.ref.aicc)

####Goodness of fit

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
  mod_CUI<- dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE)
 # V(mod_cjs_CUI)<- exp(parm[1])
  #diag(W(mod_cjs_CUI))[1:3]<- exp(parm[2:4])
  #return(mod_cjs_CUI)
  mod_CUI$FF <- mod_CUI$FF %x% diag(3)
  mod_CUI$GG <- mod_CUI$GG %x% diag(3)
  mod_CUI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_CUI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_CUI$m0 <- rep(0, 9)
  mod_CUI$C0 <- diag(9) * 0
  return(mod_CUI)
  
}

fit_cjs_CUI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_CUI, hessian=TRUE)
conv_cjs_CUI<- fit_cjs_CUI$convergence

dlmLogits_cjs_CUI<- buildFun_cjs_CUI(fit_cjs_CUI$par)

V(dlmLogits_cjs_CUI)

W(dlmLogits_cjs_CUI)

mod_cjs_CUI<- buildFun_cjs_CUI(fit_cjs_CUI$par)

#Apply Kalman filter to the model
cjs_CUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUI)

#Mean absolute percentage error
MAPE_cjs_CUI<- mean(na.omit(abs((cjs_CUIfilter$f-logit.s_cjs)/logit.s_cjs)))
rmse_cjs_CUI<- sqrt(mean(na.omit((cjs_CUIfilter$f - logit.s_cjs)^2)))


#dlmLL caculates the neg. LL (lower is better)
loglik_cjs_CUI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUI_z, addInt=FALSE))

Href_cjs_CUI<- 2*(loglik_cjs_ref - loglik_cjs_CUI)

n.coef <- 9
cjs.CUI.aic <- (2 * (loglik_cjs_CUI)) + 2 * (sum(n.coef))  
cjs.CUI.AICc <- (2 * (loglik_cjs_CUI)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1) 

#Assess Goodness of fit
resids_CUI <- residuals(cjs_CUIfilter, sd = FALSE)
plot.ts(resids_CUI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_CUI)
qqline(resids_CUI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_CUI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_CUI, lag.max = 10)

tsdiag(cjs_CUIfilter)

hist(resids_CUI)

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="Constant", Data="Mutlivariate", MAPE=MAPE_cjs_CUI, RMSE=rmse_cjs_CUI, 
                     Href=Href_cjs_CUI, AICc=cjs.CUI.AICc)
###########################################
buildFun_cjs_CUTI<- function(parm){
  mod_CUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE)
  mod_CUTI$FF <- mod_CUTI$FF %x% diag(3)
  mod_CUTI$GG <- mod_CUTI$GG %x% diag(3)
  mod_CUTI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_CUTI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_CUTI$m0 <- rep(0, 9)
  mod_CUTI$C0 <- diag(9) * 0
  return(mod_CUTI)
  
}

fit_cjs_CUTI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_CUTI, hessian=TRUE)
conv_cjs_CUTI<- fit_cjs_CUTI$convergence

dlmLogits_cjs_CUTI<- buildFun_cjs_CUTI(fit_cjs_CUTI$par)

V(dlmLogits_cjs_CUTI)

W(dlmLogits_cjs_CUTI)

mod_cjs_CUTI<- buildFun_cjs_CUTI(fit_cjs_CUTI$par)

cjs_CUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_CUTI)

#What to do about first obs??
MAPE_cjs_CUTI<- mean(na.omit(abs((cjs_CUTIfilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_CUTI<- sqrt(mean(na.omit((cjs_CUTIfilter$f - logit.s_cjs)^2)))


loglik_cjs_CUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUTI_z, addInt=FALSE))

Href_cjs_CUTI<- 2*(loglik_cjs_ref - loglik_cjs_CUTI)

n.coef<- 9
cjs.CUTI.AICc <- (2 * (loglik_cjs_CUTI)) + 2 * (sum(n.coef))  + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)   

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="CUTI", Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_CUTI, RMSE=rmse_cjs_CUTI, 
                     Href=Href_cjs_CUTI, AICc=cjs.CUTI.AICc)

#Assess Goodness of fit
resids_CUTI <- residuals(cjs_CUTIfilter, sd = FALSE)
plot.ts(resids_CUTI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_CUTI)
qqline(resids_CUTI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_CUTI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_CUTI, lag.max = 10)

tsdiag(cjs_CUTIfilter)

hist(resids_CUTI)

###########################################
buildFun_cjs_BEUTI<- function(parm){
  mod_BEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE)
  mod_BEUTI$FF <- mod_BEUTI$FF %x% diag(3)
  mod_BEUTI$GG <- mod_BEUTI$GG %x% diag(3)
  mod_BEUTI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_BEUTI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_BEUTI$m0 <- rep(0, 9)
  mod_BEUTI$C0 <- diag(9) * 0
  return(mod_BEUTI)
}

fit_cjs_BEUTI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_BEUTI, hessian=TRUE)
conv_cjs_BEUTI<- fit_cjs_BEUTI$convergence

dlmLogits_cjs_BEUTI<- buildFun_cjs_BEUTI(fit_cjs_BEUTI$par)

V(dlmLogits_cjs_BEUTI)

W(dlmLogits_cjs_BEUTI)

mod_cjs_BEUTI<- buildFun_cjs_BEUTI(fit_cjs_BEUTI$par)

cjs_BEUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_BEUTI)

#What to do about first obs??
MAPE_cjs_BEUTI<- mean(na.omit(abs((cjs_BEUTIfilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_BEUTI<- sqrt(mean(na.omit((cjs_BEUTIfilter$f - logit.s_cjs)^2)))

loglik_cjs_BEUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=BEUTI_z, addInt=FALSE))

Href_cjs_BEUTI<- 2*(loglik_cjs_ref - loglik_cjs_BEUTI)

cjs.BEUTI.AICc <- (2 * (loglik_cjs_BEUTI)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="BEUTI", Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_BEUTI, RMSE=rmse_cjs_BEUTI, 
                     Href=Href_cjs_BEUTI, AICc=cjs.BEUTI.AICc)

#Assess Goodness of fit
resids_BEUTI <- residuals(cjs_BEUTIfilter, sd = FALSE)
plot.ts(resids_BEUTI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_BEUTI)
qqline(resids_BEUTI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_BEUTI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_BEUTI, lag.max = 10)

tsdiag(cjs_BEUTIfilter)

hist(resids_BEUTI)

###########################################
buildFun_cjs_mldno3<- function(parm){
  mod_mldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE)
  mod_mldno3$FF <- mod_mldno3$FF %x% diag(3)
  mod_mldno3$GG <- mod_mldno3$GG %x% diag(3)
  mod_mldno3$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_mldno3$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_mldno3$m0 <- rep(0, 9)
  mod_mldno3$C0 <- diag(9) * 0
  return(mod_mldno3)
  
}

fit_cjs_mldno3<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_mldno3, hessian=TRUE)
conv_cjs_mldno3<- fit_cjs_mldno3$convergence

dlmLogits_cjs_mldno3<- buildFun_cjs_mldno3(fit_cjs_mldno3$par)

V(dlmLogits_cjs_mldno3)

W(dlmLogits_cjs_mldno3)

mod_cjs_mldno3<- buildFun_cjs_mldno3(fit_cjs_mldno3$par)

cjs_mldno3filter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mldno3)

#What to do about first obs??
MAPE_cjs_mldno3<- mean(na.omit(abs((cjs_mldno3filter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_mldno3<- sqrt(mean(na.omit((cjs_mldno3filter$f - logit.s_cjs)^2)))

loglik_cjs_mldno3 <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=mldno3_z, addInt=FALSE))

Href_cjs_mldno3<- 2*(loglik_cjs_ref - loglik_cjs_mldno3)

cjs.mldno3.AICc <- (2 * (loglik_cjs_mldno3)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1) 

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mldno3", Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_mldno3, RMSE=rmse_cjs_mldno3, 
                     Href=Href_cjs_mldno3, AICc=cjs.mldno3.AICc)

#Assess Goodness of fit
resids_mldno3 <- residuals(cjs_mldno3filter, sd = FALSE)
plot.ts(resids_mldno3, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_mldno3)
qqline(resids_mldno3)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_mldno3, mu = 0)$p.value

## plot ACF of innovations
acf(resids_mldno3, lag.max = 10)

tsdiag(cjs_mldno3filter)

hist(resids_mldno3)

###########################################
buildFun_cjs_mnCUI<- function(parm){
  mod_mnCUI<- dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE)
  mod_mnCUI$FF <- mod_mnCUI$FF %x% diag(3)
  mod_mnCUI$GG <- mod_mnCUI$GG %x% diag(3)
  mod_mnCUI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_mnCUI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_mnCUI$m0 <- rep(0, 9)
  mod_mnCUI$C0 <- diag(9) * 0
  return(mod_mnCUI)
  
}

fit_cjs_mnCUI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_mnCUI, hessian=TRUE)
conv_cjs_mnCUI<- fit_cjs_mnCUI$convergence

dlmLogits_cjs_mnCUI<- buildFun_cjs_mnCUI(fit_cjs_mnCUI$par)

V(dlmLogits_cjs_mnCUI)

W(dlmLogits_cjs_mnCUI)

mod_cjs_mnCUI<- buildFun_cjs_mnCUI(fit_cjs_mnCUI$par)

cjs_mnCUIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnCUI)

#What to do about first obs??
MAPE_cjs_mnCUI<- mean(na.omit(abs((cjs_mnCUIfilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_mnCUI<- sqrt(mean(na.omit((cjs_mnCUIfilter$f - logit.s_cjs)^2)))

loglik_cjs_mnCUI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=CUI_mean_z, addInt=FALSE))

Href_cjs_mnCUI<- 2*(loglik_cjs_ref - loglik_cjs_mnCUI)

cjs.mnCUI.AICc <- (2 * (loglik_cjs_mnCUI)) + 2 * (sum(n.coef))+ 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_CUI", Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_mnCUI, RMSE=rmse_cjs_mnCUI, 
                     Href=Href_cjs_mnCUI, AICc=cjs.mnCUI.AICc)
#Assess Goodness of fit
resids_mnCUI <- residuals(cjs_mnCUIfilter, sd = FALSE)
plot.ts(resids_mnCUI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_mnCUI)
qqline(resids_mnCUI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_mnCUI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_mnCUI, lag.max = 10)

tsdiag(cjs_mnCUIfilter)

hist(resids_mnCUI)

###########################################
buildFun_cjs_mnCUTI<- function(parm){
  mod_mnCUTI<- dlmModPoly(1) + dlmModReg(X=CUTI_mean_z, addInt=FALSE)
  mod_mnCUTI$FF <- mod_mnCUTI$FF %x% diag(3)
  mod_mnCUTI$GG <- mod_mnCUTI$GG %x% diag(3)
  mod_mnCUTI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_mnCUTI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_mnCUTI$m0 <- rep(0, 9)
  mod_mnCUTI$C0 <- diag(9) * 0
  return(mod_mnCUTI)
  
}

fit_cjs_mnCUTI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_mnCUTI, hessian=TRUE)
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

loglik_cjs_mnCUTI <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=CUTI_mean_z, addInt=FALSE))

Href_cjs_mnCUTI<- 2*(loglik_cjs_ref - loglik_cjs_mnCUTI)

cjs.mnCUTI.AICc <- (2 * (loglik_cjs_mnCUTI)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+ 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_CUTI", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_mnCUTI, RMSE=rmse_cjs_mnCUTI, 
                     Href=Href_cjs_mnCUTI, AICc=cjs.mnCUTI.AICc)

#Assess Goodness of fit
resids_mnCUTI <- residuals(cjs_mnCUTIfilter, sd = FALSE)
plot.ts(resids_mnCUTI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_mnCUTI)
qqline(resids_mnCUTI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_mnCUTI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_mnCUTI, lag.max = 10)

tsdiag(cjs_mnCUTIfilter)

hist(resids_mnCUTI)


###########################################
buildFun_cjs_mnBEUTI<- function(parm){
  mod_mnBEUTI<- dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE)
  mod_mnBEUTI$FF <- mod_mnBEUTI$FF %x% diag(3)
  mod_mnBEUTI$GG <- mod_mnBEUTI$GG %x% diag(3)
  mod_mnBEUTI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_mnBEUTI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_mnBEUTI$m0 <- rep(0, 9)
  mod_mnBEUTI$C0 <- diag(9) * 0
  return(mod_mnBEUTI)
  
}

fit_cjs_mnBEUTI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_mnBEUTI, hessian=TRUE)
conv_cjs_mnBEUTI<- fit_cjs_mnBEUTI$convergence

dlmLogits_cjs_mnBEUTI<- buildFun_cjs_mnBEUTI(fit_cjs_mnBEUTI$par)

V(dlmLogits_cjs_mnBEUTI)

W(dlmLogits_cjs_mnBEUTI)

mod_cjs_mnBEUTI<- buildFun_cjs_mnBEUTI(fit_cjs_mnBEUTI$par)

cjs_mnBEUTIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnBEUTI)

#What to do about first obs??
MAPE_cjs_mnBEUTI<- mean(na.omit(abs((cjs_mnBEUTIfilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_mnBEUTI<- sqrt(mean(na.omit((cjs_mnBEUTIfilter$f - logit.s_cjs)^2)))

loglik_cjs_mnBEUTI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=BEUTI_mean_z, addInt=FALSE))

Href_cjs_mnBEUTI<- 2*(loglik_cjs_ref - loglik_cjs_mnBEUTI)

cjs.mnBEUTI.AICc <- (2 * (loglik_cjs_mnBEUTI)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+ 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1) 

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_BEUTI", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_mnBEUTI, RMSE=rmse_cjs_mnBEUTI, 
                     Href=Href_cjs_mnBEUTI, AICc=cjs.mnBEUTI.AICc)

#Assess Goodness of fit
resids_mnBEUTI <- residuals(cjs_mnBEUTIfilter, sd = FALSE)
plot.ts(resids_mnBEUTI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_mnBEUTI)
qqline(resids_mnBEUTI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_mnBEUTI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_mnBEUTI, lag.max = 10)

tsdiag(cjs_mnBEUTIfilter)

hist(resids_mnBEUTI)


###########################################
buildFun_cjs_mnmldno3<- function(parm){
  mod_mnmldno3<- dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE)
  mod_mnmldno3$FF <- mod_mnmldno3$FF %x% diag(3)
  mod_mnmldno3$GG <- mod_mnmldno3$GG %x% diag(3)
  mod_mnmldno3$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_mnmldno3$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_mnmldno3$m0 <- rep(0, 9)
  mod_mnmldno3$C0 <- diag(9) * 0
  return(mod_mnmldno3)
  
}

fit_cjs_mnmldno3<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_mnmldno3, hessian=TRUE)
conv_cjs_mnmldno3<- fit_cjs_mnmldno3$convergence

dlmLogits_cjs_mnmldno3<- buildFun_cjs_mnmldno3(fit_cjs_mnmldno3$par)

V(dlmLogits_cjs_mnmldno3)

W(dlmLogits_cjs_mnmldno3)

mod_cjs_mnmldno3<- buildFun_cjs_mnmldno3(fit_cjs_mnmldno3$par)

cjs_mnmldno3filter<- dlmFilter(logit.s_cjs, mod=mod_cjs_mnmldno3)

#What to do about first obs??
MAPE_cjs_mnmldno3<- mean(na.omit(abs((cjs_mnmldno3filter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_mnmldno3<- sqrt(mean(na.omit((cjs_mnBEUTIfilter$f - logit.s_cjs)^2)))

loglik_cjs_mnmldno3 <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=mldno3_mean_z, addInt=FALSE))

Href_cjs_mnmldno3<- 2*(loglik_cjs_ref - loglik_cjs_mnmldno3)

cjs.mnmldno3.AICc <- (2 * (loglik_cjs_mnmldno3)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1) 

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_Mldno3", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_mnmldno3, RMSE=rmse_cjs_mnmldno3, 
                     Href=Href_cjs_mnmldno3, AICc=cjs.mnmldno3.AICc)

#Assess Goodness of fit
resids_mnmldno3 <- residuals(cjs_mnmldno3filter, sd = FALSE)
plot.ts(resids_mnmldno3, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_mnmldno3)
qqline(resids_mnmldno3)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_mnmldno3, mu = 0)$p.value

## plot ACF of innovations
acf(resids_mnmldno3, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_mnmldno3filter)

hist(resids_mnmldno3)

###########################################
buildFun_cjs_outflow<- function(parm){
  mod_outflow<- dlmModPoly(1) + dlmModReg(X=outflow_z, addInt=FALSE)
  mod_outflow$FF <- mod_outflow$FF %x% diag(3)
  mod_outflow$GG <- mod_outflow$GG %x% diag(3)
  mod_outflow$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_outflow$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_outflow$m0 <- rep(0, 9)
  mod_outflow$C0 <- diag(9) * 0
  return(mod_outflow)
  
}

fit_cjs_outflow<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_outflow, hessian=TRUE)
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

loglik_cjs_outflow <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=outflow_z, addInt=FALSE))

Href_cjs_outflow<- 2*(loglik_cjs_ref - loglik_cjs_outflow)

cjs.outflow.AICc <- (2 * (loglik_cjs_outflow)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_outflow", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_outflow, RMSE=rmse_cjs_outflow, 
                     Href=Href_cjs_outflow, AICc=cjs.outflow.AICc)

#Assess Goodness of fit
resids_outflow <- residuals(cjs_outflowfilter, sd = FALSE)
plot.ts(resids_outflow, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_outflow)
qqline(resids_outflow)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_outflow, mu = 0)$p.value

## plot ACF of innovations
acf(resids_outflow, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_outflowfilter)

hist(resids_outflow)

###########################################
buildFun_cjs_temp<- function(parm){
  mod_temp<- dlmModPoly(1) + dlmModReg(X=temp_z, addInt=FALSE)
  mod_temp$FF <- mod_temp$FF %x% diag(3)
  mod_temp$GG <- mod_temp$GG %x% diag(3)
  mod_temp$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_temp$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_temp$m0 <- rep(0, 9)
  mod_temp$C0 <- diag(9) * 0
  return(mod_temp)
  
}

fit_cjs_temp<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_temp, hessian=TRUE)
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

loglik_cjs_temp <- dlmLL(logit.s_cjs, dlmModPoly(2) + dlmModReg(X=temp_z, addInt=FALSE))

Href_cjs_temp<- 2*(loglik_cjs_ref - loglik_cjs_temp)

cjs.temp.AICc <- (2 * (loglik_cjs_temp)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_temperature", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_temp, RMSE=rmse_cjs_temp, 
                     Href=Href_cjs_temp, AICc=cjs.temp.AICc)

#Assess Goodness of fit
resids_temp <- residuals(cjs_tempfilter, sd = FALSE)
plot.ts(resids_temp, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_temp)
qqline(resids_temp)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_temp, mu = 0)$p.value

## plot ACF of innovations
acf(resids_temp, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_tempfilter)

hist(resids_temp)

#################################

buildFun_cjs_PNI<- function(parm){
  mod_PNI<- dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE)
  mod_PNI$FF <- mod_PNI$FF %x% diag(3)
  mod_PNI$GG <- mod_PNI$GG %x% diag(3)
  mod_PNI$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_PNI$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_PNI$m0 <- rep(0, 9)
  mod_PNI$C0 <- diag(9) * 0
  return(mod_PNI)
  
}

fit_cjs_PNI<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_PNI, hessian=TRUE)
conv_cjs_PNI<- fit_cjs_PNI$convergence

dlmLogits_cjs_PNI<- buildFun_cjs_PNI(fit_cjs_PNI$par)

V(dlmLogits_cjs_PNI)

W(dlmLogits_cjs_PNI)

mod_cjs_PNI<- buildFun_cjs_PNI(fit_cjs_PNI$par)

cjs_PNIfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_PNI)

#What to do about first obs??
MAPE_cjs_PNI<- mean(abs((cjs_PNIfilter$f-logit.s_cjs)/logit.s_cjs))

#Root mean square error
rmse_cjs_PNI<- sqrt(mean(na.omit((cjs_PNIfilter$f - logit.s_cjs)^2)))

loglik_cjs_PNI <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=PNI_z, addInt=FALSE))

Href_cjs_PNI<- 2*(loglik_cjs_ref - loglik_cjs_PNI)

cjs.PNI.AICc <- (2 * (loglik_cjs_PNI)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)    

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Mean_PNI", Model="Constant", Data="Multivariate", MAPE=MAPE_cjs_PNI, RMSE=rmse_cjs_PNI, 
                     Href=Href_cjs_PNI, AICc=cjs.PNI.AICc)

#Assess Goodness of fit
resids_PNI <- residuals(cjs_PNIfilter, sd = FALSE)
plot.ts(resids_PNI, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_PNI)
qqline(resids_PNI)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_PNI, mu = 0)$p.value

## plot ACF of innovations
acf(resids_PNI, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_PNIfilter)

hist(resids_PNI)

############################
###########################################
buildFun_cjs_SST_pre<- function(parm){
  mod_SST_pre<- dlmModPoly(1) + dlmModReg(X=SST_pre_z, addInt=FALSE)
  mod_SST_pre$FF <- mod_SST_pre$FF %x% diag(3)
  mod_SST_pre$GG <- mod_SST_pre$GG %x% diag(3)
  mod_SST_pre$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_SST_pre$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_SST_pre$m0 <- rep(0, 9)
  mod_SST_pre$C0 <- diag(9) * 0
  return(mod_SST_pre)
  
}

fit_cjs_SST_pre<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_SST_pre, hessian=TRUE)
conv_cjs_SST_pre<- fit_cjs_SST_pre$convergence

dlmLogits_cjs_SST_pre<- buildFun_cjs_SST_pre(fit_cjs_SST_pre$par)

V(dlmLogits_cjs_SST_pre)

W(dlmLogits_cjs_SST_pre)

mod_cjs_SST_pre<- buildFun_cjs_SST_pre(fit_cjs_SST_pre$par)

cjs_SST_prefilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_SST_pre)

#What to do about first obs??
MAPE_cjs_SST_pre<- mean(na.omit(abs((cjs_SST_prefilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_SST_pre<- sqrt(mean(na.omit((cjs_SST_prefilter$f - logit.s_cjs)^2)))

loglik_cjs_SST_pre <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=SST_pre_z, addInt=FALSE))

Href_cjs_SST_pre<- 2*(loglik_cjs_ref - loglik_cjs_SST_pre)

cjs.SST_pre.AICc <- (2 * (loglik_cjs_SST_pre)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Pre Entry SST", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_SST_pre, RMSE=rmse_cjs_SST_pre, 
                     Href=Href_cjs_SST_pre, AICc=cjs.SST_pre.AICc)

#Assess Goodness of fit
resids_SST_pre <- residuals(cjs_SST_prefilter, sd = FALSE)
plot.ts(resids_SST_pre, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_SST_pre)
qqline(resids_SST_pre)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_SST_pre, mu = 0)$p.value

## plot ACF of innovations
acf(resids_SST_pre, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_SST_prefilter)

hist(resids_SST_pre)

###########################################
buildFun_cjs_SST_entry<- function(parm){
  mod_SST_entry<- dlmModPoly(1) + dlmModReg(X=SST_entry_z, addInt=FALSE)
  mod_SST_entry$FF <- mod_SST_entry$FF %x% diag(3)
  mod_SST_entry$GG <- mod_SST_entry$GG %x% diag(3)
  mod_SST_entry$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_SST_entry$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_SST_entry$m0 <- rep(0, 9)
  mod_SST_entry$C0 <- diag(9) * 0
  return(mod_SST_entry)
  
}

fit_cjs_SST_entry<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_SST_entry, hessian=TRUE)
conv_cjs_SST_entry<- fit_cjs_SST_entry$convergence

dlmLogits_cjs_SST_entry<- buildFun_cjs_SST_entry(fit_cjs_SST_entry$par)

V(dlmLogits_cjs_SST_entry)

W(dlmLogits_cjs_SST_entry)

mod_cjs_SST_entry<- buildFun_cjs_SST_entry(fit_cjs_SST_entry$par)

cjs_SST_entryfilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_SST_entry)

#What to do about first obs??
MAPE_cjs_SST_entry<- mean(na.omit(abs((cjs_SST_entryfilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_SST_entry<- sqrt(mean(na.omit((cjs_SST_entryfilter$f - logit.s_cjs)^2)))

loglik_cjs_SST_entry <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=SST_entry_z, addInt=FALSE))

Href_cjs_SST_entry<- 2*(loglik_cjs_ref - loglik_cjs_SST_entry)

cjs.SST_entry.AICc <- (2 * (loglik_cjs_SST_entry)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="Entry SST", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_SST_entry, RMSE=rmse_cjs_SST_entry, 
                     Href=Href_cjs_SST_entry, AICc=cjs.SST_entry.AICc)

#Assess Goodness of fit
resids_SST_entry <- residuals(cjs_SST_entryfilter, sd = FALSE)
plot.ts(resids_SST_entry, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_SST_entry)
qqline(resids_SST_entry)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_SST_entry, mu = 0)$p.value

## plot ACF of innovations
acf(resids_SST_entry, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_SST_entryfilter)

hist(resids_SST_entry)

###########################################
buildFun_cjs_SST_6mo<- function(parm){
  mod_SST_6mo<- dlmModPoly(1) + dlmModReg(X=SST_6mo_z, addInt=FALSE)
  mod_SST_6mo$FF <- mod_SST_6mo$FF %x% diag(3)
  mod_SST_6mo$GG <- mod_SST_6mo$GG %x% diag(3)
  mod_SST_6mo$V<- diag(c(exp(parm[1]), exp(parm[2]), exp(parm[3])))
  mod_SST_6mo$W  <- diag(c(exp(parm[4]), exp(parm[5]), exp(parm[6]), exp(parm[7]), exp(parm[8]), exp(parm[9])))
  mod_SST_6mo$m0 <- rep(0, 9)
  mod_SST_6mo$C0 <- diag(9) * 0
  return(mod_SST_6mo)
  
}

fit_cjs_SST_6mo<- dlmMLE(logit.s_cjs, parm = rep(0,9), build=buildFun_cjs_SST_6mo, hessian=TRUE)
conv_cjs_SST_6mo<- fit_cjs_SST_6mo$convergence

dlmLogits_cjs_SST_6mo<- buildFun_cjs_SST_6mo(fit_cjs_SST_6mo$par)

V(dlmLogits_cjs_SST_6mo)

W(dlmLogits_cjs_SST_6mo)

mod_cjs_SST_6mo<- buildFun_cjs_SST_6mo(fit_cjs_SST_6mo$par)

cjs_SST_6mofilter<- dlmFilter(logit.s_cjs, mod=mod_cjs_SST_6mo)

#What to do about first obs??
MAPE_cjs_SST_6mo<- mean(na.omit(abs((cjs_SST_6mofilter$f-logit.s_cjs)/logit.s_cjs)))
#Root mean square error
rmse_cjs_SST_6mo<- sqrt(mean(na.omit((cjs_SST_6mofilter$f - logit.s_cjs)^2)))

loglik_cjs_SST_6mo <- dlmLL(logit.s_cjs, dlmModPoly(1) + dlmModReg(X=SST_6mo_z, addInt=FALSE))

Href_cjs_SST_6mo<- 2*(loglik_cjs_ref - loglik_cjs_SST_6mo)

cjs.SST_6mo.AICc <- (2 * (loglik_cjs_SST_6mo)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef)) + 2 * (sum(n.coef))+((2*n.coef*(n.coef+1))/length(logit.s_cjs)- n.coef-1)  

#Add to model selection table
mod_sel_df<- add_row(mod_sel_df, Covariate="6 Month SST", Model="LinTrend", Data="CJS", MAPE=MAPE_cjs_SST_6mo, RMSE=rmse_cjs_SST_6mo, 
                     Href=Href_cjs_SST_6mo, AICc=cjs.SST_6mo.AICc)

#Assess Goodness of fit
resids_SST_6mo <- residuals(cjs_SST_6mofilter, sd = FALSE)
plot.ts(resids_SST_6mo, ylab = "", xlab = "", col = "darkgrey", 
        lwd = 1.5)
abline(h = 0)
legend("topright", legend = "Residuals", lwd = 1.5, col = "darkgrey", 
       bty = "n")
qqnorm(resids_SST_6mo)
qqline(resids_SST_6mo)

## p-value for t-test of H0: E(innov) = 0
t.test(resids_SST_6mo, mu = 0)$p.value

## plot ACF of innovations
acf(resids_SST_6mo, lag.max = 10)
#could be issues with autocorrelation?
tsdiag(cjs_SST_6mofilter)

hist(resids_SST_6mo)

write.table(mod_sel_df, file='output/mod_sel_lin_uni_cjs.csv', col.names=TRUE, row.names = FALSE, sep=",")
