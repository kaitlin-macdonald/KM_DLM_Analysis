rm(list = ls())
#Load packages 
library(dlm)
library(MARSS)
library(dplyr)
library(here)
library(Metrics)
library(tidyverse)

#Load coastal upwelling and SAR data from Scheurell and Williams 2005
SAR_CUTI<- readRDS("C:/Users/krmacdon/Documents/CBR/sar_raw_data_CUI_CUTI.rds")
SAR<- SAR_CUTI[1:42,]
TT<- dim(SAR)[1]
#Put data in format to use in modeling
year_SAR<- SAR$year
logit.s_SAR<- SAR$logit.s
Apr_z<- scale(SAR$Apr.x)
Sep_z<- scale(SAR$Sep.x)
Oct_z<- scale(SAR$Oct.x)

##########################
#Run Reference model with dlm package
buildFun_SAR_ref<- function(parm){
  mod_SAR_ref<- dlmModPoly(order=1, dV=exp(parm[1]), dW=exp(parm[2])) 
  V(mod_SAR_ref)<- exp(parm[1])
  diag(W(mod_SAR_ref))[1]<- exp(parm[2])
  mod_SAR_ref$m0 <- rep(0, 2)
  mod_SAR_ref$C0 <- diag(2) * 0
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

#dlmLL caculates the neg. LL (lower is better)
loglik_SAR_ref <- dlmLL(logit.s_SAR, dlmModPoly(order=1) ) 

Compare_df<- data.frame(Covariate=NA, Model="Constant", Package="dlm", LogLik=-1*loglik_SAR_ref)

#######################
#Run same referenece model with MARSS

m=1

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_constant <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.out.constant <- MARSSkfss(dlm_constant)

#Compare log likelihood from MARSS 
dlm_constant$logLik

#Compare log likelihood from dlm (note dlm gives negative log likelihood)
-1*loglik_SAR_ref

Compare_df<- add_row(Compare_df, Covariate=NA, Model="Constant", Package="MARSS", LogLik=dlm_constant$logLik)
##########################
#Run April CUI model
buildFun_SAR_CUI<- function(parm){
  mod_SAR_CUI<- dlmModPoly(1) + dlmModReg(X=Apr_z, addInt=FALSE)
  V(mod_SAR_CUI)<- exp(parm[1])
  diag(W(mod_SAR_CUI))[1:2]<- exp(parm[2:3])
  return(mod_SAR_CUI)
  
}
#Fit the model and get convergence
fit_SAR_CUI<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_CUI, hessian=TRUE)
conv_SAR_CUI<- fit_SAR_CUI$convergence

dlmLogits_SAR_CUI<- buildFun_SAR_CUI(fit_SAR_CUI$par)
#check the estimates
V(dlmLogits_SAR_CUI)

W(dlmLogits_SAR_CUI)
#Build the model with estimates
mod_SAR_CUI<- buildFun_SAR_CUI(fit_SAR_CUI$par)

#Apply Kalman filter to the model
SAR_CUIfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_CUI)

#dlmLL caculates the neg. LL (lower is better)
loglik_SAR_CUI <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=Apr_z, addInt=FALSE))

Compare_df<- add_row(Compare_df, Covariate="April CUI", Model="Constant", Package="dlm", LogLik=-1*loglik_SAR_CUI)

###############################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- Apr_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_Apr <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.out <- MARSSkfss(dlm_Apr)

#Compare log likelihood from MARSS 
dlm_Apr$logLik

#Compare log likelihood from dlm (note dlm gives negative log likelihood)
-1*loglik_SAR_CUI

Compare_df<- add_row(Compare_df, Covariate="April CUI", Model="Constant", Package="MARSS", LogLik=dlm_Apr$logLik)

###########################################
#Run September CUI model with dlm
buildFun_SAR_Sep<- function(parm){
  mod_SAR_Sep<- dlmModPoly(1) + dlmModReg(X=Sep_z, addInt=FALSE)
  V(mod_SAR_Sep)<- exp(parm[1])
  diag(W(mod_SAR_Sep))[1:2]<- exp(parm[2:3])
  return(mod_SAR_Sep)
  
}

fit_SAR_Sep<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_Sep, hessian=TRUE)
conv_SAR_Sep<- fit_SAR_Sep$convergence

dlmLogits_SAR_Sep<- buildFun_SAR_Sep(fit_SAR_Sep$par)

V(dlmLogits_SAR_Sep)

W(dlmLogits_SAR_Sep)

mod_SAR_Sep<- buildFun_SAR_Sep(fit_SAR_Sep$par)

SAR_Sepfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_Sep)

loglik_SAR_Sep <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=Sep_z, addInt=FALSE))

Compare_df<- add_row(Compare_df, Covariate="Sept CUI", Model="Constant", Package="dlm", LogLik=-1*loglik_SAR_Sep)

###########################################
#Run September CUI model with MARSS
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- Sep_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_Sep <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.out <- MARSSkfss(dlm_Sep)

Compare_df<- add_row(Compare_df, Covariate="Sept CUI", Model="Constant", Package="MARSS", LogLik=dlm_Sep$logLik)

#########################################
#Fit dlm model for October CUI
buildFun_SAR_Oct<- function(parm){
  mod_SAR_Oct<- dlmModPoly(1) + dlmModReg(X=Oct_z, addInt=FALSE)
  V(mod_SAR_Oct)<- exp(parm[1])
  diag(W(mod_SAR_Oct))[1:2]<- exp(parm[2:3])
  return(mod_SAR_Oct)
  
}

fit_SAR_Oct<- dlmMLE(logit.s_SAR, parm = rep(0,3), build=buildFun_SAR_Oct, hessian=TRUE)
conv_SAR_Oct<- fit_SAR_Oct$convergence

dlmLogits_SAR_Oct<- buildFun_SAR_Oct(fit_SAR_Oct$par)

V(dlmLogits_SAR_Oct)

W(dlmLogits_SAR_Oct)

mod_SAR_Oct<- buildFun_SAR_Oct(fit_SAR_Oct$par)

SAR_Octfilter<- dlmFilter(logit.s_SAR, mod=mod_SAR_Oct)

loglik_SAR_Oct <- dlmLL(logit.s_SAR, dlmModPoly(1) + dlmModReg(X=Oct_z, addInt=FALSE))

Compare_df<- add_row(Compare_df, Covariate="Oct CUI", Model="Constant", Package="dlm", LogLik=-1*loglik_SAR_Oct)

###########################################
#Run October DLM model with MARSS
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- Oct_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; wo
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_Oct <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")
#kalman filter
kf.out <- MARSSkfss(dlm_Apr)

Compare_df<- add_row(Compare_df, Covariate="Oct CUI", Model="Constant", Package="MARSS", LogLik=dlm_Oct$logLik)

write.table(Compare_df, file='output/DLM_comparison.csv', col.names=TRUE, row.names = FALSE, sep=",")
