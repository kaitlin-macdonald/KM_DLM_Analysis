rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)
library(MARSS)

setwd("~/CBR/DLM_Analysis_Git")

cjs_raw<- readRDS("data/df_SAR_LGRLGA_MV_old.rds") 
marine1<- readRDS("marine1_trend.rds")
cjs_raw$M1<- marine1$estimate
#Remove NA at end of time series
cjs<- cjs_raw %>% filter(!is.na(SW_logit_phi), year<2021)
TT<- dim(cjs)[1]
################################
##Salmon Survive test model
year_cjs<- cjs$year
logit.s_SAR<- cjs$SW_logit_phi
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
marine1_z<- scale(cjs$M1)
##########################
#Run Reference model
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

##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.constant <- kf.out.constant$xtt1

## ts of E(forecasts)
fore_constant <- vector()
for (t in 1:TT) {
  fore_constant[t] <- Z[, , t] %*% eta.constant[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.constant <- kf.out.constant$Vtt1

## obs variance; 1x1 matrix
R_constant <- coef(dlm_constant, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi.constant[, , t] %*% tZ + R_constant
}

#This is the proper way to calculate likelihood!
sig.lik<- sum(log(abs(kf.out.constant$Sigma[1,1,])))
innov.lik<- (1/2)*sum((kf.out.constant$Innov^2)/kf.out.constant$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

constant_LL<- pie + sig.lik +innov.lik

#What to do about first obs??
MAPE_constant<- mean(abs((fore_constant-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_constant<- rmse(logit.s_SAR, fore_constant)

AICcConstant<- 2*constant_LL + 2*dlm_constant$num.params + ((2*dlm_constant$num.params*(dlm_constant$num.params+1))/(TT-dlm_constant$num.params-1))

mod_sel_df<- data.frame(Covariate=NA, Model="Constant", Data="SAR", MAPE=MAPE_constant, RMSE=rmse_constant, Href=0, LL=constant_LL, 
                        AICc=AICcConstant)


##########################
#Run CUI model
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- CUI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_CUI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.CUI <- MARSSkfas(dlm_CUI)

#Get forecast predictions
eta.CUI <- kf.CUI$xtt1

## ts of E(forecasts)
fore_CUI <- vector()
for (t in 1:TT) {
  fore_CUI[t] <- Z[, , t] %*% eta.CUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.CUI <- kf.CUI$Vtt1

## obs variance; 1x1 matrix
R_CUI <- coef(dlm_CUI, type = "matrix")$R

## ts of Var(forecasts)
CUI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  CUI_var[t] <- Z[, , t] %*% Phi.CUI[, , t] %*% tZ + R_CUI
}


#sig.lik.CUTI<- sum(log(abs(kf.CUI$Sigma[1,1,])))
#innov.lik.CUTI<- (1/2)*sum((kf.CUI$Innov^2)/kf.CUI$Sigma[1,1,])

sig.lik.CUI<- sum(log(abs(CUI_var)))
innov.lik.CUI<- (1/2)*sum((logit.s_SAR-fore_CUI)^2/CUI_var)
pie<- (TT/2)*log(2*pi)

#Negative log likelihood
CUI_LL<- pie + sig.lik.CUI +innov.lik.CUI

Href_cjs_CUI<- 2*(constant_LL - CUI_LL)

#What to do about first obs??
MAPE_CUI<- mean(abs((fore_CUI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_CUI<- rmse(logit.s_SAR, fore_CUI)

AICcCUI<- 2*CUI_LL + 2*dlm_CUI$num.params + ((2*dlm_CUI$num.params*(dlm_CUI$num.params+1))/(TT-dlm_CUI$num.params-1))


mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="Constant", Data="SAR", MAPE=MAPE_CUI, RMSE=rmse_CUI, 
                     Href=Href_cjs_CUI, LL=CUI_LL, AICc=AICcCUI)


###########################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- CUI_mean_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mnCUI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.mnCUI <- MARSSkfas(dlm_mnCUI)

#sig.lik.mnCUI<- sum(log(abs(kf.mnCUI$Sigma[1,1,])))
#innov.lik.mnCUI<- (1/2)*sum((kf.mnCUI$Innov^2)/kf.mnCUI$Sigma[1,1,])
#pie<- (TT/2)*log(2*pi)

#Get forecast predictions
eta.mnCUI <- kf.mnCUI$xtt1

## ts of E(forecasts)
fore_mnCUI <- vector()
for (t in 1:TT) {
  fore_mnCUI[t] <- Z[, , t] %*% eta.mnCUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.mnCUI <- kf.mnCUI$Vtt1

## obs variance; 1x1 matrix
R_mnCUI <- coef(dlm_mnCUI, type = "matrix")$R

## ts of Var(forecasts)
mnCUI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  mnCUI_var[t] <- Z[, , t] %*% Phi.mnCUI[, , t] %*% tZ + R_mnCUI
}

sig.lik.mnCUI<- sum(log(abs(mnCUI_var)))
innov.lik.mnCUI<- (1/2)*sum(((logit.s_SAR-fore_mnCUI)^2)/mnCUI_var)
pie<- (TT/2)*log(2*pi)

mnCUI_LL<- pie + sig.lik.mnCUI +innov.lik.mnCUI

Href_cjs_mnCUI<- 2*(constant_LL - mnCUI_LL)

#What to do about first obs??
MAPE_mnCUI<- mean(abs((fore_mnCUI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_mnCUI<- rmse(logit.s_SAR, fore_mnCUI)

AICcmnCUI<- 2*mnCUI_LL + 2*dlm_mnCUI$num.params + ((2*dlm_mnCUI$num.params*(dlm_mnCUI$num.params+1))/(TT-dlm_mnCUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnCUI", Model="Constant", Data="SAR", MAPE=MAPE_mnCUI, RMSE=rmse_mnCUI, 
                     Href=Href_cjs_mnCUI, LL=mnCUI_LL, AICc=AICcmnCUI)


#################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_PNI<- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.PNI <- MARSSkfas(dlm_PNI)

#sig.lik.PNI<- sum(log(abs(kf.PNI$Sigma[1,1,])))
#innov.lik.PNI<- (1/2)*sum((kf.PNI$Innov^2)/kf.PNI$Sigma[1,1,])
#pie<- (TT/2)*log(2*pi)

#PNI_LL<- pie + sig.lik.PNI +innov.lik.PNI

eta.PNI <- kf.PNI$xtt1

## ts of E(forecasts)
fore_PNI <- vector()
for (t in 1:TT) {
  fore_PNI[t] <- Z[, , t] %*% eta.PNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.PNI <- kf.PNI$Vtt1

## obs variance; 1x1 matrix
R_PNI <- coef(dlm_PNI, type = "matrix")$R

## ts of Var(forecasts)
PNI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  PNI_var[t] <- Z[, , t] %*% Phi.PNI[, , t] %*% tZ + R_PNI
}

sig.lik.PNI<- sum(log(abs(PNI_var)))
innov.lik.PNI<- (1/2)*sum(((logit.s_SAR-fore_PNI)^2)/PNI_var)
pie<- (TT/2)*log(2*pi)

PNI_LL<- pie + sig.lik.PNI +innov.lik.PNI

Href_cjs_PNI<- 2*(constant_LL - PNI_LL)

#What to do about first obs??
MAPE_PNI<- mean(abs((fore_PNI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_PNI<- rmse(logit.s_SAR, fore_PNI)

AICcPNI<- 2*PNI_LL + 2*dlm_PNI$num.params + ((2*dlm_PNI$num.params*(dlm_PNI$num.params+1))/(TT-dlm_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="Constant", Data="SAR", MAPE=MAPE_PNI, RMSE=rmse_PNI, 
                     Href=Href_cjs_PNI, LL=PNI_LL, AICc=AICcPNI)
#################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- SST_pre_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SSTpre <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.SSTpre <- MARSSkfss(dlm_SSTpre)

sig.lik.SSTpre<- sum(log(abs(kf.SSTpre$Sigma[1,1,])))
innov.lik.SSTpre<- (1/2)*sum((kf.SSTpre$Innov^2)/kf.SSTpre$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

SSTpre_LL<- pie + sig.lik.SSTpre +innov.lik.SSTpre

#Get forecasted predictions
eta.SSTpre <- kf.SSTpre$xtt1

## ts of E(forecasts)
fore_SSTpre <- vector()
for (t in 1:TT) {
  fore_SSTpre[t] <- Z[, , t] %*% eta.SSTpre[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.SSTpre <- kf.SSTpre$Vtt1

## obs variance; 1x1 matrix
R_SSTpre <- coef(dlm_SSTpre, type = "matrix")$R

## ts of Var(forecasts)
SSTpre_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  SSTpre_var[t] <- Z[, , t] %*% Phi.SSTpre[, , t] %*% tZ + R_SSTpre
}

#What to do about first obs??
MAPE_SSTpre<- mean(abs((fore_SSTpre-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SSTpre<- rmse(logit.s_SAR, fore_SSTpre)

Href_cjs_SSTpre<- 2*(constant_LL - SSTpre_LL)

AICcSSTpre<- 2*SSTpre_LL + 2*dlm_SSTpre$num.params + ((2*dlm_SSTpre$num.params*(dlm_SSTpre$num.params+1))/(TT-dlm_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTpre", Model="Constant", Data="SAR", MAPE=MAPE_SSTpre, RMSE=rmse_SSTpre, 
                     Href=Href_cjs_SSTpre, LL=SSTpre_LL, AICc=AICcSSTpre)

#################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- SST_entry_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SSTentry <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.SSTentry <- MARSSkfas(dlm_SSTentry)

#sig.lik.SSTentry<- sum(log(abs(kf.SSTentry$Sigma[1,1,])))
#innov.lik.SSTentry<- (1/2)*sum((kf.SSTentry$Innov^2)/kf.SSTentry$Sigma[1,1,])
#pie<- (TT/2)*log(2*pi)

#Get forecasted predictions
eta.SSTentry <- kf.SSTentry$xtt1

## ts of E(forecasts)
fore_SSTentry <- vector()
for (t in 1:TT) {
  fore_SSTentry[t] <- Z[, , t] %*% eta.SSTentry[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.SSTentry <- kf.SSTentry$Vtt1

## obs variance; 1x1 matrix
R_SSTentry <- coef(dlm_SSTentry, type = "matrix")$R

## ts of Var(forecasts)
SSTentry_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  SSTentry_var[t] <- Z[, , t] %*% Phi.SSTentry[, , t] %*% tZ + R_SSTentry
}

sig.lik.SSTentry<- sum(log(abs(SSTentry_var)))
innov.lik.SSTentry<- (1/2)*sum(((logit.s_SAR-fore_SSTentry)^2)/SSTentry_var)
pie<- (TT/2)*log(2*pi)


SSTentry_LL<- pie + sig.lik.SSTentry +innov.lik.SSTentry

#What to do about first obs??
MAPE_SSTentry<- mean(abs((fore_SSTentry-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SSTentry<- rmse(logit.s_SAR, fore_SSTentry)

Href_cjs_SSTentry<- 2*(constant_LL - SSTentry_LL)

AICcSSTentry<- 2*SSTentry_LL + 2*dlm_SSTentry$num.params + ((2*dlm_SSTentry$num.params*(dlm_SSTentry$num.params+1))/(TT-dlm_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTentry", Model="Constant", Data="SAR", MAPE=MAPE_SSTentry, RMSE=rmse_SSTentry, 
                     Href=Href_cjs_SSTentry, LL=SSTentry_LL, AICc=AICcSSTentry)

#################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- SST_6mo_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SST6mo <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.SST6mo <- MARSSkfas(dlm_SST6mo)

#Get forecasted predictions
eta.SST6mo <- kf.SST6mo$xtt1

## ts of E(forecasts)
fore_SST6mo <- vector()
for (t in 1:TT) {
  fore_SST6mo[t] <- Z[, , t] %*% eta.SST6mo[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.SST6mo$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_SST6mo, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.SST6mo<- sum(log(abs(fore_var)))
innov.lik.SST6mo<- (1/2)*sum(((logit.s_SAR-fore_SST6mo)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

SST6mo_LL<- pie + sig.lik.SST6mo +innov.lik.SST6mo

Href_cjs_SST6mo<- 2*(constant_LL - SST6mo_LL)


#What to do about first obs??
MAPE_SST6mo<- mean(abs((fore_SST6mo-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_SST6mo<- rmse(logit.s_SAR, fore_SST6mo)

Href_cjs_SST6mo<- 2*(constant_LL - SST6mo_LL)

AICcSST6mo<- 2*SST6mo_LL + 2*dlm_SST6mo$num.params + ((2*dlm_SST6mo$num.params*(dlm_SST6mo$num.params+1))/(TT-dlm_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST6mo", Model="Constant", Data="SAR", MAPE=MAPE_SST6mo, RMSE=rmse_SST6mo, 
                     Href=Href_cjs_SST6mo, LL=SST6mo_LL, AICc=AICcSST6mo)


#################################
m=3

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- CUI_z 
Z[1, 3, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_fullPNI<- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.fullPNI <- MARSSkfas(dlm_fullPNI)

#Get forecasted predictions
eta.fullPNI <- kf.fullPNI$xtt1

## ts of E(forecasts)
fore_fullPNI <- vector()
for (t in 1:TT) {
  fore_fullPNI[t] <- Z[, , t] %*% eta.fullPNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.fullPNI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_fullPNI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.fullPNI<- sum(log(abs(fore_var)))
innov.lik.fullPNI<- (1/2)*sum(((logit.s_SAR-fore_fullPNI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

fullPNI_LL<- pie + sig.lik.fullPNI +innov.lik.fullPNI

Href_cjs_fullPNI<- 2*(constant_LL - fullPNI_LL)

#What to do about first obs??
MAPE_fullPNI<- mean(abs((fore_fullPNI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_fullPNI<- rmse(logit.s_SAR, fore_fullPNI)

AICcfullPNI<- 2*fullPNI_LL + 2*dlm_fullPNI$num.params + ((2*dlm_fullPNI$num.params*(dlm_fullPNI$num.params+1))/(TT-dlm_fullPNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="CUI+PNI", Model="Constant", Data="SAR", MAPE=MAPE_fullPNI, RMSE=rmse_fullPNI, 
                     Href=Href_cjs_fullPNI, LL=fullPNI_LL, AICc=AICcfullPNI)
library(ggplot2)
years<- cjs$year
MR_pred<- cbind(fore_fullPNI, upper=fore_fullPNI+2*sqrt(fore_var), lower=fore_fullPNI-2*sqrt(fore_var))
fore_best<- as.data.frame(cbind(MR_pred, logit.s_SAR, years))
library(caret)
R2(fore_best$fore_fullPNI, fore_best$logit.s_SAR)
#ggplot(MR_pred, aes(x=years, y=fore.mean)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey") + geom_line() 
cjs_PNIM<- ggplot(fore_best, aes(x=years, y= plogis(logit.s_SAR))) + geom_point(color="blue") + geom_line(color="blue") + geom_ribbon(aes(ymin=plogis(lower), ymax=plogis(upper)), fill="grey", alpha=0.5)  + 
  geom_point(aes(x=years, y=plogis(fore_fullPNI)), color="red")+ geom_line(aes(x=years, y=plogis(fore_fullPNI)), color="red") +  ylab("SAR Survival") +xlab("Year")

ggsave("output/SW_CUIPNI.png", cjs_PNIM)

states<- t(dlm_fullPNI$states)
ggplot(states, aes(x=years, y= X2)) + geom_point() + geom_line() + ylab("CUI coefficient estimate")
 ggsave("Output/SW_CUI_Coef.png")
#################################
#################################
## Marine Trend
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- marine1_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_marine1 <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.marine1 <- MARSSkfas(dlm_marine1)

#Get forecasted predictions
eta.marine1 <- kf.marine1$xtt1

## ts of E(forecasts)
fore_marine1 <- vector()
for (t in 1:TT) {
  fore_marine1[t] <- Z[, , t] %*% eta.marine1[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.marine1 <- kf.marine1$Vtt1

## obs variance; 1x1 matrix
R_marine1 <- coef(dlm_marine1, type = "matrix")$R

## ts of Var(forecasts)
marine1_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  marine1_var[t] <- Z[, , t] %*% Phi.marine1[, , t] %*% tZ + R_marine1
}

sig.lik.marine1<- sum(log(abs(marine1_var)))
innov.lik.marine1<- (1/2)*sum(((logit.s_SAR-fore_marine1)^2)/marine1_var)
pie<- (TT/2)*log(2*pi)

marine1_LL<- pie + sig.lik.marine1 +innov.lik.marine1

#What to do about first obs??
MAPE_marine1<- mean(abs((fore_marine1-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_marine1<- rmse(logit.s_SAR, fore_marine1)

Href_cjs_marine1<- 2*(constant_LL - marine1_LL)

AICcmarine1<- 2*marine1_LL + 2*dlm_marine1$num.params + ((2*dlm_marine1$num.params*(dlm_marine1$num.params+1))/(TT-dlm_marine1$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="marine1", Model="Constant", Data="SAR", MAPE=MAPE_marine1, RMSE=rmse_marine1, 
                     Href=Href_cjs_marine1, LL=marine1_LL, AICc=AICcmarine1)

#################################
## Marine trend + Pacific Northwest Index
m=3

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- marine1_z  ## Nx1; predictor variable
Z[1, 3, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_PNIM <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.PNIM <- MARSSkfas(dlm_PNIM)

#Get forecasted predictions
eta.PNIM <- kf.PNIM$xtt1

## ts of E(forecasts)
fore_PNIM <- vector()
for (t in 1:TT) {
  fore_PNIM[t] <- Z[, , t] %*% eta.PNIM[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.PNIM <- kf.PNIM$Vtt1

## obs variance; 1x1 matrix
R_PNIM <- coef(dlm_PNIM, type = "matrix")$R

## ts of Var(forecasts)
PNIM_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  PNIM_var[t] <- Z[, , t] %*% Phi.PNIM[, , t] %*% tZ + R_PNIM
}

sig.lik.PNIM<- sum(log(abs(PNIM_var)))
innov.lik.PNIM<- (1/2)*sum(((logit.s_SAR-fore_PNIM)^2)/PNIM_var)
pie<- (TT/2)*log(2*pi)


PNIM_LL<- pie + sig.lik.PNIM +innov.lik.PNIM

#What to do about first obs??
MAPE_PNIM<- mean(abs((fore_PNIM-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_PNIM<- rmse(logit.s_SAR, fore_PNIM)

Href_cjs_PNIM<- 2*(constant_LL - PNIM_LL)

AICcPNIM<- 2*PNIM_LL + 2*dlm_PNIM$num.params + ((2*dlm_PNIM$num.params*(dlm_PNIM$num.params+1))/(TT-dlm_PNIM$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI_Marine", Model="Constant", Data="SAR", MAPE=MAPE_PNIM, RMSE=rmse_PNIM, 
                     Href=Href_cjs_PNIM, LL=PNIM_LL, AICc=AICcPNIM)

###################################################################################
#linear trend
m=2

## for process eqn
B <-  matrix(c(1,0,1,1), nrow=2, ncol=2) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lintrend <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lintrend <- MARSSkfss(dlm_lintrend)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lintrend <- kf.lintrend$xtt1

## ts of E(forecasts)
fore_lintrend <- vector()
for (t in 1:TT) {
  fore_lintrend[t] <- Z[, , t] %*% eta.lintrend[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lintrend$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lintrend, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.lintrend<- sum(log(abs(kf.lintrend$Sigma[1,1,])))
innov.lik.lintrend<- (1/2)*sum((kf.lintrend$Innov^2)/kf.lintrend$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

lintrend_LL<- pie + sig.lik.lintrend +innov.lik.lintrend

#What to do about first obs??
MAPE_ltconstant<- mean(abs((fore_lintrend-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltconstant<- rmse(logit.s_SAR, fore_lintrend)

Href_cjs_lt<- 2*(lintrend_LL - lintrend_LL)

AICclintrend<- 2*lintrend_LL + 2*dlm_lintrend$num.params + ((2*dlm_lintrend$num.params*(dlm_lintrend$num.params+1))/(TT-dlm_lintrend$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate=NA, Model="LinTrend", Data="SAR", MAPE=MAPE_ltconstant, RMSE=rmse_ltconstant, 
                     Href=Href_cjs_lt, LL=lintrend_LL, AICc=AICclintrend)
###########################################
#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- CUI_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_CUI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.CUI <- MARSSkfas(dlm_lt_CUI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.CUI <- kf.lt.CUI$xtt1

## ts of E(forecasts)
fore_lt_CUI <- vector()
for (t in 1:TT) {
  fore_lt_CUI[t] <- Z[, , t] %*% eta.lt.CUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.CUI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_CUI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltCUI<- sum(log(abs(fore_var)))
innov.lik.ltCUI<- (1/2)*sum(((logit.s_SAR-fore_lt_CUI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltCUI_LL<- pie + sig.lik.ltCUI +innov.lik.ltCUI

#What to do about first obs??
MAPE_ltCUI<- mean(abs((fore_lt_CUI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltCUI<- rmse(logit.s_SAR, fore_lt_CUI)

Href_cjs_ltCUI<- 2*(lintrend_LL - ltCUI_LL)

AICcltCUI<- 2*ltCUI_LL + 2*dlm_lt_CUI$num.params + ((2*dlm_lt_CUI$num.params*(dlm_lt_CUI$num.params+1))/(TT-dlm_lt_CUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="LinTrend", Data="SAR", MAPE=MAPE_ltCUI, RMSE=rmse_ltCUI, 
                     Href=Href_cjs_ltCUI, LL=ltCUI_LL, AICc=AICcltCUI)

###########################################
#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- CUI_mean_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_mnCUI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.mnCUI <- MARSSkfas(dlm_lt_mnCUI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.mnCUI <- kf.lt.mnCUI$xtt1

## ts of E(forecasts)
fore_lt_mnCUI <- vector()
for (t in 1:TT) {
  fore_lt_mnCUI[t] <- Z[, , t] %*% eta.lt.mnCUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.mnCUI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_mnCUI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltmnCUI<- sum(log(abs(fore_var)))
innov.lik.ltmnCUI<- (1/2)*sum(((logit.s_SAR-fore_lt_mnCUI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltmnCUI_LL<- pie + sig.lik.ltmnCUI +innov.lik.ltmnCUI

#What to do about first obs??
MAPE_ltmnCUI<- mean(abs((fore_lt_mnCUI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltmnCUI<- rmse(logit.s_SAR, fore_lt_mnCUI)

Href_cjs_ltmnCUI<- 2*(lintrend_LL - ltmnCUI_LL)

AICcltmnCUI<- 2*ltmnCUI_LL + 2*dlm_lt_mnCUI$num.params + ((2*dlm_lt_mnCUI$num.params*(dlm_lt_mnCUI$num.params+1))/(TT-dlm_lt_mnCUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mean CUI", Model="LinTrend", Data="SAR", MAPE=MAPE_ltmnCUI, RMSE=rmse_ltmnCUI, 
                     Href=Href_cjs_ltmnCUI, LL=ltmnCUI_LL, AICc=AICcltmnCUI)

###########################################

#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- PNI_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_PNI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.PNI <- MARSSkfas(dlm_lt_PNI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.PNI <- kf.lt.PNI$xtt1

## ts of E(forecasts)
fore_lt_PNI <- vector()
for (t in 1:TT) {
  fore_lt_PNI[t] <- Z[, , t] %*% eta.lt.PNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.PNI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_PNI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltPNI<- sum(log(abs(fore_var)))
innov.lik.ltPNI<- (1/2)*sum(((logit.s_SAR-fore_lt_PNI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltPNI_LL<- pie + sig.lik.ltPNI +innov.lik.ltPNI

#What to do about first obs??
MAPE_ltPNI<- mean(abs((fore_lt_PNI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltPNI<- rmse(logit.s_SAR, fore_lt_PNI)

Href_cjs_ltPNI<- 2*(lintrend_LL - ltPNI_LL)

AICcltPNI<- 2*ltPNI_LL + 2*dlm_lt_PNI$num.params + ((2*dlm_lt_PNI$num.params*(dlm_lt_PNI$num.params+1))/(TT-dlm_lt_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="LinTrend", Data="SAR", MAPE=MAPE_ltPNI, RMSE=rmse_ltPNI, 
                     Href=Href_cjs_ltPNI, LL=ltPNI_LL, AICc=AICcltPNI)

###########################################
#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- SST_pre_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_SSTpre <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SSTpre <- MARSSkfas(dlm_lt_SSTpre)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SSTpre <- kf.lt.SSTpre$xtt1

## ts of E(forecasts)
fore_lt_SSTpre <- vector()
for (t in 1:TT) {
  fore_lt_SSTpre[t] <- Z[, , t] %*% eta.lt.SSTpre[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.SSTpre$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_SSTpre, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltSSTpre<- sum(log(abs(fore_var)))
innov.lik.ltSSTpre<- (1/2)*sum(((logit.s_SAR-fore_lt_SSTpre)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSSTpre_LL<- pie + sig.lik.ltSSTpre +innov.lik.ltSSTpre

#What to do about first obs??
MAPE_ltSSTpre<- mean(abs((fore_lt_SSTpre-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltSSTpre<- rmse(logit.s_SAR, fore_lt_SSTpre)

Href_cjs_ltSSTpre<- 2*(lintrend_LL - ltSSTpre_LL)

AICcltSSTpre<- 2*ltSSTpre_LL + 2*dlm_lt_SSTpre$num.params + ((2*dlm_lt_SSTpre$num.params*(dlm_lt_SSTpre$num.params+1))/(TT-dlm_lt_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_pre", Model="LinTrend", Data="SAR", MAPE=MAPE_ltSSTpre, RMSE=rmse_ltSSTpre, 
                     Href=Href_cjs_ltSSTpre, LL=ltSSTpre_LL, AICc=AICcltSSTpre)

###########################################
#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- SST_entry_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_SSTentry <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SSTentry <- MARSSkfas(dlm_lt_SSTentry)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SSTentry <- kf.lt.SSTentry$xtt1

## ts of E(forecasts)
fore_lt_SSTentry <- vector()
for (t in 1:TT) {
  fore_lt_SSTentry[t] <- Z[, , t] %*% eta.lt.SSTentry[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.SSTentry$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_SSTentry, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}


#This is the proper way to calculate likelihood!
sig.lik.ltSSTentry<- sum(log(abs(fore_var)))
innov.lik.ltSSTentry<- (1/2)*sum(((logit.s_SAR-fore_lt_SSTpre)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSSTentry_LL<- pie + sig.lik.ltSSTentry +innov.lik.ltSSTentry

#What to do about first obs??
MAPE_ltSSTentry<- mean(abs((fore_lt_SSTentry-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltSSTentry<- rmse(logit.s_SAR, fore_lt_SSTentry)

Href_cjs_ltSSTentry<- 2*(lintrend_LL - ltSSTentry_LL)

AICcltSSTentry<- 2*ltSSTentry_LL + 2*dlm_lt_SSTentry$num.params + ((2*dlm_lt_SSTentry$num.params*(dlm_lt_SSTentry$num.params+1))/(TT-dlm_lt_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_entry", Model="LinTrend", Data="SAR", MAPE=MAPE_ltSSTentry, RMSE=rmse_ltSSTentry, 
                     Href=Href_cjs_ltSSTentry, LL=ltSSTentry_LL, AICc=AICcltSSTentry)

###########################################
#linear trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- SST_6mo_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_SST6mo <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SST6mo <- MARSSkfas(dlm_lt_SST6mo)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SST6mo <- kf.lt.SST6mo$xtt1

## ts of E(forecasts)
fore_lt_SST6mo <- vector()
for (t in 1:TT) {
  fore_lt_SST6mo[t] <- Z[, , t] %*% eta.lt.SST6mo[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.SST6mo$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_SST6mo, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}


#This is the proper way to calculate likelihood!
sig.lik.ltSST6mo<- sum(log(abs(fore_var)))
innov.lik.ltSST6mo<- (1/2)*sum(((logit.s_SAR-fore_lt_SST6mo)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSST6mo_LL<- pie + sig.lik.ltSST6mo +innov.lik.ltSST6mo

#What to do about first obs??
MAPE_ltSST6mo<- mean(abs((fore_lt_SST6mo-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltSST6mo<- rmse(logit.s_SAR, fore_lt_SST6mo)

Href_cjs_ltSST6mo<- 2*(lintrend_LL - ltSST6mo_LL)

AICcltSST6mo<- 2*ltSST6mo_LL + 2*dlm_lt_SST6mo$num.params + ((2*dlm_lt_SST6mo$num.params*(dlm_lt_SST6mo$num.params+1))/(TT-dlm_lt_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_6mo", Model="LinTrend", Data="SAR", MAPE=MAPE_ltSST6mo, RMSE=rmse_ltSST6mo, 
                     Href=Href_cjs_ltSST6mo, LL=ltSST6mo_LL, AICc=AICcltSST6mo)
###################################
#linear trend Marine trend
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- marine1_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_marine <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.marine <- MARSSkfas(dlm_lt_marine)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.marine <- kf.lt.marine$xtt1

## ts of E(forecasts)
fore_lt_marine <- vector()
for (t in 1:TT) {
  fore_lt_marine[t] <- Z[, , t] %*% eta.lt.marine[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi_marine <- kf.lt.marine$Vtt1

## obs variance; 1x1 matrix
R_marine <- coef(dlm_lt_marine, type = "matrix")$R

## ts of Var(forecasts)
fore_marine <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_marine[t] <- Z[, , t] %*% Phi_marine[, , t] %*% tZ + R_marine
}


#This is the proper way to calculate likelihood!
sig.lik.ltmarine<- sum(log(abs(fore_marine)))
innov.lik.ltmarine<- (1/2)*sum(((logit.s_SAR-fore_lt_marine)^2)/fore_marine)
pie<- (TT/2)*log(2*pi)

ltmarine_LL<- pie + sig.lik.ltmarine +innov.lik.ltmarine

#What to do about first obs??
MAPE_ltmarine<- mean(abs((fore_lt_marine-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltmarine<- rmse(logit.s_SAR, fore_lt_marine)

Href_cjs_ltmarine<- 2*(lintrend_LL - ltmarine_LL)


AICcltmarine<- 2*ltmarine_LL + 2*dlm_lt_marine$num.params + ((2*dlm_lt_marine$num.params*(dlm_lt_marine$num.params+1))/(TT-dlm_lt_marine$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine", Model="LinTrend", Data="SAR", MAPE=MAPE_ltmarine, RMSE=rmse_ltmarine, 
                     Href=Href_cjs_ltmarine, LL=ltmarine_LL, AICc=AICcltmarine)

###########################################
#linear trend Marine trend + pacific northwest index
m=4

## for process eqn
B <-  matrix(c(1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), nrow=4, ncol=4) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.beta1", "q.beta2") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- marine1_z
Z[1, 4, ] <- PNI_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_marinePNI <- MARSS(logit.s_SAR, inits = inits_list, model = mod_list, method="TMB")

kf.lt.marinePNI <- MARSSkfas(dlm_lt_marinePNI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.marinePNI <- kf.lt.marinePNI$xtt1

## ts of E(forecasts)
fore_lt_marinePNI <- vector()
for (t in 1:TT) {
  fore_lt_marinePNI[t] <- Z[, , t] %*% eta.lt.marinePNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi_marinePNI <- kf.lt.marinePNI$Vtt1

## obs variance; 1x1 matrix
R_marinePNI <- coef(dlm_lt_marinePNI, type = "matrix")$R

## ts of Var(forecasts)
fore_marinePNI <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_marinePNI[t] <- Z[, , t] %*% Phi_marinePNI[, , t] %*% tZ + R_marinePNI
}


#This is the proper way to calculate likelihood!
sig.lik.ltmarinePNI<- sum(log(abs(fore_marinePNI)))
innov.lik.ltmarinePNI<- (1/2)*sum(((logit.s_SAR-fore_lt_marinePNI)^2)/fore_marinePNI)
pie<- (TT/2)*log(2*pi)

ltmarinePNI_LL<- pie + sig.lik.ltmarinePNI +innov.lik.ltmarinePNI

#What to do about first obs??
MAPE_ltmarinePNI<- mean(abs((fore_lt_marinePNI-logit.s_SAR)/logit.s_SAR))
#Root mean square error
rmse_ltmarinePNI<- rmse(logit.s_SAR, fore_lt_marinePNI)

Href_cjs_ltmarinePNI<- 2*(lintrend_LL - ltmarinePNI_LL)


AICcltmarinePNI<- 2*ltmarinePNI_LL + 2*dlm_lt_marinePNI$num.params + ((2*dlm_lt_marinePNI$num.params*(dlm_lt_marinePNI$num.params+1))/(TT-dlm_lt_marinePNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine + PNI", Model="LinTrend", Data="SAR", MAPE=MAPE_ltmarinePNI, RMSE=rmse_ltmarinePNI, 
                     Href=Href_cjs_ltmarinePNI, LL=ltmarinePNI_LL, AICc=AICcltmarinePNI)

mod_selection<- mod_sel_df %>% arrange(AICc) %>% select(-MAPE)
write.csv(mod_selection, "Output/SW_Univariate_sorted.csv")

write.csv(mod_sel_df, "Output/SW_Univariate_DLM.csv")
