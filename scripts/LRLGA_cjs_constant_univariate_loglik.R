rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)
library(MARSS)

setwd("~/CBR/DLM_Analysis_Git")

#Read in data
cjs_raw<- readRDS("data/df_SAR_LGRLGA_MV_old.rds") 
marine1<- readRDS("marine1_trend.rds")

swepeak<- read.csv("data/swepeak.csv") %>% select(year, day, sum, mean)

mean_spill<- readRDS("data/mean_spill.rds") %>% mutate(year=as.double(year))

cjs_raw$M1<- marine1$estimate

#Remove NA at end of time series
cjs<- cjs_raw %>% left_join(swepeak, by="year") %>% left_join(mean_spill, by="year") %>%  filter(!is.na(cjs_logit_phi), year<2021)
corplot<- cor(cjs)
TT<- dim(cjs)[1]
################################
##Salmon Survive test model
year_cjs<- cjs$year
logit.s_cjs<- cjs$cjs_logit_phi
CUI_z<- scale(cjs$CUI)
CUTI_z<- scale(cjs$CUTI)
BEUTI_z<- scale(cjs$BEUTI)
mldno3_z<- scale(cjs$Apr.mldno3)
CUI_mean_z<- scale(cjs$mean_CUI)
CUTI_mean_z<- scale(cjs$mean_CUTI)
BEUTI_mean_z<- scale(cjs$mean_BEUTI)
mldno3_mean_z<- scale(cjs$mean_mldno3)
outflow_z<- scale(cjs$mean.out) #Mean april outflow
temp_z<- scale(cjs$mean.temp) # surface temperature (WQM) at lower granite forebay mean for april
PNI_z<- scale(cjs$Annual.PNI)
SST_pre_z<- scale(cjs$SST_pre1212)
SST_entry_z<- scale(cjs$SST_entry567)
SST_6mo_z<- scale(cjs$SST_6mo)
marine1_z<- scale(cjs$M1)
peak_day_z<- scale(cjs$day)
peak_swe_z<- scale(cjs$mean)
spill_z<- scale(cjs$per_spill)
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
dlm_constant <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
Phi <- kf.out.constant$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_constant, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik<- sum(log(abs(kf.out.constant$Sigma[1,1,])))
innov.lik<- (1/2)*sum((kf.out.constant$Innov^2)/kf.out.constant$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

constant_LL<- pie + sig.lik +innov.lik

#What to do about first obs??
MAPE_constant<- mean(abs((fore_constant-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_constant<- rmse(logit.s_cjs, fore_constant)

AICcConstant<- 2*constant_LL + 2*dlm_constant$num.params + ((2*dlm_constant$num.params*(dlm_constant$num.params+1))/(TT-dlm_constant$num.params-1))


mod_sel_df<- data.frame(Covariate=NA, Model="Constant", Data="CJS", MAPE=MAPE_constant, RMSE=rmse_constant, Href=0, LL=constant_LL, AICc=AICcConstant)


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
# Do we want to fix R in the model using the estimated observation error when no covariates are included
R <- matrix(0.151)
## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_CUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.CUI <- MARSSkfss(dlm_CUI)

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

sig.lik.CUI<- sum(log(abs(CUI_var)))
innov.lik.CUI<- (1/2)*sum(((logit.s_cjs-fore_CUI)^2)/CUI_var)
pie<- (TT/2)*log(2*pi)

CUI_LL<- pie + sig.lik.CUI +innov.lik.CUI

Href_cjs_CUI<- 2*(constant_LL - CUI_LL)

#What to do about first obs??
MAPE_CUI<- mean(abs((fore_CUI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_CUI<- rmse(logit.s_cjs, fore_CUI)

AICcCUI<- 2*CUI_LL + 2*dlm_CUI$num.params + ((2*dlm_CUI$num.params*(dlm_CUI$num.params+1))/(TT-dlm_CUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="Constant", Data="CJS", MAPE=MAPE_CUI, RMSE=rmse_CUI, 
                     Href=Href_cjs_CUI, LL=CUI_LL, AICc=AICcCUI)

###########################################
#Run CUTI models
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- CUTI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_CUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.CUTI <- MARSSkfas(dlm_CUTI)

#Get forecast predictions
eta.CUTI <- kf.CUTI$xtt1

## ts of E(forecasts)
fore_CUTI <- vector()
for (t in 1:TT) {
  fore_CUTI[t] <- Z[, , t] %*% eta.CUTI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.CUTI <- kf.CUTI$Vtt1

## obs variance; 1x1 matrix
R_CUTI <- coef(dlm_CUTI, type = "matrix")$R

## ts of Var(forecasts)
CUTI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  CUTI_var[t] <- Z[, , t] %*% Phi.CUTI[, , t] %*% tZ + R_CUTI
}


sig.lik.CUTI<- sum(log(abs(CUTI_var)))
innov.lik.CUTI<- (1/2)*sum(((logit.s_cjs-fore_CUTI)^2)/CUTI_var)
pie<- (TT/2)*log(2*pi)

CUTI_LL<- pie + sig.lik.CUTI +innov.lik.CUTI

Href_cjs_CUTI<- 2*(constant_LL - CUTI_LL)

#What to do about first obs??
MAPE_CUTI<- mean(abs((fore_CUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_CUTI<- rmse(logit.s_cjs, fore_CUTI)

AICcCUTI<- 2*CUTI_LL + 2*dlm_CUTI$num.params + ((2*dlm_CUTI$num.params*(dlm_CUTI$num.params+1))/(TT-dlm_CUTI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="CUTI", Model="Constant", Data="CJS", MAPE=MAPE_CUTI, RMSE=rmse_CUTI, 
                     Href=Href_cjs_CUTI, LL=CUTI_LL, AICc=AICcCUTI)

###########################################
## Run BEUTI
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- BEUTI_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_BEUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.BEUTI <- MARSSkfas(dlm_BEUTI) 

#Get forecast predictions
eta.BEUTI <- kf.BEUTI$xtt1

## ts of E(forecasts)
fore_BEUTI <- vector()
for (t in 1:TT) {
  fore_BEUTI[t] <- Z[, , t] %*% eta.BEUTI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.BEUTI <- kf.BEUTI$Vtt1

## obs variance; 1x1 matrix
R_BEUTI <- coef(dlm_BEUTI, type = "matrix")$R

## ts of Var(forecasts)
BEUTI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  BEUTI_var[t] <- Z[, , t] %*% Phi.BEUTI[, , t] %*% tZ + R_BEUTI
}

sig.lik.BEUTI<- sum(log(abs(BEUTI_var)))
innov.lik.BEUTI<- (1/2)*sum(((logit.s_cjs-fore_BEUTI)^2)/BEUTI_var)
pie<- (TT/2)*log(2*pi)

BEUTI_LL<- pie + sig.lik.BEUTI +innov.lik.BEUTI

Href_cjs_BEUTI<- 2*(constant_LL - BEUTI_LL)

#What to do about first obs??
MAPE_BEUTI<- mean(abs((fore_BEUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_BEUTI<- rmse(logit.s_cjs, fore_BEUTI)

AICcBEUTI<- 2*BEUTI_LL + 2*dlm_BEUTI$num.params + ((2*dlm_BEUTI$num.params*(dlm_BEUTI$num.params+1))/(TT-dlm_BEUTI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="BEUTI", Model="Constant", Data="CJS", MAPE=MAPE_BEUTI, RMSE=rmse_BEUTI, 
                     Href=Href_cjs_BEUTI, LL=BEUTI_LL, AICc=AICcBEUTI)

###########################################
## Run MLDNO3
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- mldno3_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mldno3 <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.mld <- MARSSkfss(dlm_mldno3)

sig.lik.ml<- sum(log(abs(kf.mld$Sigma[1,1,])))
innov.lik.ml<- (1/2)*sum((kf.mld$Innov^2)/kf.mld$Sigma[1,1,])
pie.ml<- (TT/2)*log(2*pi)

mld_LL<- pie.ml + sig.lik.ml +innov.lik.ml

Href_cjs_mld<- 2*(constant_LL - mld_LL)

#Get forecast predictions
eta.mld <- kf.mld$xtt1

## ts of E(forecasts)
fore_mld <- vector()
for (t in 1:TT) {
  fore_mld[t] <- Z[, , t] %*% eta.mld[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_mld<- mean(abs((fore_mld-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_mld<- rmse(logit.s_cjs, fore_mld)

AICcmld<- 2*mld_LL + 2*dlm_mldno3$num.params + ((2*dlm_mldno3$num.params*(dlm_mldno3$num.params+1))/(TT-dlm_mldno3$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="MLDNO3", Model="Constant", Data="CJS", MAPE=MAPE_mld, RMSE=rmse_mld, 
                     Href=Href_cjs_mld, LL=mld_LL, AICc=AICcmld)
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


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mnCUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.mnCUI<- (1/2)*sum(((logit.s_cjs-fore_mnCUI)^2)/mnCUI_var)
pie<- (TT/2)*log(2*pi)

mnCUI_LL<- pie + sig.lik.mnCUI +innov.lik.mnCUI

Href_cjs_mnCUI<- 2*(constant_LL - mnCUI_LL)

#What to do about first obs??
MAPE_mnCUI<- mean(abs((fore_mnCUI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_mnCUI<- rmse(logit.s_cjs, fore_mnCUI)

AICcmnCUI<- 2*mnCUI_LL + 2*dlm_mnCUI$num.params + ((2*dlm_mnCUI$num.params*(dlm_mnCUI$num.params+1))/(TT-dlm_mnCUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnCUI", Model="Constant", Data="CJS", MAPE=MAPE_mnCUI, RMSE=rmse_mnCUI, 
                     Href=Href_cjs_mnCUI, LL=mnCUI_LL, AICc=AICcmnCUI)

###########################################
## mean CUTI
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- CUTI_mean_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mnCUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.mnCUTI <- MARSSkfss(dlm_mnCUTI)

sig.lik.mnCUTI<- sum(log(abs(kf.mnCUTI$Sigma[1,1,])))
innov.lik.mnCUTI<- (1/2)*sum((kf.mnCUTI$Innov^2)/kf.mnCUTI$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

mnCUTI_LL<- pie + sig.lik.mnCUTI +innov.lik.mnCUTI

Href_cjs_mnCUTI<- 2*(constant_LL - mnCUTI_LL)

#Get forecast predictions
eta.mnCUTI <- kf.mnCUTI$xtt1

## ts of E(forecasts)
fore_mnCUTI <- vector()
for (t in 1:TT) {
  fore_mnCUTI[t] <- Z[, , t] %*% eta.mnCUTI[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_mnCUTI<- mean(abs((fore_mnCUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_mnCUTI<- rmse(logit.s_cjs, fore_mnCUTI)

AICcmnCUTI<- 2*mnCUTI_LL + 2*dlm_mnCUTI$num.params + ((2*dlm_mnCUTI$num.params*(dlm_mnCUTI$num.params+1))/(TT-dlm_mnCUTI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnCUTI", Model="Constant", Data="CJS", MAPE=MAPE_mnCUTI, RMSE=rmse_mnCUTI, 
                     Href=Href_cjs_mnCUTI, LL=mnCUTI_LL, AICc=AICcmnCUTI)

###########################################
## mean BEUTI
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- BEUTI_mean_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mnBEUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.mnBEUTI <- MARSSkfss(dlm_mnBEUTI)

sig.lik.mnBEUTI<- sum(log(abs(kf.mnBEUTI$Sigma[1,1,])))
innov.lik.mnBEUTI<- (1/2)*sum((kf.mnBEUTI$Innov^2)/kf.mnBEUTI$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

mnBEUTI_LL<- pie + sig.lik.mnBEUTI +innov.lik.mnBEUTI

#Get forecasted predictions
eta.mnBEUTI <- kf.mnBEUTI$xtt1

## ts of E(forecasts)
fore_mnBEUTI <- vector()
for (t in 1:TT) {
  fore_mnBEUTI[t] <- Z[, , t] %*% eta.mnBEUTI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.mnBEUTI <- kf.mnBEUTI$Vtt1

## obs variance; 1x1 matrix
R_mnBEUTI <- coef(dlm_mnBEUTI, type = "matrix")$R

## ts of Var(forecasts)
mnBEUTI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  mnBEUTI_var[t] <- Z[, , t] %*% Phi.mnBEUTI[, , t] %*% tZ + R_mnBEUTI
}

#What to do about first obs??
MAPE_mnBEUTI<- mean(abs((fore_mnBEUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_mnBEUTI<- rmse(logit.s_cjs, fore_mnBEUTI)

Href_cjs_BEUTI<- 2*(constant_LL - mnBEUTI_LL)

AICcmnBEUTI<- 2*mnBEUTI_LL + 2*dlm_mnBEUTI$num.params + ((2*dlm_mnBEUTI$num.params*(dlm_mnBEUTI$num.params+1))/(TT-dlm_mnBEUTI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnBEUTI", Model="Constant", Data="CJS", MAPE=MAPE_mnBEUTI, RMSE=rmse_mnBEUTI, 
                     Href=Href_cjs_BEUTI, LL=mnBEUTI_LL, AICc=AICcmnBEUTI)

###########################################
## mean MLDNO3
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- mldno3_mean_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_mnmld<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.mnmld <- MARSSkfss(dlm_mnmld)

sig.lik.mnmld<- sum(log(abs(kf.mnmld$Sigma[1,1,])))
innov.lik.mnmld<- (1/2)*sum((kf.mnmld$Innov^2)/kf.mnmld$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

mnmld_LL<- pie + sig.lik.mnmld +innov.lik.mnmld


Href_cjs_mnmld<- 2*(constant_LL - mnmld_LL)

eta.mnmld <- kf.mnmld$xtt1

## ts of E(forecasts)
fore_mnmld <- vector()
for (t in 1:TT) {
  fore_mnmld[t] <- Z[, , t] %*% eta.mnmld[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_mnmld<- mean(abs((fore_mnmld-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_mnmld<- rmse(logit.s_cjs, fore_mnmld)

AICcmnmld<- 2*mnmld_LL + 2*dlm_mnmld$num.params + ((2*dlm_mnmld$num.params*(dlm_mnmld$num.params+1))/(TT-dlm_mnmld$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnmld", Model="Constant", Data="CJS", MAPE=MAPE_mnmld, RMSE=rmse_mnmld, 
                     Href=Href_cjs_mnmld, LL=mnmld_LL, AICc=AICcmnmld)

###########################################
## outflow
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- outflow_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_outflow<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.outflow <- MARSSkfss(dlm_outflow)

sig.lik.outflow<- sum(log(abs(kf.outflow$Sigma[1,1,])))
innov.lik.outflow<- (1/2)*sum((kf.outflow$Innov^2)/kf.outflow$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

outflow_LL<- pie + sig.lik.outflow +innov.lik.outflow

Href_cjs_outflow<- 2*(constant_LL - outflow_LL)

eta.outflow <- kf.outflow$xtt1

## ts of E(forecasts)
fore_outflow <- vector()
for (t in 1:TT) {
  fore_outflow[t] <- Z[, , t] %*% eta.outflow[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_outflow<- mean(abs((fore_outflow-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_outflow<- rmse(logit.s_cjs, fore_outflow)

AICcoutflow<- 2*outflow_LL + 2*dlm_outflow$num.params + ((2*dlm_outflow$num.params*(dlm_outflow$num.params+1))/(TT-dlm_outflow$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="outflow", Model="Constant", Data="CJS", MAPE=MAPE_outflow, RMSE=rmse_outflow, 
                     Href=Href_cjs_outflow, LL=outflow_LL, AICc=AICcoutflow)
###########################################
## River temperature
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- temp_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_temp<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.temp <- MARSSkfss(dlm_temp)

sig.lik.temp<- sum(log(abs(kf.temp$Sigma[1,1,])))
innov.lik.temp<- (1/2)*sum((kf.temp$Innov^2)/kf.temp$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

temp_LL<- pie + sig.lik.temp +innov.lik.temp

Href_cjs_temp<- 2*(constant_LL - temp_LL)

eta.temp <- kf.temp$xtt1

## ts of E(forecasts)
fore_temp <- vector()
for (t in 1:TT) {
  fore_temp[t] <- Z[, , t] %*% eta.temp[, t, drop = FALSE]
}


#What to do about first obs??
MAPE_temp<- mean(abs((fore_temp-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_temp<- rmse(logit.s_cjs, fore_temp)

AICctemp<- 2*temp_LL + 2*dlm_temp$num.params + ((2*dlm_temp$num.params*(dlm_temp$num.params+1))/(TT-dlm_temp$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="temp", Model="Constant", Data="CJS", MAPE=MAPE_temp, RMSE=rmse_temp, 
                     Href=Href_cjs_temp, LL=temp_LL, AICc=AICctemp)

#################################
## SWE peak day
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- peak_day_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_peak_day<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

## fit univariate DLM
dlm_peak_day <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.pkday <- MARSSkfss(dlm_peak_day)

#Get forecasted predictions
eta.pkday <- kf.pkday$xtt1

## ts of E(forecasts)
fore_pkday <- vector()
for (t in 1:TT) {
  fore_pkday[t] <- Z[, , t] %*% eta.pkday[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.marine1 <- kf.marine1$Vtt1

## obs variance; 1x1 matrix
#R_marine1 <- coef(dlm_marine1, type = "matrix")$R

## ts of Var(forecasts)
#marine1_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  marine1_var[t] <- Z[, , t] %*% Phi.marine1[, , t] %*% tZ + R_marine1
#}

sig.lik.pkday<- sum(log(abs(kf.pkday$Sigma[1,1,])))
innov.lik.pkday<- (1/2)*sum((kf.pkday$Innov^2)/kf.pkday$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

pkday_LL<- pie + sig.lik.pkday +innov.lik.pkday

#What to do about first obs??
MAPE_pkday<- mean(abs((fore_pkday-logit.s_cjs)/logit.s_cjs), na.rm=TRUE)
#Root mean square error
rmse_pkday<- sqrt(mean((logit.s_cjs - fore_pkday)^2, na.rm = TRUE))

Href_cjs_pkday<- 2*(constant_LL - pkday_LL)

AICcpkday<- 2*pkday_LL + 2*dlm_peak_day$num.params + ((2*dlm_peak_day$num.params*(dlm_peak_day$num.params+1))/(TT-dlm_peak_day$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="peak day", Model="Constant", Data="Multi", MAPE=MAPE_pkday, RMSE=rmse_pkday, 
                     Href=Href_cjs_pkday, LL=pkday_LL, AICc=AICcpkday)

#################################
## SWE peak swe
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- peak_swe_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)


## fit univariate DLM
dlm_peak <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.pk <- MARSSkfas(dlm_peak)

#Get forecasted predictions
eta.pk <- kf.pk$xtt1

## ts of E(forecasts)
fore_pk <- vector()
for (t in 1:TT) {
  fore_pk[t] <- Z[, , t] %*% eta.pkday[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.pk <- kf.pk$Vtt1

## obs variance; 1x1 matrix
R_pk <- coef(dlm_peak, type = "matrix")$R

## ts of Var(forecasts)
peak_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  peak_var[t] <- Z[, , t] %*% Phi.pk[, , t] %*% tZ + R_pk
}

sig.lik.pk<- sum(log(abs(peak_var)))
innov.lik.pk<- (1/2)*sum(((logit.s_cjs-fore_pk)^2)/peak_var)
pie<- (TT/2)*log(2*pi)
pk_LL<- pie + sig.lik.pk +innov.lik.pk

#What to do about first obs??
MAPE_pk<- mean(abs((fore_pk-logit.s_cjs)/logit.s_cjs), na.rm=TRUE)
#Root mean square error
rmse_pk<- sqrt(mean((logit.s_cjs - fore_pk)^2, na.rm = TRUE))

Href_cjs_pk<- 2*(constant_LL - pk_LL)

AICcpk<- 2*pk_LL + 2*dlm_peak$num.params + ((2*dlm_peak$num.params*(dlm_peak$num.params+1))/(TT-dlm_peak$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="peak", Model="Constant", Data="Multi", MAPE=MAPE_pk, RMSE=rmse_pk, 
                     Href=Href_cjs_pk, LL=pk_LL, AICc=AICcpk)


#################################
## Percent spill
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- spill_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_spill <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.spill <- MARSSkfss(dlm_spill)

#Get forecasted predictions
eta.spill <- kf.spill$xtt1

## ts of E(forecasts)
fore_spill <- vector()
for (t in 1:TT) {
  fore_spill[t] <- Z[, , t] %*% eta.spill[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.marine1 <- kf.marine1$Vtt1

## obs variance; 1x1 matrix
#R_marine1 <- coef(dlm_marine1, type = "matrix")$R

## ts of Var(forecasts)
#marine1_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  marine1_var[t] <- Z[, , t] %*% Phi.marine1[, , t] %*% tZ + R_marine1
#}

sig.lik.spill<- sum(log(abs(kf.spill$Sigma[1,1,])))
innov.lik.spill<- (1/2)*sum((kf.spill$Innov^2)/kf.spill$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

spill_LL<- pie + sig.lik.spill +innov.lik.spill

#What to do about first obs??
MAPE_spill<- mean(abs((fore_spill-logit.s_cjs)/logit.s_cjs), na.rm=TRUE)
#Root mean square error
rmse_spill<- sqrt(mean((logit.s_cjs - fore_spill)^2, na.rm = TRUE))

Href_cjs_spill<- 2*(constant_LL - spill_LL)

AICcspill<- 2*spill_LL + 2*dlm_spill$num.params + ((2*dlm_spill$num.params*(dlm_spill$num.params+1))/(TT-dlm_spill$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="spill", Model="Constant", Data="Multi", MAPE=MAPE_spill, RMSE=rmse_spill, 
                     Href=Href_cjs_spill, LL=spill_LL, AICc=AICcspill)


#################################
## Pacific Northwest Index
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
dlm_PNI<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.PNI <- MARSSkfas(dlm_PNI)

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
innov.lik.PNI<- (1/2)*sum(((logit.s_cjs-fore_PNI)^2)/PNI_var)
pie<- (TT/2)*log(2*pi)

PNI_LL<- pie + sig.lik.PNI +innov.lik.PNI

Href_cjs_PNI<- 2*(constant_LL - PNI_LL)

#What to do about first obs??
MAPE_PNI<- mean(abs((fore_PNI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_PNI<- rmse(logit.s_cjs, fore_PNI)

AICcPNI<- 2*PNI_LL + 2*dlm_PNI$num.params + ((2*dlm_PNI$num.params*(dlm_PNI$num.params+1))/(TT-dlm_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="Constant", Data="CJS", MAPE=MAPE_PNI, RMSE=rmse_PNI, 
                     Href=Href_cjs_PNI, LL=PNI_LL, AICc=AICcPNI)

#################################
## Sea surface temperature pre-ocean entry
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

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SSTpre <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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

#What to do about first obs??
MAPE_SSTpre<- mean(abs((fore_SSTpre-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_SSTpre<- rmse(logit.s_cjs, fore_SSTpre)

Href_cjs_SSTpre<- 2*(constant_LL - SSTpre_LL)

AICcSSTpre<- 2*SSTpre_LL + 2*dlm_SSTpre$num.params + ((2*dlm_SSTpre$num.params*(dlm_SSTpre$num.params+1))/(TT-dlm_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTpre", Model="Constant", Data="CJS", MAPE=MAPE_SSTpre, RMSE=rmse_SSTpre, 
                     Href=Href_cjs_SSTpre, LL=SSTpre_LL, AICc=AICcSSTpre)

#################################
## Sea surface temperature at ocean entry
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

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SSTentry <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.SSTentry <- MARSSkfas(dlm_SSTentry)

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
innov.lik.SSTentry<- (1/2)*sum(((logit.s_cjs-fore_SSTentry)^2)/SSTentry_var)
pie<- (TT/2)*log(2*pi)

SSTentry_LL<- pie + sig.lik.SSTentry +innov.lik.SSTentry

#What to do about first obs??
MAPE_SSTentry<- mean(abs((fore_SSTentry-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_SSTentry<- rmse(logit.s_cjs, fore_SSTentry)

Href_cjs_SSTentry<- 2*(constant_LL - SSTentry_LL)

AICcSSTentry<- 2*SSTentry_LL + 2*dlm_SSTentry$num.params + ((2*dlm_SSTentry$num.params*(dlm_SSTentry$num.params+1))/(TT-dlm_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTentry", Model="Constant", Data="CJS", MAPE=MAPE_SSTentry, RMSE=rmse_SSTentry, 
                     Href=Href_cjs_SSTentry, LL=SSTentry_LL, AICc=AICcSSTentry)

#################################
## Sea surface temperature 6 month average
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

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_SST6mo <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.SST6mo <- MARSSkfss(dlm_SST6mo)

sig.lik.SST6mo<- sum(log(abs(kf.SST6mo$Sigma[1,1,])))
innov.lik.SST6mo<- (1/2)*sum((kf.SST6mo$Innov^2)/kf.SST6mo$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

SST6mo_LL<- pie + sig.lik.SST6mo +innov.lik.SST6mo

#Get forecasted predictions
eta.SST6mo <- kf.SST6mo$xtt1

## ts of E(forecasts)
fore_SST6mo <- vector()
for (t in 1:TT) {
  fore_SST6mo[t] <- Z[, , t] %*% eta.SST6mo[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_SST6mo<- mean(abs((fore_SST6mo-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_SST6mo<- rmse(logit.s_cjs, fore_SST6mo)

Href_cjs_SST6mo<- 2*(constant_LL - SST6mo_LL)

AICcSST6mo<- 2*SST6mo_LL + 2*dlm_SST6mo$num.params + ((2*dlm_SST6mo$num.params*(dlm_SST6mo$num.params+1))/(TT-dlm_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST6mo", Model="Constant", Data="CJS", MAPE=MAPE_SST6mo, RMSE=rmse_SST6mo, 
                     Href=Href_cjs_SST6mo, LL=SST6mo_LL, AICc=AICcSST6mo)

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
dlm_marine1 <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.marine1<- (1/2)*sum(((logit.s_cjs-fore_marine1)^2)/marine1_var)
pie<- (TT/2)*log(2*pi)

marine1_LL<- pie + sig.lik.marine1 +innov.lik.marine1

#What to do about first obs??
MAPE_marine1<- mean(abs((fore_marine1-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_marine1<- rmse(logit.s_cjs, fore_marine1)

Href_cjs_marine1<- 2*(constant_LL - marine1_LL)

AICcmarine1<- 2*marine1_LL + 2*dlm_marine1$num.params + ((2*dlm_marine1$num.params*(dlm_marine1$num.params+1))/(TT-dlm_marine1$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="marine1", Model="Constant", Data="CJS", MAPE=MAPE_marine1, RMSE=rmse_marine1, 
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
dlm_PNIM <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.PNIM<- (1/2)*sum(((logit.s_cjs-fore_PNIM)^2)/PNIM_var)
pie<- (TT/2)*log(2*pi)


PNIM_LL<- pie + sig.lik.PNIM +innov.lik.PNIM

#What to do about first obs??
MAPE_PNIM<- mean(abs((fore_PNIM-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_PNIM<- rmse(logit.s_cjs, fore_PNIM)

Href_cjs_PNIM<- 2*(constant_LL - PNIM_LL)

AICcPNIM<- 2*PNIM_LL + 2*dlm_PNIM$num.params + ((2*dlm_PNIM$num.params*(dlm_PNIM$num.params+1))/(TT-dlm_PNIM$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI_Marine", Model="Constant", Data="CJS", MAPE=MAPE_PNIM, RMSE=rmse_PNIM, 
                     Href=Href_cjs_PNIM, LL=PNIM_LL, AICc=AICcPNIM)
library(ggplot2)
years<- cjs$year
MR_pred<- cbind(fore_PNIM, upper=fore_PNIM+2*sqrt(PNIM_var), lower=fore_PNIM-2*sqrt(PNIM_var))
fore_best<- as.data.frame(cbind(MR_pred, logit.s_cjs, years))

library(caret)
R2(fore_best$fore_PNI, fore_best$logit.s_cjs)
#ggplot(MR_pred, aes(x=years, y=fore.mean)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey") + geom_line() 
cjs_PNIM<- ggplot(fore_best, aes(x=years, y= plogis(logit.s_cjs))) + geom_point(color="blue") + geom_line(color="blue") + geom_ribbon(aes(ymin=plogis(lower), ymax=plogis(upper)), fill="grey", alpha=0.5)  + 
  geom_point(aes(x=years, y=plogis(fore_PNIM)), color="red") + geom_line(aes(x=years, y=plogis(fore_PNIM)), color="red") +  ylab("SAR Survival") +xlab("Year")

ggsave("output/cjs_PNIM.png", cjs_PNIM)
##################################################
## Marine trend + outflow + river temperature

m=4

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2", "q.beta3")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- marine1_z  ## Nx1; predictor variable
Z[1, 3, ] <- outflow_z  ## Nx1; predictor variable
Z[1, 4, ] <- temp_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_Mcov <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.Mcov <- MARSSkfas(dlm_Mcov)

#Get forecasted predictions
eta.Mcov <- kf.Mcov$xtt1

## ts of E(forecasts)
fore_Mcov <- vector()
for (t in 1:TT) {
  fore_Mcov[t] <- Z[, , t] %*% eta.Mcov[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.Mcov <- kf.Mcov$Vtt1

## obs variance; 1x1 matrix
R_Mcov <- coef(dlm_Mcov, type = "matrix")$R

## ts of Var(forecasts)
Mcov_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  Mcov_var[t] <- Z[, , t] %*% Phi.Mcov[, , t] %*% tZ + R_Mcov
}

sig.lik.Mcov<- sum(log(abs(Mcov_var)))
innov.lik.Mcov<- (1/2)*sum(((logit.s_cjs-fore_Mcov)^2)/Mcov_var)
pie<- (TT/2)*log(2*pi)

Mcov_LL<- pie + sig.lik.Mcov +innov.lik.Mcov

#What to do about first obs??
MAPE_Mcov<- mean(abs((fore_Mcov-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_Mcov<- rmse(logit.s_cjs, fore_Mcov)

Href_cjs_Mcov<- 2*(constant_LL - Mcov_LL)

AICcMcov<- 2*Mcov_LL + 2*dlm_Mcov$num.params + ((2*dlm_Mcov$num.params*(dlm_Mcov$num.params+1))/(TT-dlm_Mcov$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine_FW", Model="Constant", Data="CJS", MAPE=MAPE_Mcov, RMSE=rmse_Mcov, 
                     Href=Href_cjs_Mcov, LL=Mcov_LL, AICc=AICcMcov)

#################################
## Marine trend + outflow + river temperature + Pacific Northwest Index
m=5

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2", "qbeta3", "qbeta4")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- marine1_z  ## Nx1; predictor variable
Z[1, 3, ] <- PNI_z  ## Nx1; predictor variable
Z[1, 4, ] <- outflow_z  ## Nx1; predictor variable
Z[1, 5, ] <- temp_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix
#R <- matrix(0.151)

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_FWM <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.FWM <- MARSSkfas(dlm_FWM)

#Get forecasted predictions
eta.FWM <- kf.FWM$xtt1

## ts of E(forecasts)
fore_FWM <- vector()
for (t in 1:TT) {
  fore_FWM[t] <- Z[, , t] %*% eta.FWM[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.FWM <- kf.FWM$Vtt1

## obs variance; 1x1 matrix
R_FWM <- coef(dlm_FWM, type = "matrix")$R

## ts of Var(forecasts)
FWM_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  FWM_var[t] <- Z[, , t] %*% Phi.FWM[, , t] %*% tZ + R_FWM
}

sig.lik.FWM<- sum(log(abs(FWM_var)))
innov.lik.FWM<- (1/2)*sum(((logit.s_cjs-fore_FWM)^2)/FWM_var)
pie<- (TT/2)*log(2*pi)


FWM_LL<- pie + sig.lik.FWM +innov.lik.FWM

#What to do about first obs??
MAPE_FWM<- mean(abs((fore_FWM-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_FWM<- rmse(logit.s_cjs, fore_FWM)

Href_cjs_FWM<- 2*(constant_LL - FWM_LL)

AICcFWM<- 2*FWM_LL + 2*dlm_FWM$num.params + ((2*dlm_FWM$num.params*(dlm_FWM$num.params+1))/(TT-dlm_FWM$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="FW_Marine_PNI", Model="Constant", Data="CJS", MAPE=MAPE_FWM, RMSE=rmse_FWM, 
                     Href=Href_cjs_FWM, LL=FWM_LL, AICc=AICcFWM)

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
dlm_lintrend <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
MAPE_ltconstant<- mean(abs((fore_lintrend-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltconstant<- rmse(logit.s_cjs, fore_lintrend)

Href_cjs_lt<- 2*(lintrend_LL - lintrend_LL)

AICclintrend<- 2*lintrend_LL + 2*dlm_lintrend$num.params + ((2*dlm_lintrend$num.params*(dlm_lintrend$num.params+1))/(TT-dlm_lintrend$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate=NA, Model="LinTrend", Data="CJS", MAPE=MAPE_ltconstant, RMSE=rmse_ltconstant, 
                     Href=Href_cjs_lt, LL=lintrend_LL, AICc=AICclintrend)
###########################################
#linear trend CUI
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
dlm_lt_CUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltCUI<- (1/2)*sum(((logit.s_cjs-fore_lt_CUI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltCUI_LL<- pie + sig.lik.ltCUI +innov.lik.ltCUI

#What to do about first obs??
MAPE_ltCUI<- mean(abs((fore_lt_CUI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltCUI<- rmse(logit.s_cjs, fore_lt_CUI)

Href_cjs_ltCUI<- 2*(lintrend_LL - ltCUI_LL)

AICcltCUI<- 2*ltCUI_LL + 2*dlm_lt_CUI$num.params + ((2*dlm_lt_CUI$num.params*(dlm_lt_CUI$num.params+1))/(TT-dlm_lt_CUI$num.params-1))


mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="LinTrend", Data="CJS", MAPE=MAPE_ltCUI, RMSE=rmse_ltCUI, 
                     Href=Href_cjs_ltCUI, LL=ltCUI_LL, AICc=AICcltCUI)

###########################################
#linear trend CUTI
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
Z[1, 3, ] <- CUTI_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_CUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.CUTI <- MARSSkfas(dlm_lt_CUTI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.CUTI <- kf.lt.CUTI$xtt1

## ts of E(forecasts)
fore_lt_CUTI <- vector()
for (t in 1:TT) {
  fore_lt_CUTI[t] <- Z[, , t] %*% eta.lt.CUTI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.CUTI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_CUTI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltCUTI<- sum(log(abs(fore_var)))
innov.lik.ltCUTI<- (1/2)*sum(((logit.s_cjs-fore_lt_CUTI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltCUTI_LL<- pie + sig.lik.ltCUTI +innov.lik.ltCUTI

#What to do about first obs??
MAPE_ltCUTI<- mean(abs((fore_lt_CUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltCUTI<- rmse(logit.s_cjs, fore_lt_CUTI)

Href_cjs_ltCUTI<- 2*(lintrend_LL - ltCUTI_LL)

AICcltCUTI<- 2*ltCUTI_LL + 2*dlm_lt_CUTI$num.params + ((2*dlm_lt_CUTI$num.params*(dlm_lt_CUTI$num.params+1))/(TT-dlm_lt_CUTI$num.params-1))


mod_sel_df<- add_row(mod_sel_df, Covariate="CUTI", Model="LinTrend", Data="CJS", MAPE=MAPE_ltCUTI, RMSE=rmse_ltCUTI, 
                     Href=Href_cjs_ltCUTI, LL=ltCUTI_LL, AICc=AICcltCUTI)
###########################################
#linear trend BEUTI
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
Z[1, 3, ] <- BEUTI_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_BEUTI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.BEUTI <- MARSSkfas(dlm_lt_BEUTI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.BEUTI <- kf.lt.BEUTI$xtt1

## ts of E(forecasts)
fore_lt_BEUTI <- vector()
for (t in 1:TT) {
  fore_lt_BEUTI[t] <- Z[, , t] %*% eta.lt.BEUTI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.BEUTI$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_BEUTI, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltBEUTI<- sum(log(abs(fore_var)))
innov.lik.ltBEUTI<- (1/2)*sum(((logit.s_cjs-fore_lt_BEUTI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltBEUTI_LL<- pie + sig.lik.ltBEUTI +innov.lik.ltBEUTI

#What to do about first obs??
MAPE_ltBEUTI<- mean(abs((fore_lt_BEUTI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltBEUTI<- rmse(logit.s_cjs, fore_lt_BEUTI)

Href_cjs_ltBEUTI<- 2*(lintrend_LL - ltBEUTI_LL)

AICcltBEUTI<- 2*ltBEUTI_LL + 2*dlm_lt_BEUTI$num.params + ((2*dlm_lt_BEUTI$num.params*(dlm_lt_BEUTI$num.params+1))/(TT-dlm_lt_BEUTI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="BEUTI", Model="LinTrend", Data="CJS", MAPE=MAPE_ltBEUTI, RMSE=rmse_ltBEUTI, 
                     Href=Href_cjs_ltBEUTI, LL=ltBEUTI_LL, AICc=AICcltBEUTI)
###########################################
#linear trend MLDNO3
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
Z[1, 3, ] <- mldno3_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_mldno3 <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.mldno3 <- MARSSkfss(dlm_lt_mldno3)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.mldno3 <- kf.lt.mldno3$xtt1

## ts of E(forecasts)
fore_lt_mldno3 <- vector()
for (t in 1:TT) {
  fore_lt_mldno3[t] <- Z[, , t] %*% eta.lt.mldno3[, t, drop = FALSE]
}

sig.lik.ltmldno3<- sum(log(abs(kf.lt.mldno3$Sigma[1,1,])))
innov.lik.ltmldno3<- (1/2)*sum((kf.lt.mldno3$Innov^2)/kf.lt.mldno3$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

ltmldno3_LL<- pie + sig.lik.ltmldno3 +innov.lik.ltmldno3

#What to do about first obs??
MAPE_ltmldno3<- mean(abs((fore_lt_mldno3-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltmldno3<- rmse(logit.s_cjs, fore_lt_mldno3)

Href_cjs_ltmldno3<- 2*(lintrend_LL - ltmldno3_LL)

AICcltmldno3<- 2*ltmldno3_LL + 2*dlm_lt_mldno3$num.params + ((2*dlm_lt_mldno3$num.params*(dlm_lt_mldno3$num.params+1))/(TT-dlm_lt_mldno3$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="MLDNO3", Model="LinTrend", Data="CJS", MAPE=MAPE_ltmldno3, RMSE=rmse_ltmldno3, 
                     Href=Href_cjs_ltmldno3, LL=ltmldno3_LL, AICc=AICcltmldno3)
###########################################
#linear trend outflow
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
Z[1, 3, ] <- outflow_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_outflow <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.outflow <- MARSSkfss(dlm_lt_outflow)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.outflow <- kf.lt.outflow$xtt1

## ts of E(forecasts)
fore_lt_outflow <- vector()
for (t in 1:TT) {
  fore_lt_outflow[t] <- Z[, , t] %*% eta.lt.outflow[, t, drop = FALSE]
}

sig.lik.ltoutflow<- sum(log(abs(kf.lt.outflow$Sigma[1,1,])))
innov.lik.ltoutflow<- (1/2)*sum((kf.lt.outflow$Innov^2)/kf.lt.outflow$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

ltoutflow_LL<- pie + sig.lik.ltoutflow +innov.lik.ltoutflow

#What to do about first obs??
MAPE_ltoutflow<- mean(abs((fore_lt_outflow-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltoutflow<- rmse(logit.s_cjs, fore_lt_outflow)

Href_cjs_ltoutflow<- 2*(lintrend_LL - ltoutflow_LL)

AICcltoutflow<- 2*ltoutflow_LL + 2*dlm_lt_outflow$num.params + ((2*dlm_lt_outflow$num.params*(dlm_lt_outflow$num.params+1))/(TT-dlm_lt_outflow$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Outflow", Model="LinTrend", Data="CJS", MAPE=MAPE_ltoutflow, RMSE=rmse_ltoutflow, 
                     Href=Href_cjs_ltoutflow, LL=ltoutflow_LL, AICc=AICcltoutflow)
###########################################
#linear trend river temperature
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
Z[1, 3, ] <- temp_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_temp <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.temp <- MARSSkfss(dlm_lt_temp)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.temp <- kf.lt.temp$xtt1

## ts of E(forecasts)
fore_lt_temp <- vector()
for (t in 1:TT) {
  fore_lt_temp[t] <- Z[, , t] %*% eta.lt.temp[, t, drop = FALSE]
}

sig.lik.lttemp<- sum(log(abs(kf.lt.temp$Sigma[1,1,])))
innov.lik.lttemp<- (1/2)*sum((kf.lt.temp$Innov^2)/kf.lt.temp$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

lttemp_LL<- pie + sig.lik.lttemp +innov.lik.lttemp

#What to do about first obs??
MAPE_lttemp<- mean(abs((fore_lt_temp-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_lttemp<- rmse(logit.s_cjs, fore_lt_temp)

Href_cjs_lttemp<- 2*(lintrend_LL - lttemp_LL)

AICclttemp<- 2*lttemp_LL + 2*dlm_lt_temp$num.params + ((2*dlm_lt_temp$num.params*(dlm_lt_temp$num.params+1))/(TT-dlm_lt_temp$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Temp", Model="LinTrend", Data="CJS", MAPE=MAPE_lttemp, RMSE=rmse_lttemp, 
                     Href=Href_cjs_lttemp, LL=lttemp_LL, AICc=AICclttemp)
#################################
#linear trend Pacific Northwest Index
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
dlm_lt_PNI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltPNI<- (1/2)*sum(((logit.s_cjs-fore_lt_PNI)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltPNI_LL<- pie + sig.lik.ltPNI +innov.lik.ltPNI

#What to do about first obs??
MAPE_ltPNI<- mean(abs((fore_lt_PNI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltPNI<- rmse(logit.s_cjs, fore_lt_PNI)

Href_cjs_ltPNI<- 2*(lintrend_LL - ltPNI_LL)

AICcltPNI<- 2*ltPNI_LL + 2*dlm_lt_PNI$num.params + ((2*dlm_lt_PNI$num.params*(dlm_lt_PNI$num.params+1))/(TT-dlm_lt_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="LinTrend", Data="CJS", MAPE=MAPE_ltPNI, RMSE=rmse_ltPNI, 
                     Href=Href_cjs_ltPNI, LL=ltPNI_LL, AICc=AICcltPNI)

#################################
#linear trend SWE peak day
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
Z[1, 3, ] <- peak_day_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_pkday <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.pkday <- MARSSkfss(dlm_lt_pkday)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.pkday <- kf.lt.pkday$xtt1

## ts of E(forecasts)
fore_lt_pkday <- vector()
for (t in 1:TT) {
  fore_lt_pkday[t] <- Z[, , t] %*% eta.lt.pkday[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.pkday$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_pkday, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}
sig.lik.ltpkday<- sum(log(abs(kf.lt.pkday$Sigma[1,1,])))
innov.lik.ltpkday<- (1/2)*sum((kf.lt.pkday$Innov^2)/kf.lt.pkday$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

ltpkday_LL<- pie + sig.lik.ltpkday +innov.lik.ltpkday

#What to do about first obs??
MAPE_ltpkday<- mean(abs((fore_lt_pkday-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltpkday<- rmse(logit.s_cjs, fore_lt_pkday)

Href_cjs_ltpkday<- 2*(lintrend_LL - ltpkday_LL)

AICcltpkday<- 2*ltpkday_LL + 2*dlm_lt_pkday$num.params + ((2*dlm_lt_pkday$num.params*(dlm_lt_pkday$num.params+1))/(TT-dlm_lt_pkday$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Peak day", Model="LinTrend", Data="CJS", MAPE=MAPE_ltpkday, RMSE=rmse_ltpkday, 
                     Href=Href_cjs_ltpkday, LL=ltpkday_LL, AICc=AICcltpkday)
#################################
#linear trend SWE
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
Z[1, 3, ] <- peak_swe_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_swe <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.swe <- MARSSkfas(dlm_lt_swe)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.swe <- kf.lt.swe$xtt1

## ts of E(forecasts)
fore_lt_swe <- vector()
for (t in 1:TT) {
  fore_lt_swe[t] <- Z[, , t] %*% eta.lt.swe[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.swe$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_swe, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}

#This is the proper way to calculate likelihood!
sig.lik.ltswe<- sum(log(abs(fore_var)))
innov.lik.ltswe<- (1/2)*sum(((logit.s_cjs-fore_lt_swe)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltswe_LL<- pie + sig.lik.ltswe +innov.lik.ltswe

#What to do about first obs??
MAPE_ltswe<- mean(abs((fore_lt_swe-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltswe<- rmse(logit.s_cjs, fore_lt_swe)

Href_cjs_ltswe<- 2*(lintrend_LL - ltswe_LL)

AICcltswe<- 2*ltswe_LL + 2*dlm_lt_swe$num.params + ((2*dlm_lt_swe$num.params*(dlm_lt_swe$num.params+1))/(TT-dlm_lt_swe$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Peak SWE", Model="LinTrend", Data="CJS", MAPE=MAPE_ltswe, RMSE=rmse_ltswe, 
                     Href=Href_cjs_ltswe, LL=ltswe_LL, AICc=AICcltswe)

#################################
#linear trend Spill
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
Z[1, 3, ] <- spill_z
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_lt_spill <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.spill <- MARSSkfss(dlm_lt_spill)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.spill <- kf.lt.spill$xtt1

## ts of E(forecasts)
fore_lt_spill <- vector()
for (t in 1:TT) {
  fore_lt_spill[t] <- Z[, , t] %*% eta.lt.spill[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.lt.spill$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_lt_spill, type = "matrix")$R

## ts of Var(forecasts)
fore_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
}
sig.lik.ltspill<- sum(log(abs(kf.lt.spill$Sigma[1,1,])))
innov.lik.ltspill<- (1/2)*sum((kf.lt.spill$Innov^2)/kf.lt.spill$Sigma[1,1,])
pie<- (TT/2)*log(2*pi)

ltspill_LL<- pie + sig.lik.ltspill +innov.lik.ltspill

#What to do about first obs??
MAPE_ltspill<- mean(abs((fore_lt_spill-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltspill<- rmse(logit.s_cjs, fore_lt_spill)

Href_cjs_ltspill<- 2*(lintrend_LL - ltspill_LL)

AICcltspill<- 2*ltspill_LL + 2*dlm_lt_spill$num.params + ((2*dlm_lt_spill$num.params*(dlm_lt_spill$num.params+1))/(TT-dlm_lt_spill$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Spill", Model="LinTrend", Data="CJS", MAPE=MAPE_ltspill, RMSE=rmse_ltspill, 
                     Href=Href_cjs_ltspill, LL=ltspill_LL, AICc=AICcltspill)

###########################################
#linear trend Seas Surface Temperature pre-ocean entry (Dec., Jan., Feb.)
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
dlm_lt_SSTpre <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltSSTpre<- (1/2)*sum(((logit.s_cjs-fore_lt_SSTpre)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSSTpre_LL<- pie + sig.lik.ltSSTpre +innov.lik.ltSSTpre

#What to do about first obs??
MAPE_ltSSTpre<- mean(abs((fore_lt_SSTpre-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltSSTpre<- rmse(logit.s_cjs, fore_lt_SSTpre)

Href_cjs_ltSSTpre<- 2*(lintrend_LL - ltSSTpre_LL)

AICcltSSTpre<- 2*ltSSTpre_LL + 2*dlm_lt_SSTpre$num.params + ((2*dlm_lt_SSTpre$num.params*(dlm_lt_SSTpre$num.params+1))/(TT-dlm_lt_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_pre", Model="LinTrend", Data="CJS", MAPE=MAPE_ltSSTpre, RMSE=rmse_ltSSTpre, 
                     Href=Href_cjs_ltSSTpre, LL=ltSSTpre_LL, AICc=AICcltSSTpre)

###########################################
#linear trend Sea Surface Temperature at ocean entry (May, June, July)
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
dlm_lt_SSTentry <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltSSTentry<- (1/2)*sum(((logit.s_cjs-fore_lt_SSTpre)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSSTentry_LL<- pie + sig.lik.ltSSTentry +innov.lik.ltSSTentry

#What to do about first obs??
MAPE_ltSSTentry<- mean(abs((fore_lt_SSTentry-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltSSTentry<- rmse(logit.s_cjs, fore_lt_SSTentry)

Href_cjs_ltSSTentry<- 2*(lintrend_LL - ltSSTentry_LL)

AICcltSSTentry<- 2*ltSSTentry_LL + 2*dlm_lt_SSTentry$num.params + ((2*dlm_lt_SSTentry$num.params*(dlm_lt_SSTentry$num.params+1))/(TT-dlm_lt_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_entry", Model="LinTrend", Data="CJS", MAPE=MAPE_ltSSTentry, RMSE=rmse_ltSSTentry, 
                     Href=Href_cjs_ltSSTentry, LL=ltSSTentry_LL, AICc=AICcltSSTentry)

###########################################
#linear trend Sea surface temperature 6 month average (Dec, Jan, Feb, March, April, May)
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
dlm_lt_SST6mo <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltSST6mo<- (1/2)*sum(((logit.s_cjs-fore_lt_SST6mo)^2)/fore_var)
pie<- (TT/2)*log(2*pi)

ltSST6mo_LL<- pie + sig.lik.ltSST6mo +innov.lik.ltSST6mo

#What to do about first obs??
MAPE_ltSST6mo<- mean(abs((fore_lt_SST6mo-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltSST6mo<- rmse(logit.s_cjs, fore_lt_SST6mo)

Href_cjs_ltSST6mo<- 2*(lintrend_LL - ltSST6mo_LL)


AICcltSST6mo<- 2*ltSST6mo_LL + 2*dlm_lt_SST6mo$num.params + ((2*dlm_lt_SST6mo$num.params*(dlm_lt_SST6mo$num.params+1))/(TT-dlm_lt_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_6mo", Model="LinTrend", Data="CJS", MAPE=MAPE_ltSST6mo, RMSE=rmse_ltSST6mo, 
                     Href=Href_cjs_ltSST6mo, LL=ltSST6mo_LL, AICc=AICcltSST6mo)

###########################################
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
dlm_lt_marine <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltmarine<- (1/2)*sum(((logit.s_cjs-fore_lt_marine)^2)/fore_marine)
pie<- (TT/2)*log(2*pi)

ltmarine_LL<- pie + sig.lik.ltmarine +innov.lik.ltmarine

#What to do about first obs??
MAPE_ltmarine<- mean(abs((fore_lt_marine-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltmarine<- rmse(logit.s_cjs, fore_lt_marine)

Href_cjs_ltmarine<- 2*(lintrend_LL - ltmarine_LL)


AICcltmarine<- 2*ltmarine_LL + 2*dlm_lt_marine$num.params + ((2*dlm_lt_marine$num.params*(dlm_lt_marine$num.params+1))/(TT-dlm_lt_marine$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine", Model="LinTrend", Data="CJS", MAPE=MAPE_ltmarine, RMSE=rmse_ltmarine, 
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
dlm_lt_marinePNI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

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
innov.lik.ltmarinePNI<- (1/2)*sum(((logit.s_cjs-fore_lt_marinePNI)^2)/fore_marinePNI)
pie<- (TT/2)*log(2*pi)

ltmarinePNI_LL<- pie + sig.lik.ltmarinePNI +innov.lik.ltmarinePNI

#What to do about first obs??
MAPE_ltmarinePNI<- mean(abs((fore_lt_marinePNI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltmarinePNI<- rmse(logit.s_cjs, fore_lt_marinePNI)

Href_cjs_ltmarinePNI<- 2*(lintrend_LL - ltmarinePNI_LL)


AICcltmarinePNI<- 2*ltmarinePNI_LL + 2*dlm_lt_marinePNI$num.params + ((2*dlm_lt_marinePNI$num.params*(dlm_lt_marinePNI$num.params+1))/(TT-dlm_lt_marinePNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine + PNI", Model="LinTrend", Data="CJS", MAPE=MAPE_ltmarinePNI, RMSE=rmse_ltmarinePNI, 
                     Href=Href_cjs_ltmarinePNI, LL=ltmarinePNI_LL, AICc=AICcltmarinePNI)

##################################################
## linear trend Marine trend + outflow + river temperature

m=5

## for process eqn
B <-  matrix(c(1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), nrow=5, ncol=5)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.beta1", "q.beta2", "q.beta3") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- marine1_z  ## Nx1; predictor variable
Z[1, 4, ] <- outflow_z  ## Nx1; predictor variable
Z[1, 5, ] <- temp_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_ltMFW <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.ltMFW <- MARSSkfas(dlm_ltMFW)

#Get forecasted predictions
eta.ltMFW <- kf.ltMFW$xtt1

## ts of E(forecasts)
fore_ltMFW <- vector()
for (t in 1:TT) {
  fore_ltMFW[t] <- Z[, , t] %*% eta.ltMFW[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.ltMFW <- kf.ltMFW$Vtt1

## obs variance; 1x1 matrix
R_ltMFW <- coef(dlm_ltMFW, type = "matrix")$R

## ts of Var(forecasts)
ltMFW_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  ltMFW_var[t] <- Z[, , t] %*% Phi.ltMFW[, , t] %*% tZ + R_ltMFW
}

sig.lik.ltMFW<- sum(log(abs(ltMFW_var)))
innov.lik.ltMFW<- (1/2)*sum(((logit.s_cjs-fore_ltMFW)^2)/ltMFW_var)
pie<- (TT/2)*log(2*pi)

ltMFW_LL<- pie + sig.lik.ltMFW +innov.lik.ltMFW

#What to do about first obs??
MAPE_ltMFW<- mean(abs((fore_ltMFW-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_ltMFW<- rmse(logit.s_cjs, fore_ltMFW)

Href_cjs_ltMFW<- 2*(constant_LL - ltMFW_LL)

AICcltMFW<- 2*ltMFW_LL + 2*dlm_ltMFW$num.params + ((2*dlm_ltMFW$num.params*(dlm_ltMFW$num.params+1))/(TT-dlm_ltMFW$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine + outflow + temp", Model="LinTrend", Data="CJS", MAPE=MAPE_ltMFW, RMSE=rmse_ltMFW, 
                     Href=Href_cjs_ltMFW, LL=ltMFW_LL, AICc=AICcltMFW)

#################################
## Marine trend + outflow + river temperature + Pacific Northwest Index
m=6

## for process eqn
B <-  matrix(c(1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), nrow=6, ncol=6)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.beta1", "q.beta2", "q.beta3", "q.beta4") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1, 2, ] <- rep(0, TT)
Z[1, 3, ] <- marine1_z  ## Nx1; predictor variable
Z[1, 4, ] <- PNI_z  ## Nx1; predictor variable
Z[1, 5, ] <- outflow_z  ## Nx1; predictor variable
Z[1, 6, ] <- temp_z  ## Nx1; predictor variable
A <- matrix(0)  ## 1x1; scalar = 0; 
R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)

## fit univariate DLM
dlm_FWMPNI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.FWMPNI <- MARSSkfas(dlm_FWMPNI)

#Get forecasted predictions
eta.FWMPNI <- kf.FWMPNI$xtt1

## ts of E(forecasts)
fore_FWMPNI <- vector()
for (t in 1:TT) {
  fore_FWMPNI[t] <- Z[, , t] %*% eta.FWMPNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi.FWMPNI <- kf.FWMPNI$Vtt1

## obs variance; 1x1 matrix
R_FWMPNI <- coef(dlm_FWMPNI, type = "matrix")$R

## ts of Var(forecasts)
FWMPNI_var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
  FWMPNI_var[t] <- Z[, , t] %*% Phi.FWMPNI[, , t] %*% tZ + R_FWMPNI
}

sig.lik.FWMPNI<- sum(log(abs(FWMPNI_var)))
innov.lik.FWMPNI<- (1/2)*sum(((logit.s_cjs-fore_FWMPNI)^2)/FWMPNI_var)
pie<- (TT/2)*log(2*pi)


FWMPNI_LL<- pie + sig.lik.FWMPNI +innov.lik.FWMPNI

#What to do about first obs??
MAPE_FWMPNI<- mean(abs((fore_FWMPNI-logit.s_cjs)/logit.s_cjs))
#Root mean square error
rmse_FWMPNI<- rmse(logit.s_cjs, fore_FWMPNI)

Href_cjs_FWMPNI<- 2*(constant_LL - FWMPNI_LL)

AICcFWMPNI<- 2*FWM_LL + 2*dlm_FWMPNI$num.params + ((2*dlm_FWMPNI$num.params*(dlm_FWMPNI$num.params+1))/(TT-dlm_FWMPNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="FW_Marine_PNI", Model="LinearTrend", Data="CJS", MAPE=MAPE_FWMPNI, RMSE=rmse_FWMPNI, 
                     Href=Href_cjs_FWMPNI, LL=FWMPNI_LL, AICc=AICcFWMPNI)

mod_selection<- mod_sel_df %>% arrange(AICc) %>% select(-MAPE)
write.csv(mod_selection, "Output/CJS_Univariate_sorted.csv")

write.csv(mod_sel_df, "Output/CJS_Univariate_DLM.csv")
