rm(list = ls())

library(dlm)
library(dplyr)
library(here)
library(Metrics)
library(MARSS)
library(ggplot2)

setwd("~/CBR/DLM_Analysis_Git")

#Read in data
cjs_raw<- readRDS("data/df_SAR_LGRLGA_MV_old.rds") 

#Plot the data
data_SAR<- ggplot(cjs_raw, aes(x=year, y=plogis(SW_logit_phi))) + geom_point() + geom_line() + 
  geom_point(aes(x=year, y=plogis(SAR_logit_phi)), color="blue") +geom_line(aes(x=year, y=plogis(SAR_logit_phi)), color="blue") +
  geom_point(aes(x=year, y=plogis(cjs_logit_phi)), color="purple") +geom_line(aes(x=year, y=plogis(cjs_logit_phi)), color="purple") + ylab("SAR survival") + xlab("Year")

ggsave("output/data_SAR.png", data_SAR)

marine1<- readRDS("marine1_trend.rds")

swepeak<- read.csv("data/swepeak.csv") %>% select(year, day, sum, mean)

mean_spill<- readRDS("data/mean_spill.rds") %>% mutate(year=as.double(year))

cjs_raw$M1<- marine1$estimate

cjs<- cjs_raw %>% left_join(swepeak, by="year") %>% left_join(mean_spill, by="year") %>% filter(year>=1979)


#Remove NA at end of time series
SAR<- t(as.matrix(cjs %>% dplyr::select(SW_logit_phi, SAR_logit_phi, cjs_logit_phi)))
SAR

TT<- dim(SAR)[2]
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
n=3
## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(1, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
#Z[1:n, 1:n, ] <- diag(1, n)  ## Nx1; 1's for intercept
A <- matrix(0)  ## 1x1; scalar = 0; 
#R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = c(0))

## list of model matrices & vectors
mod_list <- list(U = "zero", Q = Q, Z = Z, A = "zero", R="diagonal and unequal")

## fit univariate DLM
dlm_constant <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.out.constant <- MARSSkfss(dlm_constant)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.constant <- kf.out.constant$xtt1

## ts of E(forecasts)
fore_constant <- array(NA, dim = c(n, m, TT))
for (t in 1:TT) {
  fore_constant[, ,t] <- Z[, , t] %*% eta.constant[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
Phi <- kf.out.constant$Vtt1

## obs variance; 1x1 matrix
R_est <- coef(dlm_constant, type = "matrix")$R

ts<- n
## ts of Var(forecasts)
constant_var <- matrix(nrow=n, ncol=TT)
for (t in 1:TT) {
  for (s in 1:ts){
    tZ <- t(Z[, ,t])  ## transpose of Z
    constant_var[s, t] <- Z[s, ,t] %*% Phi[, , t] %*% tZ[,s] + R_est[s]
  }
}

#This is the proper way to calculate likelihood!
sig.lik<- sum(log(abs(kf.out.constant$Sigma[1,1,]))) +  sum(log(abs(kf.out.constant$Sigma[2,2,]))) +  sum(log(abs(kf.out.constant$Sigma[3,3,])))
innov.lik<- (1/2)*(sum((kf.out.constant$Innov[1,]^2)/kf.out.constant$Sigma[1,1,]) + sum((kf.out.constant$Innov[2,]^2)/kf.out.constant$Sigma[2,2,])+ sum((kf.out.constant$Innov[3,]^2)/kf.out.constant$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

constant_LL<- pie + sig.lik +innov.lik

fore_c<- matrix(fore_constant, nrow=n, ncol=TT)

#What to do about first obs??
MAPE_constant<- mean(abs((fore_c-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_constant<- sqrt(mean((SAR - fore_c)^2, na.rm = TRUE))

AICcConstant<- 2*constant_LL + 2*dlm_constant$num.params + ((2*dlm_constant$num.params*(dlm_constant$num.params+1))/(TT-dlm_constant$num.params-1))


mod_sel_df<- data.frame(Covariate=NA, Model="Constant", Data="Multi", MAPE=MAPE_constant, RMSE=rmse_constant, Href=0, LL=constant_LL, AICc=AICcConstant)


##########################
#Run CUI model
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- CUI_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 
#R<- matrix(0, nrow=3, ncol=3)
#diag(R) <- c("r", "r","r")  ## 1x1; scalar = r; this is the variance-covariance matrix

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R="diagonal and unequal")

## fit univariate DLM
dlm_CUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.CUI <- MARSSkfss(dlm_CUI)

#Get forecast predictions
eta.CUI <- kf.CUI$xtt1

## ts of E(forecasts)
fore_CUI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_CUI[, t] <- Z[, , t] %*% eta.CUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.CUI <- kf.CUI$Vtt1

## obs variance; 1x1 matrix
#R_CUI <- diag(coef(dlm_CUI, type = "matrix")$R)

## ts of Var(forecasts)
#ts<- n
## ts of Var(forecasts)
#CUI_var <- matrix(nrow=n, ncol=TT)
#for (t in 1:TT) {
#  for (s in 1:ts){
#    tZ <- t(Z[, ,t])  ## transpose of Z
#    CUI_var[s, t] <- Z[s, ,t] %*% Phi.CUI[, , t] %*% tZ[,s] + R_CUI[s]
#  }
#}

sig.lik.CUI<- sum(log(abs(kf.CUI$Sigma[1,1,]))) +  sum(log(abs(kf.CUI$Sigma[2,2,]))) +  sum(log(abs(kf.CUI$Sigma[3,3,])))
innov.lik.CUI<- (1/2)*(sum((kf.CUI$Innov[1,]^2)/kf.CUI$Sigma[1,1,]) + sum((kf.CUI$Innov[2,]^2)/kf.CUI$Sigma[2,2,])+ sum((kf.CUI$Innov[3,]^2)/kf.CUI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

#sig.lik.CUI<- sum(log(abs(CUI_var[n,n])))
#innov.lik.CUI<- (1/2)*sum(((logit.s_cjs-fore_CUI)^2)/CUI_var)
#pie<- (TT/2)*log(2*pi)

CUI_LL<- pie + sig.lik.CUI +innov.lik.CUI

Href_cjs_CUI<- 2*(constant_LL - CUI_LL)

#What to do about first obs??
MAPE_CUI<- mean(abs((fore_CUI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_CUI<-  sqrt(mean((SAR - fore_CUI)^2, na.rm = TRUE))

AICcCUI<- 2*CUI_LL + 2*dlm_CUI$num.params + ((2*dlm_CUI$num.params*(dlm_CUI$num.params+1))/(TT-dlm_CUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="Constant", Data="Multi", MAPE=MAPE_CUI, RMSE=rmse_CUI, 
                     Href=Href_cjs_CUI, LL=CUI_LL, AICc=AICcCUI)

###########################################
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- CUI_mean_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 
#R <- matrix("r")  ## 1x1; scalar = r; this is the variance-covariance matrix


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_mnCUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.mnCUI <- MARSSkfss(dlm_mnCUI)

#sig.lik.mnCUI<- sum(log(abs(kf.mnCUI$Sigma[1,1,])))
#innov.lik.mnCUI<- (1/2)*sum((kf.mnCUI$Innov^2)/kf.mnCUI$Sigma[1,1,])
#pie<- (TT/2)*log(2*pi)

#Get forecast predictions
eta.mnCUI <- kf.mnCUI$xtt1

## ts of E(forecasts)
fore_mnCUI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_mnCUI[,t] <- Z[, , t] %*% eta.mnCUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.mnCUI <- kf.mnCUI$Vtt1

## obs variance; 1x1 matrix
#R_mnCUI <- coef(dlm_mnCUI, type = "matrix")$R

## ts of Var(forecasts)
#mnCUI_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  mnCUI_var[t] <- Z[, , t] %*% Phi.mnCUI[, , t] %*% tZ + R_mnCUI
#}

sig.lik.mnCUI<- sum(log(abs(kf.mnCUI$Sigma[1,1,]))) +  sum(log(abs(kf.mnCUI$Sigma[2,2,]))) +  sum(log(abs(kf.mnCUI$Sigma[3,3,])))
innov.lik.mnCUI<- (1/2)*(sum((kf.mnCUI$Innov[1,]^2)/kf.mnCUI$Sigma[1,1,]) + sum((kf.mnCUI$Innov[2,]^2)/kf.mnCUI$Sigma[2,2,])+ sum((kf.mnCUI$Innov[3,]^2)/kf.mnCUI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

mnCUI_LL<- pie + sig.lik.mnCUI +innov.lik.mnCUI

Href_cjs_mnCUI<- 2*(constant_LL - mnCUI_LL)

#What to do about first obs??
MAPE_mnCUI<- mean(abs((fore_mnCUI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_mnCUI<- sqrt(mean((SAR - fore_mnCUI)^2, na.rm = TRUE))

AICcmnCUI<- 2*mnCUI_LL + 2*dlm_mnCUI$num.params + ((2*dlm_mnCUI$num.params*(dlm_mnCUI$num.params+1))/(TT-dlm_mnCUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="mnCUI", Model="Constant", Data="Multi", MAPE=MAPE_mnCUI, RMSE=rmse_mnCUI, 
                     Href=Href_cjs_mnCUI, LL=mnCUI_LL, AICc=AICcmnCUI)

#################################
## Pacific Northwest Index
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0;   ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_PNI<- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.PNI <- MARSSkfss(dlm_PNI)

eta.PNI <- kf.PNI$xtt1

## ts of E(forecasts)
fore_PNI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_PNI[,t] <- Z[, , t] %*% eta.PNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.PNI <- kf.PNI$Vtt1

## obs variance; 1x1 matrix
#R_PNI <- coef(dlm_PNI, type = "matrix")$R

## ts of Var(forecasts)
#PNI_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  PNI_var[t] <- Z[, , t] %*% Phi.PNI[, , t] %*% tZ + R_PNI
#}

sig.lik.PNI<- sum(log(abs(kf.PNI$Sigma[1,1,]))) +  sum(log(abs(kf.PNI$Sigma[2,2,]))) +  sum(log(abs(kf.PNI$Sigma[3,3,])))
innov.lik.PNI<- (1/2)*(sum((kf.PNI$Innov[1,]^2)/kf.PNI$Sigma[1,1,]) + sum((kf.PNI$Innov[2,]^2)/kf.PNI$Sigma[2,2,])+ sum((kf.PNI$Innov[3,]^2)/kf.PNI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

PNI_LL<- pie + sig.lik.PNI +innov.lik.PNI

Href_cjs_PNI<- 2*(constant_LL - PNI_LL)

#What to do about first obs??
MAPE_PNI<- mean(abs((fore_PNI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_PNI<- sqrt(mean((SAR - fore_mnCUI)^2, na.rm = TRUE))

AICcPNI<- 2*PNI_LL + 2*dlm_PNI$num.params + ((2*dlm_PNI$num.params*(dlm_PNI$num.params+1))/(TT-dlm_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="Constant", Data="Multi", MAPE=MAPE_PNI, RMSE=rmse_PNI, 
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
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- SST_pre_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0;   ## 1x1; scalar = 0; 

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_SSTpre <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.SSTpre <- MARSSkfss(dlm_SSTpre)

sig.lik.SSTpre<- sum(log(abs(kf.SSTpre$Sigma[1,1,]))) +  sum(log(abs(kf.SSTpre$Sigma[2,2,]))) +  sum(log(abs(kf.SSTpre$Sigma[3,3,])))
innov.lik.SSTpre<- (1/2)*(sum((kf.SSTpre$Innov[1,]^2)/kf.SSTpre$Sigma[1,1,]) + sum((kf.SSTpre$Innov[2,]^2)/kf.SSTpre$Sigma[2,2,])+ sum((kf.SSTpre$Innov[3,]^2)/kf.SSTpre$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

SSTpre_LL<- pie + sig.lik.SSTpre +innov.lik.SSTpre

#Get forecasted predictions
eta.SSTpre <- kf.SSTpre$xtt1

## ts of E(forecasts)
fore_SSTpre <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_SSTpre[,t] <- Z[, , t] %*% eta.SSTpre[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_SSTpre<- mean(abs((fore_SSTpre-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_SSTpre<- sqrt(mean((SAR - fore_SSTpre)^2, na.rm = TRUE))

Href_cjs_SSTpre<- 2*(constant_LL - SSTpre_LL)

AICcSSTpre<- 2*SSTpre_LL + 2*dlm_SSTpre$num.params + ((2*dlm_SSTpre$num.params*(dlm_SSTpre$num.params+1))/(TT-dlm_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTpre", Model="Constant", Data="Multi", MAPE=MAPE_SSTpre, RMSE=rmse_SSTpre, 
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
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- SST_entry_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)   ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_SSTentry <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.SSTentry <- MARSSkfss(dlm_SSTentry)

#Get forecasted predictions
eta.SSTentry <- kf.SSTentry$xtt1

## ts of E(forecasts)
fore_SSTentry <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_SSTentry[, t] <- Z[, , t] %*% eta.SSTentry[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.SSTentry <- kf.SSTentry$Vtt1

## obs variance; 1x1 matrix
#R_SSTentry <- coef(dlm_SSTentry, type = "matrix")$R

## ts of Var(forecasts)
#SSTentry_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  SSTentry_var[t] <- Z[, , t] %*% Phi.SSTentry[, , t] %*% tZ + R_SSTentry
#}

sig.lik.SSTentry<- sum(log(abs(kf.SSTentry$Sigma[1,1,]))) +  sum(log(abs(kf.SSTentry$Sigma[2,2,]))) +  sum(log(abs(kf.SSTentry$Sigma[3,3,])))
innov.lik.SSTentry<- (1/2)*(sum((kf.SSTentry$Innov[1,]^2)/kf.SSTentry$Sigma[1,1,]) + sum((kf.SSTentry$Innov[2,]^2)/kf.SSTentry$Sigma[2,2,])+ sum((kf.SSTentry$Innov[3,]^2)/kf.SSTentry$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

SSTentry_LL<- pie + sig.lik.SSTentry +innov.lik.SSTentry

#What to do about first obs??
MAPE_SSTentry<- mean(abs((fore_SSTentry-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_SSTentry<- sqrt(mean((SAR - fore_SSTentry)^2, na.rm = TRUE))

Href_cjs_SSTentry<- 2*(constant_LL - SSTentry_LL)

AICcSSTentry<- 2*SSTentry_LL + 2*dlm_SSTentry$num.params + ((2*dlm_SSTentry$num.params*(dlm_SSTentry$num.params+1))/(TT-dlm_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SSTentry", Model="Constant", Data="Multi", MAPE=MAPE_SSTentry, RMSE=rmse_SSTentry, 
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
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- SST_6mo_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_SST6mo <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.SST6mo <- MARSSkfss(dlm_SST6mo)

sig.lik.SST6mo<- sum(log(abs(kf.SST6mo$Sigma[1,1,]))) +  sum(log(abs(kf.SST6mo$Sigma[2,2,]))) +  sum(log(abs(kf.SST6mo$Sigma[3,3,])))
innov.lik.SST6mo<- (1/2)*(sum((kf.SST6mo$Innov[1,]^2)/kf.SST6mo$Sigma[1,1,]) + sum((kf.SST6mo$Innov[2,]^2)/kf.SST6mo$Sigma[2,2,])+ sum((kf.SST6mo$Innov[3,]^2)/kf.SST6mo$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

SST6mo_LL<- pie + sig.lik.SST6mo +innov.lik.SST6mo

#Get forecasted predictions
eta.SST6mo <- kf.SST6mo$xtt1

## ts of E(forecasts)
fore_SST6mo <-  matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_SST6mo[,t] <- Z[, , t] %*% eta.SST6mo[, t, drop = FALSE]
}

#What to do about first obs??
MAPE_SST6mo<- mean(abs((fore_SST6mo-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_SST6mo<- sqrt(mean((SAR - fore_SST6mo)^2, na.rm = TRUE))

Href_cjs_SST6mo<- 2*(constant_LL - SST6mo_LL)

AICcSST6mo<- 2*SST6mo_LL + 2*dlm_SST6mo$num.params + ((2*dlm_SST6mo$num.params*(dlm_SST6mo$num.params+1))/(TT-dlm_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST6mo", Model="Constant", Data="Multi", MAPE=MAPE_SST6mo, RMSE=rmse_SST6mo, 
                     Href=Href_cjs_SST6mo, LL=SST6mo_LL, AICc=AICcSST6mo)

years<- cjs_raw$year
#MR_pred<- cbind(fore_SST6mo, upper=fore_SST6mo+2*sqrt(kf.SST6mo$Sigma), lower=fore_SST6mo-2*sqrt(kf.SST6mo$Sigma))
#fore_best<- as.data.frame(cbind(t(fore_SST6mo), t(SAR), years))
#ggplot(MR_pred, aes(x=years, y=fore.mean)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey") + geom_line() 
#cjs_PNIM<- ggplot(fore_best, aes(x=years, y= plogis(logit.s_cjs))) + geom_point(color="blue") + geom_line(color="blue") + geom_ribbon(aes(ymin=plogis(lower), ymax=plogis(upper)), fill="grey", alpha=0.5)  + 
#  geom_point(aes(x=years, y=plogis(fore_PNIM)), color="red") + geom_line(aes(x=years, y=plogis(fore_PNIM)), color="red") +  ylab("SAR Survival") +xlab("Year")

#ggsave("output/Multi_SST6mo.png", cjs_PNIM)

#################################
## Percent spill
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- spill_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_spill <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.spill <- MARSSkfss(dlm_spill)

#Get forecasted predictions
eta.spill <- kf.spill$xtt1

## ts of E(forecasts)
fore_spill <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_spill[, t] <- Z[, , t] %*% eta.spill[, t, drop = FALSE]
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

sig.lik.spill<- sum(log(abs(kf.spill$Sigma[1,1,]))) +  sum(log(abs(kf.spill$Sigma[2,2,]))) +  sum(log(abs(kf.spill$Sigma[3,3,])))
innov.lik.spill<- (1/2)*(sum((kf.spill$Innov[1,]^2)/kf.spill$Sigma[1,1,]) + sum((kf.spill$Innov[2,]^2)/kf.spill$Sigma[2,2,])+ sum((kf.spill$Innov[3,]^2)/kf.spill$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

spill_LL<- pie + sig.lik.spill +innov.lik.spill

#What to do about first obs??
MAPE_spill<- mean(abs((fore_spill-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_spill<- sqrt(mean((SAR - fore_spill)^2, na.rm = TRUE))

Href_cjs_spill<- 2*(constant_LL - spill_LL)

AICcspill<- 2*spill_LL + 2*dlm_spill$num.params + ((2*dlm_spill$num.params*(dlm_spill$num.params+1))/(TT-dlm_spill$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="spill", Model="Constant", Data="Multi", MAPE=MAPE_spill, RMSE=rmse_spill, 
                     Href=Href_cjs_spill, LL=spill_LL, AICc=AICcspill)

#################################
## Percent outflow
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- outflow_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_outflow <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.outflow <- MARSSkfss(dlm_outflow)

#Get forecasted predictions
eta.outflow <- kf.outflow$xtt1

## ts of E(forecasts)
fore_outflow <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_outflow[, t] <- Z[, , t] %*% eta.outflow[, t, drop = FALSE]
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

sig.lik.outflow<- sum(log(abs(kf.outflow$Sigma[1,1,]))) +  sum(log(abs(kf.outflow$Sigma[2,2,]))) +  sum(log(abs(kf.outflow$Sigma[3,3,])))
innov.lik.outflow<- (1/2)*(sum((kf.outflow$Innov[1,]^2)/kf.outflow$Sigma[1,1,]) + sum((kf.outflow$Innov[2,]^2)/kf.outflow$Sigma[2,2,])+ sum((kf.outflow$Innov[3,]^2)/kf.outflow$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

outflow_LL<- pie + sig.lik.outflow +innov.lik.outflow

#What to do about first obs??
MAPE_outflow<- mean(abs((fore_outflow-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_outflow<- sqrt(mean((SAR - fore_outflow)^2, na.rm = TRUE))

Href_cjs_outflow<- 2*(constant_LL - outflow_LL)

AICcoutflow<- 2*outflow_LL + 2*dlm_outflow$num.params + ((2*dlm_outflow$num.params*(dlm_outflow$num.params+1))/(TT-dlm_outflow$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="outflow", Model="Constant", Data="Multi", MAPE=MAPE_outflow, RMSE=rmse_outflow, 
                     Href=Href_cjs_outflow, LL=outflow_LL, AICc=AICcoutflow)


#################################
## Marine Trend
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- marine1_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_marine1 <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.marine1 <- MARSSkfss(dlm_marine1)

#Get forecasted predictions
eta.marine1 <- kf.marine1$xtt1

## ts of E(forecasts)
fore_marine1 <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_marine1[, t] <- Z[, , t] %*% eta.marine1[, t, drop = FALSE]
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

sig.lik.marine1<- sum(log(abs(kf.marine1$Sigma[1,1,]))) +  sum(log(abs(kf.marine1$Sigma[2,2,]))) +  sum(log(abs(kf.marine1$Sigma[3,3,])))
innov.lik.marine1<- (1/2)*(sum((kf.marine1$Innov[1,]^2)/kf.marine1$Sigma[1,1,]) + sum((kf.marine1$Innov[2,]^2)/kf.marine1$Sigma[2,2,])+ sum((kf.marine1$Innov[3,]^2)/kf.marine1$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

marine1_LL<- pie + sig.lik.marine1 +innov.lik.marine1

#What to do about first obs??
MAPE_marine1<- mean(abs((fore_marine1-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_marine1<- sqrt(mean((SAR - fore_marine1)^2, na.rm = TRUE))

Href_cjs_marine1<- 2*(constant_LL - marine1_LL)

AICcmarine1<- 2*marine1_LL + 2*dlm_marine1$num.params + ((2*dlm_marine1$num.params*(dlm_marine1$num.params+1))/(TT-dlm_marine1$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="marine1", Model="Constant", Data="Multi", MAPE=MAPE_marine1, RMSE=rmse_marine1, 
                     Href=Href_cjs_marine1, LL=marine1_LL, AICc=AICcmarine1)

#################################
## SWE peak day
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- peak_day_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_peak_day <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.pkday <- MARSSkfss(dlm_peak_day)

#Get forecasted predictions
eta.pkday <- kf.pkday$xtt1

## ts of E(forecasts)
fore_pkday <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_pkday[, t] <- Z[, , t] %*% eta.pkday[, t, drop = FALSE]
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

sig.lik.pkday<- sum(log(abs(kf.pkday$Sigma[1,1,]))) +  sum(log(abs(kf.pkday$Sigma[2,2,]))) +  sum(log(abs(kf.pkday$Sigma[3,3,])))
innov.lik.pkday<- (1/2)*(sum((kf.pkday$Innov[1,]^2)/kf.pkday$Sigma[1,1,]) + sum((kf.pkday$Innov[2,]^2)/kf.pkday$Sigma[2,2,])+ sum((kf.pkday$Innov[3,]^2)/kf.pkday$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

pkday_LL<- pie + sig.lik.pkday +innov.lik.pkday

  #What to do about first obs??
MAPE_pkday<- mean(abs((fore_pkday-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_pkday<- sqrt(mean((SAR - fore_pkday)^2, na.rm = TRUE))

Href_cjs_pkday<- 2*(constant_LL - pkday_LL)

AICcpkday<- 2*pkday_LL + 2*dlm_peak_day$num.params + ((2*dlm_peak_day$num.params*(dlm_peak_day$num.params+1))/(TT-dlm_peak_day$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="peak day", Model="Constant", Data="Multi", MAPE=MAPE_pkday, RMSE=rmse_pkday, 
                     Href=Href_cjs_pkday, LL=pkday_LL, AICc=AICcpkday)

#################################
## SWE peak day
m=2

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- peak_swe_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_peak <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.pk <- MARSSkfss(dlm_peak)

#Get forecasted predictions
eta.pk <- kf.pk$xtt1

## ts of E(forecasts)
fore_pk <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_pk[, t] <- Z[, , t] %*% eta.pk[, t, drop = FALSE]
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

sig.lik.pk<- sum(log(abs(kf.pk$Sigma[1,1,]))) +  sum(log(abs(kf.pk$Sigma[2,2,]))) +  sum(log(abs(kf.pk$Sigma[3,3,])))
innov.lik.pk<- (1/2)*(sum((kf.pk$Innov[1,]^2)/kf.pk$Sigma[1,1,]) + sum((kf.pk$Innov[2,]^2)/kf.pk$Sigma[2,2,])+ sum((kf.pk$Innov[3,]^2)/kf.pk$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

pk_LL<- pie + sig.lik.pk +innov.lik.pk

#What to do about first obs??
MAPE_pk<- mean(abs((fore_pk-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_pk<- sqrt(mean((SAR - fore_pk)^2, na.rm = TRUE))

Href_cjs_pk<- 2*(constant_LL - pk_LL)

AICcpk<- 2*pk_LL + 2*dlm_peak$num.params + ((2*dlm_peak$num.params*(dlm_peak$num.params+1))/(TT-dlm_peak$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="peak", Model="Constant", Data="Multi", MAPE=MAPE_pk, RMSE=rmse_pk, 
                     Href=Href_cjs_pk, LL=pk_LL, AICc=AICcpk)


#################################
## Marine trend + Pacific Northwest Index
m=3

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- marine1_z  ## Nx1; predictor variable
Z[1:n, 3, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_PNIM <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.PNIM <- MARSSkfss(dlm_PNIM)

#Get forecasted predictions
eta.PNIM <- kf.PNIM$xtt1

## ts of E(forecasts)
fore_PNIM <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_PNIM[, t] <- Z[, , t] %*% eta.PNIM[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.PNIM <- kf.PNIM$Vtt1

## obs variance; 1x1 matrix
#R_PNIM <- coef(dlm_PNIM, type = "matrix")$R

## ts of Var(forecasts)
#PNIM_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  PNIM_var[t] <- Z[, , t] %*% Phi.PNIM[, , t] %*% tZ + R_PNIM
#}

sig.lik.PNIM<- sum(log(abs(kf.PNIM$Sigma[1,1,]))) +  sum(log(abs(kf.PNIM$Sigma[2,2,]))) +  sum(log(abs(kf.PNIM$Sigma[3,3,])))
innov.lik.PNIM<- (1/2)*(sum((kf.PNIM$Innov[1,]^2)/kf.PNIM$Sigma[1,1,]) + sum((kf.PNIM$Innov[2,]^2)/kf.PNIM$Sigma[2,2,])+ sum((kf.PNIM$Innov[3,]^2)/kf.PNIM$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)


PNIM_LL<- pie + sig.lik.PNIM +innov.lik.PNIM

#What to do about first obs??
MAPE_PNIM<- mean(abs((fore_PNIM-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_PNIM<- sqrt(mean((SAR - fore_PNIM)^2, na.rm = TRUE))

Href_cjs_PNIM<- 2*(constant_LL - PNIM_LL)

AICcPNIM<- 2*PNIM_LL + 2*dlm_PNIM$num.params + ((2*dlm_PNIM$num.params*(dlm_PNIM$num.params+1))/(TT-dlm_PNIM$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI_Marine", Model="Constant", Data="Multi", MAPE=MAPE_PNIM, RMSE=rmse_PNIM, 
                     Href=Href_cjs_PNIM, LL=PNIM_LL, AICc=AICcPNIM)

#################################
## CUI + Pacific Northwest Index
m=3

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- CUI_mean_z  ## Nx1; predictor variable
Z[1:n, 3, ] <- PNI_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_PNICUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.PNICUI <- MARSSkfss(dlm_PNICUI)

#Get forecasted predictions
eta.PNICUI <- kf.PNICUI$xtt1

## ts of E(forecasts)
fore_PNICUI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_PNICUI[, t] <- Z[, , t] %*% eta.PNICUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.PNIM <- kf.PNIM$Vtt1

## obs variance; 1x1 matrix
#R_PNIM <- coef(dlm_PNIM, type = "matrix")$R

## ts of Var(forecasts)
#PNIM_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  PNIM_var[t] <- Z[, , t] %*% Phi.PNIM[, , t] %*% tZ + R_PNIM
#}

sig.lik.PNICUI<- sum(log(abs(kf.PNICUI$Sigma[1,1,]))) +  sum(log(abs(kf.PNICUI$Sigma[2,2,]))) +  sum(log(abs(kf.PNICUI$Sigma[3,3,])))
innov.lik.PNICUI<- (1/2)*(sum((kf.PNICUI$Innov[1,]^2)/kf.PNICUI$Sigma[1,1,]) + sum((kf.PNICUI$Innov[2,]^2)/kf.PNICUI$Sigma[2,2,])+ sum((kf.PNICUI$Innov[3,]^2)/kf.PNICUI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)


PNICUI_LL<- pie + sig.lik.PNICUI +innov.lik.PNICUI

#What to do about first obs??
MAPE_PNICUI<- mean(abs((fore_PNICUI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_PNICUI<- sqrt(mean((SAR - fore_PNICUI)^2, na.rm = TRUE))

Href_cjs_PNICUI<- 2*(constant_LL - PNICUI_LL)

AICcPNICUI<- 2*PNICUI_LL + 2*dlm_PNICUI$num.params + ((2*dlm_PNICUI$num.params*(dlm_PNICUI$num.params+1))/(TT-dlm_PNICUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI_CUI", Model="Constant", Data="Multi", MAPE=MAPE_PNICUI, RMSE=rmse_PNICUI, 
                     Href=Href_cjs_PNICUI, LL=PNICUI_LL, AICc=AICcPNICUI)

####################################
## CUI + Pacific Northwest Index
m=3

## for process eqn
B <- diag(m)  ## 2x2; Identity; this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.beta1", "q.beta2")#, "q.beta3")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- CUI_mean_z  ## Nx1; predictor variable
Z[1:n, 3, ] <- outflow_z  ## Nx1; predictor variable
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_outCUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.outCUI <- MARSSkfss(dlm_outCUI)

#Get forecasted predictions
eta.outCUI <- kf.outCUI$xtt1

## ts of E(forecasts)
fore_outCUI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_outCUI[, t] <- Z[, , t] %*% eta.outCUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi.outM <- kf.outM$Vtt1

## obs variance; 1x1 matrix
#R_outM <- coef(dlm_outM, type = "matrix")$R

## ts of Var(forecasts)
#outM_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  outM_var[t] <- Z[, , t] %*% Phi.outM[, , t] %*% tZ + R_outM
#}

sig.lik.outCUI<- sum(log(abs(kf.outCUI$Sigma[1,1,]))) +  sum(log(abs(kf.outCUI$Sigma[2,2,]))) +  sum(log(abs(kf.outCUI$Sigma[3,3,])))
innov.lik.outCUI<- (1/2)*(sum((kf.outCUI$Innov[1,]^2)/kf.outCUI$Sigma[1,1,]) + sum((kf.outCUI$Innov[2,]^2)/kf.outCUI$Sigma[2,2,])+ sum((kf.outCUI$Innov[3,]^2)/kf.outCUI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)


outCUI_LL<- pie + sig.lik.outCUI +innov.lik.outCUI

#What to do about first obs??
MAPE_outCUI<- mean(abs((fore_outCUI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_outCUI<- sqrt(mean((SAR - fore_outCUI)^2, na.rm = TRUE))

Href_cjs_outCUI<- 2*(constant_LL - outCUI_LL)

AICcoutCUI<- 2*outCUI_LL + 2*dlm_outCUI$num.params + ((2*dlm_outCUI$num.params*(dlm_outCUI$num.params+1))/(TT-dlm_outCUI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="out_CUI", Model="Constant", Data="Multi", MAPE=MAPE_outCUI, RMSE=rmse_outCUI, 
                     Href=Href_cjs_outCUI, LL=outCUI_LL, AICc=AICcoutCUI)


###################################################################################
#linear trend
m=2

## for process eqn
B <-  matrix(c(1,0,1,1), nrow=2, ncol=2) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))  ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 2, ] <- rep(0, TT)
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lintrend <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lintrend <- MARSSkfss(dlm_lintrend)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lintrend <- kf.lintrend$xtt1

## ts of E(forecasts)
fore_lintrend <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lintrend[, t] <- Z[, , t] %*% eta.lintrend[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lintrend$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lintrend, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}

#This is the proper way to calculate likelihood!
sig.lik.lintrend<- sum(log(abs(kf.lintrend$Sigma[1,1,]))) +  sum(log(abs(kf.lintrend$Sigma[2,2,]))) +  sum(log(abs(kf.lintrend$Sigma[3,3,])))
innov.lik.lintrend<- (1/2)*(sum((kf.lintrend$Innov[1,]^2)/kf.lintrend$Sigma[1,1,]) + sum((kf.lintrend$Innov[2,]^2)/kf.lintrend$Sigma[2,2,])+ sum((kf.lintrend$Innov[3,]^2)/kf.lintrend$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

lintrend_LL<- pie + sig.lik.lintrend +innov.lik.lintrend

#What to do about first obs??
MAPE_ltconstant<- mean(abs((fore_lintrend-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltconstant<- sqrt(mean((SAR - fore_lintrend)^2, na.rm = TRUE))

Href_cjs_lt<- 2*(lintrend_LL - lintrend_LL)

AICclintrend<- 2*lintrend_LL + 2*dlm_lintrend$num.params + ((2*dlm_lintrend$num.params*(dlm_lintrend$num.params+1))/(TT-dlm_lintrend$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate=NA, Model="LinTrend", Data="Multi", MAPE=MAPE_ltconstant, RMSE=rmse_ltconstant, 
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
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- CUI_z
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0,0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_CUI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.CUI <- MARSSkfss(dlm_lt_CUI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.CUI <- kf.lt.CUI$xtt1

## ts of E(forecasts)
fore_lt_CUI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_CUI[,t] <- Z[, , t] %*% eta.lt.CUI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lt.CUI$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lt_CUI, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}

#This is the proper way to calculate likelihood!
sig.lik.ltCUI<- sum(log(abs(kf.lt.CUI$Sigma[1,1,]))) +  sum(log(abs(kf.lt.CUI$Sigma[2,2,]))) +  sum(log(abs(kf.lt.CUI$Sigma[3,3,])))
innov.lik.ltCUI<- (1/2)*(sum((kf.lt.CUI$Innov[1,]^2)/kf.lt.CUI$Sigma[1,1,]) + sum((kf.lt.CUI$Innov[2,]^2)/kf.lt.CUI$Sigma[2,2,])+ sum((kf.lt.CUI$Innov[3,]^2)/kf.lt.CUI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltCUI_LL<- pie + sig.lik.ltCUI +innov.lik.ltCUI

#What to do about first obs??
MAPE_ltCUI<- mean(abs((fore_lt_CUI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltCUI<- sqrt(mean((SAR - fore_lt_CUI)^2, na.rm = TRUE))

Href_cjs_ltCUI<- 2*(lintrend_LL - ltCUI_LL)

AICcltCUI<- 2*ltCUI_LL + 2*dlm_lt_CUI$num.params + ((2*dlm_lt_CUI$num.params*(dlm_lt_CUI$num.params+1))/(TT-dlm_lt_CUI$num.params-1))


mod_sel_df<- add_row(mod_sel_df, Covariate="CUI", Model="LinTrend", Data="Multi", MAPE=MAPE_ltCUI, RMSE=rmse_ltCUI, 
                     Href=Href_cjs_ltCUI, LL=ltCUI_LL, AICc=AICcltCUI)


#################################
#linear trend Pacific Northwest Index
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- PNI_z
A <- matrix(c(0,0,0), ncol=1)   ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_PNI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.PNI <- MARSSkfss(dlm_lt_PNI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.PNI <- kf.lt.PNI$xtt1

## ts of E(forecasts)
fore_lt_PNI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_PNI[,t] <- Z[, , t] %*% eta.lt.PNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lt.PNI$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lt_PNI, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}

#This is the proper way to calculate likelihood!
sig.lik.ltPNI<- sum(log(abs(kf.lt.PNI$Sigma[1,1,]))) +  sum(log(abs(kf.lt.PNI$Sigma[2,2,]))) +  sum(log(abs(kf.lt.PNI$Sigma[3,3,])))
innov.lik.ltPNI<- (1/2)*(sum((kf.lt.PNI$Innov[1,]^2)/kf.lt.PNI$Sigma[1,1,]) + sum((kf.lt.PNI$Innov[2,]^2)/kf.lt.PNI$Sigma[2,2,])+ sum((kf.lt.PNI$Innov[3,]^2)/kf.lt.PNI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltPNI_LL<- pie + sig.lik.ltPNI +innov.lik.ltPNI

#What to do about first obs??
MAPE_ltPNI<- mean(abs((fore_lt_PNI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltPNI<- sqrt(mean((SAR - fore_lt_PNI)^2, na.rm = TRUE))

Href_cjs_ltPNI<- 2*(lintrend_LL - ltPNI_LL)

AICcltPNI<- 2*ltPNI_LL + 2*dlm_lt_PNI$num.params + ((2*dlm_lt_PNI$num.params*(dlm_lt_PNI$num.params+1))/(TT-dlm_lt_PNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="PNI", Model="LinTrend", Data="Multi", MAPE=MAPE_ltPNI, RMSE=rmse_ltPNI, 
                     Href=Href_cjs_ltPNI, LL=ltPNI_LL, AICc=AICcltPNI)

###########################################
#linear trend Seas Surface Temperature pre-ocean entry (Dec., Jan., Feb.)
m=3

## for process eqn
B <-  matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 1), nrow=3, ncol=3) # this is the parameter evolution matrix (G)
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0; 
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now; this is the variance-covariance matrix
diag(Q) <- c("q.alpha", "q.slope", "q.cov") ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- SST_pre_z
A <- matrix(c(0,0,0), ncol=1)   ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_SSTpre <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SSTpre <- MARSSkfss(dlm_lt_SSTpre)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SSTpre <- kf.lt.SSTpre$xtt1

## ts of E(forecasts)
fore_lt_SSTpre <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_SSTpre[,t] <- Z[, , t] %*% eta.lt.SSTpre[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lt.SSTpre$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lt_SSTpre, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}

#This is the proper way to calculate likelihood!
sig.lik.ltSSTpre<- sum(log(abs(kf.lt.SSTpre$Sigma[1,1,]))) +  sum(log(abs(kf.lt.SSTpre$Sigma[2,2,]))) +  sum(log(abs(kf.lt.SSTpre$Sigma[3,3,])))
innov.lik.ltSSTpre<- (1/2)*(sum((kf.lt.SSTpre$Innov[1,]^2)/kf.lt.SSTpre$Sigma[1,1,]) + sum((kf.lt.SSTpre$Innov[2,]^2)/kf.lt.SSTpre$Sigma[2,2,])+ sum((kf.lt.SSTpre$Innov[3,]^2)/kf.lt.SSTpre$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltSSTpre_LL<- pie + sig.lik.ltSSTpre +innov.lik.ltSSTpre

#What to do about first obs??
MAPE_ltSSTpre<- mean(abs((fore_lt_SSTpre-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltSSTpre<-sqrt(mean((SAR - fore_lt_SSTpre)^2, na.rm = TRUE))

Href_cjs_ltSSTpre<- 2*(lintrend_LL - ltSSTpre_LL)

AICcltSSTpre<- 2*ltSSTpre_LL + 2*dlm_lt_SSTpre$num.params + ((2*dlm_lt_SSTpre$num.params*(dlm_lt_SSTpre$num.params+1))/(TT-dlm_lt_SSTpre$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_pre", Model="LinTrend", Data="Multi", MAPE=MAPE_ltSSTpre, RMSE=rmse_ltSSTpre, 
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
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- SST_entry_z
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_SSTentry <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SSTentry <- MARSSkfss(dlm_lt_SSTentry)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SSTentry <- kf.lt.SSTentry$xtt1

## ts of E(forecasts)
fore_lt_SSTentry <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_SSTentry[, t] <- Z[, , t] %*% eta.lt.SSTentry[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lt.SSTentry$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lt_SSTentry, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}


#This is the proper way to calculate likelihood!
sig.lik.ltSSTentry<- sum(log(abs(kf.lt.SSTentry$Sigma[1,1,]))) +  sum(log(abs(kf.lt.SSTentry$Sigma[2,2,]))) +  sum(log(abs(kf.lt.SSTentry$Sigma[3,3,])))
innov.lik.ltSSTentry<- (1/2)*(sum((kf.lt.SSTentry$Innov[1,]^2)/kf.lt.SSTentry$Sigma[1,1,]) + sum((kf.lt.SSTentry$Innov[2,]^2)/kf.lt.SSTentry$Sigma[2,2,])+ sum((kf.lt.SSTentry$Innov[3,]^2)/kf.lt.SSTentry$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltSSTentry_LL<- pie + sig.lik.ltSSTentry +innov.lik.ltSSTentry

#What to do about first obs??
MAPE_ltSSTentry<- mean(abs((fore_lt_SSTentry-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltSSTentry<- sqrt(mean((SAR - fore_lt_SSTentry)^2, na.rm = TRUE))

Href_cjs_ltSSTentry<- 2*(lintrend_LL - ltSSTentry_LL)

AICcltSSTentry<- 2*ltSSTentry_LL + 2*dlm_lt_SSTentry$num.params + ((2*dlm_lt_SSTentry$num.params*(dlm_lt_SSTentry$num.params+1))/(TT-dlm_lt_SSTentry$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_entry", Model="LinTrend", Data="Multi", MAPE=MAPE_ltSSTentry, RMSE=rmse_ltSSTentry, 
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
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- SST_6mo_z
A <- matrix(c(0,0,0), ncol=1)  ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_SST6mo <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.SST6mo <- MARSSkfss(dlm_lt_SST6mo)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.SST6mo <- kf.lt.SST6mo$xtt1

## ts of E(forecasts)
fore_lt_SST6mo <-  matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_SST6mo[,t] <- Z[, , t] %*% eta.lt.SST6mo[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi <- kf.lt.SST6mo$Vtt1

## obs variance; 1x1 matrix
#R_est <- coef(dlm_lt_SST6mo, type = "matrix")$R

## ts of Var(forecasts)
#fore_var <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R_est
#}


#This is the proper way to calculate likelihood!
sig.lik.ltSST6mo<- sum(log(abs(kf.lt.SST6mo$Sigma[1,1,]))) +  sum(log(abs(kf.lt.SST6mo$Sigma[2,2,]))) +  sum(log(abs(kf.lt.SST6mo$Sigma[3,3,])))
innov.lik.ltSST6mo<- (1/2)*(sum((kf.lt.SST6mo$Innov[1,]^2)/kf.lt.SST6mo$Sigma[1,1,]) + sum((kf.lt.SST6mo$Innov[2,]^2)/kf.lt.SST6mo$Sigma[2,2,])+ sum((kf.lt.SST6mo$Innov[3,]^2)/kf.lt.SST6mo$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltSST6mo_LL<- pie + sig.lik.ltSST6mo +innov.lik.ltSST6mo

#What to do about first obs??
MAPE_ltSST6mo<-  mean(abs((fore_lt_SST6mo-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltSST6mo<- sqrt(mean((SAR - fore_lt_SST6mo)^2, na.rm = TRUE))

Href_cjs_ltSST6mo<- 2*(lintrend_LL - ltSST6mo_LL)


AICcltSST6mo<- 2*ltSST6mo_LL + 2*dlm_lt_SST6mo$num.params + ((2*dlm_lt_SST6mo$num.params*(dlm_lt_SST6mo$num.params+1))/(TT-dlm_lt_SST6mo$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="SST_6mo", Model="LinTrend", Data="Multi", MAPE=MAPE_ltSST6mo, RMSE=rmse_ltSST6mo, 
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
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- marine1_z
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_marine <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.marine <- MARSSkfss(dlm_lt_marine)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.marine <- kf.lt.marine$xtt1

## ts of E(forecasts)
fore_lt_marine <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_marine[, t] <- Z[, , t] %*% eta.lt.marine[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi_marine <- kf.lt.marine$Vtt1

## obs variance; 1x1 matrix
#R_marine <- coef(dlm_lt_marine, type = "matrix")$R

## ts of Var(forecasts)
#fore_marine <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_marine[t] <- Z[, , t] %*% Phi_marine[, , t] %*% tZ + R_marine
#}


#This is the proper way to calculate likelihood!
sig.lik.ltmarine<- sum(log(abs(kf.lt.marine$Sigma[1,1,]))) +  sum(log(abs(kf.lt.marine$Sigma[2,2,]))) +  sum(log(abs(kf.lt.marine$Sigma[3,3,])))
innov.lik.ltmarine<- (1/2)*(sum((kf.lt.marine$Innov[1,]^2)/kf.lt.marine$Sigma[1,1,]) + sum((kf.lt.marine$Innov[2,]^2)/kf.lt.marine$Sigma[2,2,])+ sum((kf.lt.marine$Innov[3,]^2)/kf.lt.marine$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltmarine_LL<- pie + sig.lik.ltmarine +innov.lik.ltmarine

#What to do about first obs??
MAPE_ltmarine<- mean(abs((fore_lt_marine-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltmarine<- sqrt(mean((SAR - fore_lt_marine)^2, na.rm = TRUE))

Href_cjs_ltmarine<- 2*(lintrend_LL - ltmarine_LL)


AICcltmarine<- 2*ltmarine_LL + 2*dlm_lt_marine$num.params + ((2*dlm_lt_marine$num.params*(dlm_lt_marine$num.params+1))/(TT-dlm_lt_marine$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine", Model="LinTrend", Data="Multi", MAPE=MAPE_ltmarine, RMSE=rmse_ltmarine, 
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
Z <- array(NA, c(n, m, TT))   ## NxMxT; empty for now; matrix of predictor parameter coefficients
Z[1:n, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept
Z[1:n, 2, ] <- rep(0, TT)
Z[1:n, 3, ] <- marine1_z
Z[1:n, 4, ] <- PNI_z
A <- matrix(c(0,0,0), ncol=1) ## 1x1; scalar = 0; 


## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0, 0, 0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = "diagonal and unequal")

## fit univariate DLM
dlm_lt_marinePNI <- MARSS(logit.s_cjs, inits = inits_list, model = mod_list, method="TMB")

kf.lt.marinePNI <- MARSSkfss(dlm_lt_marinePNI)


##############From atsa forecast chapter
## forecasts of regr parameters; 2xT matrix
eta.lt.marinePNI <- kf.lt.marinePNI$xtt1

## ts of E(forecasts)
fore_lt_marinePNI <- matrix(NA, nrow=n, ncol=TT)
for (t in 1:TT) {
  fore_lt_marinePNI[, t] <- Z[, , t] %*% eta.lt.marinePNI[, t, drop = FALSE]
}

## variance of regr parameters; 1x2xT array
#Phi_marinePNI <- kf.lt.marinePNI$Vtt1

## obs variance; 1x1 matrix
#R_marinePNI <- coef(dlm_lt_marinePNI, type = "matrix")$R

## ts of Var(forecasts)
#fore_marinePNI <- vector()
#for (t in 1:TT) {
#  tZ <- matrix(Z[, , t], m, 1)  ## transpose of Z
#  fore_marinePNI[t] <- Z[, , t] %*% Phi_marinePNI[, , t] %*% tZ + R_marinePNI
#}


#This is the proper way to calculate likelihood!
sig.lik.ltmarinePNI<- sum(log(abs(kf.lt.marinePNI$Sigma[1,1,]))) +  sum(log(abs(kf.lt.marinePNI$Sigma[2,2,]))) +  sum(log(abs(kf.lt.marinePNI$Sigma[3,3,])))
innov.lik.ltmarinePNI<- (1/2)*(sum((kf.lt.marinePNI$Innov[1,]^2)/kf.lt.marinePNI$Sigma[1,1,]) + sum((kf.lt.marinePNI$Innov[2,]^2)/kf.lt.marinePNI$Sigma[2,2,])+ sum((kf.lt.marinePNI$Innov[3,]^2)/kf.lt.marinePNI$Sigma[3,3,]))
pie<- (TT/2)*log(2*pi)

ltmarinePNI_LL<- pie + sig.lik.ltmarinePNI +innov.lik.ltmarinePNI

#What to do about first obs??
MAPE_ltmarinePNI<- mean(abs((fore_lt_marinePNI-SAR)/SAR), na.rm=TRUE)
#Root mean square error
rmse_ltmarinePNI<- sqrt(mean((SAR - fore_lt_marinePNI)^2, na.rm = TRUE))

Href_cjs_ltmarinePNI<- 2*(lintrend_LL - ltmarinePNI_LL)


AICcltmarinePNI<- 2*ltmarinePNI_LL + 2*dlm_lt_marinePNI$num.params + ((2*dlm_lt_marinePNI$num.params*(dlm_lt_marinePNI$num.params+1))/(TT-dlm_lt_marinePNI$num.params-1))

mod_sel_df<- add_row(mod_sel_df, Covariate="Marine + PNI", Model="LinTrend", Data="Multi", MAPE=MAPE_ltmarinePNI, RMSE=rmse_ltmarinePNI, 
                     Href=Href_cjs_ltmarinePNI, LL=ltmarinePNI_LL, AICc=AICcltmarinePNI)

mod_selection<- mod_sel_df %>% arrange(AICc) %>% select(-MAPE)
write.csv(mod_selection, "Output/Multivariate_sorted.csv")

write.csv(mod_sel_df, "Output/SAR_Multivariate_DLM.csv")
