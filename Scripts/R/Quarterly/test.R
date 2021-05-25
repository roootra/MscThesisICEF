library(seasonal)
library(rio)
library(lubridate)
library(bvarsv)
library(VARsignR)
library(vars)
library(svars)
library(tsDyn)
library(ggplot2)
library(mFilter)
library(BVAR)

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Quarterly.xlsx", sheet=3)
dat <- na.omit(dat[dat$date > as.Date("2005-05-01"),])
dat <- dat[NROW(dat):1,]

#Deseasonalization
to_deseasonalize <- c("gdp_nominal_index", "imp_price_qoq", grep("cpi.*", colnames(dat), value=T))
dat_unseas <- dat
for(colname in to_deseasonalize){
  current_ts <- ts(dat[,colname], start=c(2005, 2), frequency=4)
  current_seas <- seas(current_ts, transform.function="none",
                       #regression.aictest=NULL, outlier=NULL,
                       #automdl=NULL,
                       seats.noadmiss="yes")
  dat_unseas[,colname] <- as.numeric(final(current_seas))
}

#Output gap
dat_unseas$nominal_gdp_gap <- hpfilter(dat_unseas$gdp_nominal_index, freq=1600)$cycle
dat_unseas$real_gdp_gap <- hpfilter(dat_unseas$gdp_real_SA_index, freq=1600)$cycle
#dat_unseas$d_real_gdp_gap <- c(NA, diff(dat_unseas$real_gdp_gap))
#dat_unseas <- na.omit(dat_unseas)
#dat_unseas$miacr_90 <- c(NA, diff(dat_unseas$miacr_90))
#dat_unseas <- na.omit(dat_unseas)
plot(dat_unseas$date, dat_unseas$real_gdp_gap)

#Cholesky decomposition
data_to_model <- dat_unseas[,c(
  "gdp_real_SA_qoq",
  "oil_USD_qoq",
  "miacr_31",
  "neer_qoq",
  "cpi_all_qoq"
)]

library(bvarr)
#chol(bayesian_model$sigma[1,,]) #upper Cholesky
bayesian_model <- bvar(data_to_model, lags=2,
                       priors = bv_priors(),
                       mh = bv_mh(),
                       n_thin=10)
bayesian_model$sigma[1,,]

Phi_bayesian <- function(bayesian_model, nahead){
  bayesian_coefs <- bayesian_model$beta[,-1,]
  nahead <- nahead + 1
  order <- bayesian_model$meta$lags
  Phi_VMA <- array(data = 0, dim=c(bayesian_model$meta$M, 
                                   bayesian_model$meta$M, 
                                   nahead, bayesian_model$meta$n_save))
  for(draw in 1:dim(Phi_VMA)[4]){
    diag(Phi_VMA[,,1,draw]) <- 1
    for(i in 1:(nahead-1)){
      sum_phi <- 0
      for(lag in 1:order){
        if(lag <= i){
          sum_phi <- sum_phi + Phi_VMA[,,i+1-lag,draw] %*% 
            bayesian_coefs[draw,
              (1+(lag-1)*bayesian_model$meta$M):((lag)*bayesian_model$meta$M),]
        }
      }
      Phi_VMA[,,(i+1),draw] <- sum_phi
    }
  }
  return(Phi_VMA)
}
Phi_bayesian(bayesian_model = bayesian_model, nahead=12)
