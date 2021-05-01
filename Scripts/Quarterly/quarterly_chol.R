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

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Quarterly.xlsx", sheet=4)
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
data_chol <- dat_unseas[,c("oil_USD_qoq", 
                           "imp_price_qoq", 
                           "miacr_31",
                           #"reserves_USD_qoq",
                           "neer_qoq",
                           "gdp_nominal_qoq", 
                           "cpi_all_qoq")]
#data_chol <- dat_unseas[,c("oil_USD_qoq", "imp_price_qoq", "reserves_USD_qoq", 
#                           "miacr_31", "neer_qoq", "d_real_gdp_gap", "cpi_all_qoq")]
data_chol$neer_qoq <- data_chol$neer_qoq * -1
data_chol$real_usd_qoq <- data_chol$real_usd_qoq * -1
data_chol$nom_usd_qoq <- data_chol$nom_usd_qoq * -1
VARselect(data_chol, lag.max=6)$selection
model_VAR <- VAR(data_chol, p = 1, type = "const")
choldec <- id.chol(model_VAR)
irf_choldec <- irf(choldec, n.ahead = 4, ortho=TRUE)
#plot(irf_choldec)
#REER gives closer results
#Reserves
sum(irf_choldec$irf$`epsilon[ reserves_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ reserves_USD_qoq ] %->% neer_qoq`)
#Oil
sum(irf_choldec$irf$`epsilon[ oil_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ oil_USD_qoq ] %->% neer_qoq`)
#MIACR 31
sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% neer_qoq`)
#NEER
sum(irf_choldec$irf$`epsilon[ neer_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ neer_qoq ] %->% neer_qoq`)
#Output
sum(irf_choldec$irf$`epsilon[ gdp_nominal_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ gdp_nominal_qoq ] %->% neer_qoq`)

#Cholesky decomposition (dollars as exchange rate)
data_chol_usd <- dat_unseas[,c("oil_USD_qoq", "imp_price_qoq", "reserves_USD_qoq",
                           "miacr_31", "nom_usd_qoq", "gdp_real_SA_qoq", "cpi_all_qoq")]
data_chol_usd$nom_usd_qoq <- data_chol_usd$nom_usd_qoq * -1
VARselect(data_chol_usd, lag.max = 4)$selection
model_VAR_usd <- VAR(data_chol_usd, p = 4, type = "const")
choldec_usd <- id.chol(model_VAR_usd)
irf_choldec_usd <- irf(choldec_usd, n.ahead = 4, ortho=TRUE)
plot(irf_choldec_usd)
#Reserves
sum(irf_choldec_usd$irf$`epsilon[ reserves_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec_usd$irf$`epsilon[ reserves_USD_qoq ] %->% nom_usd_qoq`)
#Oil
sum(irf_choldec_usd$irf$`epsilon[ oil_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec_usd$irf$`epsilon[ oil_USD_qoq ] %->% nom_usd_qoq`)
#MIACR 31
sum(irf_choldec_usd$irf$`epsilon[ miacr_31 ] %->% cpi_all_qoq`) /
  sum(irf_choldec_usd$irf$`epsilon[ miacr_31 ] %->% nom_usd_qoq`)
#NEER
sum(irf_choldec_usd$irf$`epsilon[ nom_usd_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec_usd$irf$`epsilon[ nom_usd_qoq ] %->% nom_usd_qoq`)
#Output
sum(irf_choldec_usd$irf$`epsilon[ gdp_nominal_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec_usd$irf$`epsilon[ gdp_nominal_qoq ] %->% nom_usd_qoq`)


#Smooth transition
model_cv <- id.st(model_VAR, c_lower=0.1, c_upper=0.9, c_step=3, nc=8)
irf_cv <- irf(model_cv, n.ahead = 4, ortho=TRUE)
#plot(irf_cv)
#Reserves
sum(irf_choldec$irf$`epsilon[ reserves_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ reserves_USD_qoq ] %->% neer_qoq`)
#Oil
sum(irf_choldec$irf$`epsilon[ oil_USD_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ oil_USD_qoq ] %->% neer_qoq`)
#MIACR 31
sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% neer_qoq`)
#NEER
sum(irf_choldec$irf$`epsilon[ neer_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ neer_qoq ] %->% neer_qoq`)
#Output
sum(irf_choldec$irf$`epsilon[ gdp_nominal_qoq ] %->% cpi_all_qoq`) /
  sum(irf_choldec$irf$`epsilon[ gdp_nominal_qoq ] %->% neer_qoq`)


