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

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Monthly.xlsx", sheet=3)
dat <- dat[NROW(dat):1,]

#Deseasonalization
to_deseasonalize <- c("ind_output_yoy", grep("cpi.*", colnames(dat), value=T))
dat_unseas <- dat
for(colname in to_deseasonalize){
  current_ts <- ts(dat[,colname], start=c(2005, 3), frequency=12)
  current_seas <- seas(current_ts, transform.function="none",
                       #regression.aictest=NULL, outlier=NULL,
                       #automdl=NULL,
                       seats.noadmiss="yes")
  dat_unseas[,colname] <- as.numeric(final(current_seas))
}

#Cholesky decomposition
data_chol <- dat_unseas[,c("oil_USD", "imp_price_mom", "reserves_USD_mln",
                           "miacr_31", "reer_mom", "ind_output_yoy", "cpi_all_mom")]
data_chol$reer_mom <- data_chol$reer_mom * -1
VARselect(data_chol)$selection
model_VAR <- VAR(data_chol, p = 1, type = "const")
choldec <- id.chol(model_VAR)
irf_choldec <- irf(choldec, n.ahead = 12, ortho=TRUE)
plot(irf_choldec)
#Reserves
sum(irf_choldec$irf$`epsilon[ reserves_USD_mln ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ reserves_USD_mln ] %->% reer_mom`)
#Oil
sum(irf_choldec$irf$`epsilon[ oil_USD ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ oil_USD ] %->% reer_mom`)
#MIACR 31
sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% reer_mom`)
#REER
sum(irf_choldec$irf$`epsilon[ reer_mom ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ reer_mom ] %->% reer_mom`)
#Output
sum(irf_choldec$irf$`epsilon[ ind_output_yoy ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ ind_output_yoy ] %->% reer_mom`)

#Smooth transition
model_cv <- id.st(model_VAR, c_lower=0.05, c_upper=0.95, c_step=3, nc=8)
irf_cv <- irf(model_cv, n.ahead = 12, ortho=TRUE)
plot(irf_cv)
#Reserves
sum(irf_cv$irf$`epsilon[ reserves_USD_mln ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ reserves_USD_mln ] %->% reer_mom`)
#Oil
sum(irf_cv$irf$`epsilon[ oil_USD ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ oil_USD ] %->% reer_mom`)
#MIACR 31
sum(irf_cv$irf$`epsilon[ miacr_31 ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ miacr_31 ] %->% reer_mom`)
#REER
sum(irf_cv$irf$`epsilon[ reer_mom ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ reer_mom ] %->% reer_mom`)
#Output
sum(irf_cv$irf$`epsilon[ ind_output_yoy ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ ind_output_yoy ] %->% reer_mom`)

