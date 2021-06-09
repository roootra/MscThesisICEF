library(BVAR)
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
dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Quarterly.xlsx", sheet=5)
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

dat_unseas$neer_qoq <- dat_unseas$neer_qoq*-1

#Output gap
dat_unseas$nominal_gdp_gap <- hpfilter(dat_unseas$gdp_nominal_index, freq=1600)$cycle
dat_unseas$real_gdp_gap <- hpfilter(dat_unseas$gdp_real_SA_index, freq=1600)$cycle
#dat_unseas$d_real_gdp_gap <- c(NA, diff(dat_unseas$real_gdp_gap))
#dat_unseas <- na.omit(dat_unseas)
#dat_unseas$miacr_90 <- c(NA, diff(dat_unseas$miacr_90))
#dat_unseas <- na.omit(dat_unseas)
plot(dat_unseas$date, dat_unseas$real_gdp_gap)

data_to_fit <- dat_unseas[,c("oil_USD_qoq",
                             #"imp_price_qoq",
                             "miacr_31",
                             "neer_qoq",
                             #"reserves_USD_qoq",
                             #"broad_money_SA",
                             "gdp_real_SA_qoq",
                             "cpi_all_qoq")]
data_to_fit_ts <- ts(data_to_fit)

library(dynlm)
ardl <- dynlm(cpi_all_qoq ~ L(neer_qoq, 0:3) + L(gdp_real_SA_qoq, 0:2) +
                L(oil_USD_qoq, 0:2), data=data_to_fit_ts)
summary(ardl)
cat("Short-run unconditional ERPT is: ", sum(ardl$coefficients[2:3]), ".\n", sep="")
cat("Long-run unconditional ERPT is: ", sum(ardl$coefficients[2:5]), ".\n", sep="")
library(car)
lht(ardl, "L(neer_qoq, 0:3)0 + L(neer_qoq, 0:3)1 = 0")
lht(ardl, "L(neer_qoq, 0:3)0 + L(neer_qoq, 0:3)1 + L(neer_qoq, 0:3)2 + L(neer_qoq, 0:3)3 = 0")
