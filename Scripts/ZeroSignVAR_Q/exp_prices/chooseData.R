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

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Quarterly.xlsx", sheet=3)
dat <- na.omit(dat[dat$date > as.Date("2005-01-01"),])
dat <- dat[NROW(dat):1,]

#Deseasonalization
to_deseasonalize <- c("gdp_nominal_index", "imp_price_qoq",
                      "exp_price_qoq",
                      grep("cpi.*", colnames(dat), value=T))
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

data_to_write <- dat_unseas[,c(#"date",
  "gdp_real_SA_qoq",
  "cpi_all_qoq",
  "miacr_31",
  "neer_qoq",
  "exp_price_qoq",
  "oil_USD_qoq"
  )]
#data_to_write$neer_qoq <- data_to_write$neer_qoq*-1
#data_to_write$nom_usd_qoq <- data_to_write$nom_usd_qoq*-1
VARselect(data_to_write, lag.max=6)
write.csv(data_to_write, 
          "/Users/rutra/ВШЭ/Магистратура/Thesis/Scripts/ZeroSignVAR_Q/exp_prices/data_sign_and_zero.csv",
          row.names=FALSE, col.names=NA)

