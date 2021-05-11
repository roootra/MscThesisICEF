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
library(openxlsx)

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Monthly_togap.xlsx", sheet=5)
dat <- dat[NROW(dat):1,]

#GDP gap estimation https://www.mathworks.com/help/econ/hpfilter.html#brztihs-1
dat$output_gap <- hpfilter(dat$gdp_per_cap_SA_index, freq=14400)$cycle
dat$output_gap_rate <- dat$output_gap / dat$gdp_per_cap_SA_index
dat <- dat[-1,]

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

data_to_write <- dat_unseas[,c("oil_USD", "imp_price_mom", "miacr_31", 
                             "neer_mom", "gdp_per_cap_SA", "cpi_all_mom")]
data_to_write$neer_mom <- data_to_write$neer_mom * -1
write.xlsx(data_to_write, 
          "/Users/rutra/ВШЭ/Магистратура/Thesis/Scripts/ZeroSignVAR_M/data_sign_and_zero_mon.xlsx")
