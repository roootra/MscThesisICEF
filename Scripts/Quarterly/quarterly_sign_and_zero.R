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
                               "neer_qoq",
                              "miacr_31",
                               #"reserves_USD_qoq",
                               #"broad_money_SA",
                               "gdp_real_SA_qoq", 
                               "cpi_all_qoq")]
bvar_irf_ident_mat <- matrix(nrow=NCOL(data_to_fit), ncol=NCOL(data_to_fit))
### IDENTIFICATION MATRIX
oil_order = which(colnames(data_to_fit) == "oil_USD_qoq")
impprice_order = which(colnames(data_to_fit) == "imp_price_qoq")
neer_order = which(colnames(data_to_fit) == "neer_qoq")
miacr_order = which(colnames(data_to_fit) == "miacr_31")
reserves_order = which(colnames(data_to_fit) == "reserves_USD_qoq")
money_order = which(colnames(data_to_fit) == "broad_money_SA")
#gdp_order = which(colnames(data_to_fit) == "gdp_nominal_qoq")
gdp_order = which(colnames(data_to_fit) == "gdp_real_SA_qoq")
cpi_order = which(colnames(data_to_fit) == "cpi_all_qoq")

#Oil shock
bvar_irf_ident_mat[oil_order,oil_order] <- 1 #normalization
bvar_irf_ident_mat[oil_order,impprice_order] <- 1 #imp.price, oil is costly
bvar_irf_ident_mat[oil_order,neer_order] <- -1 #oil is sold => exrate declines
bvar_irf_ident_mat[oil_order,reserves_order] <- NA #cbr buys depreciated usd??
bvar_irf_ident_mat[oil_order,money_order] <- 0 #no effect to money 
bvar_irf_ident_mat[oil_order,miacr_order] <- NA #no effect to int.rate
bvar_irf_ident_mat[oil_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[oil_order,cpi_order] <- NA #cpi?

#Import prices
bvar_irf_ident_mat[impprice_order,oil_order] <- 0 #no oil influence?
bvar_irf_ident_mat[impprice_order,impprice_order] <- 1 #normalization
bvar_irf_ident_mat[impprice_order,neer_order] <- NA #NEER?
bvar_irf_ident_mat[impprice_order,reserves_order] <- 0 #no effect to reserves
bvar_irf_ident_mat[impprice_order,money_order] <- 0 #no effect to money
bvar_irf_ident_mat[impprice_order,miacr_order] <- 0 #no effect to int.rate
bvar_irf_ident_mat[impprice_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[impprice_order,cpi_order] <- 1 #cpi grows with import prices

#NEER
bvar_irf_ident_mat[neer_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[neer_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[neer_order,neer_order] <- 1 #normalization
bvar_irf_ident_mat[neer_order,reserves_order] <- NA #but it may lower reserves, as cb compensates deprec.
bvar_irf_ident_mat[neer_order,money_order] <- NA #money
bvar_irf_ident_mat[neer_order,miacr_order] <- NA #may be positive effect to int.rate
bvar_irf_ident_mat[neer_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[neer_order,cpi_order] <- 1 #cpi grows, as imported goods become costly

#Interest rate
bvar_irf_ident_mat[miacr_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[miacr_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[miacr_order,neer_order] <- -1 #interest rate is higher => exrate declines
bvar_irf_ident_mat[miacr_order,reserves_order] <- NA #reserves?
bvar_irf_ident_mat[miacr_order,money_order] <- NA #money?
bvar_irf_ident_mat[miacr_order,miacr_order] <- 1 #normalization
bvar_irf_ident_mat[miacr_order,gdp_order] <- -1 #gdp?
bvar_irf_ident_mat[miacr_order,cpi_order] <- -1 #cpi declines

#Reserves
bvar_irf_ident_mat[reserves_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[reserves_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[reserves_order,neer_order] <- NA #reserves are bought for foreign currency
bvar_irf_ident_mat[reserves_order,reserves_order] <- 1 #normalization
bvar_irf_ident_mat[reserves_order,money_order] <- NA #effect to money?
bvar_irf_ident_mat[reserves_order,miacr_order] <- 0 #no effect to int.rate
bvar_irf_ident_mat[reserves_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[reserves_order,cpi_order] <- NA #cpi?

#Broad money
bvar_irf_ident_mat[money_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[money_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[money_order,neer_order] <- NA #more money supply => rouble's depreciation
bvar_irf_ident_mat[money_order,reserves_order] <- NA #reserves
bvar_irf_ident_mat[money_order,money_order] <- 1 #normalization
bvar_irf_ident_mat[money_order,miacr_order] <- NA #int.rate my be risen
bvar_irf_ident_mat[money_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[money_order,cpi_order] <- 1 #cpi grows with money supply

#GDP
bvar_irf_ident_mat[gdp_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[gdp_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[gdp_order,neer_order] <- NA #exrate?
bvar_irf_ident_mat[gdp_order,reserves_order] <- NA #reserves?
bvar_irf_ident_mat[gdp_order,money_order] <- NA #moneyness 
bvar_irf_ident_mat[gdp_order,miacr_order] <- NA #int.rate?
bvar_irf_ident_mat[gdp_order,gdp_order] <- 1 #normalized
bvar_irf_ident_mat[gdp_order,cpi_order] <- NA #cpi?

#CPI
bvar_irf_ident_mat[cpi_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[cpi_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[cpi_order,neer_order] <- NA #exrate?
bvar_irf_ident_mat[cpi_order,reserves_order] <- 0 #reserves?
bvar_irf_ident_mat[cpi_order,money_order] <- NA #broad money?
bvar_irf_ident_mat[cpi_order,miacr_order] <- 1 #cbr increases interest rate
bvar_irf_ident_mat[cpi_order,gdp_order] <- NA #gdp?#not positive
bvar_irf_ident_mat[cpi_order,cpi_order] <- 1 #normalization


###
bvar_model <- bvar(data=data_to_fit, lags=1)
#bvar_irf_ident <- bv_irf(horizon=4, identification = TRUE,
#                         sign_restr=bvar_irf_ident_mat, sign_lim = 5000)
bvar_irf_results <- irf(bvar_model, bv_irf(horizon=4, identification = TRUE,
                            sign_restr=bvar_irf_ident_mat, sign_lim = 1000), conf_bands=0.5)


irf_median <- apply(bvar_irf_results$irf, c(2,3,4), median)
irf_median <- bvar_irf_results$quants
irf_median_sum <- apply(irf_median, c(1,3), sum)
sum(irf_median[1,,5])/sum(irf_median[1,,2]) #oil price shock
sum(irf_median[2,,5])/sum(irf_median[2,,2]) #exrate shock
sum(irf_median[3,,5])/sum(irf_median[3,,2]) #int. rate shock
sum(irf_median[4,,5])/sum(irf_median[4,,2]) #output shock
sum(irf_median[5,,5])/sum(irf_median[5,,2]) #exog. cpi shock
test <- summary(bvar_irf_results)
