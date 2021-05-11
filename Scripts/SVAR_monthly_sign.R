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

#Sign and zero restrictions
data_to_fit <- dat_unseas[,c("oil_USD", "imp_price_mom", "miacr_31", 
                           "neer_mom", "gdp_per_cap_SA", "cpi_all_mom")]
data_to_fit$neer_mom <- data_to_fit$neer_mom * -1
VARselect(data_to_fit)$selection
model_VAR <- VAR(data_to_fit, p = 4, type = "const")

### IDENTIFICATION MATRIX
bvar_irf_ident_mat <- matrix(nrow=NCOL(data_to_fit), ncol=NCOL(data_to_fit))

oil_order = which(colnames(data_to_fit) == "oil_USD")
impprice_order = which(colnames(data_to_fit) == "imp_price_mom")
neer_order = which(colnames(data_to_fit) == "neer_mom")
miacr_order = which(colnames(data_to_fit) == "miacr_31")
gdp_order = which(colnames(data_to_fit) == "gdp_per_cap_SA")
cpi_order = which(colnames(data_to_fit) == "cpi_all_mom")

#Oil shock
bvar_irf_ident_mat[oil_order,oil_order] <- 1 #normalization
bvar_irf_ident_mat[oil_order,impprice_order] <- NA #imp.price, oil is costly
bvar_irf_ident_mat[oil_order,neer_order] <- NA #oil is sold => exrate declines
bvar_irf_ident_mat[oil_order,miacr_order] <- NA #no effect to int.rate
bvar_irf_ident_mat[oil_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[oil_order,cpi_order] <- NA #cpi?

#Import prices
bvar_irf_ident_mat[impprice_order,oil_order] <- 0 #no oil influence?
bvar_irf_ident_mat[impprice_order,impprice_order] <- 1 #normalization
bvar_irf_ident_mat[impprice_order,neer_order] <- NA #NEER?
bvar_irf_ident_mat[impprice_order,miacr_order] <- NA #no effect to int.rate
bvar_irf_ident_mat[impprice_order,gdp_order] <- 1 #gdp?
bvar_irf_ident_mat[impprice_order,cpi_order] <- 1 #cpi grows with import prices

#NEER
bvar_irf_ident_mat[neer_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[neer_order,impprice_order] <- NA #no influence to imp. price
bvar_irf_ident_mat[neer_order,neer_order] <- 1 #normalization
bvar_irf_ident_mat[neer_order,miacr_order] <- NA #may be positive effect to int.rate
bvar_irf_ident_mat[neer_order,gdp_order] <- NA #gdp?
bvar_irf_ident_mat[neer_order,cpi_order] <- 1 #cpi grows, as imported goods become costly

#Interest rate
bvar_irf_ident_mat[miacr_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[miacr_order,impprice_order] <- NA #no influence to imp. price
bvar_irf_ident_mat[miacr_order,neer_order] <- -1 #interest rate is higher => exrate declines
bvar_irf_ident_mat[miacr_order,miacr_order] <- 1 #normalization
bvar_irf_ident_mat[miacr_order,gdp_order] <- -1 #gdp?
bvar_irf_ident_mat[miacr_order,cpi_order] <- -1 #cpi declines

#GDP
bvar_irf_ident_mat[gdp_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[gdp_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[gdp_order,neer_order] <- NA #exrate?
bvar_irf_ident_mat[gdp_order,miacr_order] <- NA #int.rate?
bvar_irf_ident_mat[gdp_order,gdp_order] <- NA#normalized
bvar_irf_ident_mat[gdp_order,cpi_order] <- -1 #cpi?

#CPI
bvar_irf_ident_mat[cpi_order,oil_order] <- 0 #no oil influence
bvar_irf_ident_mat[cpi_order,impprice_order] <- 0 #no influence to imp. price
bvar_irf_ident_mat[cpi_order,neer_order] <- NA #exrate?
bvar_irf_ident_mat[cpi_order,miacr_order] <- 1 #cbr increases interest rate
bvar_irf_ident_mat[cpi_order,gdp_order] <- NA #gdp?#not positive
bvar_irf_ident_mat[cpi_order,cpi_order] <- 1 #normalization

VARselect(data_to_fit)$selection
bvar_model <- bvar(data=data_to_fit, lags=4)
bvar_irf_results <- irf(bvar_model, bv_irf(horizon=12, identification = TRUE,
                                           sign_restr=bvar_irf_ident_mat, sign_lim = 1000))
