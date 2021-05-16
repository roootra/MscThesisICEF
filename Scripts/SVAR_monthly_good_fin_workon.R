library(seasonal)
library(rio)
library(lubridate)
library(bvarsv)
library(VARsignR)
library(vars)
library(tsDyn)
library(ggplot2)
library(mFilter)
library(svars)

dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Monthly.xlsx", sheet=5)
dat <- dat[NROW(dat):1,]

#Deseasonalization
to_deseasonalize <- c("ind_output_yoy","imp_price_mom", "exp_price_mom", 
                      grep("cpi.*", colnames(dat), value=T))
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
data_chol <- dat_unseas[,c("oil_USD", "imp_price_mom", 
                           "miacr_31", "neer_mom", "gdp_per_cap_SA", "cpi_all_mom")]
#write.csv(data_chol,
#          "/Users/rutra/ВШЭ/Магистратура/Thesis/Scripts/graph_builder/data_w_impprices.csv", row.names = FALSE)
data_chol$neer_mom <- data_chol$neer_mom * -1
VARselect(data_chol)$selection
model_VAR <- VAR(data_chol, p = 4, type = "const")
choldec <- id.chol(model_VAR)
irf_choldec <- irf(choldec, n.ahead = 12, ortho=TRUE)
plot(irf_choldec)
#Import prices
sum(irf_choldec$irf$`epsilon[ imp_price_mom ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ imp_price_mom ] %->% neer_mom`)
#Oil
sum(irf_choldec$irf$`epsilon[ oil_USD ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ oil_USD ] %->% neer_mom`)
#MIACR 31
sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ miacr_31 ] %->% neer_mom`)
#NEER
sum(irf_choldec$irf$`epsilon[ neer_mom ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ neer_mom ] %->% neer_mom`)
#Output
sum(irf_choldec$irf$`epsilon[ gdp_per_cap_SA ] %->% cpi_all_mom`) /
  sum(irf_choldec$irf$`epsilon[ gdp_per_cap_SA ] %->% neer_mom`)

#Smooth transition
model_cv <- id.st(model_VAR, c_lower=0.05, c_upper=0.95, c_step=3, nc=8, c_fix=70) #70
irf_cv <- irf(model_cv, n.ahead = 12, ortho=TRUE)
plot(irf_cv)
#Import price
sum(irf_cv$irf$`epsilon[ imp_price_mom ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ imp_price_mom ] %->% neer_mom`)
#Oil
sum(irf_cv$irf$`epsilon[ oil_USD ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ oil_USD ] %->% neer_mom`)
#MIACR 31
sum(irf_cv$irf$`epsilon[ miacr_31 ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ miacr_31 ] %->% neer_mom`)
#NEER
sum(irf_cv$irf$`epsilon[ neer_mom ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ neer_mom ] %->% neer_mom`)
#Output
sum(irf_cv$irf$`epsilon[ gdp_per_cap_SA ] %->% cpi_all_mom`) /
  sum(irf_cv$irf$`epsilon[ gdp_per_cap_SA ] %->% neer_mom`)

corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}
