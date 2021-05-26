library(seasonal)
library(rio)
library(lubridate)
library(vars)
library(ggplot2)

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

data_to_model <- dat_unseas[,c(
  "gdp_real_SA_qoq",
  "cpi_all_qoq",
  "miacr_31",
  "neer_qoq",
  "oil_USD_qoq"
)]
data_to_model$neer_qoq <- data_to_model$neer_qoq*-1

require(BVAR)
require(combinat)
require(MASS)
require(expm)
bayesian_model <- bvar(data_to_model, lags=1,
                       priors = bv_priors(),
                       mh = bv_mh(),
                       n_thin=10, n_burn=95000, n_draw=100000)
horizon <- 1
nvars <- dim(irfs)[1]
signs <- array(NA, dim=c(nvars, nvars, horizon))
#supply shocks
signs[1,1,1] <- 1 #self
signs[2,1,1] <- -1 #inflation
signs[5,1,1] <- 0 #oil
#demand shocks
signs[2,2,1] <- 1 #self
signs[1,2,1] <- 1 #gdp
signs[3,2,1] <- 1 #intrate
signs[5,2,1] <- 0 #oil
#monetary shocks
signs[3,3,1] <- 1 #self
signs[1,3,1] <- -1 #gdp
signs[2,3,1] <- -1 #inflation
signs[4,3,1] <- -1 #exrate
signs[5,3,1] <- 0 #oil
#exrate shocks
signs[4,4,1] <- 1 #self
signs[2,4,1] <- 1 #inflation
signs[3,4,1] <- 1 #monetary
signs[5,4,1] <- 0 #oil
#oil shocks
signs[5,5,1] <- 1 #self
signs[4,5,1] <- -1 #exrate

#VMA representation
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
irfs <- Phi_bayesian(bayesian_model = bayesian_model, nahead=12)

#Cholesky decompositions
choldecs <- apply(bayesian_model$sigma, 1, chol)
#Sign check
sign_models <- list()
tries <- 200
succ <- 0
for(draw in 1:dim(irfs)[4]){
  P <- t(matrix(choldecs[,draw], ncol=5)) #lower Cholesky dec. of vcov mat from draw
  for(try in 1:tries){
    cat("Draw ", draw, ", try ", try, ". Accepted: ", succ,".\033[K\r", sep="")
    new_candidate <- matrix(rnorm(nvars^2), nrow=nvars)
    Q <- qr.Q(qr(new_candidate, complete=T)) #random orthogonal matrix from N(0,1)
    for(current_t in 1:horizon){
      #for(colperm in permn(1:NCOL(Q))){
      irfs_transformed <- irfs[,,current_t,draw] %*% P %*% Q#[,colperm]
      #print(diag(irfs_transformed))
      signs_check <- irfs_transformed * signs_lower[,,current_t]
      fail <- sum(na.fill(signs_check,0) < 0) #check whether restrictions are held
      if(fail == 0){
        succ = succ + 1
        sign_models[[paste0("Draw ",as.character(draw))]][[paste0("Try ", 
                                                                  as.character(try))]][["Q"]] <- Q
        sign_models[[paste0("Draw ",as.character(draw))]][[paste0("Try ", 
                                                                  as.character(try))]][["P"]] <- P
      }
      #}
    }
  }
}

signs_zeros_lower <- signs_lower
signs_zeros_lower[1:4,5,1] <- 0

for(draw in 1:dim(irfs)[4]){
  P <- t(matrix(choldecs[,draw], ncol=5)) #lower Cholesky dec. of vcov mat from draw
  for(try in 1:tries){
    cat("Draw ", draw, ", try ", try, ". Accepted: ", succ,".\033[K\r", sep="")
    new_candidate <- matrix(rnorm(nvars^2), nrow=nvars)
    Q <- qr.Q(qr(new_candidate, complete=T)) #random orthogonal matrix from N(0,1)
    
    for(current_t in 1:horizon){
      #for(colperm in permn(1:NCOL(Q))){
      irfs_transformed <- irfs[,,current_t,draw] %*% P %*% Q[,colperm]
      #print(diag(irfs_transformed))
      signs_check <- irfs_transformed * signs_lower[,,current_t]
      fail <- sum(na.fill(signs_check,0) < 0) #check whether restrictions are held
      if(fail == 0){
        succ = succ + 1
        sign_models[[paste0("Draw ",as.character(draw))]][[paste0("Try ", 
                                                                  as.character(try))]][["Q"]] <- Q
        sign_models[[paste0("Draw ",as.character(draw))]][[paste0("Try ", 
                                                                  as.character(try))]][["P"]] <- P
      }
      #}
    }
  }
}