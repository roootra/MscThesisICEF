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
                       n_thin=10, n_burn=9000, n_draw=10000)


#Pure sign approach
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

irf_ala_arias <- function(B, Sigma, p, n, horizon, LR=T){
  #B = A+ %*% A0^(-1) --- matrix of reduced parameters
  #in form B = [c, B1, ..., Bp]
  #Sigma = covariance matrix of error term
  #p = order of model
  #n = number of variables
  #horizon = number of periods to calculate IRF for
  #LR = boolean variable, whether to calculate long-run IRF
  
  A0_inv = chol(Sigma)
  A0 = solve(A0_inv)
  A_plus = B %*% A0
  
  #Finite-run IRFs
  F_mat = matrix(data=0, nrow = p*n, ncol = p*n)
  F_mat[,1:n] = t(B[-1,])
  if(p > 1){
    F_mat[1:((p-1)*n), (n+1):NCOL(F_mat)] = diag(nrow=(p-1)*n, ncol = (p-1)*n)
  }
  
  J_mat = matrix(data=0, nrow=p*n, ncol=n)
  J_mat[1:n, 1:n] = diag(nrow=n, ncol=n)
  
  irfs = matrix(nrow=(horizon+1)*n, ncol=n)
  for(time in 0:horizon){
    irfs[(1 + n*(time)):(n*(time+1)),] = 
      A0_inv %*% t(J_mat) %*% (F_mat %^%(time)) %*% J_mat 
    #delete A0_inv to get usual Phi
    #otherwise, get already Cholesky ortho IRFs
  }
  if(LR){
    A_sum = 0
    for(lag in 1:p){
      A_sum = A_sum + t(A_plus)[1:n, ((2 + n*(lag-1)):(1 + n*lag))]
    }
    irfs_lr <- solve(t(A0) - A_sum)
    irfs <- rbind(irfs, irfs_lr)
  }
  return(irfs)
}

sign_restr_ala_arias <- function(irfs, sign_matrix, tries=300, perm_Q=FALSE){
  nvars = NCOL(irfs)
  periods = NROW(irfs) / nvars
  S = na.fill(sign_matrix, 0)
  satisfying_models = list()
  for(try in 1:tries){
    #cat("Try ", try, ", successes: ", length(satisfying_models), "\r", sep="")
    X = matrix(rnorm(nvars^2), nrow=nvars)
    Q = qr.Q(qr(X))
    #each irf in period h (column-stacked) should be transformed
    irfs_transformed = matrix(nrow=nvars*periods, ncol=nvars)
    if(!perm_Q){ 
      for(period in 1:periods){ #instead of this, stack Q horisontally (columnwise) period times
        #It might be better to follow the paper here
        irfs_transformed[(1 + nvars*(period-1)):(nvars*period), 1:nvars] = 
          irfs[(1 + nvars*(period-1)):(nvars*period), 1:nvars] %*% Q
      }
    }
    else{ #TBD
      q_col_perm = permn(1:NCOL(Q))
        for(period in 1:periods){
          irfs_transformed[(1 + nvars*(period-1)):(nvars*period), 1:nvars] = 
            irfs[(1 + nvars*(period-1)):(nvars*period), 1:nvars] %*% Q
        }
    }
    #check for signs
    fail = any(S * irfs_transformed < 0) 
    if(!fail){
      satisfying_models[[length(satisfying_models) + 1]] =
        list("Q" = Q, "Transformed_irfs" = irfs_transformed)
      
    }
  }
  return(satisfying_models)
}

signs_test <- array(dim=c(5,5,2))
signs_test[,,1] <- signs[,,1]
matrix(aperm(signs_test, c(2,3,1)), nrow=10, ncol=5)

test_arias_irfs <- irf_ala_arias(B=bayesian_model$beta[1,,], Sigma=bayesian_model$sigma[1,,], 
                            p=bayesian_model$meta$lags, n=bayesian_model$meta$M, 
                            horizon=0, LR=F)
sign_restr_ala_arias(irfs=test_arias_irfs, sign_matrix = diag(nrow=5, ncol=5), tries = 10000)
sign_restr_ala_arias(irfs=test_arias_irfs, sign_matrix = signs[,,1], tries = 10000)


bayesian_sign_restr_arias <- function(B, Sigma, p, n, draws, sr_sign_matrix, 
                                      lr_sign_matrix=NULL, has_const=TRUE, tries=1000){
  sr_sign_matrix_stacked = matrix(aperm(sr_sign_matrix, 
                                         c(2,3,1)), 
                                   nrow=dim(sr_sign_matrix)[1] * dim(sr_sign_matrix)[3], 
                                   ncol=dim(sr_sign_matrix)[1])
  #sr_sign_matrix_stacked = diag(nrow=5, ncol=5) #DEBUG!
  cat("\n")
  cat("Restrictions horizon: ", (dim(sr_sign_matrix)[3]-1), "\n", sep="")
  satisfying_models = list()
  for(draw in 1:draws){
    cat("\rDraw ", draw, ". Accepted: ", length(satisfying_models), ".", sep="")
    bayesian_ortho_irfs = irf_ala_arias(B=B[draw,,], Sigma=Sigma[draw,,], 
                  p=p, n=n, horizon=(dim(sr_sign_matrix)[3]-1), LR=FALSE)
    succ_models_from_draw = sign_restr_ala_arias(irfs=bayesian_ortho_irfs, 
                         sign_matrix = sr_sign_matrix_stacked, tries = tries)
    if(length(succ_models_from_draw) != 0){
      satisfying_models[[length(satisfying_models) + 1]] = succ_models_from_draw
    }
  }
  cat("\n", sep="")
  return(satisfying_models)
}

bayesian_model_uninf <- bvar(data_to_model, lags=1,
                       priors = bv_priors(sur=bv_sur()),
                       mh = bv_mh(),
                       n_thin=10, n_burn=95000, n_draw=100000)

bayesian_sign_restr_arias(B=bayesian_model$beta,
                          Sigma=bayesian_model$sigma,
                          p=bayesian_model$meta$lags,
                          n=bayesian_model$meta$M,
                          draws=bayesian_model$meta$n_save,
                          sr_sign_matrix=signs,
                          tries=1000)

zerosign_restr_ala_arias <- function(irfs, zero_sign_matrix, tries=300, perm_Q=FALSE){
  nvars = NCOL(irfs)
  periods = NROW(irfs) / nvars
  S = na.fill(zero_sign_matrix, 0)
  Z = (zero_sign_matrix == 0)
  Z = na.fill(Z, 0)
  satisfying_models = list()
  succ = 0
  fails = 0
  cat("\n")
  for(try in 1:tries){
    #each irf in period h (column-stacked) should be transformed
    irfs_transformed = matrix(nrow=nvars*periods, ncol=nvars)
    #R-matrix recursive construction
    Q = matrix(nrow=nvars, ncol=nvars)
    for(j in 1:nvars){
      Z_j = diag(mapply(Z[,j], FUN=as.numeric))
      Z_j = Z_j[,(colSums(Z_j) != 0)]
      n_zeros = sum(Z_j)
      if(j > 1){
        if(n_zeros > 0){
          R_j = rbind(Z_j %*% irfs, t(Q[,1:(j-1)]))
        }
        else{
          R_j = t(Q[,1:(j-1)])
        }
      }
      else{
        R_j = Z_j %*% irfs
      }
      N_jminus = Null(t(R_j))
      n_j = NCOL(N_jminus)
      y_j = rnorm(n_j)
      Q[,j] = N_jminus %*% y_j / norm(y_j, type="2")
    }
    if(!perm_Q){
      for(period in 1:periods){ #instead of this, stack Q horisontally (columnwise) period times
        #It might be better to follow paper here (e_j)
        irfs_transformed[(1 + nvars*(period-1)):(nvars*period), 1:nvars] = 
          irfs[(1 + nvars*(period-1)):(nvars*period), 1:nvars] %*% Q
      }
    }
    else{ #TBD
      q_col_perm = permn(1:NCOL(Q))
      for(period in 1:periods){
        irfs_transformed[(1 + nvars*(period-1)):(nvars*period), 1:nvars] = 
          irfs[(1 + nvars*(period-1)):(nvars*period), 1:nvars] %*% Q
      }
    }
    #check for signs
    flag_fail = FALSE
    for(j in 1:nvars){
      S_j = diag(S[,j])
      S_j = t(S_j[,(colSums(S_j) != 0)])
      e_j = diag(nrow=nvars)[,j]
      if(any(S_j %*% irfs %*% e_j < 0)){
        flag_fail = TRUE
        break
      }
    }
    if(flag_fail){
      fails = fails + 1
      cat("\rTry ", try, ", successes: ", succ, ". Failed: ", fails, ".", sep="")
      }
    else{
      succ = succ + 1
      cat("\rTry ", try, ", successes: ", succ, ". Failed: ", fails, ".", sep="")
      satisfying_models[[length(satisfying_models) + 1]] =
        list("Q" = Q, "Transformed_irfs" = irfs_transformed)
    }
  }
  cat("\n")
  return(satisfying_models)
}
debug(zerosign_restr_ala_arias)
zerosign_restr_ala_arias(irfs_transformed, signs[,,1], tries = 10000)
undebug(zerosign_restr_ala_arias)

bayesian_zero_sign_restr_arias <- function(B, Sigma, p, n, draws, sr_sign_matrix, 
                                           has_const=TRUE, tries=1000){
  sr_sign_matrix_stacked = NULL
  for(period in 1:dim(sr_sign_matrix)[3]){
    sr_sign_matrix_stacked = rbind(sr_sign_matrix_stacked, sr_sign_matrix[,,period])
  }
  #sr_sign_matrix_stacked = matrix(nrow=5, ncol=5)
  #diag(sr_sign_matrix_stacked) = 1
  #sr_sign_matrix_stacked[5, 1:4] <- 0
  #DEBUG!
  cat("\n")
  cat("Restrictions horizon: ", (dim(sr_sign_matrix)[3]-1), "\n", sep="")
  satisfying_models = list()
  for(draw in 1:draws){
    cat("\rDraw ", draw, ". Accepted draws: ", length(satisfying_models), ".", sep="")
    bayesian_ortho_irfs = irf_ala_arias(B=B[draw,,], Sigma=Sigma[draw,,], 
                                        p=p, n=n, horizon=(dim(sr_sign_matrix)[3]-1), LR=FALSE)
    succ_models_from_draw = zerosign_restr_ala_arias(irfs=bayesian_ortho_irfs, 
                                                 zero_sign_matrix = sr_sign_matrix_stacked, 
                                                 tries = tries)
    if(length(succ_models_from_draw) != 0){
      satisfying_models[[length(satisfying_models) + 1]] = succ_models_from_draw
    }
  }
  cat("\n", sep="")
  return(satisfying_models)
}
debug(bayesian_zero_sign_restr_arias)
bayesian_zero_sign_restr_arias(B=bayesian_model$beta[1:100,,],
                          Sigma=bayesian_model$sigma[1:100,,],
                          p=bayesian_model$meta$lags,
                          n=bayesian_model$meta$M,
                          draws=100,#bayesian_model$meta$n_save,
                          sr_sign_matrix=signs,
                          tries=1000)
undebug(bayesian_zero_sign_restr_arias)

