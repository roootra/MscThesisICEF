library(seasonal)
library(rio)
library(lubridate)
library(vars)
library(ggplot2)
require(BVAR)

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

require(bvarr)
require(combinat)
require(MASS)
require(expm)
bayesian_model <- bvar(data_to_model, lags=2,
                       priors = bv_priors(),
                       mh = bv_mh(),
                       n_thin=10, n_burn=45000, n_draw=50000)

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
#Pure sign approach
horizon <- 1
nvars <- dim(irfs)[1]
signs <- array(NA, dim=c(nvars, nvars, horizon))
#supply shocks
signs[1,1,1] <- 1 #self
signs[2,1,1] <- -1 #inflation
#signs[5,1,1] <- 0 #oil
#demand shocks
signs[2,2,1] <- 1 #self
signs[1,2,1] <- 1 #gdp
signs[3,2,1] <- 1 #intrate
#signs[5,2,1] <- 0 #oil
#monetary shocks
signs[3,3,1] <- 1 #self
signs[1,3,1] <- -1 #gdp
signs[2,3,1] <- -1 #inflation
signs[4,3,1] <- -1 #exrate
#signs[5,3,1] <- 0 #oil
#exrate shocks
signs[4,4,1] <- 1 #self
signs[2,4,1] <- 1 #inflation
signs[3,4,1] <- 1 #monetary
#signs[5,4,1] <- 0 #oil
#oil shocks
signs[5,5,1] <- 1 #self
signs[4,5,1] <- -1 #exrate

signs_lower <- array(dim=dim(signs))
signs_lower[,,1] <- t(signs[,,1])
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
  
  irfs_fr = array(dim=c(horizon+1,n,n)) #first dim is horizon
  for(time in 0:horizon){
    irfs_fr[(time+1),,] = A0_inv %*% t(J_mat) %*% (F_mat %^%(time)) %*% J_mat 
    #delete A0_inv to get usual Phi
    #otherwise, get already Cholesky ortho IRFs
  }
  if(LR){
    A_sum = 0
    for(lag in 1:p){
      A_sum = A_sum + t(A_plus)[1:n, ((2 + n*(lag-1)):(1 + n*lag))]
    }
    irfs_lr <- solve(t(A0) - A_sum)
    return(list("Finite run IRFs"=irfs_fr, "Long run IRFs"=irfs_lr))
  }
  return(irfs_fr)
}

sign_restr_ala_arias <- function(irfs, sign_matrix, tries=300){
  nvars = NCOL(irfs)
  periods = NROW(irfs) / nvars
  S = na.fill(sign_matrix, 0)
  satisfying_models = list()
  for(try in 1:tries){
    cat("Try ", try, ", successes: ", length(satisfying_models), "\r", sep="")
    X = matrix(rnorm(nvars^2), nrow=nvars)
    Q = qr.Q(qr(X))
    #each irf in period h (column-stacked) should be transformed
    irfs_transformed = matrix(nrow=nvars*periods, ncol=nvars)
    for(period in 1:periods){
      irfs_transformed[(1 + nvars*(period-1)):(nvars*period), 1:nvars] = 
        irfs[(1 + nvars*(period-1)):(nvars*period), 1:nvars] %*% Q
    }
    #check for signs
    fail = any(S * irfs_transformed < 0) 
    if(!fail){
      satisfying_models[[length(satisfying_models) + 1]] =
        list("Q" = Q, "Transformed_irfs" = irfs_transformed)
      return(satisfying_models)
    }
  }
}

signs_test <- array(dim=c(5,5,2))
signs_test[,,1] <- signs[,,1]
matrix(aperm(signs_test, c(2,3,1)), nrow=10, ncol=5)

test_arias_irfs <- irf_ala_arias(B=bayesian_model$beta[1,,], Sigma=bayesian_model$sigma[1,,], 
                            p=bayesian_model$meta$lags, n=bayesian_model$meta$M, 
                            horizon=1, LR=F)
sign_restr_ala_arias(irfs=test_arias_irfs[1,,], sign_matrix = signs_test[,,1])


bayesian_sign_restr_arias <- function(B, Sigma, p, n, draws, sr_sign_matrix, 
                                      lr_sign_matrix=NULL, has_const=TRUE){
  sr_sign_matrix_stacked = matrix(aperm(sr_sign_matrix, 
                                         c(2,3,1)), 
                                   nrow=dim(sr_sign_matrix)[1] * dim(sr_sign_matrix)[3], 
                                   ncol=dim(sr_sign_matrix)[1])
  satisfying_models = list()
  for(draw in 1:draws){
    bayesian_ortho_irfs = irf_ala_arias(B=B[draw,,], Sigma=Sigma[draw,,], 
                  p=p, n=n, horizon=dim(sr_sign_matrix)[3], LR=FALSE)
    #stack irfs
    bayesian_ortho_irfs_stacked = matrix(aperm(bayesian_ortho_irfs, c(2,1,3)), 
                                         nrow=dim(bayesian_ortho_irfs)[2]*dim(bayesian_ortho_irfs)[1],
                                         ncol=dim(bayesian_ortho_irfs)[2])
    succ_models_from_draw = sign_restr_ala_arias(irfs=bayesian_ortho_irfs_stacked, 
                         sign_matrix = sr_sign_matrix_stacked)
    sr_sign_matrix_stacked[[length(sr_sign_matrix_stacked) + 1]] = succ_models_from_draw
  }
  return(satisfying_models)
}

bayesian_sign_restr_arias(B=bayesian_model$beta,
                          Sigma=bayesian_model$sigma,
                          p=bayesian_model$meta$lags,
                          n=bayesian_model$meta$M,
                          draws=bayesian_model$meta$n_save,
                          sr_sign_matrix=signs)
