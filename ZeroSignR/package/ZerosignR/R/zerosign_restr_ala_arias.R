zerosign_restr_ala_arias <- function(irfs, zero_sign_matrix, tries, perm_Q=FALSE){
  #' Perform simulations to identify shock under zero and sign restrictions
  #'
  #' This function simulates orthognormal random matrices Q given number of times for
  #' given IRF matrix in order to identify structural shocks. It returns a list of
  #' Q matrices and transformed IRFs that satisfy imposed restrictions.
  #'
  #' @param irfs IRF matrix of size $(nvars \times (nvars*horizons))$ ---
  #' columnwise-stacked IRFs by periods if there is more than one period of IRF.
  #' @param zero_sign_matrix Square matrix of zero and sign restrictions $(nvars \times nvars)$.
  #' Columns = shocks, structural and residual ones. Put $0$ into $(i,j)$-th cell of
  #' this matrix in order to set zero restriction to the $i$-th variable response to
  #' the shock $j$. Note that user should order
  #' variables in descending order of number of zero restrictions. If there are
  #' two or more variables with the same number of zero restrictions, their order
  #' can be arbitrary. There cannot be more than $nvars - j$ zero restrictions
  #' for $j$-th column (in total for all time periods).
  #' @param tries Number of tries of random orthonormal matrix generation
  #' @param perm_Q NOT IMPLEMENTED. Whether to account for permutation of matrix Q
  #' columns to increase method's efficiency.
  nvars = NCOL(irfs)
  periods = NROW(irfs) / nvars
  S = na.fill(zero_sign_matrix, 0)
  Z = (zero_sign_matrix == 0)
  Z = na.fill(Z, 0)
  zero_restrictions_present = ifelse(any(Z == 1), TRUE, FALSE)
  satisfying_models = list()
  succ = 0
  fails = 0
  for(try in 1:tries){
    #each irf in period h (column-stacked) should be transformed
    if(zero_restrictions_present){
      #If there are sign restrictions, impose linear constraints
      #on Q and calculate it recursively, as given by Theorem 4 in (Arias et al., 2014)
      #Recursively-constructed Q-matrix is still orthonormal!
      Q = matrix(nrow=nvars, ncol=nvars)
      for(j in 1:nvars){
        Z_j = diag(mapply(Z[,j], FUN=as.numeric))
        Z_j = t(Z_j[,(colSums(Z_j) != 0)])
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
          stopifnot("Please check that shocks are ordered in descending order of number of zero restrictions."= NCOL(Z_j) != 0)
          R_j = Z_j %*% irfs
        }
        N_jminus = Null(t(R_j))
        stopifnot("Error: check the number of restrictions for each shock: there cannot be more restrictions for one shock (in total for all IRF periods) than total number of shocks minus order of this shock."= NCOL(N_jminus) != 0)
        n_j = NCOL(N_jminus)
        y_j = rnorm(n_j)
        Q[,j] = N_jminus %*% y_j / norm(y_j, type="2")
      }
    }
    else{
      #If there are no zero restrictions, use Theorem 1
      X = matrix(rnorm(nvars^2), nrow=nvars)
      Q = qr.Q(qr(X))
    }
    #check for signs
    flag_fail = FALSE
    for(j in 1:nvars){
        S_j = diag(S[,j])
        S_j = t(S_j[,(colSums(S_j) != 0)])
        e_j = diag(nrow=nvars)[,j]
        if(any(S_j %*% irfs %*% Q[,j] < 0)){
          flag_fail = TRUE
          break
        }
      }
    if(flag_fail){
      fails = fails + 1
      #cat("\rTry ", try, ", successes: ", succ, ". Failed: ", fails, ".", sep="")
    }
    else{
      irfs_transformed = irfs %*% Q
      succ = succ + 1
      #cat("\rTry ", try, ", successes: ", succ, ". Failed: ", fails, ".", sep="")
      satisfying_models[[length(satisfying_models) + 1]] =
        list("Q" = Q, "Transformed_irfs" = irfs_transformed)
    }
  }
  #cat("\r                                                                                ")
  return(satisfying_models)
}
