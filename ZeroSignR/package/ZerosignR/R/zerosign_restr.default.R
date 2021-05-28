zerosign_restr.default <- function(B, Sigma, p, n, draws, restr_matrix, SR=TRUE, LR=FALSE,
                                   has_const=TRUE, tries=300, varnames=NULL){
  restr_matrix_stacked = NULL
  for(period in 1:dim(restr_matrix)[3]){
    restr_matrix_stacked = rbind(restr_matrix_stacked, restr_matrix[,,period])
  }
  #calculate horizon
  if(LR){
    horiz = dim(restr_matrix)[3]-2
  }
  else{
    horiz = dim(restr_matrix)[3]-1
  }
  #sr_sign_matrix_stacked = matrix(nrow=5, ncol=5)
  #diag(sr_sign_matrix_stacked) = 1
  #sr_sign_matrix_stacked[5, 1:4] <- 0
  #DEBUG!
  cat("\n")
  cat("Restrictions horizon: ",
      ifelse(LR, paste0(ifelse(horiz < 0, "Only", paste0(horiz, " +")), " long-run"), horiz), "\n", sep="")

  satisfying_models = list()
  for(draw in 1:draws){
    cat("\rDraw ", draw, " of ", draws,". Accepted draws: ", length(satisfying_models), ".", sep="")
    bayesian_ortho_irfs = irf_ala_arias(B=B[draw,,], Sigma=Sigma[draw,,],
                                        p=p, n=n, horizon=(dim(restr_matrix)[3]-1), SR=SR, LR=LR)
    succ_models_from_draw = zerosign_restr_ala_arias(irfs=bayesian_ortho_irfs,
                                                     zero_sign_matrix = restr_matrix_stacked,
                                                     tries = tries)
    if(length(succ_models_from_draw) != 0){
      satisfying_models[[length(satisfying_models) + 1]] = succ_models_from_draw
    }
  }
  cat("\n", sep="")
  if(length(satisfying_models) == 0){
    cat("No staisfying models are found!")
  }
  class(satisfying_models) <- "ZerosignR.return"
  attr(satisfying_models, "varnames") <- varnames
  return(satisfying_models)
}
