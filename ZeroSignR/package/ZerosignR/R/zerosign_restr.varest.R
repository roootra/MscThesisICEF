zerosign_restr.varest <- function(varest_model, restr_matrix, LR=FALSE, tries=300){
  require(vars)
  if(!(varest_model$type %in% c("none", "const"))){
    stop(simpleError("VARs with trend are not supported."))
  }
  B = t(Bcoef(varest_model)[,c(NCOL(Bcoef(varest_model)), (1:(NCOL(Bcoef(varest_model)) - 1)))])
  Sigma = summary(varest_model)$covres
  p = varest_model$p
  n = varest_model$K
  has_const = (varest_model$type == "const")
  varnames <- colnames(varest_test$y)
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
  cat("\n")
  cat("Restrictions horizon: ",
      ifelse(LR, paste0(ifelse(horiz < 0, "Only", paste0(horiz, " +")), " long-run"), horiz), "\n", sep="")

  satisfying_models = list()
  cat("\rFrequentist model. Accepted draws: ", length(satisfying_models), ".", sep="")
  bayesian_ortho_irfs = irf_ala_arias(B=B, Sigma=Sigma,
                                      p=p, n=n, horizon=horiz, LR=LR)
  succ_models_from_draw = zerosign_restr_ala_arias(irfs=bayesian_ortho_irfs,
                                                   zero_sign_matrix = restr_matrix_stacked,
                                                   tries = tries)
  if(length(succ_models_from_draw) != 0){
    satisfying_models[[length(satisfying_models) + 1]] = succ_models_from_draw
  }
  cat("\n", sep="")
  if(length(satisfying_models) == 0){
    cat("No staisfying models are found!")
  }
  class(satisfying_models) <- "ZerosignR.return"
  attr(satisfying_models, "varnames") <- varnames
  return(satisfying_models)
}
