#remotes::install_github("shounakch/FABLE")
library(FABLE)

library(MOFA2)
library(reticulate)
use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python312", required=TRUE)
outfile = file.path('results',"mofa_fit_1.hdf5")

source('competitors/FACTOR_ANALYSIS/FACTOR_CODE_update.R')

prepare_data_4_mofa <- function(Y){
  V <- length(Y)
  n <- nrow(Y[[1]])
  data <- list()
  samples_name <- paste0('sample_', seq(1, n))
  for(v in 1:V){
    p_v <-ncol(Y[[v]])
    features_name <- paste0('feature_', seq(1, p_v), '_view', v)
    Y_v <- data.frame(t(Y[[v]]))
    names(Y_v) <- samples_name
    row.names(Y_v) <- features_name
    data[[paste0('view_',v)]] = as.matrix(Y_v)
  }
  return(data)
}



fit_MOFA <- function(
    Y, k_0, spikeslab_factors=T, ard_factors=F, spikeslab_weights=T, ard_weights=T) {
  
  Y_mofa <- prepare_data_4_mofa(Y)
  MOFAobject <- create_mofa(Y_mofa)
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors = k_0
  model_opts$spikeslab_factors = spikeslab_factors
  model_opts$ard_factors = ard_factors # if TRUE accuracy decreases
  model_opts$spikeslab_weights = spikeslab_weights
  model_opts$ard_weights = ard_weights
  train_opts <- get_default_training_options(MOFAobject)
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  # C:\Users\mauri\anaconda3
  use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python312", required=TRUE)
  outfile = file.path('results',"mofa_fit_1.hdf5")
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  
  factors <- get_factors(MOFAobject.trained, factors = "all")
  mofa_weights_1 <- get_weights(MOFAobject.trained, views = "all", factors = "all")
  #mofa_weights_2 <- get_expectations(MOFAobject.trained, "W")
  
  return(list(
    Lambdas_hat = mofa_weights_1,
    F_hat = factors)
  )
  
}



compute_performance_mofa <- function(fit, Lambdas_0, As){
  M <- length(fit$Lambdas_hat)
  Lambda_hat_c <- do.call(rbind, fit$Lambdas_hat)
  Lambda_0_c <- Lambdas_0[[1]] %*% t(As[[1]])
  for(m in 2:M){
    Lambda_0_c <- rbind(Lambda_0_c, Lambdas_0[[m]] %*% t(As[[m]]))
  }
  rmse_all <- norm( ( tcrossprod(Lambda_hat_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F')
  print('RMSE all')
  print(rmse_all)
  rmses_view <- c()
  rmses_biview <- c()
  print('RMSE intra-view')
  for(m in 1:M){
    print(m)
    Lambda_hat_outer_s <- tcrossprod(fit$Lambdas_hat[[m]])
    Lambda_0_outer_s <- tcrossprod(Lambdas_0[[m]])
    rmses_view[m] <- norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F')
    print(rmses_view[m])
  }
  
  print('RMSE inter-view')
  for(m in 1:(M-1)){
    print(m)
    for(v in (m+1):M) {
      print(v)
      Lambda_hat_outer_sv <- fit$Lambdas_hat[[m]] %*% t(fit$Lambdas_hat[[v]])
      Lambda_0_outer_sv <- Lambdas_0[[m]] %*% t(As[[m]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
      rmses_biview <- c(rmses_biview, norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
      print(rmses_biview[length(rmses_biview)])
    }
  }
  
  
  return(list(
    rmse_all = rmse_all,
    rmses_intra = rmses_view,
    rmses_inter = rmses_biview)
  )
  
  
  
}




PseudoPosteriorSampler <- function(fit,
                                   Y,
                                   gamma0 = 1,
                                   delta0sq = 1,
                                   maxProp = 0.5,
                                   MC = 1000) {
  
  tFABLESample1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  #kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  kEst = fit$estRank
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  CovCorrectMatrix = cov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
                                        FABLEHypPars$G)
  
  varInflation = (sum(CovCorrectMatrix) / (p*(p+1)/2))^2
  
  FABLESamples = FABLESampler(Y, 
                              gamma0, 
                              delta0sq, 
                              MC,
                              U_Y,
                              V_Y,
                              svalsY,
                              kEst,
                              FABLEHypPars$tauSqEstimate,
                              FABLEHypPars$gammaDeltasq,
                              FABLEHypPars$G0,
                              varInflation)
  
  tFABLESample2 = proc.time()
  tSample = (tFABLESample2 - tFABLESample1)[3]
  
  OutputList = list("CCFABLESamples" = FABLESamples,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "varInflation" = varInflation,
                    "runTime" = tSample)
  
  return(OutputList)
  
}




compute_performance_fable <- function(
    fit, Y, Lambdas_0, As, n_MC = 500, subsample_index=1:100, UQ=T
    ){
  Y_c <- do.call(cbind, Y)
  M <- length(Lambdas_0)
  Lambda_outer_hat_fable <- fit$FABLEPostMean
  Lambda_0_c <- Lambdas_0[[1]] %*% t(As[[1]])
  for(m in 2:M){
    Lambda_0_c <- rbind(Lambda_0_c, Lambdas_0[[m]] %*% t(As[[m]]))
  }
  print(dim(tcrossprod(Lambda_0_c)))
  rmse_all <- norm( ( Lambda_outer_hat_fable - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F')
  print('RMSE all')
  print(rmse_all)
  rmses_view <- c()
  rmses_biview <- c()
  print('RMSE intra-view')
  idx_start <-1 
  p_s <- sapply(Lambdas_0, function(x) (nrow(x)))
  for(m in 1:M){
    print(m)
    Lambda_hat_outer_s <- Lambda_outer_hat_fable[idx_start:(idx_start + p_s[m]-1), idx_start:(idx_start + p_s[m]-1)]
    idx_start <- idx_start + p_s[m]
    Lambda_0_outer_s <- tcrossprod(Lambdas_0[[m]])
    rmses_view[m] <-  norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F')
    print(rmses_view[m])
  }
  
  print('RMSE inter-view')
  idx_start <-1 
  for(s in 1:(M-1)){
    print(s)
    for(v in (s+1):M) {
      print(v)
      idx_start_2 <- 1 + sum(p_s[1:(v-1)])
      Lambda_hat_outer_sv <- Lambda_outer_hat_fable[idx_start:(idx_start + p_s[s]-1), idx_start_2:(idx_start_2 + p_s[v]-1)]
      Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
      rmses_biview <- c(rmses_biview, norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
      print(rmses_biview[length(rmses_biview)])
    }
    idx_start <- idx_start + p_s[s]
  }
  
  ps <- sapply(Y, ncol)
  subsample_index_c <- subsample_index
  p_init <- ps[1]
  for(m in 2:M){
    print(m)
    subsample_index_c <- c(subsample_index_c, p_init + subsample_index)
    p_init <- p_init + ps[m]
  }
  
  out <- list(
    rmse_all = rmse_all,
    rmses_intra = rmses_view,
    rmses_inter = rmses_biview)
  #time_model_uq = time_uq[3]
  
  if(!UQ){
    return(out)
  }
  
  
  ptm <- proc.time()
  post_samples_fable <- PseudoPosteriorSampler(fable_fit_1, Y_c[,subsample_index_c], MC = n_MC) 
  
  fable_qs <- CCFABLEPostProcessing(post_samples_fable$CCFABLESamples,
                                    alpha=0.05)
  time_uq <- proc.time() - ptm; 
  
  fable_l <- fable_qs$LowerQuantileMatrix
  fable_u <- fable_qs$UpperQuantileMatrix
  
  cov_posterior_view <- c()
  cov_posterior_biview <- c()
  len_posterior_view <- c()
  len_posterior_biview <- c()
  print('coverage posterior intra-view')
  idx_start <- 0
  for(m in 1:M){
    print(m)
    Lambda_0_outer_m <- tcrossprod(Lambdas_0[[m]][subsample_index,])
    l_ci <- fable_l[(idx_start + subsample_index), (idx_start + subsample_index)]
    u_ci <- fable_u[(idx_start + subsample_index), (idx_start + subsample_index)]
    cov_posterior_view[m] <- mean((l_ci<Lambda_0_outer_m) & (u_ci>Lambda_0_outer_m))
    print(cov_posterior_view[m])
    len_posterior_view[m] <- mean(u_ci - l_ci)
    idx_start <- idx_start + length(subsample_index)
  }
  print('coverage posterior inter-view')
  idx_start_1 <- 0
  idx_start_2 <- 0
  for(m in 1:(M-1)){
    print(m)
    for(v in (m+1):M) {
      print(v)
      idx_start_2 <- (v-1) * length(subsample_index) 
      Lambda_0_outer_m <- Lambdas_0[[m]][subsample_index,] %*% t(As[[m]]) %*% As[[v]] %*% t(Lambdas_0[[v]][subsample_index,])
      l_ci <- fable_l[(idx_start_1 + subsample_index), (idx_start_2 + subsample_index)]
      u_ci <- fable_u[(idx_start_1 + subsample_index), (idx_start_2 + subsample_index)]
      cov_posterior_biview <- c(cov_posterior_biview, mean((l_ci<Lambda_0_outer_m) & (u_ci>Lambda_0_outer_m)))
      print(cov_posterior_biview[length(cov_posterior_biview)])
      len_posterior_biview <- c(len_posterior_biview, mean(u_ci - l_ci))
      
    }
    idx_start_1 <- idx_start_1 + length(subsample_index) 
    
  }
  
  
  return(list(
    rmse_all = rmse_all,
    rmses_intra = rmses_view,
    rmses_inter = rmses_biview,
    cov_posterior_intra = cov_posterior_view,
    cov_posterior_inter = cov_posterior_biview,
    len_posterior_intra = len_posterior_view,
    len_posterior_inter = len_posterior_biview
    #time_model_uq = time_uq[3]
    
  )
  )
}



fit_rotate <- function(Y_c, k=10, lambda0=5){
  p <- ncol(Y_c)
  startB <- matrix(rnorm(p*k), p, k)
  alpha <- 1/p
  lambda1 <- 0.001
  epsilon <- 0.05
  K <- k
  start <- list(B=startB, sigma=rep(1,k), theta=rep(0.5,k))
  rotate_fit <- FACTOR_ROTATE(
    Y_c, lambda0, lambda1, start, K, epsilon, alpha,TRUE,TRUE,100,TRUE
    )
  rotate_fit$FABLEPostMean <- tcrossprod(rotate_fit$B)
  return(rotate_fit)
}



