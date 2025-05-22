
library(Rcpp)
library(RcppArmadillo)

sourceCpp('helpers_functions.cpp')



estimate_latent_dimension <- function(Y, k_max){
  svd_Y <- compute_svd(Y, k_max)
  n <- nrow(Y)
  jics <- sapply(1:k_max, function(x) (compute_jic(Y, svd_Y, x)))
  plot(1:k_max, jics, type='l', xlab='k', ylab='jic', main='')
  print(paste('k_hat = ', which.min(jics)))
  return(list(k_hat = which.min(jics), jics=jics, svd_Y = svd_Y))
}

compute_jic <- function(Y, svd_Y, k){
  n <- nrow(Y); p <- ncol(Y) 
  minint <- min(n ,p)
  maxint <- max(n, p)
  M <- sqrt(n)*as.matrix(svd_Y$u[,1:k])
  Lambda <- 1/sqrt(n)* as.matrix(svd_Y$v[,1:k]) %*% diag(svd_Y$d[1:k], k, k)
  
  Y_hat <- tcrossprod(M, Lambda)
  sigma_sq_hat <- colMeans((Y - Y_hat)^2) # p * 1
  tausq_est <- (mean(colSums((Y_hat)^2) / sigma_sq_hat)) / (n * k);
  
  Lambda_est <- (sqrt(n) / (n + 1/tausq_est)) * svd_Y$v[,1:k] %*% diag(svd_Y$d[1:k], k, k);
  M_Lambda_est <- sqrt(n)* svd_Y$u[,1:k] %*% t(Lambda_est)
  res <- colSums((Y - M_Lambda_est)^2) 
  sigma_sq_hat <- res / n
  
  term1 <- -n*sum(log(sqrt(sigma_sq_hat))) - 0.5*sum(res/sigma_sq_hat)
  term1 <- term1 / (n*p)
  #print(k)
  jic <- -2*term1 + (k * maxint * log(minint) / (n*p))
  #print(jic)
  return(jic)
}



compute_svd <- function(Y, k=10, randomized_svd=F){
  n <- nrow(Y); p <- ncol(Y)
  if (n > p) {
    YtY <- t(Y) %*% Y
    if(randomized_svd){s_Y <- rsvd(YtY, k=k, nu=k, nv=k, p = 10, q = 2, sdist = "normal")}
    else{s_Y <- svd(YtY, nv=k, nu=k)}
    V <- s_Y$u[,1:k, drop=FALSE]
    D <- s_Y$d[1:k, drop=FALSE]^(1/2)
    U <- Y %*% V %*% diag(as.vector(s_Y$d[1:k, drop=FALSE]^(-1/2)))
  }
  else {
    if(randomized_svd){s_Y <- rsvd(Y, k=k, nu=k, nv=k, p = 10, q = 2, sdist = "normal")}
    else{s_Y <- svd(Y, nv=k, nu=k)}
    U <- s_Y$u[,1:k, drop=FALSE]
    D <- s_Y$d[1:k, drop=FALSE]
    V <- s_Y$v[,1:k, drop=FALSE]
  }
  return(list(u=U, d=D, v=V))
}



estimate_latent_factors_all_views <- function(Y, ks=NULL, k_max=50){
  M <- length(Y)
  n <- nrow(Y[[1]])
  F_tildes <- list()
  Ps <- list() 
  tau_sqs <- c()
  
  estimate_factors = F
  if(is.null(ks)){
    ks=c()
    estimate_factors = T
    }
  for(m in 1:M){
    print(paste0('estimating latent factors for view ', m))
    if(estimate_factors){
      est_m <- estimate_latent_dimension(Y[[m]], k_max)
      ks[m] <- est_m$k_hat
      s_Y_m <- est_m$svd_Y
      print(ks[m])
    }
    else{
      s_Y_m <- svd(Y[[m]])
    }
    print(paste0('latent factors for view ', m, ' = ', ks[m]))
    U_m = s_Y_m$u[, 1:ks[m], drop=FALSE]
    F_tildes[[m]] = U_m * sqrt(n)
    Ps[[m]] = tcrossprod(U_m)
    L_m  = sum((crossprod(U_m, Y[[m]]))^2) / n
    V_m = sum(colMeans(((diag(n) - Ps[[m]]) %*% Y[[m]] )^2))
    tau_sqs[m] <- L_m / (V_m * ks[m])
  }
  return(list(F_tildes=F_tildes, Ps=Ps, tau_sqs=tau_sqs))
}


compute_F_hat <- function(Ps, k_0=NA, tau_1 = 0.2){
  P_tilde <- Reduce("+", Ps) / length(Ps)
  s_P <- eigen(P_tilde, symmetric=T)
  if(is.na(k_0)){
    k_0 <- sum(s_P$values > tau_1)
  }
  F_hat <- sqrt(nrow(P_tilde)) * s_P$vectors[,1:k_0, drop=FALSE] 
  P_U <- tcrossprod(F_hat) / nrow(P_tilde)
  return(list(k_0 = k_0, F_hat = F_hat, P_U = P_U))
}

compute_posterior_mean_Lambda_m <- function(Y_m, tau_sq_m, F_hat){
  n <- nrow(Y_m)
  mu_m <- 1 / (n + 1/tau_sq_m) * crossprod(Y_m, F_hat)
  return(mu_m)
}

compute_clt_biview <- function(
    Y_1, Y_2, Lambda_1, Lambda_2, sigmas_sq_1, sigmas_sq_2, P_1, P_2, k_12_c = NA, tau_2=0.25){
  P_bar <- 0.5 * (P_1 + P_2)
  n <- nrow(Y_1)
  #s_P <- svd(P_bar)
  #if(is.na(k_12_c)){
  #  k_12_c <- sum( s_P$d>( 1-tau_2 ) )
  #  print(k_12_c)
  #}
  #F_hat <- sqrt(n) * s_P$u[,1:k_12_c, drop=FALSE] 
  #Lambda_1_c <- 1/n * t(Y_1) %*% F_hat
  #Lambda_2_c <- 1/n * t(Y_2) %*% F_hat
  clt_SE_12 <-  compute_sds_clt_biview(
    #Lambda_1, Lambda_2, Lambda_1_c, Lambda_2_c, sigmas_sq_1, sigmas_sq_2)
    Lambda_1, Lambda_2, Lambda_1, Lambda_2, sigmas_sq_1, sigmas_sq_2)
  return(list(clt_SE = clt_SE_12 ))
  #Lambda_1_hat = Lambda_1_c, F_hat = F_hat,
         #Lambda_2_hat = Lambda_2_c))
}


compute_posterior_mean_sigmas_sq_m <- function(
    Y_m, tau_sq_m, mu_m, v_0=1, sigma_sq_0=1){
  n <- nrow(Y_m)
  v_n <- n + v_0 
  SSTs <- apply(Y_m, 2, function(x) sum(x^2))
  SSEs <- apply(mu_m, 1, function(x) sum(x^2)) * (n + 1/tau_sq_m)
  sigmas_sq_m_est <- (v_0*sigma_sq_0 + SSTs - SSEs) / v_n
  return(sigmas_sq_m_est)
}

compute_rho <- function(Lambda_hat, sigmas_sq_hat){
  B <- compute_B(Lambda_hat, sigmas_sq_hat)
  rho <- mean(B[lower.tri(B, diag=T)])
  return(rho)
}

fit_FAMA <- function(
    Y, ks=NULL, k_0 = NA, k_max=50, sigma_sq_0 = 1, v_0 = 1, tau_1 = NA, 
    tau_2 = 0.25, clt_SE = FALSE, posterior_SE = FALSE, index_SE = 1:100){
  
  M <- length(Y)
  n <- nrow(Y[[1]])
  ps <- sapply(Y, function(x) ncol(x))
  
  ptm <- proc.time()
  factors_est = estimate_latent_factors_all_views(Y, ks, k_max)
  tau_sqs <- factors_est$tau_sqs
  if(is.na(tau_1) | tau_1 > 1 | tau_1 < 0){
    tau_1 <- 0.5 / length(Y)
  }
  overall_factors_est = compute_F_hat(factors_est$Ps, k_0, tau_1)
  k_0_hat <- overall_factors_est$k_0
  
  print(paste0('Overall latent factors = ', k_0_hat))
  
  F_hat <- overall_factors_est$F_hat
  P_U <- overall_factors_est$P_U
  Ps <- factors_est$Ps
  
  
  Lambdas_hat<-  Map(
    function(y_m, tau_m) {
      compute_posterior_mean_Lambda_m(y_m, tau_m, F_hat)
    },
    Y, tau_sqs
  )
  #sigmas_sq_hat <- mapply(
  #  compute_posterior_mean_sigmas_sq_m, Y, tau_sqs, Lambdas_hat,
  #  MoreArgs = list(v_0 = v_0, sigma_sq_0 = sigma_sq_0),
  #  SIMPLIFY = FALSE) 
  sigmas_sq_hat<-  Map(
    function(y_m, tau_m, mu_m) {
      compute_posterior_mean_sigmas_sq_m(y_m, tau_m, mu_m, v_0, sigma_sq_0)
    },
    Y, tau_sqs, Lambdas_hat
  )
  time_model_estimates <- proc.time() - ptm; 
  
  out = list(Lambdas_hat = Lambdas_hat,  sigmas_sq_hat = sigmas_sq_hat)
  out$k_0_hat = k_0_hat
  out$time_model_estimates = time_model_estimates[3]
  if((! clt_SE) & (! posterior_SE) ){
    return(out)
  }
  
  if(is.vector(index_SE)){
    index_SE <- replicate(M, index_SE, simplify = FALSE)
  }
  for(m in 1:M){
    if(max(index_SE[[m]])>ps[m]){
      index_SE[[m]] = 1:ps[m]
    }
  }
  
  
  clt_SE_intraview = c()
  clt_SE_interview = c()
  if(clt_SE){
    for(m in 1:(M-1)){
      clt_SE_intraview[[m]] = compute_sds_clt_intraview(
        Lambdas_hat[[m]][index_SE[[m]],], sigmas_sq_hat[[m]][index_SE[[m]]])
      clt_SE_interview[[m]] = list()
      for(l in (m+1):M){
        clt_SE_interview[[m]][[l]] = compute_clt_biview(
          Y[[m]][,index_SE[[m]]], Y[[l]][,index_SE[[l]]], Lambdas_hat[[m]][index_SE[[m]],], 
          Lambdas_hat[[l]][index_SE[[l]],], sigmas_sq_hat[[m]][index_SE[[m]]],
          sigmas_sq_hat[[l]][index_SE[[l]]], Ps[[m]], Ps[[l]], k_12_c=NA, tau_2=tau_2) 
      }
    }
    clt_SE_intraview[[M]] = compute_sds_clt_intraview(
      Lambdas_hat[[M]][index_SE[[m]],], sigmas_sq_hat[[M]][index_SE[[m]]])
    out$clt_SE_intraview = clt_SE_intraview
    out$clt_SE_interview = clt_SE_interview
  }
  
  
  posterior_SE_intraview = c()
  posterior_SE_interview = c()
  if(posterior_SE){
    rhos <- c()
    for(m in 1:M){
      rhos[m] = compute_rho(Lambdas_hat[[m]][index_SE[[m]],], 
                            sigmas_sq_hat[[m]][index_SE[[m]]])
    }
    out$rhos = rhos
    for(m in 1:(M-1)){
      posterior_SE_intraview[[m]] = compute_sds_posterior_distribution_intraview(
        Lambdas_hat[[m]][index_SE[[m]],], sigmas_sq_hat[[m]][index_SE[[m]]], rhos[m])
      posterior_SE_interview[[m]] = list()
      for(l in (m+1):M){
        posterior_SE_interview[[m]][[l]] =  compute_sds_posterior_distribution_biview(
          Lambdas_hat[[m]][index_SE[[m]],], Lambdas_hat[[l]][index_SE[[l]],], 
          sigmas_sq_hat[[m]][index_SE[[m]]], sigmas_sq_hat[[l]][index_SE[[l]]],
          rhos[m], rhos[l])
      }
    }
    posterior_SE_intraview[[M]] = compute_sds_posterior_distribution_intraview(
      Lambdas_hat[[M]][index_SE[[m]],], sigmas_sq_hat[[M]][index_SE[[m]]])
    out$posterior_SE_intraview = posterior_SE_intraview
    out$posterior_SE_interview = posterior_SE_interview
  }
  time_uq <- proc.time() - ptm; 
  out$time_model_estimates_and_uq = time_uq[3]
  
  
  return(out)
}

compute_CI_normal_approx <- function(Lambda_1, Lambda_2, SEs, alpha=0.05){
  Lambda_outer <- tcrossprod(Lambda_1, Lambda_2)
  if(sum(dim(Lambda_outer) == dim(SEs)) < 2){
    print('dimensions do not match')
    return()
  }
  dev <- SEs * qnorm(1-alpha/2)
  u_CI <- Lambda_outer + dev
  l_CI <- Lambda_outer - dev
  return(list(u_CI = u_CI, l_CI = l_CI))
} 

compute_performance_fama <- function(fit, Lambdas_0, As){
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
  
  cov_clt_view <- c()
  cov_clt_biview <- c()
  len_clt_view <- c()
  len_clt_biview <- c()
  print('coverage clt intra-view')
  for(m in 1:M){
    print(m)
    Lambda_0_outer_m <- tcrossprod(Lambdas_0[[m]][subsample_index,])
    ci_m <- compute_CI_normal_approx(
      fit$Lambdas_hat[[m]][subsample_index,], fit$Lambdas_hat[[m]][subsample_index,], 
      fit$clt_SE_intraview[[m]] / sqrt(n), alpha=0.05)
    cov_clt_view[m] <- mean((ci_m$l_CI<Lambda_0_outer_m) & (ci_m$u_CI>Lambda_0_outer_m))
    print(cov_clt_view[m])
    len_clt_view[m] <- mean(ci_m$u_CI - ci_m$l_CI)
  }
  print('coverage clt inter-view')
  for(m in 1:(M-1)){
    print(m)
    for(v in (m+1):M) {
      print(v)
      Lambda_0_outer_m <- Lambdas_0[[m]][subsample_index,] %*% t(As[[m]]) %*% As[[v]] %*% t(Lambdas_0[[v]][subsample_index,])
      ci_m <- compute_CI_normal_approx(
        fit$Lambdas_hat[[m]][subsample_index,], fit$Lambdas_hat[[v]][subsample_index,], 
        fit$clt_SE_interview[[m]][[v]]$clt_SE / sqrt(n), alpha=0.05)
      cov_clt_biview <- c(cov_clt_biview, mean((ci_m$l_CI<Lambda_0_outer_m) & 
                                                 (ci_m$u_CI>Lambda_0_outer_m)))
      print(cov_clt_biview[length(cov_clt_biview)])
      len_clt_biview <- c(len_clt_biview, mean(ci_m$u_CI - ci_m$l_CI))
    }
  }
  
  cov_posterior_view <- c()
  cov_posterior_biview <- c()
  len_posterior_view <- c()
  len_posterior_biview <- c()
  print('coverage posterior intra-view')
  for(m in 1:M){
    print(m)
    Lambda_0_outer_m <- tcrossprod(Lambdas_0[[m]][subsample_index,])
    ci_m <- compute_CI_normal_approx(
      fit$Lambdas_hat[[m]][subsample_index,], fit$Lambdas_hat[[m]][subsample_index,], 
      fit$posterior_SE_intraview[[m]] / sqrt(n), alpha=0.05)
    cov_posterior_view[m] <- mean((ci_m$l_CI<Lambda_0_outer_m) & (ci_m$u_CI>Lambda_0_outer_m))
    print(cov_posterior_view[m])
    len_posterior_view[m] <- mean(ci_m$u_CI - ci_m$l_CI)
    
  }
  print('coverage posterior inter-view')
  for(m in 1:(M-1)){
    print(m)
    for(v in (m+1):M) {
      print(v)
      Lambda_0_outer_m <- Lambdas_0[[m]][subsample_index,] %*% t(As[[m]]) %*% As[[v]] %*% t(Lambdas_0[[v]][subsample_index,])
      ci_m <- compute_CI_normal_approx(
        fit$Lambdas_hat[[m]][subsample_index,], fit$Lambdas_hat[[v]][subsample_index,], 
        fit$posterior_SE_interview[[m]][[v]] / sqrt(n), alpha=0.05)
      cov_posterior_biview <- c(cov_posterior_biview, mean((ci_m$l_CI<Lambda_0_outer_m) & (ci_m$u_CI>Lambda_0_outer_m)))
      print(cov_posterior_biview[length(cov_posterior_biview)])
      len_posterior_biview <- c(len_posterior_biview, mean(ci_m$u_CI - ci_m$l_CI))
      
    }
  }
  
  metrics <- list(
    rmse_all = rmse_all,
    rmses_intra = rmses_view,
    rmses_inter = rmses_biview,
    cov_clt_intra = cov_clt_view,
    cov_clt_inter = cov_clt_biview,
    cov_posterior_intra = cov_posterior_view,
    cov_posterior_inter = cov_posterior_biview,
    len_clt_intra = len_clt_view,
    len_clt_inter = len_clt_biview,
    len_posterior_intra = len_posterior_view,
    len_posterior_inter = len_posterior_biview
  )
  return(metrics)
}

