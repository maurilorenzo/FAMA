if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
#BiocManager::install("MOFA2")
setwd("C:/Users/mauri/Desktop/projects/FAMA_all/FAMA_code")
library(MOFA2)

library(LaplacesDemon)
library(MASS)
library(matrixStats)
library(readxl)
library(truncnorm)

k <- 30
k_1 <- 10
S <- 4
k_tildes <- rep(k_1, S)
p_s <- c(1000, 1000, 500, 100)
p_s <- rep(1000, S)
sum(p_s)
n <- 300

sigmas <- c(1, 1 , 0.2, 0.2)
sigmas_l <- c(0.5, 0.5, 0.5, 0.5)
sigma_u <- c(5, 5, 5, 5)

sigmas <- rep(0.5, S)
sigmas_l <- rep(0.5, S)
sigma_u <- rep(2, S)



k <- 25
k_1 <- 20
S <- 2
k_tildes <- rep(k_1, S)
p_s <- c(1000, 200)
sum(p_s)
n <- 500

sigmas <- c(0.5, 0.2)
sigmas_l <- c(0.5, 0,5)
sigma_u <- c(5, 5)


set.seed(123)
Eta_0 <- matrix(rnorm(n*k), ncol=k)
Eta_0_outer <- Eta_0 %*%t(Eta_0)


Phis_0 <- list()
Phis_0_outer <- list()
Etas_tilde_outer <- list()

Gammas_0 <- list()
Gammas_0_outer <- list()
Lambdas_0<- list()
Lambdas_0_outer<- list()
Sigmas_0 <- list()
LRs_0<-  list()
Thetas_0 <- list()
Y <- list()
Lambdas_tilde_outer <- list()
Ms_tilde_outer <- list()
Ms_tilde <-list() 
As <- list()


for(s in 1:S){
  set.seed(s)
  print(s)
  As[[s]] = matrix(0, k, k_tildes[s])
  index_ones <- sample(1:k, k_tildes[s])
  for(t in 1:k_tildes[s]){
    As[[s]][index_ones[t], t] <- 1 
  }
  Lambdas_0[[s]] <-matrix(rnorm(p_s[s]*k_tildes[s], 0, sigmas[s]), ncol=k_tildes[s])
  Lambdas_0_outer[[s]] <-  Lambdas_0[[s]]  %*% t(Lambdas_0[[s]])
  Sigmas_0[[s]] <- diag(runif(p_s[s], 0.5, 2), ncol=p_s[s], nrow=p_s[s])
  Y[[s]] <- Eta_0 %*% As[[s]] %*% t(Lambdas_0[[s]])+ #  Phis_0[[s]] %*% t(Gammas_0[[s]]) + 
    mvrnorm(n, rep(0,p_s[s]), Sigmas_0[[s]])
  s_Y_s <- svd(Y[[s]])
  k_hat <- estimate_latent_dimension(Y[[s]], 30)$k_hat
  print(k_hat)
  M_s <- s_Y_s$u[, 1:k_hat] * sqrt(n)
  Ms_tilde[[s]] <-  M_s
  M_s_outer <- M_s%*%t(M_s)
  Ms_tilde_outer[[s]] <- M_s_outer
}

M_tilde_mean <- rep(0, n,n)
for(s in 1:S){
  M_tilde_mean <- M_tilde_mean + Ms_tilde_outer[[s]]/n
}
M_tilde_mean <- M_tilde_mean/S
s_M_tilde <- svd(M_tilde_mean)
plot(s_M_tilde$d[1:60], type='l', xlim=c(0,k+k_1+10))
abline(v=k+k_1);abline(v=k)
abline(h = (1/S-0.1))
abline(h = (0.5/S), col='red')

#plot(s_M_tilde$d[-min(c(n,p))]/s_M_tilde$d[-1], type='l', xlim=c(0,k+k_1+10))
max(which(s_M_tilde$d >(0.5/S)))

k_hat <- max(which(s_M_tilde$d > (0.5/S)))
abline(h = 1/S-1/S/2)
abline(v = k_hat)
k_hat
#k_hat <- 19
M_tilde <- s_M_tilde$u[,1:k_hat] * sqrt(n)
M_tilde_outer <- M_tilde %*% t(M_tilde)
#norm((Eta_0_outer -M_tilde_outer), type='2')/norm(Eta_0_outer, type='2')
norm((Eta_0_outer -M_tilde_outer), type='F')/norm(Eta_0_outer, type='F')

compute_L_js <- function(Y, U){
  R <- solve(t(U) %*% U) %*% t(U) %*% Y 
  return(colSums(R^2))
}

compute_V_js <- function(Y, U){
  R <- Y - U%*%solve(t(U) %*% U) %*%t(U) %*% Y
  return(colMeans(R^2))
}

Lambdas_hat <- list()
for(s in 1:S){
  L_js_s <- compute_L_js(Y[[s]], M_tilde);
  V_js_s <- compute_V_js(Y[[s]], M_tilde)
  tau_2_hat_s <- 1/(k_hat*p_s[s]) * sum(L_js_s/V_js_s)
  Lambdas_hat[[s]] <- t(Y[[s]]) %*% M_tilde / (n + 1/tau_2_hat_s)
  #Lambdas_hat[[s]] <- t(Y[[s]]) %*% M_tilde / n
}

Lambda_hat <- rbind(Lambdas_hat[[1]], Lambdas_hat[[2]], Lambdas_hat[[3]],
                    Lambdas_hat[[4]])
#Lambda_hat <- rbind(Lambdas_hat[[1]], Lambdas_hat[[2]], Lambdas_hat[[3]])
dim(Lambda_hat)

Lambda_hat_outer_1 <- tcrossprod(Lambda_hat)

Lambda_0_c <- matrix(0, sum(p_s), k)
idx_start <-1 
for(s in 1:S){
  Lambda_0_c[idx_start:(idx_start + p_s[s]-1), ] <- Lambdas_0[[s]] %*% t(As[[s]])
  idx_start <- idx_start + p_s[s]
}

Lambdas_0_outer <- tcrossprod(Lambda_0_c)
norm((Lambda_hat_outer_1 - Lambdas_0_outer), type='F')/norm(Lambdas_0_outer, type='F')

for(s in 1:S){
  print(s)
  Lambda_hat_outer_s <- tcrossprod(Lambdas_hat[[s]])
  Lambda_0_outer_s <- tcrossprod(Lambdas_0[[s]])
  print( norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F'))
  for(v in (s+1):S) {
    print(v)
    Lambda_hat_outer_sv <- Lambdas_hat[[s]] %*% t(Lambdas_hat[[v]])
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}


# c-fable
Y_c <- do.call(cbind, Y)
s_Y_c <- svd(Y_c, nu = k_hat, nv=k_hat)

M_hat_fable <- s_Y_c$u[, 1:k_hat] * sqrt(n)
L_c <- compute_L_js(Y_c, M_tilde);
V_c <- compute_V_js(Y_c, M_tilde)
tau_2_hat_s <- 1/(k_hat*sum(p_s)) * sum(L_c/V_c)
Lambda_hat_fable <- t(Y_c) %*% M_hat_fable / (n + 1/tau_2_hat_s)


fama_loadings_c <- do.call(rbind, Lambdas_hat)

#print( norm( ( tcrossprod(mofa_loadings_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))
print( norm( ( tcrossprod(fama_loadings_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))
print( norm( ( tcrossprod(Lambda_hat_fable) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))


for(s in 1:S){
  print(s)
  Lambda_hat_outer_s <- tcrossprod(Lambdas_hat[[s]])
  Lambda_0_outer_s <- tcrossprod(Lambdas_0[[s]])
  print( norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F'))
}

idx_start <-1 
for(s in 1:S){
  print(s)
  Lambda_hat_fable_v <- Lambda_hat_fable[idx_start:(idx_start + p_s[s]-1), ]
  idx_start <- idx_start + p_s[s]
  Lambda_hat_outer_s <- tcrossprod(Lambda_hat_fable_v)
  Lambda_0_outer_s <- tcrossprod(Lambdas_0[[s]])
  print( norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F'))
}


for(s in 1:(S-1)){
  print(s)
  for(v in (s+1):S) {
    print(v)
    Lambda_hat_outer_sv <- Lambdas_hat[[s]] %*% t(Lambdas_hat[[v]])
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}

idx_start <- 1
for(s in 1:(S-1)){
  print(s)
  Lambda_hat_fable_s <- Lambda_hat_fable[idx_start:(idx_start + p_s[s]-1), ]
  idx_start <- idx_start + p_s[s]
  for(v in (s+1):S) {
    print(v)
    idx_start_2 <- 1 + sum(p_s[1:(v-1)])
    Lambda_hat_fable_v <- Lambda_hat_fable[idx_start_2:(idx_start_2 + p_s[v]-1), ]
    Lambda_hat_outer_sv <- Lambda_hat_fable_s %*% t(Lambda_hat_fable_v)
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}



cat("
        #########################################################
        ###           ____         __  __                     ###
        ###          | ___| /\\    |  \\/  |    /\\              ###               
        ###          | |_  /  \\   | \\  / |   /  \\             ###
        ###          | __|/ /\\ \\  | |\\/| |  / /\\ \\            ###
        ###          | | / ____ \\ | |  | | / ____ \\           ###
        ###          |_|/ /    \\_\\|_|  |_|/_/    \\_\\          ###
        ###                                                   ### 
        ######################################################### 
")



cat("
        #########################################################
        ###           __  __  ____  ______                    ### 
        ###          |  \\/  |/ __ \\|  ___ /\\    _             ### 
        ###          | \\  / | |  | | |__ /  \\ _| |_           ### 
        ###          | |\\/| | |  | |  __/ /\\ \\_   _|          ###
        ###          | |  | |  |__| | |  / ____ \\|_|            ###
        ###          |_|  |_|\\____/|_|/_/    \\_\\              ###
        ###                                                   ### 
        ####
")

################################################################################
# MOFA ####
################################################################################

data <- make_example_data(
  n_views = 2, 
  n_samples = 200, 
  n_features = 500, 
  n_factors = 10
)[[1]]

lapply(data, dim)


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

Y_mofa <- prepare_data_4_mofa(Y)
MOFAobject <- create_mofa(Y_mofa)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)
model_opts$num_factors = fit_fama_1$k_0_hat
model_opts$spikeslab_factors = T
model_opts$ard_factors = F # if TRUE accuracy decreases
model_opts$spikeslab_weights = T
model_opts$ard_weights = T

train_opts <- get_default_training_options(MOFAobject)
train_opts
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# C:\Users\mauri\anaconda3
library(reticulate)
use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python312", required=TRUE)
outfile = file.path('results',"mofa_fit_1.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

#filepath <- system.file("extdata", "model.hdf5", package = "MOFA2")
#print(filepath)
#model <- load_model(filepath)
#view_names(MOFAobject.trained)
#factors <- get_factors(MOFAobject.trained, factors = "all")
#lapply(factors,dim)
mofa_weights_1 <- get_weights(MOFAobject.trained, views = "all", factors = "all")
mofa_weights_2 <- get_expectations(MOFAobject.trained, "W")

mofa_scales <- get_scales(MOFAobject.trained)


factor_means <- get_factors(MOFAobject.trained, factors = "all", as.data.frame = FALSE)

# Extract factor variances from internal slot
factor_vars <- MOFAobject.trained@expectations$factors$var  # list: one matrix per group
factor_names <- names(factor_vars)


lapply(mofa_weights_1,dim)
lapply(mofa_weights_2,dim)

#plot_data_overview(MOFAobject.trained)

mofa_loadings_1_c <- do.call(rbind, mofa_weights_1)
mofa_loadings_2_c <- do.call(rbind, mofa_weights_2)


print( norm( ( tcrossprod(mofa_loadings_1_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))
print( norm( ( tcrossprod(mofa_loadings_2_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))
print( norm( ( tcrossprod(fama_loadings_c) - tcrossprod(Lambda_0_c)), type='F')/norm(tcrossprod(Lambda_0_c), type='F'))


for(s in 1:S){
  print(s)
  Lambda_hat_outer_s <- tcrossprod(mofa_weights_1[[paste0('view_', s)]])
  Lambda_0_outer_s <- tcrossprod(Lambdas_0[[s]] %*% t(As[[s]]))
  print( norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F'))
}

for(s in 1:S){
  print(s)
  Lambda_hat_outer_s <- tcrossprod(mofa_weights_2[[paste0('view_', s)]])
  Lambda_0_outer_s <- tcrossprod(Lambdas_0[[s]] %*% t(As[[s]] ))
  print( norm( ( Lambda_hat_outer_s - Lambda_0_outer_s), type='F')/norm(Lambda_0_outer_s, type='F'))
}



for(s in 1:(S-1)){
  print(s)
  for(v in (s+1):S) {
    print(v)
    Lambda_hat_outer_sv <- mofa_weights_1[[paste0('view_', s)]]  %*% t(mofa_weights_1[[paste0('view_', v)]])
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}

for(s in 1:(S-1)){
  print(s)
  for(v in (s+1):S) {
    print(v)
    Lambda_hat_outer_sv <- mofa_weights_2[[paste0('view_', s)]]  %*% t(mofa_weights_2[[paste0('view_', v)]])
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}

for(s in 1:(S-1)){
  print(s)
  for(v in (s+1):S) {
    print(v)
    Lambda_hat_outer_sv <- Lambdas_hat[[s]] %*% t(Lambdas_hat[[v]])
    Lambda_0_outer_sv <- Lambdas_0[[s]] %*% t(As[[s]]) %*% As[[v]] %*% t(Lambdas_0[[v]])
    print(norm( (Lambda_hat_outer_sv - Lambda_0_outer_sv), type='F')/norm(Lambda_0_outer_sv, type='F'))
  }
}

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


library(reticulate)


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



mofa_fit_1 <- fit_MOFA(Y, fit_fama_1$k_0_hat)
mofa_fit_1_metrics <- compute_performance_fable(mofa_fit_1, Lambdas_0, As)
mofa_fit_1_metrics

################################################################################
# fns ####
################################################################################

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

estimate_latent_dimension <- function(Y, k_max){
  svd_Y <- svd(Y)
  n <- nrow(Y)
  jics <- sapply(1:k_max, function(x) (compute_jic(Y, svd_Y, x)))
  plot(1:k_max, jics, type='l', xlab='k', ylab='jic', main='')
  print(paste('k_hat = ', which.min(jics)))
  return(list(k_hat = which.min(jics), jics=jics, svd_Y = svd_Y))
}

compute_L_js <- function(Y, U){
  R <- solve(t(U) %*% U) %*% t(U) %*% Y 
  return(colSums(R^2))
}

compute_V_js <- function(Y, U){
  R <- Y - U%*%solve(t(U) %*% U) %*%t(U) %*% Y
  return(colMeans(R^2))
}

