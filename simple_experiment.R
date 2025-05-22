


source('simulations/helpers_simulations.R')
source('FAMA_wrapper.R')
#source('competitors/FACTOR_ANALYSIS/FACTOR_CODE_update.R')

library(LaplacesDemon)
library(MASS)
library(matrixStats)
library(readxl)
library(truncnorm)

k <- 40
k_1 <- 20
S <- 4
k_tildes <- rep(k_1, S)
p_s <- c(1000, 1000, 1000, 1000)
p_s <- rep(5000, S)
sum(p_s)
n <- 7500

sigmas <- c(1, 1 , 0.2, 0.2)
sigmas_l <- c(0.5, 0.5, 0.5, 0.5)
sigma_u <- c(5, 5, 5, 5)

sigmas <- rep(0.5, S)
sigmas_l <- rep(0.5, S)
sigma_u <- rep(5, S)



k <- 50
k_1 <- 40
S <- 2
k_tildes <- rep(k_1, S)
p_s <- c(10000, 1000)
sum(p_s)
n <- 500

sigmas <- c(1, 0.2)
sigmas_l <- c(0.5, 0,5)
sigma_u <- c(5, 5)


k <- 70
k_1 <- 30
S <- 3
k_tildes <- rep(k_1, S)
p_s <- c(10000, 5000, 2000)
n <- 4000

sigmas <- c(1, 0.25, 0.25)
sigmas_l <- c(0.5, 0.5 ,0.5)
sigma_u <- c(5, 5, 5, 5)


k <- 30
k_1 <- 20
S <- 4
k_tildes <- rep(k_1, S)
p_s <- c(1000, 500, 200, 200)
p_s <- c(4000, 3000, 1000, 1000)

sum(p_s)
n <- 1000


sigmas <- c(1, 0.4, 0.4, 0.4)
sigmas_l <- c(0.5, 0.5 ,0.5, 0.5)*10
sigmas_u <- c(5, 5, 5, 5)*2
n <- 2000
p_s <- c(5000, 1000, 1000, 1000)

p_s <- rep(2000, S)
sigmas <- rev(sigmas)

sigmas <- rep(0.5, S)
n <- 500


sigmas <- c(0.5, 0.2)

p_s <- rep(2000, S)


S = 4

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


#As[[m]] k times k_tilde_m

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
  Sigmas_0[[s]] <- runif(p_s[s], sigmas_l[s], sigmas_u[s])
  boxplot(diag(Lambdas_0_outer[[s]]) / Sigmas_0[[s]])
  
  Y[[s]] <- Eta_0 %*% As[[s]] %*% t(Lambdas_0[[s]])+ #  Phis_0[[s]] %*% t(Gammas_0[[s]]) + 
    + matrix(rnorm(p_s[s] * n, 0 , 1), nrow=n) %*%  diag(sqrt(Sigmas_0[[s]]))
    #mvrnorm(n, rep(0,p_s[s]), Sigmas_0[[s]])
}

M <- S
mask <- rep(0, k)
for(m in 1:M){
  mask <- mask + colSums(As[[m]] %*% t(As[[m]]))
}
k_0 <- sum(mask > 0)
k_0


subsample_index <- 1:100
fit_fama_1 <- fit_FAMA(Y, clt_SE=T, posterior_SE=T, 
                       index_SE=subsample_index)
fama_fit_1_metrics <- compute_performance_fama(
  fit_fama_1, Lambdas_0, As)
fama_fit_1_metrics
fama_fit_1_metrics$rmse_all
fama_fit_1_metrics$rmses_intra
fama_fit_1_metrics$rmses_inter
fit_fama_1$time_model_estimates

Y_c <- do.call(cbind, Y)
ptm <- proc.time() 
fable_fit_1 <- PseudoPosteriorMean(Y_c, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
fable_time <- proc.time() - ptm
fable_fit_1$estRank
fable_fit_1_metrics <- compute_performance_fable(
  fable_fit_1, Y, Lambdas_0, As, n_MC=500, subsample_index=subsample_index)
fable_fit_1_metrics
fable_fit_1_metrics$rmse_all
fable_fit_1_metrics$rmses_intra
fable_fit_1_metrics$rmses_inter

fama_fit_1_metrics$rmse_all
fama_fit_1_metrics$rmses_intra
fama_fit_1_metrics$rmses_inter

K <- k + 5
p <- sum(p_s[1:S])
startB <- matrix(rnorm(p*K),p,K)
alpha <- 1/p
lambda1 <- 0.001
epsilon <- 0.05
Y_c <- do.call(cbind, Y)
start <- list(B=startB, sigma=rep(1,K), theta=rep(0.5,K))
lambda0<-5
ptm <- proc.time()
rotate_fit_5 <- FACTOR_ROTATE(Y_c,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
rotate_time <- proc.time() - ptm 
rotate_time[3]
rotate_fit_5$FABLEPostMean <- tcrossprod(rotate_fit_5$B)
rotate_fit_5_metrics <- compute_performance_fable(rotate_fit_5, Y, Lambdas_0, As, UQ=F)
rotate_fit_5_metrics
fama_fit_1_metrics$rmse_all
fama_fit_1_metrics$rmses_intra
fama_fit_1_metrics$rmses_inter

mean(fama_fit_1_metrics$rmses_inter)
mean(rotate_fit_5_metrics$rmses_inter)

n
dev.new()


ptm <- proc.time() 
mofa_fit_1 <- fit_MOFA(Y, k_0 + 5)
mofa_time <- proc.time() - ptm
mofa_fit_1_metrics <- compute_performance_mofa(
  mofa_fit_1, Lambdas_0, As)
mofa_fit_1_metrics


fama_fit_1_metrics$rmse_all
fable_fit_1_metrics$rmse_all
mofa_fit_1_metrics$rmse_all

fama_fit_1_metrics$rmses_intra
fable_fit_1_metrics$rmses_intra
mofa_fit_1_metrics$rmses_intra

fama_fit_1_metrics$rmses_inter
fable_fit_1_metrics$rmses_inter
mofa_fit_1_metrics$rmses_inter


mean(fama_fit_1_metrics$rmses_inter)
mean(fable_fit_1_metrics$rmses_inter)
mean(mofa_fit_1_metrics$rmses_inter)

dev.new()


# INITIALIZATIONS

K <- k
p <- sum(p_s)
startB <- matrix(rnorm(p*K),p,K)
alpha <- 1/p
lambda1 <- 0.001
epsilon <- 0.05
myImagePlot(abs(startB),F)
title("Initialization")

# PXEM: Dynamic Posterior Exploration (Approximate M-step)
Y_c <- do.call(cbind, Y)

start <- list(B=startB, sigma=rep(1,K), theta=rep(0.5,K))

lambda0<-5
ptm <- proc.time()
rotate_fit_5 <- FACTOR_ROTATE(Y_c,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
rotate_time <- proc.time() - ptm 
rotate_time[3]
rotate_fit_5$FABLEPostMean <- tcrossprod(rotate_fit_5$B)
rotate_fit_5_metrics <- compute_performance_fable(rotate_fit_5, Y, Lambdas_0, As, UQ=F)
rotate_fit_5_metrics
fama_fit_1_metrics$rmse_all
fama_fit_1_metrics$rmses_intra
fama_fit_1_metrics$rmses_inter

dev.new()

lambda0<-1
rotate_fit_1 <-FACTOR_ROTATE(Y_c,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
rotate_fit_1$FABLEPostMean <- tcrossprod(rotate_fit_1$B)
rotate_fit_1_metrics <- compute_performance_fable(rotate_fit_1, Y, Lambdas_0, As, UQ=F)
rotate_fit_1_metrics
fama_fit_1_metrics$rmse_all
fama_fit_1_metrics$rmses_intra
fama_fit_1_metrics$rmses_inter



lambda0<-10
rotate_fit_10 <-FACTOR_ROTATE(Y_c,lambda0,lambda1,start,K,epsilon,alpha,TRUE,TRUE,100,TRUE)
rotate_fit_10$FABLEPostMean <- tcrossprod(rotate_fit_10$B)
rotate_fit_10_metrics <- compute_performance_fable(rotate_fit_10, Y, Lambdas_0, As, UQ=F)

lambda0<-20
result_20<-FACTOR_ROTATE(Y,lambda0,lambda1,result_10,K,epsilon,alpha,TRUE,TRUE,100,TRUE)

lambda0<-30
result_30<-FACTOR_ROTATE(Y,lambda0,lambda1,result_20,K,epsilon,alpha,TRUE,TRUE,100,TRUE)


