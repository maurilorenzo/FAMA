
source('simulations/helpers_simulations.R')
source('FAMA_wrapper.R')

library(LaplacesDemon)
library(MASS)
library(matrixStats)
library(readxl)
library(truncnorm)


n_sim <- 50

# scenario 1 (unbalanced)
scenario <- 2
if(scenario == 1 | scenario == 2){
  k <- 30
  k_1 <- 15
  M <- 4
  k_tildes <- rep(k_1, M)
  p_s <- rep(2000, M)
  sigmas <- c(1, 0.5, 0.5, 0.5)
  sigmas_l <- rep(5, M)
  sigmas_u <- rep(10, M)
  n <- 500
  if(scenario == 2){
    n <- 1000
  }
}


scenario <- 11
if(scenario == 3 | scenario == 4 | scenario == 11){
  k <- 30
  k_1 <- 20
  M <- 4
  k_tildes <- rep(k_1, M)
  p_s <- c(5000, rep(1000, 3))
  sigmas <- c(1, 0.4, 0.4, 0.4)
  sigmas_l <- rep(5, M)
  sigmas_u <- rep(10, M)
  n <- 500
  if(scenario == 4){
    n <- 1000
  }
  if(scenario == 11){
    n <- 150
  }
}


if(scenario == 5 | scenario == 6){
  k <- 30
  k_1 <- 20
  M <- 4
  k_tildes <- rep(k_1, M)
  p_s <- rep(1000, M)
  sigmas <- rep(0.5, M)
  sigmas_l <- rep(5, M)
  sigmas_u <- rep(10, M)
  n <- 500
  if(scenario == 6){
    n <- 1000
  }
}

scenario <- 9
if(scenario == 7 | scenario == 8 | scenario == 9){
  k <- 30
  k_1 <- 20
  M <- 4
  k_tildes <- rep(k_1, M)
  p_s <- rep(2000, M)
  sigmas <- rep(0.5, M)
  sigmas_l <- rep(5, M)
  sigmas_u <- rep(10, M)
  n <- 500
  if(scenario == 8){
    n <- 1000
  }
  if(scenario == 9){
    n <- 150
  }
}




Gammas_0 <- list()
Gammas_0_outer <- list()
Lambdas_0<- list()
Lambdas_0_outer<- list()
Sigmas_0 <- list()
As <- list()

for(m in 1:M){
  set.seed(99 + m)
  print(m)
  As[[m]] = matrix(0, k, k_tildes[m])
  index_ones <- sample(1:k, k_tildes[m])
  for(t in 1:k_tildes[m]){
    As[[m]][index_ones[t], t] <- 1 
  }
  Lambdas_0[[m]] <-matrix(rnorm(p_s[m]*k_tildes[m], 0, sigmas[m]), 
                          ncol=k_tildes[m])
  Lambdas_0_outer[[m]] <-  Lambdas_0[[m]]  %*% t(Lambdas_0[[m]])
  Sigmas_0[[m]] <- diag(runif(p_s[m], sigmas_l[m], sigmas_u[m]), ncol=p_s[m], nrow=p_s[m])
  
}


mask <- rep(0, k)
for(m in 1:M){
  mask <- mask + colSums(As[[m]] %*% t(As[[m]]))
}
k_0 <- sum(mask > 0)
k_0

set.seed(123)

fama_results <- data.frame()
fable_results <- data.frame()
mofa_results <- data.frame()
rotate_results <- data.frame()

fama <- T; fable <- T; mofa <- T; rotate <- T; 
subsample_index<- 1:100 
for(sim in 1:n_sim){
  set.seed(sim)
  print(sim)
  Eta_0 <- matrix(rnorm(n*k), ncol=k)
  Y <- list()
  Sigmas_0 <- list()
  for(m in 1:M){
    Sigmas_0[[m]] <- runif(p_s[m], sigmas_l[m], sigmas_u[m])
    boxplot(diag(Lambdas_0_outer[[m]]) / Sigmas_0[[m]])
    
    Y[[m]] <- Eta_0 %*% As[[m]] %*% t(Lambdas_0[[m]])+  
      +matrix(rnorm(p_s[m] * n, 0 , 1), nrow=n) %*%  diag(sqrt(Sigmas_0[[m]]))
  }
  
  if(fama){
    set.seed(123)
    fit_fama_1 <- fit_FAMA(Y, clt_SE=T, posterior_SE=T, index_SE=subsample_index)
    fama_fit_1_metrics <- compute_performance_fama(fit_fama_1, Lambdas_0, As)
    fama_results <- rbind(fama_results, c(unlist(fama_fit_1_metrics), fit_fama_1$k_0_hat, fit_fama_1$time_model_estimates, fit_fama_1$time_model_estimates_and_uq))
    names(fama_results) <- c(names(unlist(fama_fit_1_metrics)), 'k_0_hat', 'time_pe', 'time_pe_uq')
  }
  if(fable){
    Y_c <- do.call(cbind, Y)
    set.seed(123)
    ptm <- proc.time() 
    fable_fit_1 <- PseudoPosteriorMean_2(Y_c, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
    fable_time <- proc.time() - ptm
    fable_fit_1_metrics <- compute_performance_fable(fable_fit_1, Y, Lambdas_0, As, n_MC = 500, subsample_index=1:100)
    fable_results <- rbind(fable_results, c(unlist(fable_fit_1_metrics),  fable_fit_1$estRank, fable_time[3]))
    names(fable_results) <- c(names(unlist(fable_fit_1_metrics)), 'k_0_hat', 'time_pe')
  }
  
  if(mofa){
    set.seed(123)
    ptm <- proc.time() 
    mofa_fit_1 <- fit_MOFA(Y, k_0 + 5)
    mofa_time <- proc.time() - ptm
    mofa_fit_1_metrics <- compute_performance_mofa(mofa_fit_1, Lambdas_0, As)
    mofa_fit_1_metrics
    mofa_results <- rbind(mofa_results, c(unlist(mofa_fit_1_metrics),  mofa_time[3]))
    names(mofa_results) <- c(names(unlist(mofa_fit_1_metrics)), 'time_pe')
  }
  
  
  if(rotate){
    Y_c <- do.call(cbind, Y)
    set.seed(123)
    ptm <- proc.time() 
    rotate_fit_1 <- fit_rotate(Y_c, k + 5)
    rotate_time <- proc.time() - ptm
    rotate_fit_metrics <- compute_performance_fable(rotate_fit_1, Y, Lambdas_0, As, UQ=F)
    rotate_results <- rbind(rotate_results, c(unlist(rotate_fit_metrics),  rotate_time[3]))
    names(rotate_results) <- c(names(unlist(rotate_fit_metrics)), 'time_pe')
  }
  
  
}


write.csv(fama_results, paste0('simulations/results/scenario_', scenario, '/fama_results.csv'))
write.csv(fable_results, paste0('simulations/results/scenario_', scenario, '/fable_results.csv'))
write.csv(mofa_results, paste0('simulations/results/scenario_', scenario, '/mofa_results.csv'))
write.csv(rotate_results, paste0('simulations/results/scenario_', scenario, '/rotate_results.csv'))
save.image(file=paste0('simulations/results/scenario_', scenario, '/scenario_', scenario,'.RData'))

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)


#scenario <- 8

fama_results <- read.csv(paste0('simulations/results/scenario_', scenario, '/fama_results.csv'))
fable_results <- read.csv(paste0('simulations/results/scenario_', scenario, '/fable_results.csv'))
rotate_results <- read.csv(paste0('simulations/results/scenario_', scenario, '/rotate_results.csv'))
mofa_results <- read.csv(paste0('simulations/results/scenario_', scenario, '/mofa_results.csv'))

fama_results <- fama_results %>% mutate(method = "FAMA")
fable_results <- fable_results %>% mutate(method = "FABLE")
mofa_results <- mofa_results %>% mutate(method = "MOFA")
rotate_results <- rotate_results %>% mutate(method = "ROTATE")

all_results <- bind_rows(fama_results, fable_results, mofa_results, rotate_results)
all_results$method <- factor(all_results$method, levels = c("FAMA", "FABLE", "MOFA", "ROTATE"))
#all_results <- bind_rows(fama_results, fable_results, rotate_results)


p1 <- ggplot(all_results, aes(x = method, y = rmse_all, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "RMSE Overall", y = "RMSE", x = "Method") +
  theme(legend.position = "none")

intra_data <- all_results %>%
  dplyr::select(method, starts_with("rmses_intra")) %>%
  pivot_longer(cols = starts_with("rmses_intra"), names_to = "metric", values_to = "value")

p2 <- ggplot(intra_data, aes(x = method, y = value, fill = method)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "RMSEs Intraview", x = "Method", y = "RMSE") +
  theme(legend.position = "none")

inter_data <- all_results %>%
  dplyr::select(method, starts_with("rmses_inter")) %>%
  pivot_longer(cols = starts_with("rmses_inter"), names_to = "metric", values_to = "value")

p3 <- ggplot(inter_data, aes(x = method, y = value, fill = method)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "RMSEs Interview", x = "Method", y = "RMSE") +
  theme(legend.position = "none")

p4 <- ggplot(all_results, aes(x = method, y = time_pe, fill = method)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Running Time", x = "Method", y = "Seconds") +
  theme(legend.position = "none")


dev.new()
(p1 / p2 / p3 / p4) + plot_layout(ncol = 4)



write.csv(all_results, paste0('simulations/results/scenario_', scenario, '/all_results.csv'))

rbind(colMeans(fama_results[ , grepl("len", names(fama_results))]),
rep(colMeans(fable_results[ , grepl("len", names(fable_results))]),2))

colMeans(fama_results)
colMeans(fable_results)


coverage_fama <- fama_results[ , grepl("cov", names(fama_results))]
coverage_fama_clt <- coverage_fama[ , grepl("clt", names(coverage_fama))]
coverage_fama_bvm <- coverage_fama[ , grepl("posterior", names(coverage_fama))]

coverage_fable <- fable_results[ , grepl("cov", names(fable_results))]

coverage_fama_clt[] <- lapply(coverage_fama_clt, as.double)
coverage_fama_bvm[] <- lapply(coverage_fama_bvm, as.double)
coverage_fable[] <- lapply(coverage_fable, as.double)

# intra 
coverage_fama_clt_intra <- coverage_fama_clt[ , grepl("intra", names(coverage_fama_clt))]
coverage_fama_bvm_intra <- coverage_fama_bvm[ , grepl("intra", names(coverage_fama_bvm))]
coverage_fable_intra <- coverage_fable[ , grepl("intra", names(coverage_fable))]

mean(colMeans(coverage_fama_clt_intra));
mean(colMeans(coverage_fama_bvm_intra));
mean(colMeans(coverage_fable_intra));


# inter 
coverage_fama_clt_inter <- coverage_fama_clt[ , grepl("inter", names(coverage_fama_clt))]
coverage_fama_bvm_inter <- coverage_fama_bvm[ , grepl("inter", names(coverage_fama_bvm))]
coverage_fable_inter <- coverage_fable[ , grepl("inter", names(coverage_fable))]

mean(colMeans(coverage_fama_clt_inter));
mean(colMeans(coverage_fama_bvm_inter));
mean(colMeans(coverage_fable_inter))


scenario
mean(colMeans(coverage_fama_clt_intra));
mean(colMeans(coverage_fama_bvm_intra));
mean(colMeans(coverage_fable_intra));
mean(colMeans(coverage_fama_clt_inter));
mean(colMeans(coverage_fama_bvm_inter));
mean(colMeans(coverage_fable_inter))

# to do

# scenario 3
# coverage + plot

# scenario 4
# run mofa
# plot

# scenario 7
# run mofa
# coverage + plot

# scenario 8
# run mofa
# plot


# scenario 5
# run mofa
# coverage + plot

# scenario 6
# all

# scenario 9
# all

# scenario 10
# all



