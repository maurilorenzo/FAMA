if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", 
                       "MultiAssayExperiment", "curatedTCGAData"))


library(curatedTCGAData)
library(SummarizedExperiment)

brca <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays = c("RNASeq2GeneNorm", "Methylation_methyl27", "miRNASeqGene"),
  version = '2.1.1',
  dry.run = FALSE
  )
class(brca)

names(experiments(brca))

rna <- assays(brca[["BRCA_RNASeq2GeneNorm-20160128"]])[[1]]
meth <- assays(brca[["BRCA_Methylation_methyl27-20160128"]])[[1]]
mirna <- assays(brca[["BRCA_miRNASeqGene-20160128"]])[[1]]

common_samples <- Reduce(intersect, list(
  colnames(rna),
  #colnames(meth),
  colnames(mirna)
))

length(common_samples)




subsample_index <- 1:100
fit_fama_1 <- fit_FAMA(Y, clt_SE=T, posterior_SE=T, 
                       index_SE=subsample_index)


Y_c <- do.call(cbind, Y)
ptm <- proc.time() 
fable_fit_1 <- PseudoPosteriorMean(Y_c, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
fable_time <- proc.time() - ptm
fable_fit_1$estRank


p <- ncol(Y_C)
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

ptm <- proc.time() 
mofa_fit_1 <- fit_MOFA(Y, k_0 + 5)
mofa_time <- proc.time() - ptm



# - preliminary oos experiments #################################################

library(mvtnorm)
fama_1_cov_intra_ests <- list()
fama_1_cov_inter_ests <- list()
for(m in 1:(M-1)){
  fama_1_cov_intra_ests[[m]] <- t(fit_fama_1$Lambdas_hat[[m]]) + diag(as.vector(fit_fama_1$sigmas_sq_hat[[m]]))
  fama_1_cov_inter_ests[[m]] <- list()
  for(v in (s+1):M) {
    print(v)
     fama_1_cov_inter_ests[[m]][[l]] = tcrossprod(fit_fama_1$Lambdas_hat[[m]], fit_fama_1$Lambdas_hat[[v]])
    }
}
fama_1_cov_intra_ests[[M]] <- t(fit_fama_1$Lambdas_hat[[M]]) + diag(as.vector(fit_fama_1$sigmas_sq_hat[[M]]))

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



fable_1_cov_est <- fable_fit_1$FABLEPostMean
compute_set <- 1:250
fama_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set],
                             sigma=fama_1_cov_intra_ests[[m]][compute_set, compute_set], log=T))
fama_1_oos_ll
fable_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=fable_1_cov_est[compute_set, compute_set], log=T))
fable_1_oos_ll
rotate_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=fama_1_cov_est[compute_set, compute_set], log=T))
fama_1_oos_ll
mofa_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=mofa_1_cov_est[compute_set, compute_set], log=T))
mofa_1_oos_ll

fama_1_oos_ll_intra_save <- c()
mofa_1_oos_ll_intra_save <- c()
rotate_1_oos_ll_intra_save <- c()
fable_1_oos_ll_intra_save <- c()

fama_1_oos_ll_inter_save <- c()
mofa_1_oos_ll_inter_save <- c()
rotate_1_oos_ll_inter_save <- c()
fable_1_oos_ll_inter_save <- c()

n_try <- 1000


for(m in 1:(M-1)){
  fama_1_oos_ll_intra_save[[m]] <- c()
  fama_1_oos_ll_intra_save[[m]] <- c()
  
  for(i in 1:n_try){
    fama_1_oos_ll_intra <- sum(dmvnorm(Y_test[[m]][, compute_set], 
                                 sigma=fama_1_cov_est[[m]][compute_set, compute_set], log=T))
    fama_1_oos_ll_intra_save[[m]][i] = fama_1_oos_ll_intra
    fama_1_oos_ll_intra <- sum(dmvnorm(Y_test[[m]][, compute_set], 
                                 sigma=fama_1_cov_est[[m]][compute_set, compute_set], log=T))
    fama_1_oos_ll_intra_save[[m]][i] = fama_1_oos_ll_intra
   
  for(v in (s+1):M) {
      print(v)
      fama_1_cov_inter_ests[[m]][[l]] = tcrossprod(fit_fama_1$Lambdas_hat[[m]], fit_fama_1$Lambdas_hat[[v]])
    }
  }
  fama_1_oos_ll_intra <- sum(dmvnorm(Y_test[[M]][, compute_set], 
                               sigma=fama_1_cov_est[[M]][compute_set, compute_set], log=T))
  fama_1_oos_ll_intra_save[[M]][i] = fama_1_oos_ll_intra  
}



for(i in 1:n_try){
  set.seed(i)
  compute_set <- sample(1:p, size=250)
  fama_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=fama_1_cov_est[compute_set, compute_set], log=T))
  fama_1_oos_ll_save[i] = fama_1_oos_ll
  mofa_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=mofa_1_cov_est[compute_set, compute_set], log=T))
  mofa_1_oos_ll_save[i] = mofa_1_oos_ll
  rotate_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=rotate_1_cov_est[compute_set, compute_set], log=T))
  rotate_1_oos_ll_save[i] = rotate_1_oos_ll
  fable_1_oos_ll <- sum(dmvnorm(Y_test[, compute_set], sigma=fable_1_cov_est[compute_set, compute_set], log=T))
  fable_1_oos_ll_save[i] = fable_1_oos_ll
}

summary(fama_1_oos_ll_save)
summary(mofa_1_oos_ll_save)
summary(rotate_1_oos_ll_save)
summary(fable_1_oos_ll_save)

lim_min <- min(c(fama_1_oos_ll_save, mofa_1_oos_ll_save, rotate_1_oos_ll_save, 
                 fable_1_oos_ll_save)) - 0.1*abs(mean(mofa_1_oos_ll_save))
lim_max <- max(c(fama_1_oos_ll_save, mofa_1_oos_ll_save, rotate_1_oos_ll_save, 
                 fable_1_oos_ll_save)) + 0.1*abs(mean(mofa_1_oos_ll_save))

par(mfrow=c(1,4))
boxplot(fama_1_oos_ll_save, ylim=c(lim_min, lim_max))
boxplot(mofa_1_oos_ll_save, ylim=c(lim_min, lim_max))
boxplot(fable_1_oos_ll_save, ylim=c(lim_min, lim_max))
boxplot(rotate_1_oos_ll_save, ylim=c(lim_min, lim_max))
