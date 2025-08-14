
#load('application/results/oos_2.RData')

file_urls <- c(
  #Clinical     = "https://linkedomics.org/data_download/TCGA-BRCA/Human__TCGA_BRCA__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
  RNAseq_ga    = "https://linkedomics.org/data_download/TCGA-BRCA/Human__TCGA_BRCA__UNC__RNAseq__GA_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz",
  CNV          = "https://linkedomics.org/data_download/TCGA-BRCA/Human__TCGA_BRCA__BI__SCNA__SNP_6.0__01_28_2016__BI__Gene__Firehose_GISTIC2.cct.gz",
  RPPA         = "https://linkedomics.org/data_download/TCGA-BRCA/Human__TCGA_BRCA__MDA__RPPA__MDA_RPPA__01_28_2016__BI__Gene__Firehose_RPPA.cct",
  Methylation  = "https://linkedomics.org/data_download/TCGA-BRCA/Human__TCGA_BRCA__JHU_USC__Methylation__Meth450__01_28_2016__BI__Gene__Firehose_Methylation_Prepocessor.cct.gz"
)

# Destination folder
dest_dir <- "application/data"
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
}


for (name in names(file_urls)) {
  url <- file_urls[name]
  # Preserve original file extension
  ext <- tools::file_ext(url)
  destfile <- file.path(dest_dir, paste0(name, ".", ext))
  download.file(url, destfile, mode = "wb")
  cat("✔ Downloaded", name, "to", destfile, "\n")
}

library(data.table)
library(matrixStats)

rna_ga_path <- "application/data/RNAseq_ga.gz"
con <- gzfile(rna_ga_path, "rt")
rna_ga_data <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
close(con)
dim(rna_ga_data)
head(colnames(rna_ga_data))
rna_ga_data[1:5, 1:5]
rownames(rna_ga_data) <- rna_ga_data[[1]]
rna_ga_data <- rna_ga_data[, -1]
dim(rna_ga_data)
sum(is.na(rna_ga_data))
head(rna_ga_data[, 1:5])
#rna_ga_data <- rna_ga_data[complete.cases(rna_ga_data), ]
dim(rna_ga_data)


cnv_path <- "application/data/CNV.gz"
con <- gzfile(cnv_path, "rt")
cnv_data <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
close(con)
dim(cnv_data)
head(colnames(cnv_data))
cnv_data[1:5, 1:5]
rownames(cnv_data) <- cnv_data[[1]]
cnv_data <- cnv_data[, -1]
dim(cnv_data)
head(cnv_data[, 1:5])



methylation_path <- "application/data/Methylation.gz"
con <- gzfile(methylation_path, "rt")
methylation_data <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
close(con)
dim(methylation_data)
head(colnames(methylation_data))
methylation_data[1:5, 1:5]
rownames(methylation_data) <- methylation_data[[1]]
methylation_data <- methylation_data[, -1]
methylation_data <- methylation_data[complete.cases(methylation_data), ]
dim(methylation_data)
head(methylation_data[, 1:5])



rppa_data <- fread('application/data/RPPA.cct', data.table = FALSE)
dim(rppa_data)
head(colnames(rppa_data))
rppa_data[1:5, 1:5]
rownames(rppa_data) <- rppa_data[[1]]
rppa_data <- rppa_data[, -1]
dim(rppa_data)
head(rppa_data[, 1:5])
rppa_data <- rppa_data[complete.cases(rppa_data), ]
dim(rppa_data)



num_genes <- 4000

#gene_variances <- apply(rna_hs_data, 1, var)
#top4000_indices <- order(gene_variances, decreasing = TRUE)[1:num_genes]
#rna_hs_top4000 <- rna_hs_data[top4000_indices,]

gene_variances <- apply(rna_ga_data, 1, var)
top4000_indices <- order(gene_variances, decreasing = TRUE)[1:num_genes]
rna_ga_top4000 <- rna_ga_data[top4000_indices,]

gene_variances <- apply(cnv_data, 1, var)
top4000_indices <- order(gene_variances, decreasing = TRUE)[1:num_genes]
cnv_top4000 <- cnv_data[top4000_indices,]

gene_variances <- apply(methylation_data, 1, var)
top4000_indices <- order(gene_variances, decreasing = TRUE)[1:num_genes]
methylation_top4000 <- methylation_data[top4000_indices,]



datasets <- list(
  RNAseq_ga = rna_ga_top4000,
  CNV = cnv_top4000,
  RPPA = rppa_data,
  methylation = methylation_top4000
)
sample_sets <- lapply(datasets, colnames)
common_samples <- Reduce(intersect, sample_sets)
cat("Number of common samples:", length(common_samples), "\n")

Y <- lapply(datasets, function(data) {
  data_common <- data[, common_samples, drop = FALSE]
  t_data <- t(data_common)
})


to_normal <- function(x){
  n <- length(x)
  std <- sd(x)
  ranks <- rank(x)
  return(qnorm(ranks/(n+1), 0 , std))
}

for(m in 1:length(Y)){
  Y_m <- apply(Y[[m]], 2, to_normal)
  Y[[m]] <- scale(Y_m[,colSds(Y_m)>0])
}


str(Y, max.level = 1)

na_check <- sapply(Y, function(mat) any(is.na(mat)))
print(na_check)

sum(is.na(Y[[m]]))
dim(Y[[m]])

# Confirm that everything is clean
if (all(!na_check)) {
  cat("✔ No missing values found in any omics matrix.\n")
} else {
  cat("⚠ Missing values found in:\n")
  print(names(na_check)[na_check])
}


set.seed(99) 
n <- nrow(Y[[1]])
M <- length(Y)
train_idx <- sample(seq_len(n), size = 0.8 * n)
Y_train <- lapply(Y, function(mat) mat[train_idx, , drop = FALSE])
Y_test  <- lapply(Y, function(mat) mat[-train_idx, , drop = FALSE])

sapply(Y_train, dim)
sapply(Y_test, dim)

k_max <- c()
for(m in 1:M){
  s_Y_m <- svd(Y_train[[m]])
  plot(cumsum(s_Y_m$d^2)/sum(s_Y_m$d^2))
  k_max[m] <- min(which(cumsum(s_Y_m$d^2)/sum(s_Y_m$d^2)>0.8))
}
k_max

#s_Y_3 <- svd(Y_train[[3]])
#plot(cumsum(s_Y_3$d^2)/sum(s_Y_3$d^2))
#which(cumsum(s_Y_3$d^2)/sum(s_Y_3$d^2)>0.9)

boxplot(colMeans(Y_train[[3]]))
#k_max[3] <- 20

subsample_index <- 1:100
set.seed(123)
ptm <- proc.time() 
#k_max <- 50
fit_fama_1 <- fit_FAMA(Y_train, clt_SE=T, posterior_SE=T, 
                       index_SE=subsample_index, k_max=k_max)
fama_time <- proc.time() - ptm
fama_time
fit_fama_1$k_0_hat

Y_train_c <- do.call(cbind, Y_train)
ptm <- proc.time() 
fable_fit_1 <- PseudoPosteriorMean(Y_train_c, gamma0 = 1, delta0sq = 1, maxProp = 0.95)
fable_time <- proc.time() - ptm
fable_fit_1$estRank

Y_train_c <- do.call(cbind, Y_train)
set.seed(123)
ptm <- proc.time() 
rotate_fit_1 <- fit_rotate(Y_train_c, fit_fama_1$k_0_hat + 5)
rotate_time <- proc.time() - ptm
rotate_time[3]
rotate_fit_1$cov_est <- tcrossprod(rotate_fit_1$B)

set.seed(123)
ptm <- proc.time() 
rotate_fit_2 <- fit_rotate(Y_train_c, fable_fit_1$estRank + 5)
rotate_time_2 <- proc.time() - ptm
rotate_time_2[3]
rotate_fit_2$cov_est <- tcrossprod(rotate_fit_2$B)


set.seed(123)
ptm <- proc.time() 
mofa_fit_1 <- fit_MOFA(Y_train, fit_fama_1$k_0_hat + 5)
mofa_time <- proc.time() - ptm
mofa_time

set.seed(123)
ptm <- proc.time() 
mofa_fit_2 <- fit_MOFA(Y_train,  fable_fit_1$estRank + 5) 
mofa_time_2 <- proc.time() - ptm
mofa_time_2
#save.image(file=paste0('application/results/oos_3.RData'))

#load('application/results/oos_2.RData')
# - preliminary oos experiments #################################################
fama_cov_all <- tcrossprod(do.call(rbind, fit_fama_1$Lambdas_hat))+
  diag(as.vector(do.call(c, fit_fama_1$sigmas_sq_hat)))
fable_cov_all <- fable_fit_1$FABLEPostMean
rotate_cov_all <- tcrossprod(rotate_fit_1$B) + diag(as.vector(rotate_fit_1$sigma^2))
rotate_cov_all_2 <- tcrossprod(rotate_fit_2$B) + diag(as.vector(rotate_fit_2$sigma^2))

mofa_fit_1$sigmas_sq_hat <- list()
for(m in 1:length(Y)){
  est <- mofa_fit_1$F_hat$group1 %*% t(mofa_fit_1$Lambdas_hat[[m]])
  res <- Y_train[[m]] - est
  sigmas_sq <- colMeans(res^2) 
  mofa_fit_1$sigmas_sq_hat[[m]] <- sigmas_sq
}

mofa_cov_all <- tcrossprod(do.call(rbind, mofa_fit_1$Lambdas_hat)) +
  diag(as.vector(do.call(c, mofa_fit_1$sigmas_sq_hat)))
dim(mofa_cov_all)

mofa_fit_2$sigmas_sq_hat <- list()
for(m in 1:length(Y)){
  est <- mofa_fit_2$F_hat$group1 %*% t(mofa_fit_2$Lambdas_hat[[m]])
  res <- Y_train[[m]] - est
  sigmas_sq <- colMeans(res^2) 
  mofa_fit_2$sigmas_sq_hat[[m]] <- sigmas_sq
}

mofa_cov_all_2 <- tcrossprod(do.call(rbind, mofa_fit_2$Lambdas_hat))+
  diag(as.vector(do.call(c, mofa_fit_2$sigmas_sq_hat)))

library(mvtnorm)

Y_test_c <- do.call(cbind, Y_test)
dim(Y_test_c)
fama_oos_lls <- dmvnorm(Y_test_c, sigma=fama_cov_all, log=T)
rotate_oos_lls <- dmvnorm(Y_test_c, sigma=rotate_cov_all, log=T)
sum(fama_oos_lls); sum(rotate_oos_lls)

mofa_oos_lls <- dmvnorm(Y_test_c, sigma=mofa_cov_all, log=T)
sum(fama_oos_lls); sum(fable_oos_lls);  sum(mofa_oos_lls);  sum(rotate_oos_lls)
sum(fama_oos_lls); sum(mofa_oos_lls);


rotate_oos_lls_2 <- dmvnorm(Y_test_c, sigma=rotate_cov_all_2, log=T)
sum(fama_oos_lls); sum(rotate_oos_lls);  sum(rotate_oos_lls_2)

fable_oos_lls <- dmvnorm(Y_test_c, sigma=fable_cov_all, log=T)
sum(fama_oos_lls); sum(fable_oos_lls); sum(rotate_oos_lls)

mofa_oos_lls <- dmvnorm(Y_test_c, sigma=mofa_cov_all, log=T)
sum(fama_oos_lls); sum(fable_oos_lls);  sum(mofa_oos_lls);  sum(rotate_oos_lls)
sum(fama_oos_lls); sum(mofa_oos_lls);


mofa_oos_lls_2 <- dmvnorm(Y_test_c, sigma=mofa_cov_all_2, log=T)
sum(fama_oos_lls);sum(fable_oos_lls);  sum(mofa_oos_lls); sum(mofa_oos_lls_2); sum(rotate_oos_lls); sum(rotate_oos_lls_2)

t1 <- t.test(fama_oos_lls, fable_oos_lls, alternative='greater', paired=TRUE)
t1

t2 <- t.test(fama_oos_lls, rotate_oos_lls, alternative='greater', paired=TRUE)
t2


t2 <- t.test(fama_oos_lls, rotate_oos_lls, alternative='greater', paired=TRUE)
t2

t2.2 <- t.test(fama_oos_lls, rotate_oos_lls_2, alternative='greater', paired=TRUE)
t2.2

t3 <- t.test(fama_oos_lls, mofa_oos_lls, alternative='greater', paired=TRUE)
t3

t3.2 <- t.test(fama_oos_lls, mofa_oos_lls_2, alternative='greater', paired=TRUE)
t3.2


fama_oos_lls_intra <- list()
fama_oos_lls_intra <- list()
fable_oos_lls_intra <- list()
rotate_oos_lls_intra <- list()
rotate_oos_lls_2_intra <- list()
mofa_oos_lls_intra <- list()
mofa_oos_lls_2_intra <- list()

idx_start <- 1
for(m in 1:M){
  print(m)
  idx_m <- idx_start:(idx_start + ncol(Y_test[[m]])-1)
  idx_start <- idx_start + ncol(Y_test[[m]])
  fama_cov_m <- fama_cov_all[idx_m, idx_m]
  fable_cov_m <- fable_cov_all[idx_m, idx_m]
  rotate_cov_m <- rotate_cov_all[idx_m, idx_m]
  rotate_cov_m_2 <- rotate_cov_all_2[idx_m, idx_m]
  mofa_cov_m <- mofa_cov_all[idx_m, idx_m]
  mofa_cov_m_2 <- mofa_cov_all_2[idx_m, idx_m]
  fama_oos_lls_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=fama_cov_m, log=T)
  fable_oos_lls_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=fable_cov_m, log=T)
  rotate_oos_lls_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=rotate_cov_m, log=T)
  rotate_oos_lls_2_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=rotate_cov_m_2, log=T)
  mofa_oos_lls_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=mofa_cov_m, log=T)
  mofa_oos_lls_2_intra[[m]] <- dmvnorm(Y_test[[m]], sigma=mofa_cov_m_2, log=T)
  print(sum(fama_oos_lls_intra[[m]]))
  print(sum(mofa_oos_lls_intra[[m]]))
  print(sum(rotate_oos_lls_intra[[m]]))
  print(sum(mofa_oos_lls_2_intra[[m]]))
  print(sum(rotate_oos_lls_2_intra[[m]]))
}

test <- T
for(m in 1:M){
  print(m)
  print(sum(fama_oos_lls_intra[[m]]))
  print(sum(fable_oos_lls_intra[[m]]))
  print(sum(rotate_oos_lls_intra[[m]]))
  print(sum(rotate_oos_lls_2_intra[[m]]))
  print(sum(mofa_oos_lls_intra[[m]]))
  print(sum(mofa_oos_lls_2_intra[[m]]))
  
  if(test){
    test.1 <- t.test(fama_oos_lls_intra[[m]], fable_oos_lls_intra[[m]], alternative='greater', paired=TRUE)
    print(test.1$p.value)
    test.1 <- t.test(fama_oos_lls_intra[[m]], rotate_oos_lls_intra[[m]], alternative='greater', paired=TRUE)
    print(test.1$p.value)
    test.1 <- t.test(fama_oos_lls_intra[[m]], rotate_oos_lls_2_intra[[m]], alternative='greater', paired=TRUE)
    print(test.1$p.value)
    test.1 <- t.test(fama_oos_lls_intra[[m]], mofa_oos_lls_intra[[m]], alternative='greater', paired=TRUE)
    print(test.1$p.value)
    test.1 <- t.test(fama_oos_lls_intra[[m]], mofa_oos_lls_2_intra[[m]], alternative='greater', paired=TRUE)
    print(test.1$p.value)
  }
}

fama_oos_lls_inter <- list()
fama_oos_lls_inter <- list()
fable_oos_lls_inter <- list()
rotate_oos_lls_inter <- list()
rotate_oos_lls_2_inter <- list()
mofa_oos_lls_inter <- list()
mofa_oos_lls_2_inter <- list()


p_s <- sapply(Y, ncol)
idx_start <- 1
for(m in 1:(M-1)){
  print(m)
  idx_m <- idx_start:(idx_start + ncol(Y_test[[m]])-1)
  idx_start <- idx_start + ncol(Y_test[[m]])
  fama_oos_lls_inter[[m]] <- list()
  fable_oos_lls_inter[[m]] <- list()
  rotate_oos_lls_inter[[m]] <- list()
  rotate_oos_lls_2_inter[[m]] <- list() 
  mofa_oos_lls_inter[[m]] <- list()
  mofa_oos_lls_2_inter[[m]] <- list()

  for(v in (m+1):M){
    print(v)
    idx_v <- sum(p_s[1:(v-1)]):sum(p_s[1:(v)])
    idx_tot <- c(idx_m, idx_v) 
    fama_cov_m <- fama_cov_all[idx_tot, idx_tot]
    fable_cov_m <- fable_cov_all[idx_tot, idx_tot]
    rotate_cov_m <- rotate_cov_all[idx_tot, idx_tot]
    rotate_cov_m_2 <- rotate_cov_all_2[idx_tot, idx_tot]
    mofa_cov_m <- mofa_cov_all[idx_tot, idx_tot]
    mofa_cov_m_2 <- mofa_cov_all_2[idx_tot, idx_tot]
    fama_oos_lls_inter[[m]][[v]] <- dmvnorm(Y_test_c[, idx_tot], sigma=fama_cov_m, log=T)
    fable_oos_lls_inter[[m]][[v]] <- dmvnorm(Y_test_c[,idx_tot], sigma=fable_cov_m, log=T)
    rotate_oos_lls_inter[[m]][[v]] <- dmvnorm(Y_test_c[,idx_tot], sigma=rotate_cov_m, log=T)
    rotate_oos_lls_2_inter[[m]][[v]] <- dmvnorm(Y_test_c[,idx_tot], sigma=rotate_cov_m_2, log=T)
    mofa_oos_lls_inter[[m]][[v]] <- dmvnorm(Y_test_c[,idx_tot], sigma=mofa_cov_m, log=T)
    mofa_oos_lls_2_inter[[m]][[v]] <- dmvnorm(Y_test_c[,idx_tot], sigma=mofa_cov_m_2, log=T)
    print(sum(fama_oos_lls_inter[[m]][[v]]))
    print(sum(fable_oos_lls_inter[[m]][[v]]))
    print(sum(rotate_oos_lls_inter[[m]][[v]]))
    print(sum(rotate_oos_lls_2_inter[[m]][[v]]))
    print(sum(mofa_oos_lls_inter[[m]][[v]]))
    print(sum(mofa_oos_lls_2_inter[[m]][[v]]))

  }
  
}

for(m in 1:(M-1)){
  print(m)
  
  
  for(v in (m+1):M){
    print(v)
    print(sum(fama_oos_lls_inter[[m]][[v]]))
    print(sum(fable_oos_lls_inter[[m]][[v]]))
    print(sum(rotate_oos_lls_inter[[m]][[v]]))
    print(sum(rotate_oos_lls_2_inter[[m]][[v]]))
    print(sum(mofa_oos_lls_inter[[m]][[v]]))
    print(sum(mofa_oos_lls_2_inter[[m]][[v]]))
    
    if(test){
      if(any(fama_oos_lls_inter[[m]][[v]] == -Inf)){
        print('fama - inf')
        
        } else {
        if(any(fable_oos_lls_inter[[m]][[v]] == -Inf)) {
          print('fable - inf')
        }
        else{
          test.1 <- t.test(fama_oos_lls_inter[[m]][[v]], fable_oos_lls_inter[[m]][[v]], alternative='greater', paired=TRUE)
          print(test.1$p.value)
        }
        
        if(any(rotate_oos_lls_inter[[m]][[v]] == -Inf)) {
          print('rotate - inf')
        }
        else{
          test.1 <- t.test(fama_oos_lls_inter[[m]][[v]], rotate_oos_lls_inter[[m]][[v]], alternative='greater', paired=TRUE)
          print(test.1$p.value)
        }
        if(any(rotate_oos_lls_2_inter[[m]][[v]] == -Inf)) {
          print('rotate 2 - inf')
        }
        else{
          test.1 <- t.test(fama_oos_lls_inter[[m]][[v]], rotate_oos_lls_2_inter[[m]][[v]], alternative='greater', paired=TRUE)
          print(test.1$p.value)
        }
        if(any(mofa_oos_lls_inter[[m]][[v]] == -Inf)) {
          print('mofa - inf')
        }
        else{
          test.1 <- t.test(fama_oos_lls_inter[[m]][[v]], mofa_oos_lls_inter[[m]][[v]], alternative='greater', paired=TRUE)
          print(test.1$p.value)
        }
        if(any(mofa_oos_lls_2_inter[[m]][[v]] == -Inf)) {
          print('mofa 2 - inf')
        }
        else{
          test.1 <- t.test(fama_oos_lls_inter[[m]][[v]], mofa_oos_lls_2_inter[[m]][[v]], alternative='greater', paired=TRUE)
          print(test.1$p.value)
        }
        
      }
    }
  }
  
}




