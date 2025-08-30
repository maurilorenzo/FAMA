
library(FABLE)

# slight modification of code from https://github.com/shounakch/FABLE/tree/main
# to implement the FABLE methodology from
# Shounak Chattopadhyay, Anru R. Zhang, and David B. Dunson, (2024) 
# "Blessing of dimension in Bayesian inference on covariance matrices"
# arXiv preprint arXiv:2404.03805


PseudoPosteriorMean_2 <- function(Y,
                                gamma0 = 1,
                                delta0sq = 1,
                                maxProp = 0.5) {
  
  tFABLEPostMean1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  kEst = RankEstimator(Y, 
                       U_Y,
                       V_Y,
                       svalsY,
                       kMax)
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      gamma0,
                                      delta0sq)
  
  Part1 = FABLEHypPars$G
  Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  OutputList = list("FABLEPostMean" = CovEst,
                    "Lambda_outer" = Part1,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean)
  
  return(OutputList)
  
}


PseudoPosteriorSampler_2 <- function(fit,
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
