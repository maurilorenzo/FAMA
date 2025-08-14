#include <C:/Users/mauri/AppData/Local/R/win-library/4.3/RcppArmadillo/include/RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::mat compute_sds_clt_intraview(
    const arma::mat& Lambda, const arma::vec& Sigma_2s
) {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat sds = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      sds(j, l) = sqrt(Lambda_outer(j, j)*Sigma_2s(l) + Lambda_outer(l, l)*Sigma_2s(j) +
        std::pow(Lambda_outer(j, l),2) + Lambda_outer(l, l) *  Lambda_outer(j, j));
      sds(l, j) = sds(j, l); 
    }
    sds(j, j) =sqrt(4*Lambda_outer(j, j)*Sigma_2s(j) + 2*std::pow(Lambda_outer(j, j),2) );
  }
  sds(p - 1, p - 1) = sqrt(4*Lambda_outer(p-1, p-1)*Sigma_2s(p-1) + 
    2*std::pow(Lambda_outer(p-1, p-1),2));
  
  return sds;
}

// [[Rcpp::export]]
arma::mat compute_sds_clt_biview(
    const arma::mat& Lambda_1, const  arma::mat& Lambda_2, 
    const arma::mat& D_Lambda_1, const  arma::mat& D_Lambda_2, 
    const arma::vec& Sigma_2s_1, const arma::vec& Sigma_2s_2
) {
  
  arma::mat Lambda_outer = Lambda_1 * Lambda_2.t(); 
  arma::vec Lambda_outer_1 = arma::sum(arma::square(Lambda_1), 1);; 
  arma::vec Lambda_outer_2 = arma::sum(arma::square(Lambda_2), 1);; 
  arma::vec D_Lambda_1_outer = arma::sum(arma::square(D_Lambda_1), 1);; 
  arma::vec D_Lambda_2_outer = arma::sum(arma::square(D_Lambda_2), 1);; 
  
  int p_1 = Lambda_1.n_rows; 
  int p_2 = Lambda_2.n_rows; 
  
  arma::mat sds = arma::zeros<arma::mat>(p_1, p_2);
  for (int j = 0; j < p_1; ++j) {
    for (int l = 0; l < p_2; ++l) {
      sds(j, l) = sqrt(Lambda_outer_1(j)*Sigma_2s_2(l) + Lambda_outer_2(l)*Sigma_2s_1(j) +
        std::pow(Lambda_outer(j, l),2) + D_Lambda_1_outer(j) *  D_Lambda_2_outer(l));
    }
  }
 return sds;
}

// [[Rcpp::export]]
arma::mat compute_sds_posterior_distribution_intraview(
    const arma::mat& Lambda, const arma::vec& Sigma_2s, const double rho=1
){
  arma::vec Lambda_inner_prods = arma::sum(arma::square(Lambda), 1); 
  int p = Lambda.n_rows; 
  arma::mat sds = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      sds(j, l) = rho * ( sqrt( Sigma_2s(l) * Lambda_inner_prods(j)  +
          Sigma_2s(j) * Lambda_inner_prods(l) ) ); 
      sds(l, j) = sds(j, l); 
    }
    sds(j, j) = 2 * rho * sqrt( Lambda_inner_prods(j)*Sigma_2s(j) ) ;
  }
  sds(p - 1, p - 1) = 2 * rho * sqrt( Lambda_inner_prods(p-1)*Sigma_2s(p-1) );
  
  return sds;
}

// [[Rcpp::export]]
arma::mat compute_sds_posterior_distribution_biview(
    const arma::mat& Lambda_1, const  arma::mat& Lambda_2, 
    const arma::vec& Sigma_2s_1, const arma::vec& Sigma_2s_2,
    const double rho_1=1, const double rho_2=1 
){
  arma::vec Lambda_1_inner_prods = arma::sum(arma::square(Lambda_1), 1); 
  arma::vec Lambda_2_inner_prods = arma::sum(arma::square(Lambda_2), 1); 
  
  int p_1 = Lambda_1.n_rows; 
  int p_2 = Lambda_2.n_rows; 
  
  arma::mat sds = arma::zeros<arma::mat>(p_1, p_2);
  
  for (int j = 0; j < p_1; ++j) {
    for (int l = 0; l < p_2; ++l) {
      sds(j, l) = sqrt(std::pow(rho_1, 2) * ( Sigma_2s_1(j) * Lambda_2_inner_prods(l) ) +
        std::pow(rho_2, 2) * ( Sigma_2s_2(l) * Lambda_1_inner_prods(j) ) ); 
    }
  }
  return sds;
}




// [[Rcpp::export]]
arma::mat compute_B(
    const arma::mat& Lambda, const arma::vec& Sigma_2s) 
  {
  arma::mat Lambda_outer = Lambda * Lambda.t(); 
  int p = Lambda.n_rows; 
  arma::mat B = arma::zeros<arma::mat>(p, p);
  
  for (int j = 0; j < p - 1; ++j) {
    for (int l = j + 1; l < p; ++l) {
      B(j, l) = sqrt(1 + (Lambda_outer(j, j) * Lambda_outer(l, l) + Lambda_outer(j, l) * Lambda_outer(j, l)) / 
        (Sigma_2s(j) * Lambda_outer(l, l) + Sigma_2s(l) * Lambda_outer(j, j)));
      B(l, j) = B(j, l); 
    }
    B(j, j) = sqrt(1 + Lambda_outer(j, j) / (2 * Sigma_2s(j)));
  }
  B(p - 1, p - 1) = sqrt(1 + Lambda_outer(p - 1, p - 1) / (2 * Sigma_2s(p - 1)));
  
  return B;
}


