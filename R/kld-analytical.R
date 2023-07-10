# Analytical Kullback-Leibler divergence for specific classes of functions

#' Analytical KL divergence for two uni- or multivariate Gaussian distributions
#'
#' This function computes \eqn{D_{KL}(p||q)}, where \eqn{p\sim \mathcal{N}(\mu_1,\Sigma_1)}
#' and \eqn{q\sim \mathcal{N}(\mu_2,\Sigma_2)}.
#' @inherit
#' @param mu1 A numeric vector (mean of true Gaussian)
#' @param sigma1 A s.p.d. matrix (Covariance matrix of true Gaussian)
#' @param mu2 A numeric vector (mean of approximate Gaussian)
#' @param sigma2 A s.p.d. matrix  (Covariance matrix of approximate Gaussian)
#' @return A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kl_div_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
#'                 mu2 = rep(1,2), sigma2 = matrix(c(1,0.5,0.5,1), nrow = 2))
kl_div_gaussian <- function(mu1, sigma1, mu2, sigma2) {

    sigma1 <- as.matrix(sigma1)
    sigma2 <- as.matrix(sigma2)

    d <- length(mu1)
    stopifnot(length(mu2) == d, dim(sigma1) == d, dim(sigma2) == d)

    ldet1   <- log(det(sigma1))
    ldet2   <- log(det(sigma2))
    isig2   <- solve(sigma2)
    diff_mu <- mu1 - mu2

    0.5*(ldet2 - ldet1 + tr(isig2 %*% sigma1) + as.vector(t(diff_mu) %*% isig2 %*% diff_mu) - d)

}

#' Analytical KL divergence for two univariate exponential distributions
#'
#' This function computes \eqn{D_{KL}(p||q)}, where \eqn{p\sim \text{Exp}(\lambda_1)}
#' and \eqn{q\sim \text{Exp}(\lambda_2)}, in rate parametrization.
#'
#' @param lambda1 A scalar (rate parameter of true exponential distribution)
#' @param lambda2 A scalar (rate parameter of approximate exponential distribution)
#' @returns A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' kl_div_exponential(lambda1 = 1, lambda2 = 2)
kl_div_exponential <- function(lambda1, lambda2) {

    log(lambda1) - log(lambda2) + lambda2 / lambda1 - 1

}

