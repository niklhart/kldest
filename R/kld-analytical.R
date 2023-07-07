# Analytical Kullback-Leibler divergence for specific classes of functions

#' Analytical KL divergence for two uni- or multivariate Gaussian distributions
#'
#' This function computes $D_{KL}(p||q)$, where $p\sim \mathcal{N}(\mu_1,\Sigma_1)$
#' and $q\sim \mathcal{N}(\mu_2,\Sigma_2)$.
#'
#' @param mu1 A numeric vector (mean of true Gaussian)
#' @param sigma1 A s.p.d. matrix (Covariance matrix of true Gaussian)
#' @param mu2 A numeric vector (mean of approximate Gaussian)
#' @param sigma2 A s.p.d. matrix  (Covariance matrix of approximate Gaussian)
#' @return A scalar (the Kullback-Leibler divergence)
kl_divergence_gaussian <- function(mu1, sigma1, mu2, sigma2) {

    sigma1 <- as.matrix(sigma1)
    sigma2 <- as.matrix(sigma2)

    D     <- length(mu1)
    ldet1 <- log(det(sigma1))
    ldet2 <- log(det(sigma2))
    isig2 <- solve(sigma2)
    dmu   <- mu1 - mu2

    0.5*(ldet2 - ldet1 + tr(isig2 %*% sigma1) + t(dmu) %*% isig2 %*% dmu - D)

}

#' Analytical KL divergence for two exponential distributions (1D)
#'
#' This function computes $D_{KL}(p||q)$, where $p\sim \text{Exp}(\lambda_1)$
#' and $q\sim \text{Exp}(\lambda_2)$, in rate parametrization.
#'
#' @param lambda1 A scalar (rate parameter of true exponential distribution)
#' @param lambda2 A scalar (rate parameter of approximate exponential distribution)
#' @returns A scalar (the Kullback-Leibler divergence)
kl_divergence_exponential <- function(lambda1, lambda2) {

    log(lambda1) - log(lambda2) + lambda2 / lambda1 - 1
}

