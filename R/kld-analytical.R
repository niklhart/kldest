# Analytical Kullback-Leibler divergence for specific classes of functions

#' Analytical KL divergence for two uni- or multivariate Gaussian distributions
#'
#' This function computes \eqn{D_{KL}(p||q)}, where \eqn{p\sim \mathcal{N}(\mu_1,\Sigma_1)}
#' and \eqn{q\sim \mathcal{N}(\mu_2,\Sigma_2)}.
#'
#' @param mu1 A numeric vector (mean of true Gaussian)
#' @param sigma1 A s.p.d. matrix (Covariance matrix of true Gaussian)
#' @param mu2 A numeric vector (mean of approximate Gaussian)
#' @param sigma2 A s.p.d. matrix  (Covariance matrix of approximate Gaussian)
#' @return A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' kld_gaussian(mu1 = 1, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kld_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
#'                 mu2 = rep(1,2), sigma2 = matrix(c(1,0.5,0.5,1), nrow = 2))
kld_gaussian <- function(mu1, sigma1, mu2, sigma2) {

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
#' kld_exponential(lambda1 = 1, lambda2 = 2)
kld_exponential <- function(lambda1, lambda2) {

    log(lambda1) - log(lambda2) + lambda2 / lambda1 - 1

}


#' Analytical KL divergence for two uniform distributions
#'
#' This function computes \eqn{D_{KL}(p||q)}, where \eqn{p\sim \text{U}(a_1,b_1)}
#' and \eqn{q\sim \text{U}(a_2,b_2)}, with \eqn{a_2<a_1<b_1<b_2}.
#'
#' @param a1,b1 Range of true uniform distribution
#' @param a2,b2 Range of approximate uniform distribution
#' @returns A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' kld_uniform(a1 = 0, b1 = 1, a2 = 0, b2 = 2)
kld_uniform <- function(a1, b1, a2, b2) {

    stopifnot(a1 < b1, a2 <= a1, b2 >= b1)

    log((b2-a2)/(b1-a1))
}

#' Analytical KL divergence between a uniform and a Gaussian distribution
#'
#' This function computes \eqn{D_{KL}(p||q)}, where \eqn{p\sim \text{U}(a,b)}
#' and \eqn{q\sim \mathcal{N}(\mu,\sigma^2)}.
#'
#' @param a,b Parameters of uniform (true) distribution
#' @param mu,sigma2 Parameters of Gaussian (approximate) distribution
#' @returns A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' kld_uniform_gaussian(a = 0, b = 1, mu = 0, sigma2 = 1)
kld_uniform_gaussian <- function(a = 0, b = 1, mu = 0, sigma2 = 1) {

    stopifnot(a < b, sigma2 > 0)

    log(sqrt(2*pi*sigma2)/(b-a)) + ((b-mu)^3 - (a-mu)^3) / (6*(b-a)*sigma2)

}

#' Analytical KL divergence for two discrete distributions
#'
#' @param P,Q Numerical arrays with the same dimensions, representing discrete
#'     probability distributions
#' @returns A scalar (the Kullback-Leibler divergence)
#' @export
#' @examples
#' # 1-D example
#' P <- 1:4/10
#' Q <- rep(0.25,4)
#' kld_discrete(P,Q)
#'
#' # The above example in 2-D
#' P <- matrix(1:4/10,nrow=2)
#' Q <- matrix(0.25,nrow=2,ncol=2)
#' kld_discrete(P,Q)
#'
kld_discrete <- function(P,Q) {

    # Input checking
    P <- as.array(P)
    Q <- as.array(Q)
    if (any(dim(P) != dim(Q))) stop("Inputs must have the same dimensions.")
    if (any(P < 0 | Q < 0)) stop("Input arrays must be nonnegative.")
    if (any(abs(c(sum(P),sum(Q)) - 1) > sqrt(.Machine$double.eps))) stop("Input arrays must sum up to 1.")

    # Calculation, taking care of edge case P[i] == 0:
    posP <- P > 0
    sum(P[posP] * log(P[posP]/Q[posP]))
}



