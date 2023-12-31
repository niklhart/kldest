% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kld-estimation-nearest-neighbours.R
\name{kld_est_nn}
\alias{kld_est_nn}
\title{k-nearest neighbour KL divergence estimator}
\usage{
kld_est_nn(X, Y = NULL, q = NULL, k = 1L, eps = 0, log.q = FALSE)
}
\arguments{
\item{X, Y}{\code{n}-by-\code{d} and \code{m}-by-\code{d} matrices, representing \code{n} samples from
the true distribution \eqn{P} and \code{m} samples from the approximate distribution
\eqn{Q}, both in \code{d} dimensions. Vector input is treated as a column matrix.
\code{Y} can be left blank if \code{q} is specified (see below).}

\item{q}{The density function of the approximate distribution \eqn{Q}. Either
\code{Y} or \code{q} must be specified.}

\item{k}{The number of nearest neighbours to consider for NN density estimation.
Larger values for \code{k} generally increase bias, but decrease variance of the
estimator. Defaults to \code{k = 1}.}

\item{eps}{Error bound in the nearest neighbour search. A value of \code{eps = 0}
(the default) implies an exact nearest neighbour search, for \code{eps > 0}
approximate nearest neighbours are sought, which may be somewhat faster for
high-dimensional problems.}

\item{log.q}{If \code{TRUE}, function \code{q} is the log-density rather than the density
of the approximate distribution \eqn{Q} (default: \code{log.q = FALSE}).}
}
\value{
A scalar, the estimated Kullback-Leibler divergence \eqn{\hat D_{KL}(P||Q)}.
}
\description{
This function estimates Kullback-Leibler divergence \eqn{D_{KL}(P||Q)} between
two continuous distributions \eqn{P} and \eqn{Q} using nearest-neighbour (NN)
density estimation in a Monte Carlo approximation of \eqn{D_{KL}(P||Q)}.
}
\details{
Input for estimation is a sample \code{X} from \eqn{P} and either the density
function \code{q} of \eqn{Q} (one-sample problem) or a sample \code{Y} of \eqn{Q}
(two-sample problem). In the two-sample problem, it is the estimator in Eq.(5)
of Wang et al. (2009). In the one-sample problem, the asymptotic bias (the
expectation of a Gamma distribution) is substracted, see Pérez-Cruz (2008),
Eq.(18).

References:

Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
Theory, Vol. 55, No. 5 (2009).

Pérez-Cruz, "Kullback-Leibler Divergence Estimation of Continuous
Distributions", IEEE International Symposium on Information Theory (2008).
}
\examples{
# KL-D between one or two samples from 1-D Gaussians:
set.seed(0)
X <- rnorm(100)
Y <- rnorm(100, mean = 1, sd = 2)
q <- function(x) dnorm(x, mean = 1, sd =2)
kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
kld_est_nn(X, Y)
kld_est_nn(X, q = q)
kld_est_nn(X, Y, k = 5)
kld_est_nn(X, q = q, k = 5)
kld_est_brnn(X, Y)


# KL-D between two samples from 2-D Gaussians:
set.seed(0)
X1 <- rnorm(100)
X2 <- rnorm(100)
Y1 <- rnorm(100)
Y2 <- Y1 + rnorm(100)
X <- cbind(X1,X2)
Y <- cbind(Y1,Y2)
kld_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
             mu2 = rep(0,2), sigma2 = matrix(c(1,1,1,2),nrow=2))
kld_est_nn(X, Y)
kld_est_nn(X, Y, k = 5)
kld_est_brnn(X, Y)
}
