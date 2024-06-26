% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kld-uncertainty.R
\name{kld_ci_subsampling}
\alias{kld_ci_subsampling}
\title{Uncertainty of KL divergence estimate using Politis/Romano's subsampling bootstrap.}
\usage{
kld_ci_subsampling(
  X,
  Y = NULL,
  q = NULL,
  estimator = kld_est_nn,
  B = 500L,
  alpha = 0.05,
  subsample.size = function(x) x^(2/3),
  convergence.rate = sqrt,
  method = c("quantile", "se"),
  include.boot = FALSE,
  n.cores = 1L,
  ...
)
}
\arguments{
\item{X, Y}{\code{n}-by-\code{d} and \code{m}-by-\code{d} data frames or matrices (multivariate
samples), or numeric/character vectors (univariate samples, i.e. \code{d = 1}),
representing \code{n} samples from the true distribution \eqn{P} and \code{m}
samples from the approximate distribution \eqn{Q} in \code{d} dimensions.
\code{Y} can be left blank if \code{q} is specified (see below).}

\item{q}{The density function of the approximate distribution \eqn{Q}. Either
\code{Y} or \code{q} must be specified. If the distributions are all continuous or
all discrete, \code{q} can be directly specified as the probability density/mass
function. However, for mixed continuous/discrete distributions, \code{q} must
be given in decomposed form, \eqn{q(y_c,y_d)=q_{c|d}(y_c|y_d)q_d(y_d)},
specified as a named list with field \code{cond} for the conditional density
\eqn{q_{c|d}(y_c|y_d)} (a function that expects two arguments \code{y_c} and
\code{y_d}) and \code{disc} for the discrete marginal density \eqn{q_d(y_d)} (a
function that expects one argument \code{y_d}). If such a decomposition is not
available, it may be preferable to instead simulate a large sample from
\eqn{Q} and use the two-sample syntax.}

\item{estimator}{The Kullback-Leibler divergence estimation method; a
function expecting two inputs (\code{X} and \code{Y} or \code{q}, depending on arguments
provided). Defaults to \code{kld_est_nn}.}

\item{B}{Number of bootstrap replicates (default: \code{500}), the larger, the
more accurate, but also more computationally expensive.}

\item{alpha}{Error level, defaults to \code{0.05}.}

\item{subsample.size}{A function specifying the size of the subsamples,
defaults to \eqn{f(x) = x^{2/3}}.}

\item{convergence.rate}{A function computing the convergence rate of the
estimator as a function of sample sizes. Defaults to \eqn{f(x) = x^{1/2}}.
If \code{convergence.rate} is \code{NULL}, it is estimated empirically from the
sample(s) using \code{kldest::convergence_rate()}.}

\item{method}{Either \code{"quantile"} (the default), also known as the reverse
percentile method, or \code{"se"} for a normal approximation of the KL
divergence estimator using the standard error of the subsamples.}

\item{include.boot}{Boolean, \code{TRUE} means KL divergence estimates on subsamples
are included in the returned list. Defaults to \code{FALSE}.}

\item{n.cores}{Number of cores to use in parallel computing (defaults to \code{1},
which means that no parallel computing is used).
To use this option, the \code{parallel} package must be installed and the OS
must be of UNIX type (i.e., not Windows). Otherwise, \code{n.cores} will be
reset to \code{1}, with a message.}

\item{...}{Arguments passed on to \code{estimator}, i.e. via the call
\code{estimator(X, Y = Y, ...)} or \code{estimator(X, q = q, ...)}.}
}
\value{
A list with the following fields:
\itemize{
\item \code{"est"} (the estimated KL divergence),
\item \code{"ci"} (a length \code{2} vector containing the lower and upper limits of the
estimated confidence interval).
\item \code{"boot"} (a length \code{B} numeric vector with KL divergence estimates on
the bootstrap subsamples), only included if \code{include.boot = TRUE},
}
}
\description{
This function computes a confidence interval for KL divergence based on the
subsampling bootstrap introduced by Politis and Romano. See \strong{Details} for
theoretical properties of this method.
}
\details{
In general terms, tetting \eqn{b_n} be the subsample size for a sample of
size \eqn{n}, and \eqn{\tau_n} the convergence rate of the estimator, a
confidence interval calculated by subsampling has asymptotic coverage
\eqn{1 - \alpha} as long as \eqn{b_n/n\rightarrow 0},
\eqn{b_n\rightarrow\infty} and \eqn{\frac{\tau_{b_n}}{\tau_n}\rightarrow 0}.

In many cases, the convergence rate of the nearest-neighbour based KL
divergence estimator is \eqn{\tau_n = \sqrt{n}} and the condition on the
subsample size reduces to \eqn{b_n/n\rightarrow 0} and \eqn{b_n\rightarrow\infty}.
By default, \eqn{b_n = n^{2/3}}. In a two-sample problem, \eqn{n} and \eqn{b_n}
are replaced by effective sample sizes \eqn{n_\text{eff} = \min(n,m)} and
\eqn{b_{n,\text{eff}} = \min(b_n,b_m)}.

Reference:

Politis and Romano, "Large sample confidence regions based on subsamples under
minimal assumptions", The Annals of Statistics, Vol. 22, No. 4 (1994).
}
\examples{
# 1D Gaussian (one- and two-sample problems)
set.seed(0)
X <- rnorm(100)
Y <- rnorm(100, mean = 1, sd = 2)
q <- function(x) dnorm(x, mean =1, sd = 2)
kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
kld_est_nn(X, Y = Y)
kld_est_nn(X, q = q)
kld_ci_subsampling(X, Y)$ci
kld_ci_subsampling(X, q = q)$ci

}
