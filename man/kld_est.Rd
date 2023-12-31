% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kld-estimation-interfaces.R
\name{kld_est}
\alias{kld_est}
\title{Kullback-Leibler divergence estimator for discrete, continuous or mixed data.}
\usage{
kld_est(
  X,
  Y = NULL,
  q = NULL,
  estimator.continuous = kld_est_nn,
  estimator.discrete = kld_est_discrete,
  vartype = NULL
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

\item{estimator.continuous, estimator.discrete}{KL divergence estimators for
continuous and discrete data, respectively. Both are functions with two
arguments \code{X} and \code{Y} or \code{X} and \code{q}, depending on whether a two-sample or
one-sample problem is considered. Defaults are \code{kld_est_nn} and
\code{kld_est_discrete}, respectively.}

\item{vartype}{A length \code{d} character vector, with \code{vartype[i] = "c"} meaning
the \code{i}-th variable is continuous, and \code{vartype[i] = "d"} meaning it is
discrete. If unspecified, \code{vartype} is \code{"c"} for numeric columns and \code{"d"}
for character or factor columns. This default will mostly work, except if
levels of discrete variables are encoded using numbers (e.g., \code{0} for
females and \code{1} for males) or for count data.}
}
\value{
A scalar, the estimated Kullback-Leibler divergence \eqn{\hat D_{KL}(P||Q)}.
}
\description{
For two mixed continuous/discrete distributions with densities \eqn{p} and
\eqn{q}, and denoting \eqn{x = (x_\text{c},x_\text{d})}, the Kullback-Leibler
divergence \eqn{D_{KL}(p||q)} is given as
\deqn{D_{KL}(p||q) = \sum_{x_d} \int p(x_c,x_d) \log\left(\frac{p(x_c,x_d)}{q(x_c,x_d)}\right)dx_c.}
Conditioning on the discrete variables \eqn{x_d}, this can be re-written as
\deqn{D_{KL}(p||q) = \sum_{x_d} p(x_d) D_{KL}\big(p(\cdot|x_d)||q(\cdot|x_d)\big) +
D_{KL}\big(p_{x_d}||q_{x_d}\big).}
Here, the terms
\deqn{D_{KL}\big(p(\cdot|x_d)||q(\cdot|x_d)\big)}
are approximated via nearest neighbour- or kernel-based density estimates on
the datasets \code{X} and \code{Y} stratified by the discrete variables, and
\deqn{D_{KL}\big(p_{x_d}||q_{x_d}\big)}
is approximated using relative frequencies.
}
\examples{
# 2D example, two samples
set.seed(0)
X <- data.frame(cont  = rnorm(10),
                discr = c(rep('a',4),rep('b',6)))
Y <- data.frame(cont  = c(rnorm(5), rnorm(5, sd = 2)),
                discr = c(rep('a',5),rep('b',5)))
kld_est(X, Y)

# 2D example, one sample
set.seed(0)
X <- data.frame(cont  = rnorm(10),
                discr = c(rep(0,4),rep(1,6)))
q <- list(cond = function(xc,xd) dnorm(xc, mean = xd, sd = 1),
          disc = function(xd) dbinom(xd, size = 1, prob = 0.5))
kld_est(X, q = q, vartype = c("c","d"))
}
