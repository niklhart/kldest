
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kldest: Kullback-Leibler divergence estimation

<!-- badges: start -->

[![R-CMD-check](https://github.com/niklhart/kldest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/niklhart/kldest/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/niklhart/kldest/branch/master/graph/badge.svg)](https://app.codecov.io/gh/niklhart/kldest?branch=master)
<!-- badges: end -->

The goal of `kldest` is to estimate Kullback-Leibler (KL) divergence
$D_{KL}(P||Q)$ between two probability distributions $P$ and $Q$ based
on:

- a sample $x_1,...,x_n$ from $P$ and the probability density $q$ of
  $Q$, or
- samples $x_1,...,x_n$ from $P$ and $y_1,...,y_m$ from $Q$.

The distributions $P$ and $Q$ may be uni- or multivariate, and they may
be discrete, continuous or mixed discrete/continuous.

Different estimation algorithms are provided for continuous
distributions, either based on nearest neighbour density estimation or
kernel density estimation. Confidence intervals for KL divergence can
also be computed, either via subsampling (preferred) or bootstrapping.

## Installation

You can install kldest from CRAN:

``` r
install.packages("kldest")
```

Alternatively, can install the development version of kldest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("niklhart/kldest")
```

## A minimal example for KL divergence estimation

KL divergence estimation based on nearest neighbour density estimates is
the most flexible approach.

``` r
library(kldest)
```

Set a seed for reproducibility

``` r
set.seed(0)
```

### KL divergence between 1-D Gaussians

Analytical KL divergence:

``` r
kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#> [1] 0.4431472
```

Estimate based on two samples from these Gaussians:

``` r
X <- rnorm(100)
Y <- rnorm(100, mean = 1, sd = 2)
kld_est_nn(X, Y)
#> [1] 0.2169136
```

Estimate based on a sample from the first Gaussian and the density of
the second:

``` r
q <- function(x) dnorm(x, mean = 1, sd =2)
kld_est_nn(X, q = q)
#> [1] 0.6374628
```

Uncertainty quantification via subsampling:

``` r
kld_ci_subsampling(X, q = q)
#> $est
#> [1] 0.6374628
#> 
#> $ci
#>      2.5%     97.5% 
#> 0.2601375 0.9008446
```

### KL divergence between 2-D Gaussians

Analytical KL divergence between an uncorrelated and a correlated
Gaussian:

``` r
kld_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
             mu2 = rep(0,2), sigma2 = matrix(c(1,1,1,2),nrow=2))
#> [1] 0.5
```

Estimate based on two samples from these Gaussians:

``` r
X1 <- rnorm(100)
X2 <- rnorm(100)
Y1 <- rnorm(100)
Y2 <- Y1 + rnorm(100)
X <- cbind(X1,X2)
Y <- cbind(Y1,Y2)

kld_est_nn(X, Y)
#> [1] 0.3358918
```
