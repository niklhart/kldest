# Kernel density-based Kullback-Leibler divergence estimators


#' 1-D kernel density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by a kernel density estimate using function 'density' from
#' package 'stats'. Only the two-sample, not the one-sample problem is implemented.
#'
#' @inherit kld_est_nn return
#' @param X,Y  Numeric vectors or single-column matrices, representing samples
#'    from the true distribution \eqn{P} and the approximate distribution
#'    \eqn{Q}, respectively.
#' @param MC A boolean: use a Monte Carlo approximation instead of numerical
#'    integration via the trapezoidal rule (default: `FALSE`)?
#' @param ... Further parameters to passed on to `stats::density` (e.g.,
#'    argument `bw`)
#' @examples
#' # KL-D between two samples from 1D Gaussians:
#' set.seed(0)
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kld_est_kde1(X,Y)
#' kld_est_kde1(X,Y, MC = TRUE)
#' @export
kld_est_kde1 <- function(X, Y, MC = FALSE, ...) {

    # input processing
    X <- as.vector(X)
    Y <- as.vector(Y)

    # using base R's density function
    dX <- density(X, ...)
    dY <- density(Y, from = head(dX$x,1), to = tail(dX$x,1), ...)

    if (MC) {
        pX <- approxfun(dX)
        pY <- approxfun(dY)

        # MC average
        mean(log(pX(X)/pY(X)))

    } else {
        # trapezoidal integration
        trapz(h  = dX$x[2] - dX$x[1],
              fx = dX$y*log(dX$y/dY$y))
    }
}



#' 2-D kernel density-based estimation of Kullback-Leibler divergence
#'
#' @description
#' This estimation method approximates the densities of the unknown bivariate
#' distributions \eqn{P} and \eqn{Q} by kernel density estimates using function
#' 'bkde' from package 'KernSmooth'. If 'KernSmooth' is not installed, a message
#' is issued and the (much) slower function 'kld_est_kde' is used instead.
#'
#' @inherit kld_est_nn return
#' @param X,Y  `n`-by-`2` and `m`-by-`2` matrices, representing `n` samples from
#'    the bivariate true distribution \eqn{P} and `m` samples from the approximate
#'    distribution \eqn{Q}, respectively.
#' @param MC A boolean: use a Monte Carlo approximation instead of numerical
#'    integration via the trapezoidal rule (default: `FALSE`)? Currently, this
#'    option is not implemented, i.e. a value of `TRUE` results in an error.
#' @param hX,hY Bandwidths for the kernel density estimates of \eqn{P} and \eqn{Q},
#'    respectively. The default `NULL` means they are determined by argument `rule`.
#' @param rule A heuristic to derive parameters `hX` and `hY`, default is
#'    `"Silverman", which means that
#'    `\deqn{h_i = \sigma_i\left(\frac{4}{(2+d)n}\right)^{1/(d+4)}.}
#' @param eps A nonnegative scalar; if `eps > 0`, \eqn{Q} is estimated as a mixture
#'    between the kernel density estimate and a uniform distribution on the computational
#'    grid. The weight of the uniform component is `eps` times the maximum density
#'    estimate of \eqn{Q}. This increases the robustness of the estimator at the
#'    expense of an additional bias. Defaults to `eps = 1e-5`.
#' @examples
#' # KL-D between two samples from 2-D Gaussians:
#' set.seed(0)
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' Y1 <- rnorm(1000)
#' Y2 <- Y1 + rnorm(1000)
#' X <- cbind(X1,X2)
#' Y <- cbind(Y1,Y2)
#' kld_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
#'              mu2 = rep(0,2), sigma2 = matrix(c(1,1,1,2),nrow=2))
#' kld_est_kde2(X,Y)
#' @export
kld_est_kde2 <- function(X, Y, MC = FALSE, hX = NULL, hY = NULL,
                        rule = c("Silverman","Scott"), eps = 1e-5) {


    # fallback if package 'KernSmooth' is not installed
    if (!requireNamespace("KernSmooth", quietly = TRUE)) {
        msg <- paste("Using the slower function 'kldest:::kld_est_kde' since",
                     "package 'KernSmooth' is not installed.")
        message(msg)

        return(kld_est_kde(X = X, Y = Y, hX = hX, hY = hY, rule = rule))
    }

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)

    # check that we are indeed in the 2D case
    d <- 2L
    stopifnot(ncol(X) == d, ncol(Y) == d)

    n <- nrow(X)
    m <- nrow(Y)

    # Use heuristics rule for bandwidth parameters, if unspecified
    rule <- match.arg(rule)
    if (is.null(hX)) {
        sd_X <- apply(X, MARGIN = 2, FUN = sd)
        hX <- switch(rule,
                     Silverman = sd_X*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_X/n^(1/(d+4))
        )
    }
    if (is.null(hY)) {
        sd_Y <- apply(Y, MARGIN = 2, FUN = sd)
        hY <- switch(rule,
                     Silverman = sd_Y*(4/((d+2)*m))^(1/(d+4)),
                     Scott     =  sd_Y/m^(1/(d+4))
        )
    }

    # using 2D kernel density estimation implemented in package 'KernSmooth'
    range.x <- list(range(X[,1]),range(X[,2]))
    dX <- KernSmooth::bkde2D(x = X, bandwidth = hX, range.x = range.x)
    dY <- KernSmooth::bkde2D(x = Y, bandwidth = hY, range.x = range.x, truncate = FALSE)

    if (MC) {
        stop("Not implemented yet...")
#        pX <- approxfun(dX)
#        pY <- approxfun(dY)
#        mean(log(pX(X)/pY(X)))

    } else {
        # numerical integration
        h <- c(dX$x1[2]-dX$x1[1], dX$x2[2] - dX$x2[1])
        if (eps > 0) {
            dY$fhat[dY$fhat == 0] <- eps * max(dY$fhat)
            normY <- trapz(h = h, fx = dY$fhat)
            dY$fhat <- dY$fhat/normY
        }
        fX <- dX$fhat * log(dX$fhat/dY$fhat)
        fX[dX$fhat == 0] <- 0
        trapz(h = h, fx = fX)
    }

}


#' Kernel density-based Kullback-Leibler divergence estimation in any dimension
#'
#' Disclaimer: this function doesn't use binning and/or the fast Fourier transform
#' and hence, it is extremely slow even for moderate datasets. For this reason,
#' it is not exported currently.
#'
#' This estimation method approximates the densities of the unknown distributions
#' \eqn{P} and \eqn{Q} by kernel density estimates, using a sample size- and
#' dimension-dependent bandwidth parameter and a Gaussian kernel. It works for
#' any number of dimensions but is very slow.
#'
#' @inherit kld_est_nn return
#' @param X,Y  `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#' @param hX,hY Positive scalars or length `d` vectors, representing bandwidth
#'    parameters (possibly different in each component) for the density estimates
#'    of \eqn{P} and \eqn{Q}, respectively. If unspecified, a heurestic specified
#'    via the `rule` argument is used.
#' @param rule A heuristic for computing arguments `hX` and/or `hY`. The default
#'    `"silverman"` is Silverman's rule
#'    \deqn{h_i = \sigma_i\left(\frac{4}{(2+d)n}\right)^{1/(d+4)}.}
#'    As an alternative, Scott's rule `"scott"` can be used,
#'    \deqn{h_i = \frac{\sigma_i}{n^{1/(d+4)}}.}
#' @example examples/continuous-estimators.R
kld_est_kde <- function(X, Y, hX = NULL, hY = NULL, rule = c("Silverman","Scott")) {

    # input processing
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    d <- ncol(X)
    m <- nrow(Y)
    stopifnot(ncol(Y) == d)

    # Use heuristics rule for bandwidth parameters, if unspecified
    rule <- match.arg(rule)
    if (is.null(hX)) {
        sd_X <- apply(X, MARGIN = 2, FUN = sd)
        hX <- switch(rule,
                     Silverman = sd_X*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_X/n^(1/(d+4))
        )
    }
    if (is.null(hY)) {
        sd_Y <- apply(Y, MARGIN = 2, FUN = sd)
        hY <- switch(rule,
                     Silverman = sd_Y*(4/((d+2)*n))^(1/(d+4)),
                     Scott     =  sd_Y/n^(1/(d+4))
        )
    }

    kXX <- vapply(1:d,
                  FUN = function(k) dnorm(outer(X[ ,k], X[ ,k], FUN = "-"),
                                          sd = hX[k]),
                  FUN.VALUE = matrix(0,nrow=n,ncol=n))
    kXY <- vapply(1:d,
                  FUN = function(k) dnorm(outer(X[ ,k], Y[ ,k], FUN = "-"),
                                          sd = hY[k]),
                  FUN.VALUE = matrix(0,nrow=n,ncol=m))

    p_X <- rowMeans(apply(kXX, MARGIN = 1:2, FUN = prod))
    q_X <- rowMeans(apply(kXY, MARGIN = 1:2, FUN = prod))

    # KL-divergence using estimated densities
    mean(log(p_X/q_X))

}
