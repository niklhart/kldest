# Assessing uncertainty of KL-divergence estimation

#' Uncertainty of KL divergence estimate using Efron's bootstrap.
#'
#' This function computes a confidence interval for KL divergence based on Efron's
#' bootstrap. The approach only works for kernel density-based estimators since
#' nearest neighbour-based estimators cannot deal with the ties produced when
#' sampling with replacement.
#'
#' Reference:
#'
#' Efron, "Bootstrap Methods: Another Look at the Jackknife", The Annals of
#' Statistics, Vol. 7, No. 1 (1979).
#'
#' @param X,Y `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#' @param estimator A function expecting two inputs `X` and `Y`, the
#'     Kullback-Leibler divergence estimation method. Defaults to `kld_est_kde1`,
#'     which can only deal with one-dimensional two-sample problems (i.e.,
#'     `d = 1` and `q = NULL`).
#' @param B Number of bootstrap replicates (default: `500`), the larger, the
#'     more accurate, but also more computationally expensive.
#' @param alpha Error level, defaults to `0.05`.
#' @returns A list with the fields `"est"` (the estimated KL divergence),
#'    `"boot"` (a length `B` numeric vector with KL divergence estimates on
#'    the bootstrap samples), and `"ci"` (a length `2` vector containing the
#'    lower and upper limits of the estimated confidence interval).
#' @examples
#' # 1D Gaussian, two samples
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kld_est_kde1(X, Y)
#' kld_ci_bootstrap(X, Y)
#'
#' @export
kld_ci_bootstrap <- function(X, Y, estimator = kld_est_kde1, B = 500L, alpha = 0.05) {

    X <- as.matrix(X)
    Y <- as.matrix(Y)

    n <- nrow(X)
    m <- nrow(Y)
    d <- ncol(X)

    kld_hat <- estimator(X,Y)
    kld_boot <- numeric(B)

    for (b in 1:B) {
        iX <- sample.int(n, replace = TRUE)
        iY <- sample.int(m, replace = TRUE)

        kld_boot[b] <- estimator(X[iX, ], Y[iY, ])
    }

    # computation of confidence levels
    q_boot <- quantile(kld_boot, probs = c(1-alpha/2, alpha/2), na.rm = TRUE)
    ci_boot <- 2*kld_hat - q_boot

    list(
        est  = kld_hat,
        boot = kld_boot,
        ci   = ci_boot
    )
}



#' Uncertainty of KL divergence estimate using Politis/Romano's subsampling bootstrap.
#'
#' This function computes a confidence interval for KL divergence based on the
#' subsampling bootstrap by Politis and Romano. The calculated interval has
#' asymptotic coverage \eqn{1 - \alpha} as long as \eqn{b_n/n\rightarrow 0},
#' \eqn{b_n\rightarrow\infty} and \eqn{\frac{\tau_b}{\tau_n}\rightarrow 0}.
#'
#' Reference:
#'
#' Politis and Romano, "Large sample confidence regions based on subsamples under
#' minimal assumptions", The Annals of Statistics, Vol. 22, No. 4 (1994).
#'
#' @param X,Y `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#'    `Y` can be left blank if `q` is specified (see below).
#' @param q The density function of the approximate distribution \eqn{Q}. Either
#'    `Y` or `q` must be specified.
#' @param estimator The Kullback-Leibler divergence estimation method; a
#'    function expecting two inputs (`X` and `Y` or `q`, depending on arguments
#'    provided). Defaults to `kld_est_nn`.
#' @param B Number of bootstrap replicates (default: `500`), the larger, the
#'     more accurate, but also more computationally expensive.
#' @param alpha Error level, defaults to `0.05`.
#' @param subsample.size A function specifying the size of the subsamples,
#'     defaults to \eqn{f(x) = x^{2/3}}.
#' @param convergence.rate A function computing the convergence rate of the
#'     estimator as a function of sample sizes. Defaults to \eqn{f(x) = x^{1/2}}.
#' @returns A list with the fields `"est"` (the estimated KL divergence),
#'    `"boot"` (a length `B` numeric vector with KL divergence estimates on
#'    the bootstrap subsamples), and `"ci"` (a length `2` vector containing the
#'    lower and upper limits of the estimated confidence interval).
#' @examples
#' # 1D Gaussian (one- and two-sample problems)
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' q <- function(x) dnorm(x, mean =1, sd = 2)
#' kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kld_est_nn(X, Y = Y)
#' kld_est_nn(X, q = q)
#' kld_ci_subsampling(X, Y)$ci
#' kld_ci_subsampling(X, q = q)$ci
#'
#' @export
kld_ci_subsampling <- function(X, Y = NULL, q = NULL, estimator = kld_est_nn,
                               B = 500L, alpha = 0.05,
                               subsample.size = function(x) x^(2/3),
                               convergence.rate = sqrt) {

    # Important dimensions for X
    X <- as.matrix(X)
    n <- nrow(X)
    sn <- subsample.size(n)

    # check validity of input: one- or two-sample problem?
    two.sample <- is_two_sample(Y, q)

    if (two.sample) {
        Y <- as.matrix(Y)
        m <- nrow(Y)
        kld_hat  <- estimator(X, Y = Y)

        sm <- subsample.size(m)
        seff <- min(sn,sm)
        neff <- min(n,m)
    } else {
        kld_hat <- estimator(X, q = q)

        seff <- sn
        neff <- n
    }

    kld_boot <- numeric(B)

    for (b in 1:B) {
        iX <- sample.int(n, size = sn)
        if (two.sample) {
            iY <- sample.int(m, size = sm)
            kld_boot[b] <- estimator(X[iX, ], Y = Y[iY, ])
        } else {
            kld_boot[b] <- estimator(X[iX, ], q = q)
        }
    }

    # computation of confidence levels
    z_star <- convergence.rate(seff) * (kld_boot - kld_hat)
    crit_val <- quantile(z_star, probs = c(1-alpha/2, alpha/2), na.rm = TRUE)
    ci_boot <- kld_hat - crit_val/convergence.rate(neff)
    names(ci_boot) <- names(ci_boot)[2:1]

    list(
        est  = kld_hat,
        boot = kld_boot,
        ci   = ci_boot
    )
}