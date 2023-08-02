# Assessing uncertainty of KL-divergence estimation

# .632 bootstrapping?
# sampling from distributions (in particular q)?
# combination of both?

#' Uncertainty of KL divergence estimate using the subsampling bootstrap.
#'
#' This function computes a confidence interval for KL divergence based on the
#' subsampling bootstrap by Politis and Romano. The calculated interval has
#' asymptotic coverage \eqn{1 - \alpha} as long as \eqn{b_n/n\rightarrow 0},
#' \eqn{b_n\rightarrow\infty} and \eqn{\frac{\tau_b}{\tau_n}\rightarrow 0}.
#'
#' Reference:
#' Politis and Romano, "Large sample confidence regions based on subsamples under
#' minimal assumptions", The Annals of Statistics, Vol. 22, No. 4 (1994).
#'
#' @param X,Y `n`-by-`d` and `m`-by-`d` matrices, representing `n` samples from
#'    the true distribution \eqn{P} and `m` samples from the approximate distribution
#'    \eqn{Q}, both in `d` dimensions. Vector input is treated as a column matrix.
#' @param estimator A function handle expecting two inputs `X` and `Y`, the
#'     Kullback-Leibler divergence estimation method. Defaults to `kld_est_1nn`.
#' @param B Number of bootstrap replicates (default: `100`), the larger, the
#'     more accurate, but also more computationally expensive.
#' @param alpha Error level, defaults to `0.05`.
#' @param size A function specifying the size of the subsamples, defaults to
#'     \eqn{f(x) = x^{2/3}}.
#' @param rate A function computing the convergence rate of the estimator as a
#'     function of sample sizes. Defaults to `sqrt`.
#' @returns A list with the fields `"kld"` (the estimated KL divergence),
#'    `"boot"` (a length `B` numeric vector with KL divergence estimates on
#'    the bootstrap subsamples), and `"ci"` (a length `2` vector containing the
#'    lower and upper limits of the estimated confidence interval).
#' @examples
#' # 1D Gaussian
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1, sd = 2)
#' kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
#' kld_est_1nn(X, Y)
#' kld_ci_subboot(X, Y)
#'
#' @export
kld_ci_subboot <- function(X, Y, estimator = kld_est_1nn, B = 100L, alpha = 0.05,
                           size = function(x) x^(2/3), rate = sqrt) {

    X <- as.matrix(X)
    Y <- as.matrix(Y)

    n <- nrow(X)
    m <- nrow(Y)

    stopifnot(n >= 5, m >= 5)

    kld_hat  <- estimator(X, Y)
    kld_boot <- numeric(B)

    sn <- size(n)
    sm <- size(m)
    seff <- min(sn,sm)
    neff <- min(n,m)

    for (b in 1:B) {
        iX <- sample.int(n, size = sn)
        iY <- sample.int(m, size = sm)
        kld_boot[b] <- estimator(X[iX, ], Y[iY, ])
    }

    # computation of confidence levels
    z_star <- rate(seff) * (kld_boot - kld_hat)
    crit_val <- quantile(z_star, probs = c(1-alpha/2, alpha/2), na.rm = TRUE)
    ci_boot <- kld_hat - crit_val/rate(neff)

    list(
        dir  = kld_hat,
        boot = kld_boot,
        ci   = ci_boot
    )
}
