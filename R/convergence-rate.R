# convergence rate of estimators


#' Empirical convergence rate of a KL divergence estimator
#'
#' Subsampling-based confidence intervals computed by `kld_ci_subsampling()`
#' require the convergence rate of the KL divergence estimator as an input. The
#' default rate of `0.5` assumes that the variance term dominates the bias term.
#' For high-dimensional problems, depending on the data, the convergence rate
#' might be lower. This function allows to empirically derive the convergence
#' rate.
#'
#' @inheritParams kld_est
#' @param estimator A KL divergence estimator.
#' @param n.sizes Number of different subsample sizes to use (default: `4`).
#' @param spacing.factor Multiplicative factor controlling the spacing of sample
#'     sizes (default: `1.5`).
#' @param typical.subsample A function that produces a typical subsample size,
#'     used as the geometric mean of subsample sizes (default: `sqrt(n)`).
#' @param B Number of subsamples to draw per subsample size.
#' @param plot A boolean (default: `FALSE`) controlling whether to produce a
#'     diagnostic plot visualizing the fit.
#' @returns A scalar, the parameter \eqn{\beta} in the empirical convergence
#'     rate \eqn{n^-\beta} of the `estimator` to the true KL divergence.
#'     It can be used in the `convergence.rate` argument of `kld_ci_subsampling()`
#'     as `convergence.rate = function(n) n^beta`.
#' @examples
#'     # NN method usually has a convergence rate around 0.5:
#'     set.seed(0)
#'     convergence_rate(kld_est_nn, X = rnorm(1000), Y = rnorm(1000, mean = 1, sd = 2))
#'
#' @details
#' References:
#'
#' Politis, Romano and Wolf, "Subsampling", Chapter 8 (1999), for theory.
#'
#' The implementation has been adapted from lecture notes by C. J. Geyer,
#' https://www.stat.umn.edu/geyer/5601/notes/sub.pdf
#'
#' @export
convergence_rate <- function(estimator, X, Y = NULL, q = NULL,
                             n.sizes = 4, spacing.factor = 1.5,
                             typical.subsample = function(n) sqrt(n),
                             B = 500L, plot = FALSE) {

    two.sample <- is_two_sample(Y, q)

    # important dimensions for X and Y
    if (is.vector(X)) X <- as.matrix(X)
    n <- nrow(X)

    # subsample rule
    subsample.sizes <- function(n) {
        b_min <- typical.subsample(n) / sqrt((n.sizes-1)*spacing.factor)
        if (b_min < 2) stop("Subsample size too small to continue.")
        floor(b_min * spacing.factor^(0:(n.sizes-1)))
    }

    if (two.sample) {

        Y <- as.matrix(Y) # number of dimensions must be the same in X and Y
        m <- nrow(Y)      # number of samples in Y

        # determine subsample sizes from input parameters
        bn <- subsample.sizes(n)
        bm <- subsample.sizes(m)
        b  <- pmin(bn,bm)

        theta.hat <- estimator(X, Y = Y)

        theta.star <- matrix(NA, B, n.sizes)
        for (i in 1:B) {
            X.star <- X
            Y.star <- Y
            # backwards nested subsampling
            for (j in n.sizes:1) {
                X.star <- sample(X.star, bn[j], replace = FALSE)
                Y.star <- sample(Y.star, bm[j], replace = FALSE)
                theta.star[i, j] <- estimator(X.star, Y = Y.star)
            }
        }

    } else { # one sample

        b <- subsample.sizes(n)

        theta.hat <- estimator(X, q = q)

        theta.star <- matrix(NA, B, n.sizes)
        for (i in 1:B) {
            X.star <- X
            # backwards nested subsampling
            for (j in n.sizes:1) {
                X.star <- sample(X.star, b[j], replace = FALSE)
                theta.star[i, j] <- estimator(X.star, q = q)
            }
        }
    }

    zmat <- theta.star - theta.hat

    # Catch issues with the estimator
    if (!is.finite(theta.hat)) return(NA_real_)

    # calculate quantile differences
    l_probs <- seq(0.05, 0.45, by = 0.05)
    u_probs <- seq(0.55, 0.95, by = 0.05)
    lqmat <- log(apply(zmat, MARGIN = 2, FUN = function(x) quantile(x, u_probs))
                 - apply(zmat, MARGIN = 2, FUN = function(x) quantile(x, l_probs)))
    dimnames(lqmat) <- list(NULL,b)

    y <- colMeans(lqmat)
    beta <- -cov(y, log(b)) / var(log(b))

    if (plot) {
        inter <- mean(y) + beta * mean(log(b))
        plot(rep(b, each = nrow(lqmat)), as.vector(lqmat),
             xlab = if (two.sample) "Effective subsample size" else "Subsample size",
             ylab = "log(high quantile - low quantile)",
             main = paste0("Empirical convergence rate (beta = ",signif(beta,3),")"),
             log = "x")
        lines(b, inter-beta*log(b), col = "red")
    }
    beta
}
