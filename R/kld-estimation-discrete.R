# Kullback-Leibler divergence estimators between discrete distributions

#' Plug-in KL divergence estimator for samples from discrete distributions
#'
#' @inherit kld_est_kde return
#' @param X,Y  `n`-by-`d` and `m`-by-`d` matrices or data frames, representing
#'    `n` samples from the true discrete distribution \eqn{P} and `m` samples from
#'    the approximate discrete distribution \eqn{Q}, both in `d` dimensions.
#'    Vector input is treated as a column matrix. Argument `Y` can be omitted if
#'    argument `q` is given (see below).
#' @param q The probability mass function of the approximate distribution
#'    \eqn{Q}. Currently, the one-sample problem is only implemented for `d=1`.
#' @examples
#' # 1D example, two samples
#' X <- c(rep('M',5),rep('F',5))
#' Y <- c(rep('M',6),rep('F',4))
#' kld_est_discrete(X, Y)
#'
#' # 1D example, one sample
#' X <- c(rep(0,4),rep(1,6))
#' q <- function(x) dbinom(x, size = 1, prob = 0.5)
#' kld_est_discrete(X, q = q)
#'
#' @export
kld_est_discrete <- function(X, Y = NULL, q = NULL) {

    # get important dimensions
    X <- as.data.frame(X)
    d <- ncol(X) # number of dimensions
    n <- nrow(X) # number of samples in X

    # one- or two-sample problem?
    two.sample <- is_two_sample(Y,q)

    if (two.sample) {
        Y <- as.data.frame(Y)
        m <- nrow(Y) # number of samples in Y

        # check consistency of dimensions
        stopifnot(ncol(Y) == d)

        # ensure consistency of factor levels
        uLvl <- if (d == 1) {
            list(unique(c(X[[1]],Y[[1]])))
        } else {
            mapply(function(x,y) unique(c(x,y)),
                   lapply(X, unique),
                   lapply(Y, unique),
                   SIMPLIFY = FALSE)
        }
    } else {
        uLvl <- if (d == 1) {
            list(unique(X[[1]]))
        } else {
            lapply(X, unique)
        }
    }

    # convert all columns to factor to correctly account for missing levels
    X[] <- lapply(1:d, function(k) factor(X[[k]], levels = uLvl[[k]]))

    # frequency counts via tables
    tX <- table(X)

    if (two.sample) {

        # repeat the above for Y
        Y[] <- lapply(1:d, function(k) factor(Y[[k]], levels = uLvl[[k]]))
        tY <- table(Y)

        # estimate KL divergence via relative frequencies
        kld_discrete(tX/n, tY/m)

    } else {

        if (d > 1) stop("One-sample version currently only implemented in 1D.")
        # The code below just works for a single discrete variable. For multiple
        # discrete variables, compute one representative per interaction.
        # I still have to think about the expected format.
        xcat <- uLvl[[1]]

        # estimate KL divergence via relative frequencies
        kld_discrete(tX/n, q(xcat))
    }

}

