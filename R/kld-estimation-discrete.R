# Kullback-Leibler divergence estimators between discrete distributions

#' Plug-in KL divergence estimator for samples from discrete distributions
#'
#' @inherit kld_est_kde return
#' @param X,Y Two samples from discrete distributions, specified as vectors,
#'    matrices or data frames. Argument `Y` can be omitted if argument `q` is
#'    given (see below).
#' @param q The probability mass function of the approximate distribution \eqn{Q}.
#' @examples
#' # 1D example
#' X <- c(rep('M',5),rep('F',5))
#' Y <- c(rep('M',6),rep('F',4))
#' kld_est_discrete(X, Y)
#' @export
kld_est_discrete <- function(X, Y = NULL, q = NULL) {

    # get important dimensions
    X <- as.data.frame(X)
    d <- ncol(X) # number of dimensions
    n <- nrow(X) # number of samples in X

    # one- or two-sample problem?
    two.sample <- !is.null(Y)
    if (!xor(is.null(Y),is.null(q))) stop("Either input `Y` or `q` must be specified.")
    if (!is.null(q) && d > 1) stop("One-sample version currently only implemented in 1D.")

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

        # TODO: the code below just works for a single discrete variable.
        #       for multiple discrete variables, compute one representative per
        #       interaction. I still have to think about the expected format.
        xcat <- uLvl[[1]]

        # estimate KL divergence via relative frequencies
        kld_discrete(tX/n, q(xcat))
    }

}


#' Zhang/Grabchak KL divergence estimator for samples from discrete distributions
#'
#' CAVE: not implemented yet!
#'
#' The estimator is that from Eq. (1.3) from Zhang and Grabchak (2014).
#'
#' Reference:
#' Zhang and Grabchak, "Nonparametric Estimation of Kullback-Leibler Divergence",
#' Neural Computation 26, 2570–2593 (2014).
#'
#' @inherit kld_est_kde return
#' @param X,Y Two samples from discrete distributions, specified as vectors,
#'    matrices or data frames.
#' @examples
#' # 1D example
#' X <- c(rep('M',5),rep('F',5))
#' Y <- c(rep('M',6),rep('F',4))
#' kld_est_discrete(X, Y)
#' @export
kld_est_zg2014 <- function(X, Y) {

    # get important dimensions
    X <- as.data.frame(X)
    Y <- as.data.frame(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d)

    # ensure consistency of factor levels
    uXY <- if (d == 1) {
        list(unique(c(X[[1]],Y[[1]])))
    } else {
        mapply(function(x,y) unique(c(x,y)),
               lapply(X, unique),
               lapply(Y, unique))
    }

    # convert all columns to factors to correctly account for missing levels
    X[] <- lapply(1:d, function(k) factor(X[[k]], levels = uXY[[k]]))
    Y[] <- lapply(1:d, function(k) factor(Y[[k]], levels = uXY[[k]]))

    # frequency counts via tables
    tX <- table(X)
    tY <- table(Y)

    pX <- stop("Implement this!")

    # estimate KL divergencce via relative frequencies
    kl_div_discrete(tX/n, tY/m)
}

