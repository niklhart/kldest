# Interfaces for Kullback-Leibler divergence estimation


#' Kullback-Leibler divergence estimator for discrete, continuous or mixed data.
#'
#' For two mixed continuous/discrete distributions with densities \eqn{p} and
#' \eqn{q}, and denoting \eqn{x = (x_\text{c},x_\text{d})}, the Kullback-Leibler
#' divergence \eqn{D_{KL}(p||q)} is given as
#' \deqn{D_{KL}(p||q) = \sum_{x_d} \int p(x_c,x_d) \log\left(\frac{p(x_c,x_d)}{q(x_c,x_d)}\right)dx_c.}
#' Conditioning on the discrete variables \eqn{x_d}, this can be re-written as
#' \deqn{D_{KL}(p||q) = \sum_{x_d} p(x_d) D_{KL}\big(p(\cdot|x_d)||q(\cdot|x_d)\big) +
#' D_{KL}\big(p_{x_d}||q_{x_d}\big).}
#' Here, the terms
#' \deqn{D_{KL}\big(p(\cdot|x_d)||q(\cdot|x_d)\big)}
#' are approximated via nearest neighbour- or kernel-based density estimates on
#' the datasets `X` and `Y` stratified by the discrete variables, and
#' \deqn{D_{KL}\big(p_{x_d}||q_{x_d}\big)}
#' is approximated using relative frequencies.
#'
#' @inherit kld_est_kde return
#' @param X,Y `n`-by-`d` and `m`-by-`d` data frames or matrices (multivariate
#'    samples), or numeric/character vectors (univariate samples, i.e. `d = 1`),
#'    representing `n` samples from the true distribution \eqn{P} and `m`
#'    samples from the approximate distribution \eqn{Q} in `d` dimensions.
#'    `Y` can be left blank if `q` is specified (see below).
#' @param q The density function of the approximate distribution \eqn{Q}. Either
#'    `Y` or `q` must be specified. If the distributions are all continuous or
#'    all discrete, `q` can be directly specified as the probability density/mass
#'    function. However, for mixed continuous/discrete distributions, `q` must
#'    be given in decomposed form, \eqn{q(y_c,y_d)=q_{c|d}(y_c|y_d)q_d(y_d)},
#'    specified as a named list with field `cond` for the conditional density
#'    \eqn{q_{c|d}(y_c|y_d)} (a function that expects two arguments `y_c` and
#'    `y_d`) and `disc` for the discrete marginal density \eqn{q_d(y_d)} (a
#'    function that expects one argument `y_d`). If such a decomposition is not
#'    available, it may be preferable to instead simulate a large sample from
#'    \eqn{Q} and use the two-sample syntax.
#' @param estimator.continuous,estimator.discrete KL divergence estimators for
#'    continuous and discrete data, respectively. Both are functions with two
#'    arguments `X` and `Y` or `X` and `q`, depending on whether a two-sample or
#'    one-sample problem is considered. Defaults are `kld_est_nn` and
#'    `kld_est_discrete`, respectively.
#' @param vartype A length `d` character vector, with `vartype[i] = "c"` meaning
#'    the `i`-th variable is continuous, and `vartype[i] = "d"` meaning it is
#'    discrete. If unspecified, `vartype` is `"c"` for numeric columns and `"d"`
#'    for character or factor columns. This default will mostly work, except if
#'    levels of discrete variables are encoded using numbers (e.g., `0` for
#'    females and `1` for males) or for count data.
#' @examples
#' # 2D example, two samples
#' set.seed(0)
#' X <- data.frame(cont  = rnorm(10),
#'                 discr = c(rep('a',4),rep('b',6)))
#' Y <- data.frame(cont  = c(rnorm(5), rnorm(5, sd = 2)),
#'                 discr = c(rep('a',5),rep('b',5)))
#' kld_est(X, Y)
#'
#' # 2D example, one sample
#' set.seed(0)
#' X <- data.frame(cont  = rnorm(10),
#'                 discr = c(rep(0,4),rep(1,6)))
#' q <- list(cond = function(xc,xd) dnorm(xc, mean = xd, sd = 1),
#'           disc = function(xd) dbinom(xd, size = 1, prob = 0.5))
#' kld_est(X, q = q, vartype = c("c","d"))
#' @export
kld_est <- function(X, Y = NULL, q = NULL, estimator.continuous = kld_est_nn,
                    estimator.discrete = kld_est_discrete, vartype = NULL) {

    # guess vartype from X
    if (is.null(vartype)) {
        vartype <- if (is.numeric(X)) {
            "c"
        } else if (is.character(X) || is.factor(X)) {
            "d"
        } else if (is.data.frame(X)) {
            ifelse(vapply(X, is.numeric, logical(1)),
                   yes = "c",
                   no  = "d")
        } else stop("Invalid class of input X.")
    }

    # one- or two-sample problem?
    two.sample <- is_two_sample(Y, q)

    # handle all-continuous or all-discrete variable case
    if (all(vartype == "c")) {
        if (two.sample) {
            return(estimator.continuous(X, Y = Y))
        } else {
            return(estimator.continuous(X, q = if (is.list(q)) q$cond else q))
        }
    } else if (all(vartype == "d")) {
        if (two.sample) {
            return(estimator.discrete(X, Y = Y))
        } else {
            return(estimator.discrete(X, q = if (is.list(q)) q$disc else q))
        }
    }

    # now we know it's a mixed discrete/continuous dataset
    Xcont  <- as.data.frame(X[vartype == "c"])
    Xdisc  <- X[vartype == "d"]
    iXdisc <- interaction(Xdisc, drop = TRUE)
    sXcont <- lapply(X = split(Xcont, f = iXdisc), FUN = as.matrix)

    if (two.sample) {
        Ycont <- as.data.frame(Y[vartype == "c"])
        Ydisc <- Y[vartype == "d"]
        # deal with possibility of missing levels
        iYdisc <- interaction(Ydisc, drop = TRUE)
        allLevels <- union(levels(iXdisc),levels(iYdisc))
        levels(iXdisc) <- allLevels
        levels(iYdisc) <- allLevels

        sXcont <- lapply(X = split(Xcont, f = iXdisc), FUN = as.matrix)
        sYcont <- lapply(X = split(Ycont, f = iYdisc), FUN = as.matrix)

        KLcont <- mapply(FUN = estimator.continuous, X = sXcont, Y = sYcont)
        KLdisc <- estimator.discrete(X = Xdisc, Y = Ydisc)

    } else {
        sXcont <- lapply(X = split(Xcont, f = iXdisc), FUN = as.matrix)
        qXcond <- lapply(Xdisc[match(levels(iXdisc),iXdisc), ],
                         function(d) {force(d); function(c) q$cond(c,d)})
        KLcont <- mapply(FUN = estimator.continuous, X = sXcont, q = qXcond)
        KLdisc <- estimator.discrete(X = Xdisc, q = q$disc)
    }

    # return compound KL divergence
    tdisc <- as.vector(table(iXdisc))
    KLcont[tdisc == 0] <- 0           # ensure "0 * NA = 0" (missing levels in X are ok)
    n <- sum(tdisc)
    sum(KLcont * tdisc/n) + KLdisc

}
