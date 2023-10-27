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
#' @param X,Y Data frames or matrices with the same number of columns `d`
#'    (multivariate samples), or numeric/character vectors (univariate samples,
#'    i.e. `d=1`), representing `n` samples from the true distribution \eqn{P}
#'    and `m` samples from the approximate distribution \eqn{Q}, both in `d`
#'    dimensions. `Y` can be left blank if `q` is specified (see below).
#' @param q The density function of the approximate distribution \eqn{Q}. Either
#'    `Y` or `q` must be specified. In general, `q` must be given in decomposed
#'    form, \eqn{q(y_c|y_d)q(y_d)}, specified as a named list with field `cond`
#'    for the conditional density \eqn{q(y_c|y_d)} (a function that expects two
#'    arguments `y_c` and `y_d`) and `disc` for the discrete marginal density
#'    \eqn{q(y_d)} (a function that expects one argument `y_d`). If this is not
#'    possible, it may be preferable to simulate a large sample from \eqn{Q} and
#'    use the two-sample syntax instead. For compatibility with the continuous
#'    and discrete one-sample estimators, if the sample(s) is/are all continuous
#'    or all discrete, instead of specifying `q` as a length 2 list, it may also
#'    be given as a function handle computing the continous or discrete density.
#' @param estimator.continuous,estimator.discrete KL divergence estimators for
#'    continuous and discrete data, respectively. Both are function with two
#'    arguments `X` and `Y` or `X` and `q`, depending on whether a two-sample or
#'    one-sample problem is considered. Defaults are `kld_est_nn` and
#'    `kld_est_discrete`, respectively.
#' @param vartype A length `d` character vector, with `vartype[i] = "c"` meaning
#'    the `i`-th variable is continuous, and `vartype[i] = "d"` meaning it is
#'    discrete. If unspecified, `vartype` is `"c"` for numeric columns and `"d"`
#'    for character or factor columns. This default will not work if levels of
#'    discrete variables are encoded using numbers (e.g., `0` for females and
#'    `1` for males) or for count data.
#' @examples
#' # 2D example, two samples
#' X <- data.frame(cont  = rnorm(10),
#'                 discr = c(rep('a',4),rep('b',6)))
#' Y <- data.frame(cont  = c(rnorm(5), rnorm(5, sd = 2)),
#'                 discr = c(rep('a',5),rep('b',5)))
#' kld_est(X, Y)
#'
#' # 2D example, one sample
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
        } else if (is.character(X)) {
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
    Xcont  <- as.matrix(X[vartype == "c"])
    Xdisc  <- X[vartype == "d"]
    iXdisc <- interaction(Xdisc, drop = TRUE)
    sXcont <- split(Xcont, f = iXdisc)

    if (two.sample) {
        Ycont <- as.matrix(Y[vartype == "c"])
        Ydisc <- Y[vartype == "d"]
        sYcont <- split(Ycont, f = interaction(Ydisc, drop = TRUE))

        KLcont <- mapply(FUN = estimator.continuous, X = sXcont, Y = sYcont)
        KLdisc <- estimator.discrete(Xdisc,Ydisc)

    } else {

        qXcond <- lapply(Xdisc[match(levels(iXdisc),iXdisc), ],
                         function(d) {force(d); function(c) q$cond(c,d)})
        KLcont <- mapply(FUN = estimator.continuous, X = sXcont, q = qXcond)
        KLdisc <- estimator.discrete(Xdisc, q = q$disc)
    }


    # return compound KL divergence
    tdisc <- table(iXdisc)
    n <- sum(tdisc)
    sum(KLcont * tdisc/n) + KLdisc

}