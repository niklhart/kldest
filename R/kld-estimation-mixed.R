# Kullback-Leibler divergence estimators for mixed continuous/discrete distributions

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
#'    i.e. `d=1`).
#' @param method The estimation method to be used for the continuous variables,
#'    either `"1nn"` (default) for 1-nearest-neighbour density estimation,
#'    `"gknn"` for generalized k-nearest-neighbour density estimation
#'    or `"kde"` for kernel density estimation using a Gaussian kernel.
#' @param vartype A length `d` character vector, with `vartype[i] = "c"` meaning
#'    the `i`-th variable is continuous, and `vartype[i] = "d"` meaning it is
#'    discrete. If unspecified, `vartype` is `"c"` for numeric columns and `"d"`
#'    for character or factor columns. This default will not work if levels of
#'    discrete variables are encoded using numbers (e.g., `0` for females and
#'    `1` for males) or for count data.
#' @param ... Inputs passed on to density estimation functions.
#' @examples
#' # 2D example
#' X <- data.frame(cont  = rnorm(10),
#'                 discr = c(rep('a',4),rep('b',6)))
#' Y <- data.frame(cont  = c(rnorm(5), rnorm(5, sd = 2)),
#'                 discr = c(rep('a',5),rep('b',5)))
#' kld_est(X,Y)
#' @export
kld_est <- function(X, Y, method = c("1nn","gknn","kde"), vartype = NULL, ...) {

    # guess vartype from X and Y
    if (is.null(vartype)) {
        vartype <- if (is.numeric(X)) {
            "c"
        } else if (is.character(X)) {
            "d"
        } else if (is.data.frame(X)) {
            ifelse(vapply(X, is.numeric,logical(1)),
                   yes = "c",
                   no  = "d")
        } else stop("Invalid class of input X.")
    }

    # process method argument (not required for discrete variables)
    if (!all(vartype == "d")) {
        kld_est_continuous <- switch(match.arg(method),
                                    "1nn"  = kld_est_1nn,
                                    "gknn" = kld_est_gknn,
                                    "kde"  = kld_est_kde)
    }

    # handle all-continuous or all-discrete variable case
    if (all(vartype == "c")) {
        return(kldest_continuous(X,Y, ...))
    } else if (all(vartype == "d")) {
        return(kldest_discrete(X, Y))
    }

    # now we know it's a mixed discrete/continuous dataset
    Xcont  <- as.matrix(X[vartype == "c"])
    Xdisc  <- X[vartype == "d"]
    iXdisc <- interaction(Xdisc, drop = TRUE)
    sXcont <- split(Xcont, f = iXdisc)

    Ycont <- as.matrix(Y[vartype == "c"])
    Ydisc <- Y[vartype == "d"]
    sYcont <- split(Ycont, f = interaction(Ydisc, drop = TRUE))

    # compute KL divergence
    KLcont <- mapply(FUN = kld_est_continuous, X = sXcont, Y = sYcont, MoreArgs = list(...))
    KLdisc <- kld_est_discrete(Xdisc,Ydisc)

    # return compound KL divergence
    sum(KLcont * table(iXdisc)/n) + KLdisc

}
