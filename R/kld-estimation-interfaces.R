# Interfaces for Kullback-Leibler divergence estimation

#' Kullback-Leibler divergence estimation between continuous distributions
#'
#' This function is an interface to several different estimation methods for
#' Kullback-Leibler divergence \eqn{D_{KL}(P||Q)} between two continuous
#' distributions \eqn{P} and \eqn{Q} based on a sample `X` from \eqn{P} and
#' either the density function `q` of \eqn{Q} (one-sample problem) or a sample
#' `Y` of \eqn{Q} (two-sample problem).
#'
#' @param method The estimation method to be used, either `nn` for nearest-neighbour
#'    based estimation (the default), `gnn` for generalized nearest neightbour
#'    estimation, `brnn` for bias-reduced generalized nearest neighbour estimation,
#'    `kde1` for kernel density based estimation in 1D or `kde2` for kernel density
#'    based estimation in 2D. The last option requires package `KernSmooth`
#'    to be installed.
#' @param ... Arguments to be passed to the respective method (see links below).
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{D_{KL}(P||Q)}.
#' @example examples/continuous-estimators.R
#' @seealso [kld_est_nn()], [kld_est_gnn()], [kld_est_brnn()], [kld_est_kde1()],
#'     [kld_est_kde2()] for possible arguments passed to these methods.
#' @export
kld_cont <- function(..., method = c("nn","gnn","brnn","kde1","kde2")) {

    switch(match.arg(method),
           nn  = kld_est_nn(...),
           gnn = kld_est_gnn(...),
           brnn = kld_est_brnn(...),
           kde1 = kld_est_kde1(...),
           kde2 = kld_est_kde2(...))
}

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
#' @param ... Inputs passed on to function [kld_cont()].
#' @param vartype A length `d` character vector, with `vartype[i] = "c"` meaning
#'    the `i`-th variable is continuous, and `vartype[i] = "d"` meaning it is
#'    discrete. If unspecified, `vartype` is `"c"` for numeric columns and `"d"`
#'    for character or factor columns. This default will not work if levels of
#'    discrete variables are encoded using numbers (e.g., `0` for females and
#'    `1` for males) or for count data.
#' @examples
#' # 2D example
#' X <- data.frame(cont  = rnorm(10),
#'                 discr = c(rep('a',4),rep('b',6)))
#' Y <- data.frame(cont  = c(rnorm(5), rnorm(5, sd = 2)),
#'                 discr = c(rep('a',5),rep('b',5)))
#' kld_est(X,Y)
#' @export
kld_est <- function(X, Y, ..., vartype = NULL) {

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

    # handle all-continuous or all-discrete variable case
    if (all(vartype == "c")) {
        return(kld_cont(X,Y, ...))
    } else if (all(vartype == "d")) {
        return(kld_est_discrete(X, Y))
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
    KLcont <- mapply(FUN = kld_cont, X = sXcont, Y = sYcont, MoreArgs = list(...))
    KLdisc <- kld_est_discrete(Xdisc,Ydisc)

    # return compound KL divergence
    tdisc <- table(iXdisc)
    n <- sum(tdisc)
    sum(KLcont * tdisc/n) + KLdisc

}
