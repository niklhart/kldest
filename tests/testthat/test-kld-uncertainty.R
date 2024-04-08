
test_that("kld_ci_subsampling has the correct coverage on an easy example", {

    # input parameters
    set.seed(123456)
    n <- 100
    alpha <- 0.4
    KL_true <- kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
    nrep <- 20

    # 1-/2-sample cases
    q <- function(x) dnorm(x, mean =1, sd = 2)
    cov_1s <- numeric(nrep)
    cov_2s <- numeric(nrep)
    cov_1s_se <- numeric(nrep)
    cov_2s_se <- numeric(nrep)
    for (i in 1:nrep) {
        X <- rnorm(n)
        Y <- rnorm(n, mean = 1, sd = 2)
        ci_1s <- kld_ci_subsampling(X, q = q, B = 100, alpha = alpha)$ci
        ci_2s <- kld_ci_subsampling(X, Y = Y, B = 100, alpha = alpha)$ci
        ci_1s_se <- kld_ci_subsampling(X, q = q, B = 100, alpha = alpha, method = "se")$ci
        ci_2s_se <- kld_ci_subsampling(X, Y = Y, B = 100, alpha = alpha, method = "se")$ci
        cov_1s[i] <- ci_1s[1] <= KL_true & ci_1s[2] >= KL_true
        cov_2s[i] <- ci_2s[1] <= KL_true & ci_2s[2] >= KL_true
        cov_1s_se[i] <- ci_1s_se[1] <= KL_true & ci_1s_se[2] >= KL_true
        cov_2s_se[i] <- ci_2s_se[1] <= KL_true & ci_2s_se[2] >= KL_true
    }
    # 20% tolerance since nrep and B are very small
    expect_equal(mean(cov_2s), 1-alpha, tolerance = 0.2)
    expect_equal(mean(cov_1s), 1-alpha, tolerance = 0.2)
    expect_equal(mean(cov_2s_se), 1-alpha, tolerance = 0.2)
    expect_equal(mean(cov_1s_se), 1-alpha, tolerance = 0.2)

})


test_that("kld_ci_subsampling parallelization works", {

    # NOTE: there is still randomness in the test due to different seeds used in
    # the parallel processes. I use a huge value of B and a large tolerance to
    # compensate for this. But it would be better to have a truly reproducible
    # way of handling this

    # input parameters
    set.seed(123456)
    n <- 100
    alpha <- 0.4
    B <- 5000
    X <- rnorm(n)

    # 1-sample case only
    q <- function(x) dnorm(x, mean =1, sd = 2)

    ci_1s_seq <- kld_ci_subsampling(X, q = q, B = B, alpha = alpha, n.cores = 1)$ci
    ci_1s_par <- kld_ci_subsampling(X, q = q, B = B, alpha = alpha, n.cores = 2)$ci

    expect_equal(ci_1s_par,ci_1s_seq, tolerance = 0.05)

})


test_that("kld_ci_subsampling works with mixed data", {

    # input parameters
    set.seed(123456)
    n <- 100
    alpha <- 0.4

    # 2-sample case only
    X <- data.frame(c = rnorm(n),
                    d = rbinom(n, size = 1, prob = 0.5))
    Y <- data.frame(c = rnorm(n, mean = 1, sd = 2),
                    d = rbinom(n, size = 1, prob = 0.3))

    # only check that subsampling runs without errors, no reference value
    expect_no_error(kld_ci_subsampling(X, Y = Y, B = 50, alpha = alpha, estimator = kld_est,
                                       vartype = c("c","d")))

})


test_that("kld_ci_bootstrap has the correct coverage on an easy example", {

    # input parameters
    set.seed(123456)
    n <- 100
    alpha <- 0.4
    KL_true <- kld_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
    nrep <- 20

    # 2-sample case only
    cov_2s <- numeric(nrep)
    for (i in 1:nrep) {
        X <- rnorm(n)
        Y <- rnorm(n, mean = 1, sd = 2)
        ci_2s <- kld_ci_bootstrap(X, Y, B = 100, alpha = alpha)$ci
        cov_2s[i] <- ci_2s[1] <= KL_true & ci_2s[2] >= KL_true
    }
    expect_equal(mean(cov_2s), 1-alpha, tolerance = 0.1)

})

test_that("include_boot option works", {

    set.seed(123456)

    X <- 1:10
    Y <- 11:20
    B <- 20
    KL_boot <- kld_ci_bootstrap(X, Y, B = B, include.boot = TRUE)
    KL_subs <- kld_ci_subsampling(X, Y, B = B, include.boot = TRUE)

    expect_length(KL_boot$boot, n = B)
    expect_length(KL_subs$boot, n = B)
})


test_that("kld_ci_bootstrap can deal with duplicates", {

    # input parameters
    set.seed(123456)

    X <- rep(0:10,2)
    Y <- X
    KL_se <- kld_ci_bootstrap(X, Y, method = "se")
    KL_q  <- kld_ci_bootstrap(X, Y, method = "quantile")

    expect_equal(KL_se$est, 0)
    expect_equal(KL_q$est, 0)
    expect_false(any(is.na(KL_se$ci)))
    expect_false(any(is.na(KL_q$ci)))

})


test_that("kld_ci_subsampling handles convergence.rate = NULL correctly", {

    set.seed(123456)

    X <- 1:10
    Y <- 11:20
    B <- 20
    KL_rate <- kld_ci_subsampling(X, Y, B = B, convergence.rate = NULL)
    KL_half <- kld_ci_subsampling(X, Y, B = B)

    expect_equal(KL_rate$est, KL_half$est)
    expect_equal(KL_rate$ci, KL_half$ci, tolerance = 0.1)
})


