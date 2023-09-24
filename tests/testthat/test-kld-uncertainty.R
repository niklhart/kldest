
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
    for (i in 1:nrep) {
        X <- rnorm(n)
        Y <- rnorm(n, mean = 1, sd = 2)
        ci_1s <- kld_ci_subsampling(X, q = q, B = 100, alpha = alpha)$ci
        ci_2s <- kld_ci_subsampling(X, Y, B = 100, alpha = alpha)$ci
        cov_1s[i] <- ci_1s[1] <= KL_true & ci_1s[2] >= KL_true
        cov_2s[i] <- ci_2s[1] <= KL_true & ci_2s[2] >= KL_true
    }
    expect_equal(mean(cov_2s), 1-alpha, tolerance = 0.1)
    expect_equal(mean(cov_1s), 1-alpha, tolerance = 0.17)

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
