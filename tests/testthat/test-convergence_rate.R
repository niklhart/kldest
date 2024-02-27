test_that("Function convergence_rate works in a simple 1D problem", {

    set.seed(0)

    n <- 1000
    X <- rnorm(n)
    Y <- rnorm(n, mean = 1, sd = 2)
    q <- function(x) dnorm(x, mean =1, sd = 2)
    estimator <- kld_est_nn

    # 1-/2-sample problems
    beta1s <- convergence_rate(estimator = kld_est_nn, X = X, q = q, B = 500)
    beta2s <- convergence_rate(estimator = kld_est_nn, X = X, Y = Y, B = 500)

    # check against theoretical rate (with tolerance)
    expect_equal(beta1s, 0.5, tolerance = 0.03)
    expect_equal(beta2s, 0.5, tolerance = 0.04)

})
