context("Correctness of estimation methods")

test_that("all estimation algorithms work well in an easy 1D example", {

    set.seed(123456)

    n <- 10000
    m <- 10001
    mu1 <- 0
    mu2 <- 5
    sd1 <- 1
    sd2 <- 2
    X <- rnorm(n, mean = mu1, sd = sd1)
    Y <- rnorm(m, mean = mu2, sd = sd2)

    kld_ref  <- kl_div_gaussian(mu1 = mu1, sigma1 = sd1^2, mu2 = mu2, sigma2 = sd2^2)
    kld_est1 <- kldest_density(X,Y)
    kld_est2 <- kl_universal_1nn(X,Y)
    kld_est3 <- kl_generalized_knn_eps(X,Y, warn.max.k = FALSE)

    # large
    tol <- 0.05
    expect_equal(kld_ref, kld_est1, tolerance = tol)
    expect_equal(kld_ref, kld_est2, tolerance = tol)
    expect_equal(kld_ref, kld_est3, tolerance = tol)
})
