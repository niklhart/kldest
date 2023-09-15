test_that("all estimation algorithms work well in an easy 1D example", {

    set.seed(123456)

    n <- 10000
    m <- 10001
    mu1 <- 0
    mu2 <- 5
    sd1 <- 1
    sd2 <- 5
    X <- rnorm(n, mean = mu1, sd = sd1)
    Y <- rnorm(m, mean = mu2, sd = sd2)

    kld_ref  <- kld_gaussian(mu1 = mu1, sigma1 = sd1^2, mu2 = mu2, sigma2 = sd2^2)
    kld_est1 <- kld_est_kde1(X,Y)
    kld_est2 <- kld_est_nn(X,Y)
    kld_est3 <- kld_est_brnn(X,Y, warn.max.k = FALSE)

    # large tol, to account for estimation variance
    tol <- 0.05
    expect_equal(kld_est1, kld_ref, tolerance = tol)
    expect_equal(kld_est2, kld_ref, tolerance = tol)
    expect_equal(kld_est3, kld_ref, tolerance = tol)
})


test_that("The density based estimator agrees with a hardcoded result", {
    X <- 1:2
    Y <- X+1      # same SD as X, and m=n

    h <- abs(X[1]-X[2])/sqrt(2)*(2/3)^0.2  # Silverman formula for h

    d0 <- dnorm(0/h)
    d1 <- dnorm(1/h)
    d2 <- dnorm(2/h)

    KL_ref <- 0.5 * log((d0+d1)/(d1+d2))
    KL_num <- kld_est_kde(X, Y)

    expect_equal(KL_ref,KL_num)
})

test_that("The nearest neighbour based estimator agrees with hardcoded results", {

    X <- 1:2
    Y <- X+0.5      # same SD as X, and m=n

    KL_ref <- 0
    KL_num <- kld_est_nn(X, Y)

    expect_equal(KL_ref,KL_num)

    X <- c(1,2,4)
    Y <- c(1.5,3)

    KL_ref <- -log(2)
    KL_num <- kld_est_nn(X, Y)

    expect_equal(KL_ref,KL_num)

})

