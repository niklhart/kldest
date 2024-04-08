test_that("All estimators for continuous data work well in an easy 1D example", {

    set.seed(123456)

    n <- 10000      # slight difference in sample size to detect
    m <- 10001      # coding problems for unequal sample sizes

    mu1 <- 0
    mu2 <- 5
    sd1 <- 1
    sd2 <- 5
    X <- rnorm(n, mean = mu1, sd = sd1)
    Y <- rnorm(m, mean = mu2, sd = sd2)
    q <- function(x) dnorm(x, mean = mu2, sd = sd2)

    kld_ref  <- kld_gaussian(mu1 = mu1, sigma1 = sd1^2, mu2 = mu2, sigma2 = sd2^2)

    kld_kde1a <- kld_est_kde1(X,Y)
    kld_kde1b <- kld_est_kde1(X,Y, MC = TRUE)

    kld_nn1a <- kld_est_nn(X,Y, eps = 0.01) # approximate 1-NN
    kld_nn1b <- kld_est_nn(X,Y, eps = 0)    # exact 1-NN
    kld_nn1c <- kld_est_nn(X,q = q)
    kld_nn1d <- kld_est_nn(X,q = function(x) log(q(x)), log.q = TRUE)
    kld_nn2  <- kld_est_nn(X,Y, k = 2)

    kld_brnn <- kld_est_brnn(X,Y, warn.max.k = FALSE)

    # large tol, to account for estimation variance
    tol <- 0.02
    expect_equal(kld_kde1a, kld_ref, tolerance = tol)
    expect_equal(kld_kde1b, kld_ref, tolerance = tol)
    expect_equal(kld_nn1a,  kld_ref, tolerance = tol)
    expect_equal(kld_nn1b,  kld_ref, tolerance = tol)
    expect_equal(kld_nn1c,  kld_ref, tolerance = tol)
    expect_equal(kld_nn1d,  kld_ref, tolerance = tol)
    expect_equal(kld_nn2,   kld_ref, tolerance = tol)
    expect_equal(kld_brnn,  kld_ref, tolerance = tol)

})

test_that("All estimators for continuous data work well in an easy 2D example", {

    set.seed(123456)

    n <- 10000      # slight difference in sample size to detect
    m <- 10001      # coding problems for unequal sample sizes

    X1 <- rnorm(n)
    X2 <- rnorm(n)
    X <- cbind(X1,X2)

    Y1 <- rnorm(m)
    Y2 <- Y1 + rnorm(m)
    Y <- cbind(Y1,Y2)

    mu     <- rep(0,2)
    Sigma1 <- diag(2)
    Sigma2 <- matrix(c(1,1,1,2),nrow=2)

    q <- function(x) mvdnorm(x, mu = mu, Sigma = Sigma2)

    kld_ref <- kld_gaussian(mu1 = mu, sigma1 = Sigma1,
                            mu2 = mu, sigma2 = Sigma2)
    kld_kde2 <- kld_est_kde2(X, Y)

    kld_nnXY <- kld_est_nn(X, Y)
    kld_nnXq <- kld_est_nn(X, q = q)
    kld_brnn <- kld_est_brnn(X,Y, warn.max.k = FALSE)

    # large tol, to account for estimation variance
    tol <- 0.1
    expect_equal(kld_kde2, kld_ref, tolerance = tol)
    expect_equal(kld_nnXY, kld_ref, tolerance = tol)
    expect_equal(kld_nnXq, kld_ref, tolerance = tol)
    expect_equal(kld_brnn, kld_ref, tolerance = tol)

})


test_that("Kernel density based estimator agrees with a hardcoded result in 1D", {
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

test_that("Nearest neighbour based estimator agrees with hardcoded results", {

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

test_that("Bias-reduced NN behaves as expected", {

    set.seed(123456)

    # if max.k == 1, bias reduced NN is identical to plain NN
    X <- rnorm(10, mean = 0, sd = 1)
    Y <- rnorm(11, mean = 5, sd = 5)

    KL_nn1 <- kld_est_nn(X, Y)
    KL_brnn <- kld_est_brnn(X, Y, max.k = 1, warn.max.k = FALSE)

    expect_equal(KL_nn1,KL_brnn)

    # if max.k is too small, BRNN will trigger a warning
    expect_warning(kld_est_brnn(X, Y, max.k = 1, warn.max.k = TRUE))


    # if n or m are too small, NA is returned
    KL_nn1_NA <- kld_est_nn(1, 1:10)
    KL_nn2_NA <- kld_est_nn(1:10, 1, k = 2)
    KL_brnn1_NA <- kld_est_nn(1, 1:10)
    KL_brnn2_NA <- kld_est_nn(1:10, numeric(0))

    expect_equal(KL_nn1_NA, NA_real_)
    expect_equal(KL_nn2_NA, NA_real_)
    expect_equal(KL_brnn1_NA, NA_real_)
    expect_equal(KL_brnn2_NA, NA_real_)

})


test_that("KL-D estimation for discrete variables works", {

    # 1D example: invariant to type of input
    Xn <- c(1,1,2,2)
    Yn <- c(1,2,2,2)
    Xc <- as.character(Xn)
    Yc <- as.character(Yn)
    Xd <- data.frame(Xn)
    Yd <- data.frame(Yn)

    KL_num <- kld_est_discrete(Xn, Yn)
    KL_chr <- kld_est_discrete(Xc, Yc)
    KL_df  <- kld_est_discrete(Xd, Yd)
    KL_ref <- kld_discrete(c(0.5,0.5), c(0.25,0.75))

    expect_equal(KL_num,KL_ref)
    expect_equal(KL_chr,KL_ref)
    expect_equal(KL_df, KL_ref)

    # 2D example
    X2 <- matrix(c(1,1,2,1,2,2),ncol=2)
    Y2 <- matrix(c(1,1,2,2,1,2,1,2),ncol=2)

    KL2_num <- kld_est_discrete(X2, Y2)
    KL2_ref <- kld_discrete(matrix(c(1,0,1,1)/3,nrow=2),
                           matrix(0.25,nrow=2,ncol=2))

    expect_equal(KL2_num,KL2_ref)

    # 1D example with unobserved factor levels
    Xu <- factor(rep(2,4), levels = 1:2)
    Yu <- factor(rep(1:2,2), levels = 1:2)

    KLu_num <- kld_est_discrete(Xu,Yu)
    KLu_ref <- log(2)

    expect_equal(KLu_num,KLu_ref)

    # 1D example with one sample
    X <- c(0,0,1,1,1)
    q <- function(x) dbinom(x, size = 1, prob = 0.5)

    KL_q   <- kld_est_discrete(X, q = q)
    KL_ref <- kld_discrete(c(0.4,0.6), c(0.5,0.5))

    expect_equal(KL_q,KL_ref)

    # 2D example with one sample errors, but check that it runs up to this point
    # (once the behaviour is implemented, expand this test!)
    expect_error(kld_est_discrete(X = X2, q = force),
                 "One-sample version currently only implemented in 1D.")

})


