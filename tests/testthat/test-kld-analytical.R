test_that("Discrete KL-D calculation works", {
    # 1D example
    P <- c(1,0)
    Q <- c(exp(-1),1-exp(-1))
    D_KL <- kld_discrete(P,Q)
    expect_equal(D_KL, 1)

    # 2D example
    P2 <- matrix(c(0.5,0,0,0.5),nrow=2)
    Q2 <- matrix(c(0.5*exp(-1),1-exp(-1),0,0.5*exp(-1)),nrow=2)
    D_KL2 <- kld_discrete(P2,Q2)
    expect_equal(D_KL2, 1)

    # Infinite KL-D
    P3 <- c(0.5,0.5)
    Q3 <- c(0,1)
    D_KL3 <- kld_discrete(P3,Q3)
    expect_equal(D_KL3, Inf)

    # Invalid probabilities
    P4 <- c(1.2,-0.2)
    Q4 <- c(0.5,0.5)
    expect_error(kld_discrete(P4,Q4), "Input arrays must be nonnegative.")
    expect_error(kld_discrete(1, Q4), "Inputs must have the same dimensions.")
    expect_error(kld_discrete(1, 2),  "Input arrays must sum up to 1.")

})

test_that("KL-D of independent Gaussians is additive", {

    sigma1 <- matrix(
        c(2,1,0,0,
          1,2,0,0,
          0,0,2,1,
          0,0,1,2),
        nrow=4)
    sigma2 <- diag(4)
    mu <- rep(0,4)

    KL_4D <- kld_gaussian(mu1 = mu, sigma1 = sigma1,
                          mu2 = mu, sigma2 = sigma2)
    KL_2D <- kld_gaussian(mu1 = rep(0,2), sigma1 = constDiagMatrix(dim=2, diag=2,offDiag=1),
                          mu2 = rep(0,2), sigma2 = diag(2))
    expect_equal(KL_4D,2*KL_2D)

})


test_that("KL-D is 0 for identical distribution in analytical formulas", {

    # uniform distribution
    a <- -1
    b <-  2
    KL_num <- kld_uniform(a1=a,b1=b,
                             a2=a,b2=b)
    expect_equal(KL_num, 0)

    # exponential distribution
    lambda <- 2
    KL_num <- kld_exponential(lambda1 = lambda, lambda2 = lambda)
    expect_equal(KL_num, 0)

    # 3D normal distribution
    mu <- 1:3
    sigma <- constDiagMatrix(dim=3, diag = 7, offDiag = 2)

    KL_num <- kld_gaussian(mu1 = mu, sigma1 = sigma,
                           mu2 = mu, sigma2 = sigma)
    expect_equal(KL_num, 0)

})

test_that("KL-D between uniform and Gaussian behaves as expected", {

    # sd scaling
    KL_1 <- kld_uniform_gaussian(a = -1,
                                 b = 1,
                                 mu = 0,
                                 sigma2 = 1)
    KL_10 <- kld_uniform_gaussian(a = -10,
                                 b = 10,
                                 mu = 0,
                                 sigma2 = 100)
    expect_equal(KL_1, KL_10)

    # erroneous input
    expect_error(kld_uniform_gaussian(a = 1))
    expect_error(kld_uniform_gaussian(sigma2 = -1))

})

