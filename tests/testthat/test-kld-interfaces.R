test_that("kld_est works as expected for discrete data", {

    # 1D example (two samples): invariant to type of input
    Xn <- c(1,1,2,2)
    Yn <- c(1,2,2,2)
    Xc <- as.character(Xn)
    Yc <- as.character(Yn)
    Xf <- as.factor(Xn)
    Yf <- as.factor(Yn)
    Xnd <- data.frame(Xn)
    Ynd <- data.frame(Yn)
    Xcd <- data.frame(Xc)
    Ycd <- data.frame(Yc)

    KL_num <- kld_est(Xn, Yn, vartype = "d")
    KL_chr <- kld_est(Xc, Yc)
    KL_fct <- kld_est(Xf, Yf)
    KL_ndf <- kld_est(Xnd, Ynd, vartype = "d")
    KL_cdf <- kld_est(Xcd, Ycd)
    KL_ref <- kld_discrete(c(0.5,0.5), c(0.25,0.75))

    expect_equal(KL_num, KL_ref)
    expect_equal(KL_chr, KL_ref)
    expect_equal(KL_fct, KL_ref)
    expect_equal(KL_ndf, KL_ref)
    expect_equal(KL_cdf, KL_ref)

    # 2D example (two samples)
    X2n <- matrix(c(1,1,2,1,2,2),    ncol=2)
    Y2n <- matrix(c(1,1,2,2,1,2,1,2),ncol=2)
    X2c <- matrix(as.character(X2n),ncol=2)
    Y2c <- matrix(as.character(Y2n),ncol=2)
    X2dnn <- as.data.frame(X2n)
    Y2dnn <- as.data.frame(Y2n)
    X2dcc <- as.data.frame(X2c)
    Y2dcc <- as.data.frame(Y2c)
    X2dnc <- X2dnn; X2dnc$V1 <- as.character(X2dnn$V1)
    Y2dnc <- Y2dnn; Y2dnc$V1 <- as.character(Y2dnn$V1)

    KL2_num <- kld_est(X2n, Y2n, vartype = c("d","d"))
    KL2_chr <- kld_est(X2c, Y2c)
    KL2_dnn <- kld_est(X2dnn, Y2dnn, vartype = c("d","d"))
    KL2_dcc <- kld_est(X2dcc, Y2dcc)
    KL2_dnc <- kld_est(X2dnc, Y2dnc, vartype = c("d","d"))

    KL2_ref <- kld_discrete(matrix(c(1,0,1,1)/3,nrow=2),
                            matrix(0.25,nrow=2,ncol=2))

    expect_equal(KL2_num, KL2_ref)
    expect_equal(KL2_chr, KL2_ref)
    expect_equal(KL2_dnn, KL2_ref)
    expect_equal(KL2_dcc, KL2_ref)
    expect_equal(KL2_dnc, KL2_ref)

    # 1D example (one sample)
    Xd <- c(0,0,1,1,1)
    qd <- function(x) dbinom(x, size = 1, prob = 0.5)

    KLd_est <- kld_est(Xd, q = qd, vartype = "d")
    KLd_ref <- kld_est_discrete(Xd, q = qd)

    expect_equal(KLd_est, KLd_ref)

})


test_that("kld_est works as expected for numeric data", {

    # 1-D example (two samples)
    X1 <- c(1,1,2,2)
    Y1 <- c(1,2,2,2)

    KL1_est <- kld_est(X1, Y1)
    KL1_ref <- kld_est_nn(X1, Y1)

    expect_equal(KL1_est, KL1_ref)

    # 2-D example (two samples)
    X2m <- matrix(c(1,1,2,1,2,2),    ncol=2)
    Y2m <- matrix(c(1,1,2,2,1,2,1,2),ncol=2)
    X2d <- as.data.frame(X2m)
    Y2d <- as.data.frame(Y2m)

    KL2m_est <- kld_est(X2m, Y2m)
    KL2d_est <- kld_est(X2d, Y2d)
    KL2_ref  <- kld_est_nn(X2m, Y2m)

    expect_equal(KL2m_est,KL2_ref)
    expect_equal(KL2d_est,KL2_ref)

    # 1D example (one sample)
    Xc <- rnorm(10)
    qc <- dnorm

    KLc_est <- kld_est(Xc, q = qc)
    KLc_ref <- kld_est_nn(Xc, q = qc)

    expect_equal(KLc_est, KLc_ref)

})


test_that("kld_est works as expected for mixed data", {

    # check that heuristic for detecting column type works
    Xnnn <- data.frame(A = 1:5, B = 6:10, C = c(1,1,2,2,2))
    Ynnn <- data.frame(A = 11:14, B = 15:18, C = c(1,2,1,2))
    Xnnc <- Xnnn; Xnnc$C <- as.character(Xnnc$C)
    Ynnc <- Ynnn; Ynnc$C <- as.character(Ynnc$C)

    KLnnn_est <- kld_est(Xnnn, Ynnn, vartype = c("c","c","d"))
    KLnnc_est <- kld_est(Xnnc, Ynnc)

    expect_equal(KLnnn_est, KLnnc_est)

    # check that computed KL-D agrees with hardcoded mixed KL-D
    X1 <- Xnnn[Xnnn$C == 1, c("A","B")]
    X2 <- Xnnn[Xnnn$C == 2, c("A","B")]
    Y1 <- Ynnn[Ynnn$C == 1, c("A","B")]
    Y2 <- Ynnn[Ynnn$C == 2, c("A","B")]

    p1 <- mean(Xnnn$C == 1); p2 <- 1 - p1
    q1 <- mean(Ynnn$C == 1); q2 <- 1 - q1

    KL_ref <- p1*kld_est_nn(X1, Y1) + p2*kld_est_nn(X2, Y2) + kld_discrete(c(p1,p2),c(q1,q2))

    expect_equal(KLnnn_est,KL_ref)

    # 2D example, one sample
    X <- data.frame(A = rnorm(5),
                    B = c(0,0,1,1,1))
    q <- list(cond = function(xc,xd) dnorm(xc, mean = xd, sd = 1),
              disc = function(xd) dbinom(xd, size = 1, prob = 0.5))

    KL_Xq_est <- kld_est(X, q = q, vartype = c("c","d"))

    p0 <- mean(X$B == 0); p1 <- 1 - p0

    KL_Xq_ref <- p0 * kld_est_nn(X$A[X$B == 0], q = function(x) dnorm(x, mean = 0)) +
                 p1 * kld_est_nn(X$A[X$B == 1], q = function(x) dnorm(x, mean = 1)) +
                      kld_discrete(c(p0,p1), c(0.5,0.5))

    expect_equal(KL_Xq_est,KL_Xq_ref)

})


test_that("kld_est works for mixed data with missing level in X", {

    # check that heuristic for detecting column type works
    Xnnn <- data.frame(A = 1:5, B = 6:10, C = c(1,1,2,2,2))
    Ynnn <- data.frame(A = 11:16, B = 15:20, C = c(1,2,1,2,3,3))

    KLnnn_est <- kld_est(Xnnn, Ynnn, vartype = c("c","c","d"))

    # check that computed KL-D agrees with hardcoded mixed KL-D
    X1 <- Xnnn[Xnnn$C == 1, c("A","B")]
    X2 <- Xnnn[Xnnn$C == 2, c("A","B")]
    Y1 <- Ynnn[Ynnn$C == 1, c("A","B")]
    Y2 <- Ynnn[Ynnn$C == 2, c("A","B")]

    p1 <- mean(Xnnn$C == 1); p2 <- 1 - p1; p3 <- 0
    q1 <- mean(Ynnn$C == 1)
    q2 <- mean(Ynnn$C == 2)
    q3 <- mean(Ynnn$C == 2)

    KL_ref <- p1*kld_est_nn(X1, Y1) + p2*kld_est_nn(X2, Y2) + kld_discrete(c(p1,p2,p3),c(q1,q2,q3))

    expect_equal(KLnnn_est,KL_ref)

})


