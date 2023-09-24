test_that("kld_cont works as expected", {

    X <- 1:2
    Y <- X+1      # same SD as X, and m=n

    KL_nn_ref  <- kld_est_nn(X, Y)
    KL_nn_test <- kld_cont(X, Y, method = "nn")

    KL_brnn_ref  <- kld_est_brnn(X, Y)
    KL_brnn_test <- kld_cont(X, Y, method = "brnn")

    KL_kde1_ref  <- kld_est_kde1(X, Y)
    KL_kde1_test <- kld_cont(X, Y, method = "kde1")

    expect_equal(KL_nn_test,  KL_nn_ref)
    expect_equal(KL_brnn_test,KL_brnn_ref)
    expect_equal(KL_kde1_test,KL_kde1_ref)

})


test_that("kld_est works as expected for discrete data", {

    # 1D example: invariant to type of input
    Xn <- c(1,1,2,2)
    Yn <- c(1,2,2,2)
    Xc <- as.character(Xn)
    Yc <- as.character(Yn)
    Xnd <- data.frame(Xn)
    Ynd <- data.frame(Yn)
    Xcd <- data.frame(Xc)
    Ycd <- data.frame(Yc)

    KL_num <- kld_est(Xn, Yn, vartype = "d")
    KL_chr <- kld_est(Xc, Yc)
    KL_ndf  <- kld_est(Xnd, Ynd, vartype = "d")
    KL_cdf  <- kld_est(Xcd, Ycd)
    KL_ref <- kld_discrete(c(0.5,0.5), c(0.25,0.75))

    expect_equal(KL_num,KL_ref)
    expect_equal(KL_chr,KL_ref)
    expect_equal(KL_ndf, KL_ref)
    expect_equal(KL_cdf, KL_ref)

    # 2D example
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

    expect_equal(KL2_num,KL2_ref)
    expect_equal(KL2_chr,KL2_ref)
    expect_equal(KL2_dnn,KL2_ref)
    expect_equal(KL2_dcc,KL2_ref)
    expect_equal(KL2_dnc,KL2_ref)

})


test_that("kld_est works as expected for numeric data", {

    # 1-D test
    X1 <- c(1,1,2,2)
    Y1 <- c(1,2,2,2)

    KL1_est  <- kld_est(X1, Y1)
    KL1_cont <- kld_cont(X1, Y1)

    expect_equal(KL1_est,KL1_cont)

    # 2-D test
    X2m <- matrix(c(1,1,2,1,2,2),    ncol=2)
    Y2m <- matrix(c(1,1,2,2,1,2,1,2),ncol=2)
    X2d <- as.data.frame(X2m)
    Y2d <- as.data.frame(Y2m)

    KL2m_est  <- kld_est(X2m, Y2m)
    KL2d_est  <- kld_est(X2d, Y2d)
    KL2_cont <- kld_cont(X2m, Y2m)

    expect_equal(KL2m_est,KL2_cont)
    expect_equal(KL2d_est,KL2_cont)

})


test_that("kld_est works as expected for mixed data", {

    # check that heuristic for detecting column type works
    Xnn <- data.frame(A = c(1,1,1,2,2),B = c(1,1,2,2,2))
    Ynn <- data.frame(A = c(1,1,2,2), B = c(1,2,1,2))
    Xnc <- Xnn; Xnc$B <- as.character(Xnc$B)
    Ync <- Ynn; Ync$B <- as.character(Ync$B)

    KLnn_est <- kld_est(Xnn, Ynn, vartype = c("c","d"))
    KLnc_est <- kld_est(Xnc, Ync)

    expect_equal(KLnn_est, KLnc_est)

    # check that computed KL-D agrees with hardcoded mixed KL-D
    X1 <- Xnn$A[Xnn$B == 1]
    X2 <- Xnn$A[Xnn$B == 2]
    Y1 <- Ynn$A[Ynn$B == 1]
    Y2 <- Ynn$A[Ynn$B == 2]

    p1 <- mean(Xnn$B == 1); p2 <- 1 - p1
    q1 <- mean(Ynn$B == 1); q2 <- 1 - q1

    KL_ref <- p1*kld_cont(X1, Y1) + p2*kld_cont(X2, Y2) + kld_discrete(c(p1,p2),c(q1,q2))

    expect_equal(KLnn_est,KL_ref)

})

