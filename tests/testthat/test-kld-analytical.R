test_that("Discrete KL-D calculation works", {
    # 1D example
    P <- c(1,0)
    Q <- c(exp(-1),1-exp(-1))
    D_KL <- kl_div_discrete(P,Q)
    expect_equal(D_KL, 1)

    # 2D example
    P2 <- matrix(c(0.5,0,0,0.5),nrow=2)
    Q2 <- matrix(c(0.5*exp(-1),1-exp(-1),0,0.5*exp(-1)),nrow=2)
    D_KL2 <- kl_div_discrete(P2,Q2)
    expect_equal(D_KL2, 1)

    # Infinite KL-D
    P3 <- c(0.5,0.5)
    Q3 <- c(0,1)
    D_KL3 <- kl_div_discrete(P3,Q3)
    expect_equal(D_KL3, Inf)

})


