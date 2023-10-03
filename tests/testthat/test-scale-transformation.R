test_that("Transformation to uniform scale works on a hardcoded example", {

    # 1D example
    X <- matrix(c(1,2,4))
    Y <- matrix(c(3,6))
    L <- to_uniform_scale(X, Y)

    uXY <- do.call(rbind,L)
    su <- 1:5 / 5

    expect_equal(sort(uXY), su)

    # 2D example
    X2 <- matrix(c(1,2,3,4),nrow=2)
    Y2 <- matrix(c(3,4,1,2),nrow=2)
    L2 <- to_uniform_scale(X2, Y2)

    uXY2  <- do.call(rbind,L2)
    uXY2s <- apply(uXY2, MARGIN = 2, FUN = sort)
    su2 <- 1:4 / 4

    expect_equal(uXY2s, matrix(su2,nrow=4,ncol=2))

})

