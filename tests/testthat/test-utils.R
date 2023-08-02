test_that("Trapezoidal integration works", {
    # 1D example (linear)
    int_1d <- trapz(h=1, fx = 1:4)
    expect_equal(int_1d, 7.5)

    # 2D example (linear)
    f <- function(x,y) x + 2*y
    x <- 0:3
    y <- 0:3
    fx <- outer(x,y,f)
    int_2d <- trapz(h=c(1,1), fx = fx)
    expect_equal(int_2d, 40.5)

})

test_that("constDiagMatrix works as expected", {

    # 1D example (edge case)
    M <- constDiagMatrix(diag = 3)
    R <- matrix(3)
    expect_equal(M, R)

    # 3D example
    M <- constDiagMatrix(dim = 3, diag = 3, offDiag = 1)
    R <- matrix(c(3,1,1,
                  1,3,1,
                  1,1,3), nrow = 3)
    expect_equal(M, R)
})


test_that("Matrix trace works", {

    M <- matrix(c(1,2,3,
                  1,2,3,
                  1,2,3), nrow = 3)
    expect_equal(kldest:::tr(M), 6)
})

test_that("Function 'combinations' works as expected", {

    comb <- combinations(a = c('a','b'), b = c('a','c','d'), c = 0:1)
    ref  <- data.frame(
        a = c('a','b','a','b','a','b','a','b','a','b','a','b'),
        b = c('a','a','c','c','d','d','a','a','c','c','d','d'),
        c = rep(0:1, each = 6)
    )
    expect_equal(comb,ref)
})

