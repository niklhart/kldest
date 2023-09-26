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

test_that("Function 'mvdnorm' works as expected", {

    # 1D example
    tst1 <- mvdnorm(x = 2, mu = 1, Sigma = 2)
    ref1 <-   dnorm(x = 2, mean = 1, sd = sqrt(2))
    expect_equal(tst1,ref1)

    # Independent 2D example
    tst2 <- mvdnorm(x = c(2,2), mu = c(1,1), Sigma = diag(1:2))
    ref2 <- prod(dnorm(x = c(2,2), mean = c(1,1), sd = sqrt(1:2)))
    expect_equal(tst2,ref2)

    # Correlated 2D example (determinant 1, evaluated at center)
    tst3 <- mvdnorm(x = c(0,0), mu = c(0,0), Sigma = matrix(c(2,1,1,1),nrow=2))
    ref3 <- 1/(2*pi)
    expect_equal(tst3,ref3)

    # Error for singular covariance matrix
    expect_error(mvdnorm(x = c(2,2), mu = c(1,1), Sigma = matrix(1,nrow=2,ncol=2)))

})


test_that("function `is_two_sample` works as expected", {

    expect_true(is_two_sample(Y = 1, q = NULL))
    expect_true(is_two_sample(Y = matrix(1), q = NULL))
    expect_true(is_two_sample(Y = data.frame(1), q = NULL))

    expect_false(is_two_sample(Y = NULL, q = force))

    expect_error(is_two_sample(Y = NULL,  q = NULL))
    expect_error(is_two_sample(Y = 1,     q = force))
    expect_error(is_two_sample(Y = force, q = NULL))

})
