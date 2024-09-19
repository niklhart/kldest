# experimental features not included into the current version.


#' Nearest-neighbour density ratio based KL divergence estimator.
#'
#' Currently, this estimator is not part of the kldest package. Its performance
#' for small samples is poor, since a single point in `Y` for which no point in
#' `X` is closer than the `k`-th nearest neighbour in `Y` suffices for the KL
#' divergence estimate to be `-Inf`.
#'
#' @inherit kld_est_kde return
#' @inheritParams kld_est_nn
#' @param k Number of nearest neighbours to consider for NN density estimation.
#'    Defaults to `k = neff^(1/(d+1))`, where `neff = min(n,m)`. The choice of
#'    `k` is quite important for accuracy and the default may not be optimal.
#' @examples
#' # 1D example
#' set.seed(0)
#' X <- rnorm(100)
#' Y <- rnorm(100, mean = 1)
#' kld_est_nndr(X,Y)
kld_est_nndr <- function(X, Y, k = NULL, eps = 0) {

    # get important dimensions
    X <- as.matrix(X)
    Y <- as.matrix(Y)  # number of dimensions must be the same in X and Y
    d <- ncol(X)       # number of dimensions
    n <- nrow(X)       # number of samples in X
    m <- nrow(Y)       # number of samples in Y

    # default number of nearest neighbours (Noshad et al., 2017)
    if (is.null(k)) k <- ceiling(min(n,m)^(1/(d+1)))

    # combined dataset
    Z <- rbind(X,Y)

    idx.ZY <- RANN::nn2(Z, Y, k = k, eps = eps)$nn.idx
    Mi <- rowSums(idx.ZY > n)
    Ni <- k - Mi

    # return
    -1/n * sum(log(Ni / (Mi + 1)))

}




#' Zhang/Grabchak KL divergence estimator for samples from discrete distributions
#'
#' CAVE: not implemented yet!
#'
#' The estimator is that from Eq. (1.3) from Zhang and Grabchak (2014).
#'
#' Reference:
#' Zhang and Grabchak, "Nonparametric Estimation of Kullback-Leibler Divergence",
#' Neural Computation 26, 2570–2593 (2014).
#'
#' @inherit kld_est_kde return
#' @param X,Y Two samples from discrete distributions, specified as vectors,
#'    matrices or data frames.
#' @examples
#' # 1D example
#' X <- c(rep('M',5),rep('F',5))
#' Y <- c(rep('M',6),rep('F',4))
#' kld_est_discrete(X, Y)
#' @export
kld_est_zg2014 <- function(X, Y) {

    # get important dimensions
    X <- as.data.frame(X)
    Y <- as.data.frame(Y)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d)

    # ensure consistency of factor levels
    uXY <- if (d == 1) {
        list(unique(c(X[[1]],Y[[1]])))
    } else {
        mapply(function(x,y) unique(c(x,y)),
               lapply(X, unique),
               lapply(Y, unique))
    }

    # convert all columns to factors to correctly account for missing levels
    X[] <- lapply(1:d, function(k) factor(X[[k]], levels = uXY[[k]]))
    Y[] <- lapply(1:d, function(k) factor(Y[[k]], levels = uXY[[k]]))

    # frequency counts via tables
    tX <- table(X)
    tY <- table(Y)

    pX <- stop("Implement this!")

    # estimate KL divergencce via relative frequencies
    kl_div_discrete(tX/n, tY/m)
}


#' Generalized k-nearest neighbour KL divergence estimator.
#'
#' This function implements the generalized k-nearest neighbour estimator in
#' Wang et al. (2009), Eq.(17). In this estimator, the number of nearest
#' neighbours to consider may differ between samples and sample points.
#'
#' The method is fully functional, but since it is not really useful in practice,
#' I haven't included it in the package.
#'
#' Reference:
#' Wang, Kulkarni and Verdú, "Divergence Estimation for Multidimensional
#' Densities Via k-Nearest-Neighbor Distances", IEEE Transactions on Information
#' Theory, Vol. 55, No. 5 (2009).
#'
#' @inherit kld_est_nn params return
#' @param l,k Scalars or numeric vectors of length `n` and `m`, representing the
#'   number of nearest neighbours to use for nearest neighbour density estimation
#'   of P and Q, respectively, for each of the data points \code{X[i,]}, with
#'   \code{i} ranging from \code{1} to \code{n}, and \code{Y[j,]}, with
#'   \code{j} ranging from \code{1} to \code{m}. The default is `l = k = 1`.
#'   In the special case that `l = k` and `k` is scalar, the estimator coincides
#'   with `kld_est_nn(X, ..., k = k).
#'
kld_est_gnn <- function(X, Y = NULL, q = NULL, l = k, k = 1, eps = 0, log.q = FALSE) {

    # get important dimensions
    X <- as.matrix(X)
    d <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X

    # check validity of input l
    if (length(l) == 1) l <- rep(l, n)
    stopifnot(length(l) == n, max(l) < n)

    # nearest neighbours within X (delete the NN in X to X, which is X itself)
    nnXX <- RANN::nn2(X, X, k = max(l)+1, eps = eps)$nn.dists[ ,-1, drop = FALSE]

    # distance to the l-th nearest neighbour
    rho_l <- vapply(1:n, function(i) nnXX[i,l[i]], 1)

    # check validity of input: one- or two-sample problem?
    two.sample <- is_two_sample(Y,q)

    if (two.sample) {
        # two-sample problem
        Y <- as.matrix(Y)
        m <- nrow(Y) # number of samples in Y

        # check validity of input k
        if (length(k) == 1) k <- rep(k, n)
        stopifnot(ncol(Y) == d, length(k) == n, max(k) <= m)

        # nearest neighbours to X within Y
        nnYX <- RANN::nn2(Y, X, k = max(k), eps = eps)$nn.dists

        # distance to the k-th nearest neighbour
        nu_k <- vapply(1:n, function(i) nnYX[i,k[i]], 1)

        # equation (17) from Wang et al.
        return(d*mean(log(nu_k/rho_l)) + mean(digamma(l) - digamma(k)) + log(m/(n-1)))

    } else {
        # one-sample problem
        log_phat_X <- -log(n-1) + lgamma(0.5*d+1) - 0.5*d*log(pi) - d*log(rho_l) + digamma(l)
        log_q_X    <- if (log.q) {
            apply(X, MARGIN = 1, FUN = q)
        } else {
            log(apply(X, MARGIN = 1, FUN = q))
        }
        return(mean(log_phat_X-log_q_X))
    }

}





#' Neural KL divergence estimation (Donsker-Varadhan representation) using `torch`
#'
#' Disclaimer: this is a simple test implementation which is not optimized by
#' any means. In particular:
#' - it only has a single hidden layer
#' - it uses standard gradient descient on the full dataset
#'
#' Estimation is done as described for mutual information in Belghazi et al.
#' (see ref. below), except that standard gradient descent is used on the full
#' samples X and Y instead of using batches. Indeed, in the case where X and Y
#' have a different length, batch sampling is not that straightforward. Network
#' architecture is a fully connected network with a single hidden layer.
#'
#' Reference: Belghazi et al., Mutual Information Neural Estimation,
#'            PMLR 80:531-540, 2018.
#'
#' @param d_hidden Number of nodes in hidden layer (default: `32`)
#' @param learning_rate Learning rate during gradient descent (default: `1e-4`)
#' @param epochs Number of training epochs (default: `200`)
#' @param device Calculation device, either `"cpu"` (default), `"cuda"` or `"mps"`.
#'
#' @examples
#' # 2D example
#' # analytical solution
#' kld_gaussian(mu1 = rep(0,2), sigma1 = diag(2),
#'              mu2 = rep(0,2), sigma2 = matrix(c(1,1,1,2),nrow=2))
#' # sample generation
#' set.seed(0)
#' nxy <- 1000
#' X1 <- rnorm(nxy)
#' X2 <- rnorm(nxy)
#' Y1 <- rnorm(nxy)
#' Y2 <- Y1 + rnorm(nxy)
#' X <- cbind(X1,X2)
#' Y <- cbind(Y1,Y2)
#' # Estimation
#' kld_est_nn(X, Y)
#' kld_est_neural(X, Y)
#' @export
kld_est_neural <- function(X, Y, d_hidden = 1024, learning_rate = 1e-4,
                           epochs = 5000, device = "cpu") {

    # Input dimensionality
    d_in <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d_in)

    # turn X and Y into torch tensors
    X <- torch::torch_tensor(X)
    Y <- torch::torch_tensor(Y)

    # output dimensionality (number of predicted features)
    d_out <- 1

    # weights connecting input to hidden layer
    w1 <- torch::torch_randn(d_in, d_hidden, requires_grad = TRUE)
    # weights connecting hidden to output layer
    w2 <- torch::torch_randn(d_hidden, d_out, requires_grad = TRUE)

    # hidden layer bias
    b1 <- torch::torch_zeros(1, d_hidden, requires_grad = TRUE)
    # output layer bias
    b2 <- torch::torch_zeros(1, d_out, requires_grad = TRUE)

    ### training loop ----------------------------------------

    for (t in 1:epochs) {

        ### -------- Forward pass --------

        x_pred <- X$mm(w1)$add(b1)$relu()$mm(w2)$add(b2)
        y_pred <- Y$mm(w1)$add(b1)$relu()$mm(w2)$add(b2)

        ### -------- Compute loss --------
        # L(theta) = log(mean(exp(f(Y,theta)))) - mean(f(X,theta))
        loss <- y_pred$exp()$mean()$log()$subtract(x_pred$mean())

        if (t %% 10 == 0)
            cat("Epoch: ", t, "   Loss: ", loss$item(), "\n")

        ### -------- Backpropagation --------

        # compute gradient of loss w.r.t. all tensors with
        # requires_grad = TRUE
        loss$backward()

        ### -------- Update weights --------

        # Wrap in with_no_grad() because this is a part we don't
        # want to record for automatic gradient computation
        with_no_grad({
            w1 <- w1$sub_(learning_rate * w1$grad)
            w2 <- w2$sub_(learning_rate * w2$grad)
            b1 <- b1$sub_(learning_rate * b1$grad)
            b2 <- b2$sub_(learning_rate * b2$grad)

            # Zero gradients after every pass, as they'd
            # accumulate otherwise
            w1$grad$zero_()
            w2$grad$zero_()
            b1$grad$zero_()
            b2$grad$zero_()
        })

    }


    # ds <- torch::dataset(
    #     name = "samples",
    #
    #     initialize = function(idx) {
    #         self$x <- torch::torch_tensor(X)
    #         self$y <- torch::torch_tensor(Y)
    #     },
    #
    #     .getbatch = function(idx) {
    #         list(x = self$x[idx, , drop = FALSE],
    #              y = self$y[idx])
    #     },
    #
    #     .getitem = function(i) {
    #         list(x = self$x[i, ],
    #              y = self$y[i])
    #     },
    #
    #     .length = function() {
    #         self$y$size()[[1]]
    #     }
    # )
    #
    # dl <- torch::dataloader(ds, batch_size = 64, shuffle = TRUE)
    #
    # net <- torch::nn_module(
    #
    #     "easy mlp",
    #
    #     initialize = function() {
    #
    #         self$fc1 <- torch::nn_linear(in_features = length(X) + length(Y),
    #                                      out_features = 128)
    #         self$fc2 <- torch::nn_linear(in_features = 128,
    #                                      out_features = 1)
    #
    #     },
    #
    #     forward = function(x) {
    #
    #         x |>
    #             self$fc1() |>
    #             torch::nnf_relu() |>
    #             self$fc2() |>
    #             torch::torch_flatten()
    #
    #     }
    # )
    #
    # loss <- stop("TODO!!")
    #
    # fitted <- net |>
    #     setup(
    #         loss = nnf_mse_loss,
    #         optimizer = optim_adam
    #     ) |>
    #     set_opt_hparams(lr = 0.0001) |>
    #     fit(train_dl, epochs = 50, valid_data = valid_dl)


}


#' Neural Estimation of KL divergence (modularized)
#'
#' @param d_hidden Number of nodes in hidden layer (default: `32`)
#' @param learning_rate Learning rate during gradient descent (default: `1e-4`)
#' @param epochs Number of training epochs (default: `200`)
#' @param device Calculation device, either `"cpu"` (default), `"cuda"` or `"mps"`.
#' Reference:
#' Belghazi et al., Mutual Information Neural Estimation, PMLR 80:531-540, 2018.
kld_est_neural2 <- function(X, Y, d_hidden = 32, learning_rate = 1e-4,
                           epochs = 200, device = "cpu") {


    # FOR TESTING ONLY, REMOVE LATER
    set.seed(0)
    nxy <- 1000
    X1 <- rnorm(nxy)
    X2 <- rnorm(nxy)
    Y1 <- rnorm(nxy)
    Y2 <- Y1 + rnorm(nxy)
    X <- cbind(X1,X2)
    Y <- cbind(Y1,Y2)

    d_hidden <- 1024
    learning_rate <- 1e-4
    epochs <- 5000

    # Input dimensionality
    d_in <- ncol(X) # number of dimensions, must be the same in X and Y
    n <- nrow(X) # number of samples in X
    m <- nrow(Y) # number of samples in Y

    # check consistency of dimensions
    stopifnot(ncol(Y) == d_in)

    # turn X and Y into torch tensors
    X <- torch::torch_tensor(X)
    Y <- torch::torch_tensor(Y)

    # neural network structure
    net <- nn_module(
        initialize = function(d_in,d_hidden) {

            self$fc1 <- nn_linear(in_features = d_in,     out_features = d_hidden)
            self$fc2 <- nn_linear(in_features = d_hidden, out_features = 1)

        },
        forward = function(X,Y) {
            Xp <- X |>
                self$fc1() |>
                torch::nnf_relu() |>
                self$fc2() |>
                torch::torch_flatten()
            Yp <- Y |>
                self$fc1() |>
                torch::nnf_relu() |>
                self$fc2() |>
                torch::torch_flatten()
            torch::torch_cat(Xp,Yp)
        }
    )

    # loss function
    loss <- function(input, target) {
        fx <- input[1:n,]
        fy <- input[n+(1:m),]
        fy$exp()$mean()$log()$subtract(fx$mean())
    }

    # optimizer

    # Here, I'd like to use the ADAM optimizer. However, since X and Y are
    # concatenated, I cannot simply batch sample from c(X,Y)! What to do in this
    # case?


}
