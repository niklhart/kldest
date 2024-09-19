# Neural KL divergence estimation

#' Neural KL divergence estimation (Donsker-Varadhan representation) using `torch`
#'
#' Estimation of KL divergence between continuous distributions based on the
#' Donsker-Varadhan representation
#' \deqn{D_{KL}(P||Q) = \sup_{f} E_P[f(X)] - \log\left(E[e^{f(X)}]\right)}
#' using Monte Carlo averages to approximate the expectations, and optimizing
#' over a class of neural networks. The `torch` package is required to use this
#' function.
#'
#' Disclaimer: this is a simple test implementation which is not optimized by
#' any means. In particular:
#' - it uses a fully connected network with (only) a single hidden layer
#' - it uses standard gradient descient on the full dataset and not more advanced
#'     estimators
#' Also, he syntax is likely to change in the future.
#'
#' Estimation is done as described for mutual information in Belghazi et al.
#' (see ref. below), except that standard gradient descent is used on the full
#' samples X and Y instead of using batches. Indeed, in the case where X and Y
#' have a different length, batch sampling is not that straightforward.
#'
#' Reference: Belghazi et al., Mutual Information Neural Estimation,
#'            PMLR 80:531-540, 2018.
#'
#' @param X,Y  `n`-by-`d` and `m`-by-`d` numeric matrices, representing `n`
#'    samples from the true distribution \eqn{P} and `m` samples from the
#'    approximate distribution \eqn{Q}, both in `d` dimensions. Vector input is
#'    treated as a column matrix.
#' @param d_hidden Number of nodes in hidden layer (default: `32`)
#' @param learning_rate Learning rate during gradient descent (default: `1e-4`)
#' @param epochs Number of training epochs (default: `200`)
#' @param device Calculation device, either `"cpu"` (default), `"cuda"` or `"mps"`.
#' @param verbose Generate progress report to consolue during training of the
#'     neutral network (default: `FALSE`)?
#' @returns A scalar, the estimated Kullback-Leibler divergence \eqn{\hat D_{KL}(P||Q)}.
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
                           epochs = 5000, device = c("cpu","cuda","mps"),
                           verbose = FALSE) {

    # error if package torch is not installed
    if (!requireNamespace("torch", quietly = TRUE)) {
        stop("Package 'torch' is required to use function 'kldest::kld_est_neural'.")
    }

    # check device, fallback to "cpu"
    device <- switch(
        match.arg(device),
        cuda = if (torch::cuda_is_available()) "cuda" else {
            message("cuda unavailable, using cpu instead."); "cpu"
        },
        mps  = if (torch::backends_mps_is_available()) "mps" else {
            message("mps unavailable, using cpu instead."); "cpu"
        },
        cpu  = "cpu"
    )

    # transform input to matrix format and check that it is numeric
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    stopifnot(is.numeric(X), is.numeric(Y))

    # Dimensions
    d_in <- ncol(X) # number of dimensions, must be the same in X and Y
    n    <- nrow(X) # number of samples in X
    m    <- nrow(Y) # number of samples in Y
    d_out <- 1      # single output

    # check consistency of dimensions
    stopifnot(ncol(Y) == d_in)

    # turn X and Y into torch tensors
    X <- torch::torch_tensor(X)#$to(device = device)
    Y <- torch::torch_tensor(Y)#$to(device = device)

    # weights connecting input to hidden / hidden to output layer
    w1 <- torch::torch_randn(d_in, d_hidden, requires_grad = TRUE)#$to(device = device)
    w2 <- torch::torch_randn(d_hidden, d_out, requires_grad = TRUE)#$to(device = device)

    # hidden / output layer bias
    b1 <- torch::torch_zeros(1, d_hidden, requires_grad = TRUE)#$to(device = device)
    b2 <- torch::torch_zeros(1, d_out, requires_grad = TRUE)#$to(device = device)

    ### training loop ----------------------------------------

    for (t in 1:epochs) {

        ### -------- Forward pass --------
        x_pred <- X$mm(w1)$add(b1)$relu()$mm(w2)$add(b2)
        y_pred <- Y$mm(w1)$add(b1)$relu()$mm(w2)$add(b2)

        ### -------- Compute loss --------
        # L(theta) = log(mean(exp(f(Y,theta)))) - mean(f(X,theta))
        loss <- y_pred$exp()$mean()$log()$subtract(x_pred$mean())


        if (t == 1 && is.nan(loss$item())) {
            message("Function cannot be evaluated at the initial point, please try again.")
            break()
        }

        if (verbose && (t %% 10 == 0) )
            cat("Epoch: ", t, "   Loss: ", loss$item(), "\n")

        ### -------- Backpropagation --------

        # compute gradient of loss w.r.t. all tensors with
        # requires_grad = TRUE
        loss$backward()

        ### -------- Update weights --------

        # Wrap in with_no_grad() because this is a part we don't
        # want to record for automatic gradient computation
        torch::with_no_grad({
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

    # Estimated KL divergence = negative loss of optimized network
    -loss$item()
}
