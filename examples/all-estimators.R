# KL-D between two samples from 1D Gaussians:
X <- rnorm(1000)
Y <- rnorm(1000, mean = 1, sd = 2)
kl_div_gaussian(mu1 = 0, sigma1 = 1, mu2 = 1, sigma2 = 2^2)
kldest_density(X, Y)
kl_universal_1nn(X, Y)
kl_generalized_knn_eps(X, Y)
