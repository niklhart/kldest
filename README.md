
# Package `kldest`

This package contains tools for estimating Kullback-Leibler divergence based on univariate or multivariate samples. It is based on package `mimicvine` developed by Aleskandra Khatova in her master thesis.

## Development plan

### Implementation

* [ ] Density estimator in higher dimensions based on `ks::kde`; runtime comparison?
* [x] ~~k-NN estimator (replacing the 1-NN implementation, with default `k=l=1`)~~
* [ ] bootstrap for uncertainty of estimators
	* subsampling bootstrap by Politis and Romano (1994), described [here](https://www.stat.umn.edu/geyer/5601/notes/sub.pdf)
	* `TODO`: theoretical investigation of rate of convergence based on sample size `n` (first in the balanced case, then in a more general setting: maybe take `min(n,m)`?). For now, I use the square root law. From my understanding, this has proven to be the rate of convergence of different variants of nearest-neighbour based estimation.
* [ ] Implement discrete case (Zhang / Grabchak, 2014)? 


### Vignettes

* [x] ~~Vignette for 1-D examples (Gaussian, Uniform, Exponential), comparing different algorithms (density-based with/without MC, density-ratio based, nearest-neighbour-based), dependency of runtime and bias/variance on sample size~~
* [ ] ~~Vignette for 2-D examples (Gaussian, uncorrelated and correlated), comparing different algorithms (nearest-neighbour-based and `KernSmooth::bkde2D`), dependency of runtime and bias/variance on sample size~~
* [ ] Vignette for increasing dimensionality (nearest-neighbour-based)
* [ ] Vignette for evaluating of uncertainty quantification of estimators on different examples.

### Tests



### Documentation

* [ ] More detailed documentation of algorithms, including formulas and 


### To remember

* [ ] Think about the interface and possible scope of the package (e.g, other divergences?)
