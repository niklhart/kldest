
# Package `kldest`

This package contains tools for estimating Kullback-Leibler divergence based on univariate or multivariate samples. 

## Development plan

### Implementation

* [ ] theoretical investigation of rate of convergence based on sample size `n` (first in the balanced case, then in a more general setting: maybe take `min(n,m)`?). Subsampling bootstrap by Politis and Romano (1994), described [here](https://www.stat.umn.edu/geyer/5601/notes/sub.pdf)
* [ ] implement bias-reduced variant in discrete case (Zhang / Grabchak, 2014)
* [ ] migrate scripts to website

