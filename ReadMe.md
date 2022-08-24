<div > 
  <img src="vignettes/sgs_icon.png" width="35%" height="35%">
</div>

Inference in Bayesian Networks
-----------
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`SGS` is an R package for inference in Bayesian networks. It allows for efficient exact and approximate inference that works both in low- and high-dimensional settings. Efficient marginalization is reached by splitting the calculation into sub-calculations of lower dimensionality. 
This code is an implementation of the paper [High-Dimensional Inference in Bayesian Networks](https://arxiv.org/abs/2112.09217), [arXiv preprint](https://arxiv.org/abs/2112.09217).

Implemented exact inference methods:
- SubGroupSeparation (fastest)
- Junction-tree algorithm
- Complete enumeration

Implemented approximate inference methods:
- SubGroupSeparation (highest accuracy)
- Loopy belief propagation
- Markov chain Monte Carlo (MCMC) sampling


Installation
-----------

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/SGS`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is possible to use the `install_github`
tool from an R session:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz", "RBGL"))

library("devtools")
install_github("cbg-ethz/SGS")
```

The packages "graph", "Rgraphviz" and "RBGL" need to be installed from BioConductor, as they are not hosted on CRAN.

`SGS` requires R `>= 3.5`, and depends on
`bitops` and
`methods`. Other packages are requested in
order to plot graphs, but are not mandatory.


Examples
-------

```
library(SGS)

# create BN and label variables 
set.seed(6)
myBayesNet <- randomBN(3)
myBayesNet@variables <- c("rain", "sprinkler", "wet grass")
plot_bn(myBayesNet)

# what's the probability of having rain and wet grass at the same time?
# define observed variables and calculate marginal probability
myObserved <- list(observed.vars=c("rain", "wet grass"), observed.vals=c(2,2))
exactInference(myBayesNet,myObserved)

# another example: 
# let's learn the Bayesian network from the "Asia dataset"
asia_bn <- learn_bn(Asia)
plot_bn(asia_bn)

# now we can do the inference on the learned Bayesian network
myObserved <- list(observed.vars=c("X", "D"), observed.vals=c(1,1))
exactInference(asia_bn, myObserved)
```

Benchmark Results 
-------

We benchmarked the performance of our SGS method against standard inference schemes (Gibbs sampling and loopy belief propagation) over a broad range of different Bayesian networks. The results are summarized in the Figure below (lower is better), displaying the normalized root mean squared error (NRMSE). To reproduce the results, run the scripts in the benchmark folder. 

# ![SGS](https://github.com/cbg-ethz/SGS/blob/master/vignettes/benchmark.png)

Reference
---------

If you find this code useful, please consider citing:

Fritz M. Bayer, Giusi Moffa, Niko Beerenwinkel, Jack Kuipers. [High-Dimensional Inference in Bayesian Networks](https://arxiv.org/abs/2112.09217), [arXiv preprint](https://arxiv.org/abs/2112.09217), 2021

```
@article{bayer2021marginalization,
  title={High-Dimensional Inference in Bayesian Networks},
  author={Bayer, Fritz M and Moffa, Giusi and Beerenwinkel, Niko and Kuipers, Jack},
  journal={arXiv preprint arXiv:2112.09217},
  year={2021}
}
```
