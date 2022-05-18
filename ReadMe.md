Efficient Exact and Approximate Inference in Bayesian Networks
-----------
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`SubGroupSeparation` is an R package for marginalization in Bayesian networks. It allows for efficient exact and approximate inference that works both in low- and high-dimensional settings. As illustrated below, efficient marginalization is reached by splitting the calculation into sub-calculations of lower dimensionality.

Installation
-----------

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/SubGroupSeparation`
from a terminal, or `make install` from within the package source folder.

`SubGroupSeparation` requires R `>= 3.5`, and depends on
`bitops` and
`methods`. Other packages are requested in
order to plot graphs, but are not mandatory.


Example 
-------

```{r eval=FALSE}
library(SubGroupSeparation)

# create BN and label variables 
set.seed(6)
myBayesNet <- randomBN(3)
myBayesNet@variables <- c("rain", "sprinkler", "wet grass")
plot(myBayesNet)

# what's the probability of having rain and wet grass at the same time?
# define observed variables and calculate marginal probability
myObserved <- list(observed.vars=c("rain", "wet grass"), observed.vals=c(2,2))
exactInference(myBayesNet,myObserved)
```

Benchmark Results 
-------

# ![SGS](https://github.com/cbg-ethz/SubGroupSeparation/blob/master/vignettes/figures/benchmark.png)
