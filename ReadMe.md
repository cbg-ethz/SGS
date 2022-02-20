SubGroupSeparation
========

R package for marginalization in Bayesian networks: exact and approximate inference

Introduction
-----------

Bayesian Networks are probabilistic graphical models that can compactly represent dependencies among random variables. This package provides an efficient method for marginalization in Bayesian networks, including both exact and approximate inference options. It can be used to deal with missing data or answer probabilistic queries.

`SubGroupSeparation` is efficient tool for calculating marginal probabilities in Bayesian networks that works in both low- and high-dimensional settings.

Installation
-----------

The latest development version of `SubGroupSeparation` can be found on GitHub
[here](https://github.com/cbg-ethz/SubGroupSeparation).

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/SubGroupSeparation`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is possible to use the `install_github`
tool from an R session:

```{r eval=FALSE}
library("devtools")
load_all()
```

`SubGroupSeparation` requires R `>= 3.5`, and depends on
`bitops`, `igraph`, `Matrix`, `graph` and
`methods`. Package `Rgraphviz` is requested in
order to plot graphs, but is not mandatory.


