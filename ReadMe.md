SubGroupSeparation
========

R package for marginalization in Bayesian networks: Integrating exact and approximate inference

Introduction
-----------

Bayesian Networks are probabilistic graphical models that can compactly represent dependencies among random variables. This package allows for efficient margainlization in Bayesian networks with exact and approximate inference options. It can be used to handle missing data or to answer probabilistic quieries. 

`SubGroupSeparation` is an efficient tool for calculating marginal probabilities in Bayesian networks that works in both low- and high-dimensional settings.

Installation
-----------

The latest development version of `SubGroupSeparation` can be found on GitHub
[here](https://github.com/cbg-ethz/SubGroupSeparation).

In order to install the package, it suffices to launch
`R CMD INSTALL path/to/SubGroupSeparation`
from a terminal, or `make install` from within the package source folder.

Being hosted on GitHub, it is also possible to use the `install_github`
tool from an R session:

```{r eval=FALSE}
library("devtools")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install_github("cbg-ethz/SubGroupSeparation")
```

`SubGroupSeparation` requires R `>= 3.5`, and depends on
`bitops`, `igraph`, `Matrix`, `graph` and
`methods`. Package `Rgraphviz` is requested in
order to plot graphs, but is not mandatory.

Reference
---------
If you use `SubGroupSeparation` in your work, please cite it as:
```
Fritz M. Bayer, Giusi Moffa, Niko Beerenwinkel, Jack Kuipers. "Marginalization in Bayesian Networks: Integrating Exact and Approximate Inference" arXiv preprint, 2021
```
