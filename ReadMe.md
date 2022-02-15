SubGroupSeparation
========

R package for marginalization in Bayesian networks: Integrating exact and approximate inference

Introduction
-----------

Bayesian networks have shown to be a valuable metric for modelling the underlying dependencies between variables of interest and to facilitate causal discovery. Encoding the dependencies amongst all variables, Bayesian networks are particularly suitable for dealing with missing data and prior knowledge. However, such computations require the calculation of the marginal probability distribution, which remains a fundamental problem for many statistical and scientific studies. 

`SubGroupSeparation` provides an efficient tool for calculating marginal probabilities in Bayesian networks that works in both low- and high-dimensional settings.

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
