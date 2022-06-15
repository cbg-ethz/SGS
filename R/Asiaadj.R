#' Asiamat
#'
#' An adjacency matrix representing the ground truth DAG used to generate a synthetic dataset from 
#' Lauritzen and Spiegelhalter (1988) about lung
#' diseases (tuberculosis, lung cancer or bronchitis) and visits to Asia.
#' 
#' @source \url{https://www.bnlearn.com/bnrepository/}
#' @format A binary matrix with 8 rows and 8 columns representing an adjacency matrix of a DAG with 8 nodes:
#' \itemize{
#' \item D (dyspnoea), binary 1/0 corresponding to "yes" and "no"
#' \item T (tuberculosis), binary 1/0 corresponding to "yes" and "no"
#' \item L (lung cancer), binary 1/0 corresponding to "yes" and "no"
#' \item B (bronchitis), binary 1/0 corresponding to "yes" and "no"
#' \item A (visit to Asia), binary 1/0 corresponding to "yes" and "no"
#' \item S (smoking), binary 1/0 corresponding to "yes" and "no"
#' \item X (chest X-ray), binary 1/0 corresponding to "yes" and "no"
#' \item E (tuberculosis versus lung cancer/bronchitis), binary 1/0 corresponding to "yes" and "no"
#' }
#'@references Lauritzen S, Spiegelhalter D (1988). `Local Computation with Probabilities on Graphical Structures and their Application to Expert Systems (with discussion)'.
#'Journal of the Royal Statistical Society: Series B 50, 157-224.
#'
"Asiamat"

#' Asia dataset
#'
#' A synthetic dataset from Lauritzen and Spiegelhalter (1988) about lung
#' diseases (tuberculosis, lung cancer or bronchitis) and visits to Asia.
#'
#' @source \url{https://www.bnlearn.com/bnrepository/}
#' @format A data frame with 5000 rows and 8 binary variables:
#' \itemize{
#' \item D (dyspnoea), binary 1/0 corresponding to "yes" and "no"
#' \item T (tuberculosis), binary 1/0 corresponding to "yes" and "no"
#' \item L (lung cancer), binary 1/0 corresponding to "yes" and "no"
#' \item B (bronchitis), binary 1/0 corresponding to "yes" and "no"
#' \item A (visit to Asia), binary 1/0 corresponding to "yes" and "no"
#' \item S (smoking), binary 1/0 corresponding to "yes" and "no"
#' \item X (chest X-ray), binary 1/0 corresponding to "yes" and "no"
#' \item E (tuberculosis versus lung cancer/bronchitis), binary 1/0 corresponding to "yes" and "no"
#' }
#'@references Lauritzen S, Spiegelhalter D (1988). `Local Computation with Probabilities on Graphical Structures and their Application to Expert Systems (with discussion)'.
#'Journal of the Royal Statistical Society: Series B 50, 157-224.
#'
"Asia"