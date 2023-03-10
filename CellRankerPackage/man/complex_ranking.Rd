% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CellRanker Functions.R
\name{complex_ranking}
\alias{complex_ranking}
\title{Helper function that ranks by the following:
Sensitivity - probability that a gene is expressed in given celltype
Specificity - probaility that gene is not expressed in other celltypes
Ubiquitousness index - probability that a gene is expressed at random,
 i.e a housekeeping gene}
\usage{
complex_ranking(
  celltype.list,
  pangao.list,
  ranking.list,
  weights = c(1, 1, 1),
  NA.penalty = 1e-04
)
}
\arguments{
\item{ranking.list}{A list vectors containing numeric rankings for each cluster.}

\item{weights}{Wieghts to change the importance of sensitivity, specificity, and ui}

\item{NA.penalty}{A penalty to incur on the ranking if the database does not contain
information on a given value. Necessary to prevent genes with little documentation
from being wieghted signficantly}

\item{pangao.listA}{list the same length as the number of clusters. Each list
entry ontains a sliced dataframe with information regarding each relevant marker
gene to the respective cluster.}
}
\description{
#'Helper function for celltype_ranking that intially ranks cell types
by number of occurances
#'@param celltype.list A list the same length as the number of clusters
that contains a vector of all possible celltyeps for each respective
cluster
}
