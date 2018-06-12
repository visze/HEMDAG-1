\name{HEMDAG-package}
\alias{HEMDAG-package}
\alias{HEMDAG}
\docType{package}
\title{HEMDAG: Hierarchical Ensemble Methods for Directed Acyclic Graphs}

\description{
The HEMDAG package provides an implementation of several Hierarchical Ensemble Methods for DAGs. HEMDAG can be used to enhance the predictions of virtually any flat learning methods, by taking into account the hierarchical nature of the classes of a bio-ontology. HEMDAG is specifically designed for exploiting the hierarchical relationships of DAG-structured taxonomies, such as the Human Phenotype Ontology (HPO) or the Gene Ontology (GO), but it can be also safely applied to tree-structured taxonomies (as FunCat), since trees are DAGs. HEMDAG scale nicely both in terms of the complexity of the taxonomy and in the cardinality of the examples.
}

\details{
The HEMDAG package provides many utility functions to handle graph data structures and implements several Hierarchical Ensemble Methods for DAGs:
\enumerate{
\item \bold{HTD-DAG}: Hierarchical Top Down (\code{\link{HTD-DAG}});
\item \bold{TPR-DAG}: True-Path Rule (\code{\link{TPR-DAG}});
\item \bold{DESCENS}: Descendants Ensemble Classifier (\code{\link{DESCENS}});
\item \bold{MAX}, \bold{AND}, \bold{OR}: Heuristic Methods, \cite{Obozinski et al.} (\code{\link{Heuristic-Methods}});
}}

\author{
\emph{Marco Notaro} and \emph{Giorgio Valentini},
\href{https://sites.google.com/site/anacletolaboratory/}{AnacletoLab},
DI, Dipartimento di Informatica,
Universita' degli Studi di Milano

Maintainer: 
\emph{Marco Notaro} 
\email{marco.notaro@unimi.it}
}

\references{
Marco Notaro, Max Schubach, Peter N. Robinson and Giorgio Valentini, 
\emph{Prediction of Human Phenotype Ontology terms by means of Hierarchical Ensemble methods},
BMC Bioinformatics 2017, 18(1):449, doi:\href{https://doi.org/10.1186/s12859-017-1854-y}{10.1186/s12859-017-1854-y}
}

\keyword{package}
